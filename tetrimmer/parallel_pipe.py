import multiprocessing
import multiprocessing.connection

def worker_wrapper(worker_id, worker_task, args_queue, pipe):
    """Module-level wrapper so nothing captures `self` (avoids pickling the manager)."""
    def send_and_receive(result):
        pipe.send(result)
        return pipe.recv()
    try:
        worker_task(worker_id, args_queue, send_and_receive)
    finally:
        try:
            pipe.close()
        except Exception:
            pass

class ChattyParallelProcessor:
    """
    Chatty multiprocessing manager with per-worker queues + a Pipe for request/response.
    - worker_task MUST be a module-level callable: worker_task(worker_id, queue, send_and_receive)
    - args is a list of picklable task-argument lists/tuples. Heavy/non-picklable objects must
      be created inside the worker, not passed in args.
    """
    def __init__(self, num_workers=1, task_func_dict=None, worker_task=None, args=None):
        self.num_workers = int(num_workers)
        self.worker_task = worker_task
        self.process_result_funcs = task_func_dict or {}
        self.pipes = []
        self.processes = []
        self.queues = []
        self.args = list(args) if args is not None else []
        self.cur_arg_ct = 0

    def start_workers(self):
        for i in range(self.num_workers):
            parent_conn, child_conn = multiprocessing.Pipe()
            q = multiprocessing.Queue()
            self.pipes.append(parent_conn)
            self.queues.append(q)
            p = multiprocessing.Process(
                target=worker_wrapper,                 # module-level, NOT bound to self
                args=(i, self.worker_task, q, child_conn),
            )
            p.start()
            self.processes.append(p)
            child_conn.close()  # main never uses child end

    def add_args(self):
        # ensure one STOP per worker
        self.args.extend(["STOP"] * self.num_workers)
        # seed each worker with one task (or STOP)
        for i in range(self.num_workers):
            self.queues[i].put(self.args[i] if i < len(self.args) else "STOP")
        self.cur_arg_ct += min(self.num_workers, len(self.args))

    def _schedule_next(self, worker_id):
        q = self.queues[worker_id]
        if self.cur_arg_ct < len(self.args):
            q.put(self.args[self.cur_arg_ct])
            self.cur_arg_ct += 1
        else:
            q.put("STOP")

    def run(self):
        if not self.args:
            return
        self.start_workers()
        self.add_args()

        active = self.num_workers
        while active > 0 and self.pipes:
            ready = multiprocessing.connection.wait(self.pipes)
            for conn in ready:
                try:
                    worker_id, request_args, request_type = conn.recv()
                except EOFError:
                    try:
                        conn.close()
                    except Exception:
                        pass
                    if conn in self.pipes:
                        self.pipes.remove(conn)
                    active -= 1
                    continue

                resp = ()
                if request_type == "TASK_COMPLETE":
                    self._schedule_next(worker_id)
                elif request_type == "TASK_ERROR":
                    # Let caller handle logging if they provided a handler
                    handler = self.process_result_funcs.get("TASK_ERROR")
                    if handler:
                        try:
                            handler(request_args)  # e.g., (name, err, traceback)
                        except Exception:
                            pass
                    # Still schedule the next task so the worker keeps going
                    self._schedule_next(worker_id)
                else:
                    handler = self.process_result_funcs.get(request_type)
                    resp = handler(request_args) if handler else ()

                try:
                    conn.send(resp)
                except (BrokenPipeError, EOFError):
                    try:
                        conn.close()
                    except Exception:
                        pass
                    if conn in self.pipes:
                        self.pipes.remove(conn)
                    active -= 1

        for p in self.processes:
            try:
                p.join()
            except Exception:
                pass

    def shutdown(self):
        for q in self.queues:
            try:
                q.put("STOP")
            except Exception:
                pass
        for conn in self.pipes:
            try:
                conn.close()
            except Exception:
                pass
        for p in self.processes:
            try:
                if p.is_alive():
                    p.terminate()
            except Exception:
                pass
