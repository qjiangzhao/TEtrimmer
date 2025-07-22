import multiprocessing
import multiprocessing.connection


'''
Method for asynchronously processing tasks from beginning to end, but with communication to and waiting on results from main

args must contain all of the information required to run the tasks

piped results from workers are strictly of format (worker_id, (*request_args,), request_type,)

Workers may handle sends as lists.
'''
class ChattyParallelProcessor:
	def __init__(self, num_workers = 1, task_func_dict = None, worker_task = None, args = None):
		self.num_workers = num_workers
		self.worker_task = worker_task
		self.process_result_funcs = task_func_dict
		self.pipes = []
		self.processes = []
		self.queues = []
		
		self.args = args
		
		self.cur_arg_ct = 0
		
	def start_workers(self):
		for i in range(self.num_workers):
			parent_conn, child_conn = multiprocessing.Pipe()
			this_queue = multiprocessing.Queue()
			self.pipes.append(parent_conn)
			p = multiprocessing.Process(target=self._worker_wrapper, args=(i, self.worker_task, this_queue, child_conn))
			p.start()
			self.processes.append(p)
			child_conn.close()  # Close the child end in the main process to avoid resource leaks
			self.queues.append(this_queue)

	#For load balancing, passing args in order from biggest -> smallest is wise.
	#Adds the first num workers args to initialize the workers
	def add_args(self):		
		self.args.extend(['STOP']*self.num_workers) #Add required number of stop signals
		#Add initial arguments to the queue, add more as tasks finish
		for i in range(0, self.num_workers):
			q = self.queues[i]
			q.put(self.args[i])
			
		self.cur_arg_ct += self.num_workers #Update place
		
	def _worker_wrapper(self, worker_id, worker_task, args, pipe):
		def send_and_receive(result):
			pipe.send(result) #Send result to main
			from_main = pipe.recv()  # Wait for acknowledgment from main that the result has been processed and give results
			return from_main
		try:
			self.worker_task(worker_id, args, send_and_receive)
		finally:
			pipe.close()  # Ensure pipe is closed when done

	def run(self):
		if len(self.args) > 0:
			self.start_workers()
			self.add_args()
			
			active_workers = self.num_workers
			while active_workers > 0:
				ready = multiprocessing.connection.wait(self.pipes)
				for conn in ready:
					try:
						result = conn.recv()
						worker_id = result[0]
						request_args = result[1]
						request_type = result[2]
						
						res_set = ()
						
						#print(worker_id, request_args, request_type)
						#print("")
						
						#Add next task or stop;
						if request_type == "TASK_COMPLETE":
							q = self.queues[result[0]]
							q.put(self.args[self.cur_arg_ct])
							self.cur_arg_ct += 1
						else:
							res_set = self.process_result_funcs[request_type](request_args)
						
						conn.send(res_set)  # Send acknowledgment
					except EOFError:
						conn.close()
						self.pipes.remove(conn)
						active_workers -= 1
			for p in self.processes:
				p.join()

	#Elegant shutdown
	def shutdown(self):
		#Put stop signal into each queue
		for q in self.queues:
			q.put("STOP")
		#Close communication
		for conn in self.pipes:
			conn.close()
		#Kill child processes
		for p in self.processes:
			if p.is_alive():
				p.terminate()
	
'''
def process_result_func(args):
	print(f"Main received: {args[1]} from worker {args[0]}")
	# Simulate processing time
	return args
	
def end_whopple(args):
	print(f"Main finally received: {args[1]} from worker {args[0]}")
	# Simulate processing time
	return args
	
	
#needs to be extended to 
def worker_task(my_id, queue, send_and_wait):
	while True:
		next_input = queue.get()
		#The worker is done.
		if next_input == "STOP":
			break
		
		
		
		main_ct = send_and_wait((my_id, (my_id, "start_sleep",), "MID_PROC",))
		# Simulate further processing with random duration
		sleepy = random.randint(1, 10)
		time.sleep(sleepy/10)
		send_and_wait((my_id, (my_id, "end_sleep",), "END_PROC",))
		
		
		send_and_wait((my_id, (my_id, "done",), "TASK_COMPLETE",))

if __name__ == "__main__":
	#Just in cast multiple functions are needed.
	task_names = ["MID_PROC", "END_PROC"]
	task_funcs = [process_result_func, end_whopple]
	tfd = dict(zip(task_names, task_funcs))
	
	long_args = list(range(0, 1000))
	processor = ChattyParallelProcessor(num_workers=20, task_func_dict = tfd, worker_task=worker_task, args = long_args)
	try:
		processor.run()
	except KeyboardInterrupt:
		processor.shutdown()
'''
	