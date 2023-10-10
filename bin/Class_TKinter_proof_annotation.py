import os
import shutil
import subprocess
from tkinter import Tk, Frame, Button, messagebox, Scrollbar, Canvas, Label, Menu, BooleanVar
import click

copy_history = []


@click.command()
@click.option('--te_trimmer_proof_annotation_dir', '-i', required=True, type=str,
              help='Supply the TE Trimmer output directory path')
@click.option('--output_dir', '-o', default=None, type=str,
              help='Define the output directory for TE proof annotation. Default: parent directory of input file')
def proof_annotation(te_trimmer_proof_annotation_dir, output_dir):
    """
    This tool can help do quick proof annotation

    python ./path_to_TE_Trimmer_bin/Class_TKinter_proof_annotation.py -i "TE_Trimmer_output_folder"
    """

    if output_dir is None:
        output_dir = os.path.dirname(te_trimmer_proof_annotation_dir)

    # so.path.abspath(__file__) will return the current executable python file
    bin_py_path = os.path.dirname(os.path.abspath(__file__))

    # Path to the Aliview program
    aliview_path = os.path.join(bin_py_path, "aliview/aliview")

    # Define output folders
    consensus_folder = os.path.join(output_dir, "Proof_annotation_consensus_folder")
    need_more_extension = os.path.join(output_dir, "Proof_annotation_need_more_extension")
    trash = os.path.join(output_dir, "Proof_annotation_trash")
    low_copy_elements = os.path.join(output_dir, "Proof_annotation_low_copy_elements")

    # Check if output folder exists, otherwise create it.
    for dir_path in [consensus_folder, need_more_extension, trash, low_copy_elements]:
        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)

    consensus_folder = os.path.abspath(consensus_folder)
    need_more_extension = os.path.abspath(need_more_extension)
    trash = os.path.abspath(trash)
    low_copy_elements = os.path.abspath(low_copy_elements)

    root = Tk()
    root.title("TE Trimmer")
    root.geometry('1200x800')

    # Initialize the BooleanVar here, after the root Tk instance is created
    show_confirmation = BooleanVar(value=True)

    canvas = Canvas(root, bg='white')
    scrollbar = Scrollbar(root, orient='vertical', command=canvas.yview)
    canvas.configure(yscrollcommand=scrollbar.set)
    scrollbar.pack(side='right', fill='y')
    canvas.pack(side='left', fill='both', expand=True)

    def scroll(event):
        if event.num == 4 or event.delta == 120:
            canvas.yview_scroll(-1, "units")
        elif event.num == 5 or event.delta == -120:
            canvas.yview_scroll(1, "units")

    root.bind("<MouseWheel>", scroll)
    root.bind("<Button-4>", scroll)
    root.bind("<Button-5>", scroll)

    frame = Frame(canvas, bg='white')
    canvas_frame = canvas.create_window((0, 0), window=frame, anchor='nw')

    def update_scrollregion(event):
        canvas.configure(scrollregion=canvas.bbox('all'))

    def on_configure(event):
        canvas.itemconfig(canvas_frame, width=event.width)

    canvas.bind('<Configure>', on_configure)
    frame.bind('<Configure>', update_scrollregion)

    def open_file(filename, button, source_dir=te_trimmer_proof_annotation_dir):
        def _open_file(event):
            if filename.lower().endswith(('.fa', '.fasta')):
                subprocess.run([aliview_path, os.path.join(source_dir, filename)])
            elif filename.lower().endswith('.pdf'):
                subprocess.run(['gio', 'open', os.path.join(source_dir, filename)])
            button.config(bg='yellow')

        return _open_file

    def disable_confirmation():
        """Disable the show_confirmation variable."""
        show_confirmation.set(False)

    def copy_file(filename, button, target_directory, is_low_copy=False):
        def _copy_file(event):

            if is_low_copy:
                file_path = os.path.join(te_trimmer_proof_annotation_dir, "Low_copy_TE", filename)
                print(file_path)
            else:
                file_path = os.path.join(te_trimmer_proof_annotation_dir, filename)

            # Check if the confirmation should be shown
            if not show_confirmation.get() or messagebox.askokcancel("Confirmation",
                                                                     f"Do you want to copy '{filename}' to '{target_directory}'?"):
                os.makedirs(target_directory, exist_ok=True)
                try:
                    shutil.copy(file_path, target_directory)
                    button.config(bg='green')

                    # Check if copy_history length exceeds 100
                    if len(copy_history) >= 100:
                        copy_history.pop(0)  # Remove the oldest entry

                    copy_history.append(
                        (os.path.join(target_directory, filename), target_directory, button))  # Append to history
                except Exception as e:
                    messagebox.showerror("Error", f"An error occurred while copying the file: {str(e)}")

        return _copy_file

    def undo_last_copy():
        if not copy_history:
            messagebox.showinfo("Info", "No actions to undo.")
            return

        last_copied_file, target_directory, last_button = copy_history[-1]

        if not os.path.exists(last_copied_file):
            messagebox.showerror("Error", f"File {last_copied_file} doesn't exist. Can't undo.")
            return

        if messagebox.askokcancel("Confirmation", f"Do you want to undo the copy of '{last_copied_file}'?"):
            try:
                os.remove(last_copied_file)
                copy_history.pop()  # Remove the last action from the history
                last_button.config(bg='white')  # Reset the button color
                messagebox.showinfo("Info", f"Successfully removed '{last_copied_file}'.")
            except Exception as e:
                messagebox.showerror("Error", f"An error occurred while deleting the file: {str(e)}")

    def clear_frame():
        for widget in frame.winfo_children():
            widget.destroy()

    def load_files(start, end, source_dir=te_trimmer_proof_annotation_dir):
        clear_frame()
        sorted_files = [f for f in sorted(os.listdir(source_dir)) if f != "Low_copy_TE"]

        for i, filename in enumerate(sorted_files[start:end], start=start):
            line_number = Label(frame, text=str(i + 1), bg='white')
            line_number.grid(row=i - start, column=0)

            file_button = Button(frame, text=filename, anchor='w', bg='white')
            file_button.grid(row=i - start, column=1, sticky='ew')
            file_button.bind('<Double-Button-1>', open_file(filename, file_button, source_dir))

            button_frame = Frame(frame, bg='white')
            button_frame.grid(row=i - start, column=2, sticky='e')

            if not source_dir.endswith("Low_copy_TE"):
                copy_button = Button(button_frame, text="Consensus", bg='white', fg='black')
                copy_button.grid(row=0, column=0, padx=5)
                copy_button.bind('<Button-1>', copy_file(filename, copy_button, consensus_folder))

                more_extend_button = Button(button_frame, text="Extension", bg='white', fg='black')
                more_extend_button.grid(row=0, column=1, padx=5)
                more_extend_button.bind('<Button-1>', copy_file(filename, more_extend_button, need_more_extension))

                discard_button = Button(button_frame, text="Trash", bg='white', fg='black')
                discard_button.grid(row=0, column=2, padx=5)
                discard_button.bind('<Button-1>', copy_file(filename, discard_button, trash))

            else:
                low_copy_button = Button(button_frame, text="Low_copy", bg='white', fg='black')
                low_copy_button.grid(row=0, column=3, padx=5)
                low_copy_button.bind('<Button-1>', copy_file(filename, low_copy_button, low_copy_elements, True))

            button_frame.grid_columnconfigure(0, weight=1)
            button_frame.grid_rowconfigure(0, weight=1)

            frame.grid_columnconfigure(1, weight=1)
            frame.grid_columnconfigure(2, weight=0)


    # Set menu bar
    menubar = Menu(root)
    root.config(menu=menubar)

    fileMenu = Menu(menubar, tearoff=0)
    menubar.add_cascade(label="High copy TE", menu=fileMenu)

    sorted_files = sorted(os.listdir(te_trimmer_proof_annotation_dir))
    for i in range(0, len(sorted_files), 1000):
        end = min(i + 1000, len(sorted_files))
        fileMenu.add_command(label=f"{i + 1}-{end}", command=lambda i=i, end=end: load_files(i, end))

    low_copy_menu = Menu(menubar, tearoff=0)
    menubar.add_cascade(label="Low copy", menu=low_copy_menu)
    low_copy_sorted_files = sorted(os.listdir(os.path.join(te_trimmer_proof_annotation_dir, "Low_copy_TE")))
    for i in range(0, len(low_copy_sorted_files), 1000):
        end = min(i + 1000, len(low_copy_sorted_files))
        low_copy_menu.add_command(label=f"{i + 1}-{end}",
                                  command=lambda i=i, end=end: load_files(i, end,
                                                                          source_dir=os.path.join(
                                                                              te_trimmer_proof_annotation_dir,
                                                                              "Low_copy_TE")))

    # Add the Settings Menu to the main menubar
    confirm_menu = Menu(menubar, tearoff=0)
    menubar.add_cascade(label="Settings", menu=confirm_menu)

    # Add a checkbutton to the settings menu
    confirm_menu.add_checkbutton(label="Show Confirmation Window", onvalue=True, offvalue=False,
                                 variable=show_confirmation)

    # Add a command to the settings menu to disable the confirmation window
    confirm_menu.add_command(label="Don't Show Confirmation Window", command=disable_confirmation)

    menubar.add_command(label="Undo", command=undo_last_copy)

    load_files(0, min(1000, len(sorted_files)))

    root.mainloop()


if __name__ == '__main__':
    proof_annotation()

