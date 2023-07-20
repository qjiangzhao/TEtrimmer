import os
import shutil
import subprocess
from tkinter import Tk, Frame, Button, messagebox, Scrollbar, Canvas, Label, Menu

# Specify the directories you want to use
source_directory = r"/work/ur376715/TE_analysis/Software_develpment_for_TE_curation/version_18/output_dir/Multiple_sequence_alignment"
target_directory_copy = r"/work/ur376715/TE_analysis/Software_develpment_for_TE_curation/version_14/output_dir/Sequences_for_consensus_generation"
target_directory_more_extend = r"/work/ur376715/TE_analysis/Software_develpment_for_TE_curation/version_14/output_dir/Sequences_need_more_extention"
target_directory_discard = r"/work/ur376715/TE_analysis/Software_develpment_for_TE_curation/version_14/output_dir/Trash"

# Path to the AliView program
aliview_path = "/work/ur376715/TE_analysis/Software_develpment_for_TE_curation/version_14/aliview/aliview"

def open_file(filename, button):
    def _open_file(event):
        if filename.lower().endswith(('.fa', '.fasta')):
            subprocess.run([aliview_path, os.path.join(source_directory, filename)])
        elif filename.lower().endswith('.pdf'):
            subprocess.run(['gio', 'open', os.path.join(source_directory, filename)])
        button.config(bg='yellow')
    return _open_file

def copy_file(filename, button, target_directory):
    def _copy_file(event):
        if messagebox.askokcancel("Confirmation", f"Do you want to copy '{filename}' to '{target_directory}'?"):
            if not os.path.exists(target_directory):
                os.makedirs(target_directory)
            try:
                shutil.copy(os.path.join(source_directory, filename), target_directory)
                button.config(bg='green')
            except Exception as e:
                messagebox.showerror("Error", f"An error occurred while copying the file: {str(e)}")
    return _copy_file

root = Tk()
root.title("TE Trimmer")
root.geometry('1200x800')

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

# Sort files in the source_directory
sorted_files = sorted(os.listdir(source_directory))

def clear_frame():
    for widget in frame.winfo_children():
        widget.destroy()

def load_files(start, end):
    clear_frame()
    for i, filename in enumerate(sorted_files[start:end], start=start):
        line_number = Label(frame, text=str(i + 1), bg='white')
        line_number.grid(row=i-start, column=0)

        file_button = Button(frame, text=filename, anchor='w', bg='white')
        file_button.grid(row=i-start, column=1, sticky='ew')
        file_button.bind('<Double-Button-1>', open_file(filename, file_button))

        button_frame = Frame(frame, bg='white')  # Create a new frame for the buttons
        button_frame.grid(row=i-start, column=2, sticky='e')  # Place the button frame on the right side

        copy_button = Button(button_frame, text="Consensus", bg='white', fg='black')
        copy_button.grid(row=0, column=0, padx=5)  # Use grid layout manager inside the button frame
        copy_button.bind('<Button-1>', copy_file(filename, copy_button, target_directory_copy))  # Update binding

        more_extend_button = Button(button_frame, text="Extension", bg='white', fg='black')
        more_extend_button.grid(row=0, column=1, padx=5)  # Use grid layout manager inside the button frame
        more_extend_button.bind('<Button-1>', copy_file(filename, more_extend_button, target_directory_more_extend))  # Update binding

        discard_button = Button(button_frame, text="Trash", bg='white', fg='black')
        discard_button.grid(row=0, column=2, padx=5)  # Use grid layout manager inside the button frame
        discard_button.bind('<Button-1>', copy_file(filename, discard_button, target_directory_discard))  # Update binding

        # Configure column and row weights to expand the button_frame
        button_frame.grid_columnconfigure(0, weight=1)
        button_frame.grid_rowconfigure(0, weight=1)

        # Configure column weights to expand the file button and align buttons to the right
        frame.grid_columnconfigure(1, weight=1)
        frame.grid_columnconfigure(2, weight=0)

menubar = Menu(root)
root.config(menu=menubar)

# Create a pull-down menu and add it to the menu bar
fileMenu = Menu(menubar, tearoff=0)
menubar.add_cascade(label="File Range", menu=fileMenu)

for i in range(0, len(sorted_files), 1000):
    end = min(i+1000, len(sorted_files))
    fileMenu.add_command(label=f"{i+1}-{end}", command=lambda i=i, end=end: load_files(i, end))

load_files(0, min(1000, len(sorted_files)))

root.mainloop()
