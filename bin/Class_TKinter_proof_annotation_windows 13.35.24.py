import os
import shutil
from tkinter import Tk, Frame, Button, messagebox, Scrollbar, Canvas, Label

# Specify the directories you want to use
source_directory = r"G:\TE_manual_curation\Software_develop\test_for_html_index"
target_directory_copy = r"G:\TE_manual_curation\Software_develop\destination_folder_copy"
target_directory_more_extend = r"G:\TE_manual_curation\Software_develop\destination_folder_more_extend"
target_directory_discard = r"G:\TE_manual_curation\Software_develop\destination_folder_discard"


def open_file(filename, button):
    def _open_file(event):
        # Open the file with the default application
        os.startfile(os.path.join(source_directory, filename))
        # Change the button color after the file is opened
        button.config(bg='yellow')

    return _open_file


def copy_file(filename, button, target_directory):
    def _copy_file(event):
        # Show a confirmation window before copying
        if messagebox.askokcancel("Confirmation", f"Do you want to copy '{filename}' to '{target_directory}'?"):
            # Check if the target directory exists. If not, create it.
            if not os.path.exists(target_directory):
                os.makedirs(target_directory)

            # Copy the file to the target directory
            try:
                shutil.copy(os.path.join(source_directory, filename), target_directory)
                button.config(bg='green')
            except Exception as e:
                messagebox.showerror("Error", f"An error occurred while copying the file: {str(e)}")

    return _copy_file


# Create a GUI window
root = Tk()
root.title("TE Trimmer")
root.geometry('1200x800')

# Create a Canvas and a Scrollbar
canvas = Canvas(root, bg='white')
scrollbar = Scrollbar(root, orient='vertical', command=canvas.yview)
canvas.configure(yscrollcommand=scrollbar.set)
scrollbar.pack(side='right', fill='y')
canvas.pack(side='left', fill='both', expand=True)

# Enable mouse wheel scrolling
root.bind("<MouseWheel>", lambda event: canvas.yview_scroll(int(-1 * (event.delta / 120)), "units"))

# Create a Frame inside the Canvas
frame = Frame(canvas, bg='white')
canvas_frame = canvas.create_window((0, 0), window=frame, anchor='nw')


# Update the scrollregion after widgets are placed on the Canvas
def update_scrollregion(event):
    canvas.configure(scrollregion=canvas.bbox('all'))


# Update the size of the Canvas window
def on_configure(event):
    canvas.itemconfig(canvas_frame, width=event.width)


canvas.bind('<Configure>', on_configure)
frame.bind('<Configure>', update_scrollregion)

# Add files and buttons to the Frame
i = 0
for filename in os.listdir(source_directory):
    line_number = Label(frame, text=str(i + 1), bg='white')
    line_number.grid(row=i, column=0)

    file_button = Button(frame, text=filename, anchor='w', bg='white')
    file_button.grid(row=i, column=1, sticky='ew')
    file_button.bind('<Double-Button-1>', open_file(filename, file_button))

    copy_button = Button(frame, text="Consensus", bg='white', fg='black')
    copy_button.grid(row=i, column=2)
    copy_button.bind('<Button-1>', copy_file(filename, copy_button, target_directory_copy))

    more_extend_button = Button(frame, text="Extension", bg='white', fg='black')
    more_extend_button.grid(row=i, column=3)
    more_extend_button.bind('<Button-1>', copy_file(filename, more_extend_button, target_directory_more_extend))

    discard_button = Button(frame, text="Trash", bg='white', fg='black')
    discard_button.grid(row=i, column=4)
    discard_button.bind('<Button-1>', copy_file(filename, discard_button, target_directory_discard))

    i += 1

# Run the GUI
root.mainloop()
