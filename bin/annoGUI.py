import os
import re
import shutil
import subprocess
from tkinter import Tk, Frame, Button, messagebox, Scrollbar, Canvas, Label, Menu, BooleanVar, Toplevel
import click
from functools import partial
import platform


#####################################################################################################
# Code block: make Aliveiw available to be used
#####################################################################################################
# Change Aliview permission
def change_permissions_recursive(input_dir, mode):
    try:
        for dirpath, dirnames, filenames in os.walk(input_dir):
            os.chmod(dirpath, mode)
            for filename in filenames:
                os.chmod(os.path.join(dirpath, filename), mode)
    except PermissionError:
        click.echo("TETrimmer don't have right to change permissions. Pleas use sudo to run TETrimmer")
        return False
    return True


# Define Aliview software path and change permission
bin_py_path = os.path.dirname(os.path.abspath(__file__))
aliview_path = os.path.join(bin_py_path, "aliview/aliview")
change_permissions_recursive(aliview_path, 0o755)

#####################################################################################################
# Code block: set click command
#####################################################################################################

@click.command()
@click.option('--te_trimmer_proof_annotation_dir', '-i', required=True, type=str,
              help='Supply the TETrimmer output directory path')
@click.option('--output_dir', '-o', default=None, type=str,
              help='Define the output directory for TE proof annotation. Default: input folder directory')
def proof_annotation(te_trimmer_proof_annotation_dir, output_dir):
    """
    This tool can help do quick proof annotation

    python ./path_to_TETrimmer_bin/Class_TKinter_proof_annotation.py -i "TETrimmer_output_folder"
    """

    # Define empty list to store copy history, which enable undo button
    copy_history = []

    # Detect system OS type
    os_type = platform.system()

    # If the -o option is not given, use the parent directory of -i as output directory.
    if output_dir is None:
        output_dir = te_trimmer_proof_annotation_dir

    # Define Aliview software path and change permission
    bin_py_path = os.path.dirname(os.path.abspath(__file__))
    aliview_path = os.path.join(bin_py_path, "aliview/aliview")
    if os_type == "Windows":
        aliview_path = os.path.join(bin_py_path, r"aliview\aliview.jar")

    # Define output folders, create them when they are not found
    consensus_folder = os.path.abspath(os.path.join(output_dir, "Proof_annotation_consensus_folder"))
    need_more_extension = os.path.abspath(os.path.join(output_dir, "Proof_annotation_need_more_extension"))
    others_dir = os.path.abspath(os.path.join(output_dir, "Proof_annotation_others"))
    low_copy_elements = os.path.abspath(os.path.join(output_dir, "Proof_annotation_low_copy_elements"))
    rescue_skip_elements = os.path.abspath(os.path.join(output_dir, "Proof_annotation_rescued_skip_elements"))

    for dir_path in [consensus_folder, need_more_extension, low_copy_elements]:
        print(dir_path)
        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)

    #####################################################################################################
    # Code block: build TKinter window
    #####################################################################################################

    # Initialize Tk window
    root = Tk()
    root.title("TETrimmer proof annotation tool")
    root.geometry('1200x800')

    # Create canvas on root
    canvas = Canvas(root, bg='white')

    # fill='both'  This tells the canvas to expand and fill both the X-axis (horizontally) and Y-axis (vertically)
    # in its parent widget. The canvas will take up as much space as possible in both directions.
    # expand=True  If the window is resized, the canvas will grow or shrink accordingly.
    canvas.pack(side='left', fill='both', expand=True)

    #####################################################################################################
    # Code block: set a vertical scroll bar
    #####################################################################################################
    # Set scrollbar. command=canvas.yview links the scrollbar to the vertical view of the canvas.
    # This means that scrolling the scrollbar will move the canvas vertically.
    scrollbar = Scrollbar(root, orient='vertical', command=canvas.yview)
    canvas.configure(yscrollcommand=scrollbar.set)
    scrollbar.pack(side='right', fill='y')

    def scroll(event, canvas):
        if event.num == 4 or event.delta > 0:
            canvas.yview_scroll(-1, "units")
        elif event.num == 5 or event.delta < 0:
            canvas.yview_scroll(1, "units")

    # <MouseWheel> is typically used on Windows and macOS to handel mouse wheel movements
    # <Button-4> and <Button-5> are for Linux/Unix systems
    root.bind("<MouseWheel>", lambda event: scroll(event, canvas))
    root.bind("<Button-4>", lambda event: scroll(event, canvas))
    root.bind("<Button-5>", lambda event: scroll(event, canvas))

    #####################################################################################################
    # Code block: re-size canvas and frame size when the windows size is modified
    #####################################################################################################

    # Frame is a container widget, it is used to organize and group other widgets together
    frame = Frame(canvas, bg='white')  # canvas specifies the parent widget in which this frame will be placed
    canvas_frame = canvas.create_window((0, 0), window=frame, anchor='nw')

    def on_configure(event, canvas, canvas_frame):
        canvas.itemconfig(canvas_frame, width=event.width)

    def update_scrollregion(event, canvas):
        canvas.configure(scrollregion=canvas.bbox('all'))

    canvas.bind('<Configure>', lambda event: on_configure(event, canvas, canvas_frame))
    frame.bind('<Configure>', lambda event: update_scrollregion(event, canvas))

    #####################################################################################################
    # Code block: show the start page
    #####################################################################################################
    global logo_label
    global text_label

    log_text = "████████╗███████╗████████╗██████╗ ██╗███╗   ███╗███╗   ███╗███████╗██████╗\n"\
               "╚══██╔══╝██╔════╝╚══██╔══╝██╔══██╗██║████╗ ████║████╗ ████║██╔════╝██╔══██╗\n"\
               "   ██║   █████╗     ██║   ██████╔╝██║██╔████╔██║██╔████╔██║█████╗  ██████╔╝\n"\
               "   ██║   ██╔══╝     ██║   ██╔══██╗██║██║╚██╔╝██║██║╚██╔╝██║██╔══╝  ██╔══██╗\n"\
               "   ██║   ███████╗   ██║   ██║  ██║██║██║ ╚═╝ ██║██║ ╚═╝ ██║███████╗██║  ██║\n"\
               "   ╚═╝   ╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝╚═╝     ╚═╝╚═╝     ╚═╝╚══════╝╚═╝  ╚═╝\n"\


    initial_text = "TETrimmer proof annotation assistant tool\n\n" \
                   "Introduction:\n\n" \
                   "We highly recommend to do proof annotation to increase your TE annotation quality.\n\n"\
                   "1, Click the buttons in the menu bar, which corresponds to the different annotation status.\n\n" \
                   "2, All files in the chose folder (menu button) will be displayed.\n\n" \
                   "3, For each TE output, you can find four files including <seq_name.anno.fa>, <seq_name.fa>, <seq_name.bed>, and " \
                   "<seq_name.pdf>.\n\n" \
                   "   <seq_name.anno.fa> is the multiple sequence alignment (MSA) file before cleaning\n" \
                   "   <seq_name.fa> is the multiple sequence alignment (MSA) file after cleaning\n" \
                   "   <seq_name.pdf> contains four plots used to evaluate annotation quality\n" \
                   "   <seq_name.bed> contains sequence position information at the genome of MSA. Used for further extension.\n\n" \
                   "4, Double click <seq_name.pdf> and evaluate annotation quality.\n\n" \
                   "5, If you are satisfied with the result, click 'Consensus' button behind <seq_name.fa>. This MSA " \
                   "file will go to <Proof_annotation_consensus_folder>.\n\n" \
                   "6, If you are not satisfied, double click <seq_name.fa> or <seq_name.anno.fa> to modify MSA and " \
                   "save it to consensus folder.\n\n" \
                   "7, If you want more extension for the MSA, click 'Extension' behind <seq_name.bed> and this " \
                   "bed file will be salved to <Proof_annotation_need_more_extension> folder.\n\n" \
                   "8, For low copy element, check the pdf file and decide if to include it into the final consensus " \
                   "library. Note: low copy element do not have multiple sequence alignment file"

    # Display ASCII logo with 'Courier' font
    if os_type == "Linux":
        logo_font = ('DejaVu Sans Mono', 10)
    elif os_type == "Darwin":  # macOS
        logo_font = ('Courier', 13)
    else:
        logo_font = ('Courier', 5)
    logo_label = Label(canvas, text=log_text, bg='white', font=logo_font, justify='left', wraplength=1100)
    logo_label.pack(pady=10)

    # Display the explanatory text with 'Arial' font
    text_label = Label(canvas, text=initial_text, bg='white', font=('Arial', 15), justify='left', wraplength=1100)
    text_label.pack(pady=10)

    #####################################################################################################
    # Code block: Define functions
    #####################################################################################################

    def numerical_sort_key(filename):
        numbers = re.findall(r'\d+', filename)
        return [int(num) for num in numbers]

    def clear_frame():
        for widget in frame.winfo_children():
            widget.destroy()

    def destroy_initial_label():
        global logo_label
        if logo_label:
            logo_label.destroy()
            logo_label = None
        global text_label
        if text_label:
            text_label.destroy()
            text_label = None

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
                copy_history.pop()
                last_button.config(bg='white')
                messagebox.showinfo("Info", f"Successfully removed '{last_copied_file}'.")
            except Exception as e:
                messagebox.showerror("Error", f"An error occurred while deleting the file: {str(e)}")

    def copy_file(filename, button, target_directory, source_dir, parent_win):
        def _copy_file(event):

            file_path = os.path.join(source_dir, filename)
            if not show_confirmation.get() or messagebox.askokcancel(
                    "Confirmation", f"Do you want to copy '{filename}' to '{target_directory}'?", parent=parent_win):
                os.makedirs(target_directory, exist_ok=True)
                try:
                    shutil.copy(file_path, target_directory)
                    button.configure(bg='green')
                    button.update_idletasks()
                    if len(copy_history) >= 100:
                        copy_history.pop(0)
                    copy_history.append((os.path.join(target_directory, filename), target_directory, button))
                except Exception as e:
                    messagebox.showerror("Error", f"An error occurred while copying the file: {str(e)}",
                                         parent=parent_win)
        return _copy_file

    def load_files(start, end, frame, canvas, source_dir=te_trimmer_proof_annotation_dir):
        clear_frame()
        canvas.yview_moveto(0)  # Reset scrollbar to top
        if not os.path.exists(source_dir) or not os.listdir(source_dir):
            label = Label(frame, text="No files found here, try another folder.", bg='white')
            label.pack(pady=20)
            return
        # Sort files
        sorted_files = sorted(os.listdir(source_dir), key=numerical_sort_key)
        for i, filename in enumerate(sorted_files[start:end], start=start):
            # Add line number into canvas frame
            line_number = Label(frame, text=str(i + 1), bg='white')
            line_number.grid(row=i - start, column=0)

            # Add file name button into canvas frame
            file_button = Button(frame, text=filename, anchor='w', bg='white')
            file_button.grid(row=i - start, column=1, sticky='ew')
            # Bind with open_file function to open file
            file_button.bind('<Double-Button-1>', open_file(filename, file_button, source_dir))

            # Build a child button_frame inside frame
            button_frame = Frame(frame, bg='white')
            button_frame.grid(row=i - start, column=2, sticky='e')

            # If it's a folder in "Clustered_proof_annotation", add a label with the file count
            if source_dir.endswith("Clustered_proof_annotation") and os.path.isdir(os.path.join(source_dir, filename)):
                file_count = len(os.listdir(os.path.join(source_dir, filename)))
                file_count_label = Label(frame, text=f"({file_count} files)", bg='white')
                file_count_label.grid(row=i - start, column=3)

            # Don't show "Consensus" and "Extension" button for low copy and clustered proof annotation
            if not source_dir.endswith("Low_copy_TE") and not source_dir.endswith("Clustered_proof_annotation")\
                    and not source_dir.endswith("Skipped_TE"):
                # Create "Consensus" button inside button_frame
                copy_button = Button(button_frame, text="Consensus", bg='white', fg='black')
                copy_button.grid(row=0, column=0, padx=5)
                # Bind "Consensus" button with copy_file function with specific source and destination folder
                copy_button.bind('<Button-1>', copy_file(filename, copy_button, consensus_folder,
                                                         source_dir, root))

                # Define "Extension" button
                more_extend_button = Button(button_frame, text="Extension", bg='white', fg='black')
                more_extend_button.grid(row=0, column=1, padx=5)
                # Bind "Extension" button with copy_file function with different destination folder
                more_extend_button.bind('<Button-1>', copy_file(filename, more_extend_button, need_more_extension,
                                                                source_dir, root))
            # Add rescue button for skipped TE
            if os.path.basename(source_dir) == "Skipped_TE":
                rescue_skip_button = Button(button_frame, text="Rescue skip", bg='white', fg='black')
                rescue_skip_button.grid(row=0, column=1, padx=5)
                rescue_skip_button.bind('<Button-1>', copy_file(filename, rescue_skip_button, rescue_skip_elements,
                                                                source_dir, root))

            # Add "Save to low copy" button for low copy TE
            if os.path.basename(source_dir) == "Low_copy_TE":
                # Define "Low_copy" button with specific source and destination folder
                low_copy_button = Button(button_frame, text="Save to low copy", bg='white', fg='black')
                low_copy_button.grid(row=0, column=1, padx=5)
                low_copy_button.bind('<Button-1>', copy_file(filename, low_copy_button, low_copy_elements,
                                                             source_dir, root))
            button_frame.grid_columnconfigure(0, weight=1)
            button_frame.grid_rowconfigure(0, weight=1)
            frame.grid_columnconfigure(1, weight=1)
            frame.grid_columnconfigure(2, weight=0)

    def load_files_with_destroy(start, end, frame, canvas, path):
        destroy_initial_label()
        load_files(start, end, frame, canvas, source_dir=path)

    # Set nested function to enable this function containing parameters to be used by button
    def open_file(filename, button, source_dir):
        def _open_file(event):
            filepath = os.path.join(source_dir, filename)
            if os_type == "Darwin":
                button.config(fg='red')  # Change button text color under macOS system
                button.update_idletasks()  # Update UI immediately
            else:
                button.config(bg='yellow')  # Change button color
                button.update_idletasks()  # Update UI immediately
            if os.path.isdir(filepath):
                open_folder(filename, source_dir)
            else:
                if filename.lower().endswith(('.fa', '.fasta')):
                    if os_type == "Windows":
                        subprocess.run(["java", "-jar", aliview_path, filepath])
                    else:
                        subprocess.run([aliview_path, filepath])
                elif filename.lower().endswith('.pdf'):
                    if os_type == "Linux":
                        subprocess.run(['xdg-open', filepath])
                    elif os_type == "Darwin":  # macOS
                        subprocess.run(['open', filepath])
                    elif os_type == "Windows":
                        os.startfile(filepath)
                    else:
                        subprocess.run(['xdg-open', filepath])
                elif filename.lower().endswith(('.txt', '.py', '.csv', '.md', '.bed')):
                    if os_type == "Linux":
                        text_editor = 'gedit'  # Replace 'gedit' with your preferred text editor
                        subprocess.run([text_editor, filepath])
                    elif os_type == "Darwin":  # macOS
                        subprocess.run(['open', '-a', 'TextEdit', filepath])
                    elif os_type == "Windows":
                        notepad_path = 'notepad.exe'  # or path to another text editor if preferred
                        subprocess.run([notepad_path, filepath])
                    else:
                        text_editor = 'gedit'  # Fallback for other systems
                        subprocess.run([text_editor, filepath])
                else:  # Fallback for other file types
                    if os_type == "Linux":
                        subprocess.run(['xdg-open', filepath])
                    elif os_type == "Darwin":  # macOS
                        subprocess.run(['open', filepath])
                    elif os_type == "Windows":
                        os.startfile(filepath)
                    else:
                        subprocess.run(['xdg-open', filepath])

        return _open_file

    #####################################################################################################
    # Code block: Define child canvas
    #####################################################################################################

    def child_load_files(start, end, frame, canvas, source_dir, current_win):
        canvas.yview_moveto(0)  # Reset scrollbar to top
        if not os.path.exists(source_dir) or not os.listdir(source_dir):
            label = Label(frame, text="No files found here, try another folder.", bg='white')
            label.pack(pady=20)
            return
        # Sort files
        sorted_files = [f for f in sorted(os.listdir(source_dir))]
        for i, filename in enumerate(sorted_files[start:end], start=start):
            # Add line number into canvas frame
            line_number = Label(frame, text=str(i + 1), bg='white')
            line_number.grid(row=i - start, column=0)

            # Add file name button into canvas frame
            file_button = Button(frame, text=filename, anchor='w', bg='white')
            file_button.grid(row=i - start, column=1, sticky='ew')
            # Bind with child_open_file function to open file
            file_button.bind('<Double-Button-1>', child_open_file(filename, file_button, source_dir))

            # Build a child button_frame inside frame
            button_frame = Frame(frame, bg='white')
            button_frame.grid(row=i - start, column=2, sticky='e')

            # Create "Consensus" button inside button_frame
            copy_button = Button(button_frame, text="Consensus", bg='white', fg='black')
            copy_button.grid(row=0, column=0, padx=5)
            # Bind "Consensus" button with copy_file function with specific source and destination folder
            copy_button.bind('<Button-1>', copy_file(filename, copy_button, consensus_folder,
                                                     source_dir, current_win))

            # Define "Extension" button
            more_extend_button = Button(button_frame, text="Extension", bg='white', fg='black')
            more_extend_button.grid(row=0, column=1, padx=5)
            # Bind "Extension" button with copy_file function with different destination folder
            more_extend_button.bind('<Button-1>', copy_file(filename, more_extend_button, need_more_extension,
                                                            source_dir, current_win))

            # Define "Others" button
            others_button = Button(button_frame, text="Others", bg='white', fg='black')
            others_button.grid(row=0, column=2, padx=5)
            # Bind "Extension" button with copy_file function with different destination folder
            others_button.bind('<Button-1>', copy_file(filename, more_extend_button, others_dir,
                                                       source_dir, current_win))

            button_frame.grid_columnconfigure(0, weight=1)
            button_frame.grid_rowconfigure(0, weight=1)
            frame.grid_columnconfigure(1, weight=1)
            frame.grid_columnconfigure(2, weight=0)

    def child_open_file(filename, button, source_dir):
        def _child_open_file(event):
            filepath = os.path.join(source_dir, filename)
            if os_type == "Darwin":
                button.config(fg='red')  # Change button text color under macOS system
                button.update_idletasks()  # Update UI immediately
            else:
                button.config(bg='yellow')  # Change button color
                button.update_idletasks()  # Update UI immediately

            if filename.lower().endswith(('.fa', '.fasta')):
                if os_type == "Windows":
                    subprocess.run(["java", "-jar", aliview_path, filepath])
                else:
                    subprocess.run([aliview_path, filepath])
            elif filename.lower().endswith('.pdf'):
                if os_type == "Linux":
                    subprocess.run(['xdg-open', filepath])
                elif os_type == "Darwin":  # macOS
                    subprocess.run(['open', filepath])
                elif os_type == "Windows":
                    os.startfile(filepath)
                else:
                    subprocess.run(['xdg-open', filepath])
            elif filename.lower().endswith(('.txt', '.py', '.csv', '.md', '.bed')):
                if os_type == "Linux":
                    text_editor = 'gedit'  # Replace 'gedit' with your preferred text editor
                    subprocess.run([text_editor, filepath])
                elif os_type == "Darwin":  # macOS
                    subprocess.run(['open', '-a', 'TextEdit', filepath])
                elif os_type == "Windows":
                    notepad_path = 'notepad.exe'  # or path to another text editor if preferred
                    subprocess.run([notepad_path, filepath])
                else:
                    text_editor = 'gedit'  # Fallback for other systems
                    subprocess.run([text_editor, filepath])
            else:  # Fallback for other file types
                if os_type == "Linux":
                    subprocess.run(['xdg-open', filepath])
                elif os_type == "Darwin":  # macOS
                    subprocess.run(['open', filepath])
                elif os_type == "Windows":
                    os.startfile(filepath)
                else:
                    subprocess.run(['xdg-open', filepath])

        return _child_open_file

    def open_folder(folder_n, source_dir):

        # Create a new top-level window
        folder_window = Toplevel()
        folder_window.title(folder_n)
        folder_window.geometry('800x600')

        # Create canvas on the new window
        folder_canvas = Canvas(folder_window, bg='white')
        folder_canvas.pack(side='left', fill='both', expand=True)

        # Set scrollbar for the canvas
        folder_scrollbar = Scrollbar(folder_window, orient='vertical', command=folder_canvas.yview)
        folder_canvas.configure(yscrollcommand=folder_scrollbar.set)
        folder_scrollbar.pack(side='right', fill='y')

        # Set scrollbar control
        folder_window.bind("<MouseWheel>", lambda event: scroll(event, folder_canvas))
        folder_window.bind("<Button-4>", lambda event: scroll(event, folder_canvas))
        folder_window.bind("<Button-5>", lambda event: scroll(event, folder_canvas))

        # Create frame for the canvas
        folder_frame = Frame(folder_canvas, bg='white')
        folder_canvas_frame = folder_canvas.create_window((0, 0), window=folder_frame, anchor='nw')

        # Load file, show maximum 1000 file
        child_load_files(0, 1000, folder_frame, folder_canvas, os.path.join(source_dir, folder_n), folder_window)

        # Bind events to the new window's canvas and frame
        folder_canvas.bind('<Configure>', lambda event: on_configure(event, folder_canvas, folder_canvas_frame))
        folder_frame.bind('<Configure>', lambda event: update_scrollregion(event, folder_canvas))

        folder_window.mainloop()

    #####################################################################################################
    # Code block: Add menu bar for mother window
    #####################################################################################################

    # Create mother menu
    menubar = Menu(root)

    # Show menu on the window
    root.configure(menu=menubar)

    annotation_folders = ["Clustered_proof_annotation", "Low_copy_TE", "Skipped_TE"]

    # Create sub-menu
    for annotation in annotation_folders:
        annotationMenu = Menu(menubar, tearoff=0)
        menubar.add_cascade(label=annotation, menu=annotationMenu)
        annotation_path = os.path.join(te_trimmer_proof_annotation_dir, annotation)

        # Give hits when folder isn't found
        if not os.path.exists(annotation_path):
            annotationMenu.add_command(label="Folder not detected! Use correct input folder path",
                                       command=partial(messagebox.showerror,
                                                       "Error", "Please use the correct input directory."
                                                                " three folder should be contained in your input path, "
                                                                "including 'Clustered_proof_annotation', "
                                                                "'Skipped_TE', and 'Low_copy_TE'"))
            continue

        # Sort files inside each folder
        sorted_files_annotation = sorted(os.listdir(annotation_path))

        # Give No files found when the folder is empty
        if not sorted_files_annotation:
            annotationMenu.add_command(label="No files found here, try another folder.")
        else:
            # Show maximum 1000 files on the window
            for i in range(0, len(sorted_files_annotation), 1000):
                end = min(i + 1000, len(sorted_files_annotation))
                annotationMenu.add_command(label=f"{i + 1}-{end}",
                                           command=partial(load_files_with_destroy,
                                                           i, end, frame, canvas,annotation_path))

    # Add confirm menu button
    confirm_menu = Menu(menubar, tearoff=0)
    # Set confirm_menu child menu
    # Initialize the BooleanVar to store "Show Confirmation Window" status
    show_confirmation = BooleanVar(value=True)
    confirm_menu.add_checkbutton(label="Show Confirmation Window", onvalue=True, offvalue=False,
                                 variable=show_confirmation)
    menubar.add_cascade(label="Settings", menu=confirm_menu)

    # Add Undo button
    undo_menu = Menu(menubar, tearoff=0)
    # Add Undo child menu
    undo_menu.add_command(label="undo last copy", command=undo_last_copy)
    menubar.add_cascade(label="Undo", menu=undo_menu)

    root.mainloop()


if __name__ == '__main__':
    proof_annotation()
