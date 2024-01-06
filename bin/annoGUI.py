import os
import shutil
import subprocess
from tkinter import Tk, Frame, Button, messagebox, Scrollbar, Canvas, Label, Menu, BooleanVar
import click
from functools import partial
import platform


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

copy_history = []

# Detect OS
os_type = platform.system()


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
    low_copy_elements = os.path.abspath(os.path.join(output_dir, "Proof_annotation_low_copy_elements"))

    for dir_path in [consensus_folder, need_more_extension, low_copy_elements]:
        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)

    # Initialize Tk window
    root = Tk()
    root.title("TETrimmer proof annotation tool")
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

    # Use Aliview to open fasta file
    def open_file(filename, button, source_dir=te_trimmer_proof_annotation_dir):
        def _open_file(event):
            filepath = os.path.join(source_dir, filename)
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
            button.config(bg='yellow')

        return _open_file

    def disable_confirmation():
        show_confirmation.set(False)

    def copy_file(filename, button, target_directory, source_dir=te_trimmer_proof_annotation_dir):
        def _copy_file(event):
            file_path = os.path.join(source_dir, filename)
            if not show_confirmation.get() or messagebox.askokcancel("Confirmation", f"Do you want to copy '{filename}' to '{target_directory}'?"):
                os.makedirs(target_directory, exist_ok=True)
                try:
                    shutil.copy(file_path, target_directory)
                    button.config(bg='green')
                    button.update_idletasks()
                    if len(copy_history) >= 100:
                        copy_history.pop(0)
                    copy_history.append((os.path.join(target_directory, filename), target_directory, button))
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
                copy_history.pop()
                last_button.config(bg='white')
                messagebox.showinfo("Info", f"Successfully removed '{last_copied_file}'.")
            except Exception as e:
                messagebox.showerror("Error", f"An error occurred while deleting the file: {str(e)}")

    def clear_frame():
        for widget in frame.winfo_children():
            widget.destroy()

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
                   "We highly recommend to do proof annotation for 'Recommend_check_annotation' and 'Need_check_annotation" \
                   " ,which can dramatically increase your TE annotation quality.\n\n"\
                   "1, Click the buttons in the menu bar, which corresponds to the different annotation status.\n\n" \
                   "2, All files in the chose folder (button) will be displayed.\n\n" \
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
        logo_font = ('Courier', 10)
    else:
        logo_font = ('Courier', 5)
    logo_label = Label(canvas, text=log_text, bg='white', font=logo_font, justify='left', wraplength=1100)
    logo_label.pack(pady=10)

    # Display the explanatory text with 'Arial' font
    text_label = Label(canvas, text=initial_text, bg='white', font=('Arial', 15), justify='left', wraplength=1100)
    text_label.pack(pady=10)

    def destroy_initial_label():
        global logo_label
        if logo_label:
            logo_label.destroy()
            logo_label = None
        global text_label
        if text_label:
            text_label.destroy()
            text_label = None

    def load_files_with_destroy(start, end, path):
        destroy_initial_label()
        load_files(start, end, source_dir=path)

    def load_files(start, end, source_dir=te_trimmer_proof_annotation_dir):
        clear_frame()
        canvas.yview_moveto(0)  # Reset scrollbar to top
        if not os.path.exists(source_dir) or not os.listdir(source_dir):
            label = Label(frame, text="No files found here, try another folder.", bg='white')
            label.pack(pady=20)
            return
        sorted_files = [f for f in sorted(os.listdir(source_dir))]
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
                copy_button.bind('<Button-1>', copy_file(filename, copy_button, consensus_folder, source_dir=source_dir))
                more_extend_button = Button(button_frame, text="Extension", bg='white', fg='black')
                more_extend_button.grid(row=0, column=1, padx=5)
                more_extend_button.bind('<Button-1>', copy_file(filename, more_extend_button, need_more_extension, source_dir=source_dir))
            else:
                low_copy_button = Button(button_frame, text="Low_copy", bg='white', fg='black')
                low_copy_button.grid(row=0, column=3, padx=5)
                low_copy_button.bind('<Button-1>', copy_file(filename, low_copy_button, low_copy_elements, source_dir=source_dir))
            button_frame.grid_columnconfigure(0, weight=1)
            button_frame.grid_rowconfigure(0, weight=1)
            frame.grid_columnconfigure(1, weight=1)
            frame.grid_columnconfigure(2, weight=0)

    menubar = Menu(root)
    root.config(menu=menubar)

    annotation_folders = ["Perfect_annotation", "Good_annotation", "Recommend_check_annotation", "Need_check_annotation", "Low_copy_TE"]

    for annotation in annotation_folders:
        annotationMenu = Menu(menubar, tearoff=0)
        menubar.add_cascade(label=annotation, menu=annotationMenu)
        annotation_path = os.path.join(te_trimmer_proof_annotation_dir, annotation)
        if not os.path.exists(annotation_path):
            annotationMenu.add_command(label="Folder not detected! Use correct input folder path",
                                       command=partial(messagebox.showerror,
                                                       "Error", "Please use the correct input directory."
                                                       " Five folder should be contained in your input path, "
                                                       "including 'Perfect_annotation', 'Good_annotation', "
                                                       "'Recommend_check_annotation', 'need_check_annotation', "
                                                       "and 'Low_copy_TE'"))
            continue
        sorted_files_annotation = sorted(os.listdir(annotation_path))
        if not sorted_files_annotation:
            annotationMenu.add_command(label="No files found here, try another folder.")
        else:
            for i in range(0, len(sorted_files_annotation), 1000):
                end = min(i + 1000, len(sorted_files_annotation))
                annotationMenu.add_command(label=f"{i + 1}-{end}", command=partial(load_files_with_destroy, i, end, annotation_path))

    confirm_menu = Menu(menubar, tearoff=0)
    menubar.add_cascade(label="Settings", menu=confirm_menu)
    confirm_menu.add_checkbutton(label="Show Confirmation Window", onvalue=True, offvalue=False, variable=show_confirmation)
    confirm_menu.add_command(label="Don't Show Confirmation Window", command=disable_confirmation)
    menubar.add_command(label="Undo", command=undo_last_copy)

    canvas_frame = canvas.create_window((0, 0), window=frame, anchor='nw')

    root.mainloop()


if __name__ == '__main__':
    proof_annotation()
