import subprocess
import sys


def install_and_import(required_packages_dict):
    for package in required_packages_dict:
        try:
            __import__(package)
        except ImportError:
            try:
                print(f"{package} was not found. Installing it automatically.")
                subprocess.check_call([sys.executable, "-m", "pip", "install", required_packages_dict[package]])
                print(f"{package} was successfully installed.")
            except subprocess.CalledProcessError as e:
                print(
                    f"\nRequired Python packages are missing and cannot be installed automatically. Installation failed with error {e.stderr}"
                    "\nPlease install 'click' and 'biopython' using 'pip install'.\n")
                return


required_packages = {'click': 'click', 'Bio': 'biopython', 'numpy': 'numpy', 'pandas': 'pandas', 'plotly': 'plotly'}
install_and_import(required_packages)


import os
import re
import shutil
import subprocess
from tkinter import Tk, Frame, Button, messagebox, Scrollbar, Canvas, Label, Menu, BooleanVar, Toplevel, simpledialog, Text, Entry
import click
import traceback
from functools import partial
import platform
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Plotter import GUI_plotter
# Import cleaning module
from crop_end_divergence import crop_end_div
from crop_end_gap import crop_end_gap
from remove_gap import remove_gaps_with_similarity_check

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
        click.echo("TEtrimmer don't have right to change permissions. Pleas use sudo to run TEtrimmer")
        return False
    return True


# Detect system OS type
os_type = platform.system()

# Define Aliview software path and change permission
bin_py_path = os.path.dirname(os.path.abspath(__file__))
aliview_folder = os.path.join(bin_py_path, "aliview")
change_permissions_recursive(aliview_folder, 0o755)

aliview_path = os.path.join(bin_py_path, "aliview/aliview")
if os_type == "Windows":
    aliview_path = os.path.join(bin_py_path, r"aliview\aliview.jar")

#####################################################################################################
# Code block: set click command
#####################################################################################################

@click.command()
@click.option('--te_trimmer_proof_curation_dir', '-i', default=None, type=str,
              help='Define TEtrimmer proof curation output path. '
                   'Like <TEtrimmer_output_path>/<TEtrimmer_for_proof_curation>')
@click.option('--output_dir', '-o', default=None, type=str,
              help='Output directory. Default: input directory')
@click.option('--genome_file', '-g', required=True, type=str,
              help='Genome fasta file path.')
def proof_curation(te_trimmer_proof_curation_dir, output_dir, genome_file):
    """
    This tool can help do quick proof curation

    python ./annoGUI.py -i <TEtrimmer_output_folder> -g <genome_file>
    """

    # Make directory for temporary files
    temp_folder = os.path.join(bin_py_path, "temp_folder")
    os.makedirs(temp_folder, exist_ok=True)

    # Define empty list to store copy history, which enable undo button
    copy_history = []

    # Define cleaning module global parameters
    # Define cleaning module global parameters
    global crop_div_thr_g, crop_div_win_g, crop_gap_thr_g, crop_gap_win_g
    global column_gap_thr_g, simi_check_gap_thr_g, similarity_thr_g, min_nucleotide_g

    crop_div_thr_g = 0.65
    crop_div_win_g = 40
    crop_gap_thr_g = 0.05
    crop_gap_win_g = 150
    column_gap_thr_g = 0.8
    simi_check_gap_thr_g = 0.4
    similarity_thr_g = 0.7
    min_nucleotide_g = 5

    # If the -i option is None define the default input directory
    if te_trimmer_proof_curation_dir is None:
        te_trimmer_proof_curation_dir = os.path.abspath(os.path.join(bin_py_path, os.pardir))

    # If the -o option is not given, use the parent directory of -i as output directory.
    if output_dir is None:
        output_dir = os.path.join(te_trimmer_proof_curation_dir, "TEtrimmer_proof_anno_results")
    # Define output folders, create them when they are not found
    consensus_folder = os.path.abspath(os.path.join(output_dir, "Proof_curation_consensus_folder"))
    need_more_extension = os.path.abspath(os.path.join(output_dir, "Proof_curation_need_more_extension"))
    others_dir = os.path.abspath(os.path.join(output_dir, "Proof_curation_discard"))
    low_copy_elements = os.path.abspath(os.path.join(output_dir, "Proof_curation_low_copy_elements"))
    rescue_skip_elements = os.path.abspath(os.path.join(output_dir, "Proof_curation_rescued_skip_elements"))

    for dir_path in [consensus_folder, need_more_extension, others_dir, low_copy_elements, rescue_skip_elements]:
        os.makedirs(dir_path, exist_ok=True)

    #####################################################################################################
    # Code block: build TKinter window
    #####################################################################################################

    # Initialize Tk window
    root = Tk()
    root.title("TEtrimmer Proof Curation Tool")
    root.geometry('900x750')

    # Create canvas on root
    canvas = Canvas(root, bg='white')

    # fill='both'  This tells the canvas to expand and fill both the X-axis (horizontally) and Y-axis (vertically)
    # in its parent widget. The canvas will take up as much space as possible in both directions.
    # expand=True  If the window is resized, the canvas will grow or shrink accordingly.
    canvas.pack(side='left', fill='both', expand=True)

    #####################################################################################################
    # Code block: get genome length file
    #####################################################################################################
    # Check if genome length file available, otherwise create it
    def calculate_genome_length(genome_file, genome_length_output):
        """
        Calculate the length of each sequence in a genome file in FASTA format
        and write the lengths to an output file.

        :param genome_file: str, path to genome file in FASTA format
        :return: str, path to the output file containing sequence names and lengths
        """
        genome_lengths = {}
        with open(genome_file, "r") as f:
            current_seq = None
            current_length = 0
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_seq is not None:
                        genome_lengths[current_seq] = current_length
                    # If chromosome header contains empty spaces, only consider the content before the first space
                    current_seq = line[1:].split(" ")[0]
                    current_length = 0
                else:
                    current_length += len(line)
            # Add length of last sequence to dictionary
            if current_seq is not None:
                genome_lengths[current_seq] = current_length

        # Write lengths to output file
        with open(genome_length_output, "w") as out:
            for seq_name, length in genome_lengths.items():
                out.write(f"{seq_name}\t{length}\n")

        return genome_length_output

    def read_genome_lengths(genome_length_file):

        chrom_sizes = {}
        with open(genome_length_file, 'r') as file:
            for line in file:
                parts = line.strip().split()
                chrom_sizes[parts[0]] = int(parts[1])
        return chrom_sizes

    # Check is genome length file exist otherwise create it
    # Define genome length file name
    genome_length_f = os.path.join(bin_py_path, f"{os.path.basename(genome_file)}.length")
    try:
        if os.path.isfile(genome_length_f):
            chrom_size = read_genome_lengths(genome_length_f)
        else:
            # Generate genome length file when it can't be found
            genome_length_f = calculate_genome_length(genome_file, genome_length_f)
            chrom_size = read_genome_lengths(genome_length_f)
    except Exception as e:
        click.echo(f"\nError while generate genome length file. Error: {str(e)}\n")
        return

    #####################################################################################################
    # Code block: extension module
    #####################################################################################################
    def fresh_child_canvas(child_frame, child_canvas, child_source_dir, current_win):
        # Record the current scroll position and button states
        scroll_position = child_canvas.yview()[0]

        button_states = {}

        for row in range(len(child_frame.grid_slaves(column=0))):
            row_widgets = child_frame.grid_slaves(row=row)

            # Before sort looks like
            # [<tkinter.Frame object .!toplevel.!canvas.!frame.!frame>, <tkinter.Button object .!toplevel.!canvas.!frame.!button>, <tkinter.Label object .!toplevel.!canvas.!frame.!label>]
            # After sort looks like
            # [<tkinter.Label object .!toplevel.!canvas.!frame.!label>, <tkinter.Button object .!toplevel.!canvas.!frame.!button>, <tkinter.Frame object .!toplevel.!canvas.!frame.!frame>]
            row_widgets.sort(key=lambda widget: widget.grid_info()["column"])

            # This corresponds to file name button
            filename = row_widgets[1].cget("text")
            file_bg = row_widgets[1].cget("bg")

            # Update background dictionary
            button_states[filename] = [file_bg]

            # Get the button frame which is the third widget in the row
            button_frame = row_widgets[2]

            # Iterate over the children of the button frame
            for button in button_frame.winfo_children():
                # Get the button text and background color
                button_bg = button.cget("bg")

                # Update button_states to include other button background
                button_states[filename].append(button_bg)

        # Destroy all widgets in child_frame before reload it
        for widget in child_frame.winfo_children():
            widget.destroy()

        # Reload the child canvas to show the new file
        child_load_files(0, 1000, child_frame, child_canvas, child_source_dir, current_win,
                         scroll_position=scroll_position, button_states=button_states)

    # Generate bed file based on fasta header
    # The fasta header could look like
    """
    > scaffold_1:343622-349068(+)
    > scaffold_1:346171-347467(+)
    > scaffold_1:346385-347665(+)
    > scaffold_1:346200-347196(+)
    """
    def fasta_header_to_bed(input_file, output_file):
        with open(input_file, 'r') as fasta, open(output_file, 'w') as bed:
            for line in fasta:
                if line.startswith('>'):
                    header = line.strip().lstrip('>')  # Remove '>'
                    parts = header.rsplit(':', 2)  # Split from right side
                    scaffold = parts[0]
                    range_strand = parts[1].split('(')
                    range_part = range_strand[0]
                    strand = range_strand[1].split(')')[0]
                    start, end = range_part.split('-')
                    bed.write(f'{scaffold}\t{start}\t{end}\tTEtrimmer\t0\t{strand}\n')
        return output_file

    # Avoid to use bedtools slop and getfasta to make it compatible with Windows system
    def extend_bed_regions(bed_file, left_extension, right_extension, chrom_sizes, output_bed):

        with open(bed_file, 'r') as infile, open(output_bed, 'w') as outfile:
            for line in infile:
                parts = line.strip().split()
                chrom, start, end, strand = parts[0], int(parts[1]), int(parts[2]), parts[5]

                # Adjust extensions based on strand
                if strand == '+':
                    new_start = max(0, start - left_extension)
                    new_end = min(chrom_sizes[chrom], end + right_extension)
                else:  # For anti sense strand
                    new_start = max(0, start - right_extension)  # Extend "start" less, as it's the "end"
                    new_end = min(chrom_sizes[chrom], end + left_extension)  # Extend "end" more, as it's the "start"

                # Write the extended region to the output BED file
                outfile.write(f"{chrom}\t{new_start}\t{new_end}\tTEtrimmer\t0\t{strand}\n")
        return output_bed

    def read_bed(bed_file):

        with open(bed_file, 'r') as file:
            for line in file:
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue
                yield {
                    'chrom': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'strand': parts[5]
                }

    def extract_fasta_from_bed(genome_fasta, bed_file, output_fasta):

        genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, 'fasta'))
        with open(output_fasta, 'w') as out_fasta:
            for region in read_bed(bed_file):
                chrom = region['chrom']
                start = region['start']
                end = region['end']
                strand = region['strand']
                if chrom in genome:
                    sequence = genome[chrom].seq[start:end]
                    if strand == '-':
                        sequence = sequence.reverse_complement()  # Reverse complement if on negative strand
                    seq_record = SeqRecord(sequence, id=f"{chrom}:{start}-{end}({strand})", description="")
                    SeqIO.write(seq_record, out_fasta, 'fasta')
        return output_fasta

    # Combined extension module
    def extension_function(input_fasta_n, ext_button, source_dir, output_dir, parent_win, chrom_s, child_frame,
                           child_canvas, child_source_dir, current_win):
        def _extension_function(event):
            left_ex = simpledialog.askinteger("Input", "Enter left extension length (bp):",
                                              parent=parent_win, minvalue=0, initialvalue=1000)
            right_ex = simpledialog.askinteger("Input", "Enter right extension length (bp):",
                                               parent=parent_win, minvalue=0, initialvalue=1000)

            if input_fasta_n.lower().endswith(('.fa', '.fasta')):
                input_fasta_f = os.path.join(source_dir, input_fasta_n)
                base_name = os.path.splitext(input_fasta_n)[0]
                output_fasta = os.path.join(source_dir, f"{input_fasta_n}_{left_ex}_{right_ex}.fa")

                try:
                    # Generate bed file based on fasta header
                    input_fasta_bed = fasta_header_to_bed(input_fasta_f, os.path.join(output_dir, f"{base_name}.bed"))

                    # Do bed file extension
                    input_fasta_after_ex_bed = extend_bed_regions(input_fasta_bed, left_ex, right_ex, chrom_s,
                                                                  os.path.join(output_dir,
                                                                               f"{base_name}_{left_ex}_{right_ex}.bed"))

                    # Get fasta file based on the extended bed file
                    extract_fasta_from_bed(genome_file, input_fasta_after_ex_bed, output_fasta)

                    if os_type == "Darwin":
                        ext_button.config(fg='red')  # Change button text color under macOS system
                    else:
                        ext_button.config(bg='light green')  # Change button color
                    ext_button.update_idletasks()

                    # Fresh child canvas
                    fresh_child_canvas(child_frame, child_canvas, child_source_dir, current_win)

                except Exception as e:
                    click.echo(f"An error occurred during extension: \n {traceback.format_exc()}")
                    messagebox.showerror("Error", f"An error occurred during extension: {str(e)}", parent=parent_win)

        return _extension_function

    #####################################################################################################
    # Code block: set plotter function
    #####################################################################################################
    def plotter_function(input_fasta_n, button, source_dir, output_dir, genome_file, parent_win):
        def _plotter_function(event):
            if input_fasta_n.lower().endswith(('.fa', '.fasta')):
                input_file = os.path.join(source_dir, input_fasta_n)

                try:
                    GUI_plotter(input_file, output_dir, genome_file)

                    if os_type == "Darwin":
                        button.config(fg='red')  # Change button text color under macOS system
                        button.update_idletasks()  # Update UI immediately
                    else:
                        button.config(bg='light green')  # Change button color
                        button.update_idletasks()  # Update UI immediately
                    button.update_idletasks()

                except Exception as e:
                    click.echo(f"An error occurred during plotting: \n {traceback.format_exc()}")
                    messagebox.showerror("Error", f"Plotting failed. Make sure BLAST is correctly installed: {str(e)}",
                                         parent=parent_win)
        return _plotter_function

    #####################################################################################################
    # Code block: define cleaning functions
    #####################################################################################################

    # Define cleaning functions using global parameters
    def crop_end_div_gui(input_fasta_n, button, source_dir, output_dir, parent_win,
                         child_frame, child_canvas, child_source_dir, current_win):
        def _crop_end_div_gui(event):
            try:
                input_file = os.path.join(source_dir, input_fasta_n)
                output_file = os.path.join(output_dir, f"{input_fasta_n}_cld.fa")

                crop_end_div(input_file, output_file, threshold=crop_div_thr_g, window_size=crop_div_win_g)

                if os_type == "Darwin":
                    button.config(fg='red')  # Change button text color under macOS system
                    button.update_idletasks()  # Update UI immediately
                else:
                    button.config(bg='light green')  # Change button color
                    button.update_idletasks()  # Update UI immediately
                button.update_idletasks()

                # Fresh child canvas
                fresh_child_canvas(child_frame, child_canvas, child_source_dir, current_win)

            except Exception as e:
                click.echo(f"An error occurred for crop end by divergence: \n {traceback.format_exc()}")
                messagebox.showerror("Error", f"MSA cleaning crop end by divergence failed: {str(e)}",
                                     parent=parent_win)

        return _crop_end_div_gui

    def crop_end_gap_gui(input_fasta_n, button, source_dir, output_dir, parent_win,
                         child_frame, child_canvas, child_source_dir, current_win):
        def _crop_end_gap_gui(event):
            try:
                input_file = os.path.join(source_dir, input_fasta_n)
                output_file = os.path.join(output_dir, f"{input_fasta_n}_clg.fa")

                crop_end_gap(input_file, output_file, gap_threshold=crop_gap_thr_g, window_size=crop_gap_win_g)

                if os_type == "Darwin":
                    button.config(fg='red')  # Change button text color under macOS system
                    button.update_idletasks()  # Update UI immediately
                else:
                    button.config(bg='light green')  # Change button color
                    button.update_idletasks()  # Update UI immediately
                button.update_idletasks()

                # Fresh child canvas
                fresh_child_canvas(child_frame, child_canvas, child_source_dir, current_win)

            except Exception as e:
                click.echo(f"An error occurred for crop end by gap: \n {traceback.format_exc()}")
                messagebox.showerror("Error", f"MSA cleaning crop end by gap failed: {str(e)}",
                                     parent=parent_win)

        return _crop_end_gap_gui

    def remove_gaps_with_similarity_check_gui(input_fasta_n, button, source_dir, output_dir, parent_win,
                                              child_frame, child_canvas, child_source_dir, current_win):
        def _remove_gaps_with_similarity_check_gui(event):
            try:
                input_file = os.path.join(source_dir, input_fasta_n)
                output_file = os.path.join(output_dir, f"{input_fasta_n}_gr.fa")

                remove_gaps_with_similarity_check(input_file, output_file, gap_threshold=column_gap_thr_g,
                                                  simi_check_gap_thr=simi_check_gap_thr_g,
                                                  similarity_thr=similarity_thr_g,
                                                  min_nucleotide=min_nucleotide_g)
                if os_type == "Darwin":
                    button.config(fg='red')  # Change button text color under macOS system
                    button.update_idletasks()  # Update UI immediately
                else:
                    button.config(bg='light green')  # Change button color
                    button.update_idletasks()  # Update UI immediately
                button.update_idletasks()

                # Fresh child canvas
                fresh_child_canvas(child_frame, child_canvas, child_source_dir, current_win)

            except Exception as e:
                click.echo(f"An error occurred for cleaning gap columns: \n {traceback.format_exc()}")
                messagebox.showerror("Error", f"MSA cleaning remove gap column failed: {str(e)}",
                                     parent=parent_win)

        return _remove_gaps_with_similarity_check_gui

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


    initial_text = "Manual proof curation is highly recommended to improve the quality of TE annotations.\n\n" \
                   "1, Click the <Clustered_proof_curation> button in the menu bar.\n\n" \
                   "   TEs with more than 90% identity are grouped into one cluster.\n\n" \
                   "2, Click each <Cluster> button.\n\n" \
                   "   For each TE, you can find four files: \n\n" \
                   "     <seq_name.cluster.fa>          multiple sequence alignment (MSA) before clustering\n" \
                   "     <seq_name.raw.fa>               MSA file before cleaning\n" \
                   "     <seq_name.fa>                      MSA file corresponding to TE consensus sequence\n" \
                   "     <seq_name.pdf>                    report plots file for evaluating annotation quality\n\n" \
                   "3, Double click <seq_name.pdf> to evaluate annotation quality.\n\n" \
                   "4, If satisfied, click <Consensus> button next to <seq_name.fa>.\n" \
                   "   This MSA file will be copied to <Proof_curation_consensus_folder>.\n\n" \
                   "5, If not satisfied, check <seq_name.fa> or <seq_name.raw.fa> \n\n" \
                   "6, For more extension, click <Extension> button next to the fasta file\n" \
                   "   Use <TEAid> button to generate interactive report plots. \n\n" \
                   "7, If still not satisfied, check <seq_name.cluster.fa>.\n" \
                   "   Select the sequences you want to use and click <Extension> to find the boundary\n\n" \
                   "8, For skipped and low copy elements, evaluate by checking the PDF file.\n"

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
    # Code block: Define GUI functions
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

    # Set function to allow to modify MSA cleaning parameters
    def show_settings_dialog():
        settings_window = Toplevel(root)
        settings_window.title("Modify Cleaning Parameters")
        settings_window.geometry('300x400')

        def save_settings():
            global crop_div_thr_g, crop_div_win_g, crop_gap_thr_g, crop_gap_win_g
            global column_gap_thr_g, simi_check_gap_thr_g, similarity_thr_g, min_nucleotide_g

            try:
                crop_div_thr_g = float(crop_div_thr_entry.get())
                crop_div_win_g = int(crop_div_win_entry.get())
                crop_gap_thr_g = float(crop_gap_thr_entry.get())
                crop_gap_win_g = int(crop_gap_win_entry.get())
                column_gap_thr_g = float(column_gap_thr_entry.get())
                simi_check_gap_thr_g = float(simi_check_gap_thr_entry.get())
                similarity_thr_g = float(similarity_thr_entry.get())
                min_nucleotide_g = int(min_nucleotide_entry.get())
                settings_window.destroy()

            except ValueError as e:
                messagebox.showerror("Invalid input", f"Please enter valid values. Error: {str(e)}")

        Label(settings_window, text="Crop End Divergence Threshold:").pack(pady=2)
        crop_div_thr_entry = Entry(settings_window, bg='white')
        crop_div_thr_entry.insert(0, str(crop_div_thr_g))
        crop_div_thr_entry.pack(pady=2)

        Label(settings_window, text="Crop End Divergence Window Size:").pack(pady=2)
        crop_div_win_entry = Entry(settings_window, bg='white')
        crop_div_win_entry.insert(0, str(crop_div_win_g))
        crop_div_win_entry.pack(pady=2)

        Label(settings_window, text="Crop End Gap Threshold:").pack(pady=2)
        crop_gap_thr_entry = Entry(settings_window, bg='white')
        crop_gap_thr_entry.insert(0, str(crop_gap_thr_g))
        crop_gap_thr_entry.pack(pady=2)

        Label(settings_window, text="Crop End Gap Window Size:").pack(pady=2)
        crop_gap_win_entry = Entry(settings_window, bg='white')
        crop_gap_win_entry.insert(0, str(crop_gap_win_g))
        crop_gap_win_entry.pack(pady=2)

        Label(settings_window, text="Column Gap Threshold:").pack(pady=2)
        column_gap_thr_entry = Entry(settings_window, bg='white')
        column_gap_thr_entry.insert(0, str(column_gap_thr_g))
        column_gap_thr_entry.pack(pady=2)

        Label(settings_window, text="Similarity Check Gap Threshold:").pack(pady=2)
        simi_check_gap_thr_entry = Entry(settings_window, bg='white')
        simi_check_gap_thr_entry.insert(0, str(simi_check_gap_thr_g))
        simi_check_gap_thr_entry.pack(pady=2)

        Label(settings_window, text="Similarity Threshold:").pack(pady=2)
        similarity_thr_entry = Entry(settings_window, bg='white')
        similarity_thr_entry.insert(0, str(similarity_thr_g))
        similarity_thr_entry.pack(pady=2)

        Label(settings_window, text="Min Nucleotide number:").pack(pady=2)
        min_nucleotide_entry = Entry(settings_window, bg='white')
        min_nucleotide_entry.insert(0, str(min_nucleotide_g))
        min_nucleotide_entry.pack(pady=2)

        Button(settings_window, text="Save", command=save_settings).pack(pady=2)

    # Define function to allow retrieve copy operation
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

                if os_type == "Darwin":
                    last_button.config(fg='black')  # Change button text color under macOS system
                    last_button.update_idletasks()  # Update UI immediately
                else:
                    last_button.config(bg='white')  # Change button color
                    last_button.update_idletasks()  # Update UI immediately

                messagebox.showinfo("Info", f"Successfully removed '{last_copied_file}'.")
            except Exception as e:
                messagebox.showerror("Error", f"An error occurred while deleting the file: {str(e)}")

    # copy_file function is bundled with copy buttons like Consensus
    def copy_file(filename, button, target_directory, source_dir, parent_win):
        def _copy_file(event):

            file_path = os.path.join(source_dir, filename)
            if not show_confirmation.get() or messagebox.askokcancel(
                    "Confirmation", f"Do you want to copy '{filename}' to '{target_directory}'?", parent=parent_win):
                os.makedirs(target_directory, exist_ok=True)
                try:
                    shutil.copy(file_path, target_directory)
                    if os_type == "Darwin":
                        button.config(fg='red')  # Change button text color under macOS system
                        button.update_idletasks()  # Update UI immediately
                    else:
                        button.config(bg='light green')  # Change button color
                        button.update_idletasks()  # Update UI immediately
                    button.update_idletasks()
                    if len(copy_history) >= 100:
                        copy_history.pop(0)
                    copy_history.append((os.path.join(target_directory, filename), target_directory, button))
                except Exception as e:
                    messagebox.showerror("Error", f"An error occurred while copying the file: {str(e)}",
                                         parent=parent_win)
        return _copy_file

    def load_files(start, end, frame, canvas, source_dir=te_trimmer_proof_curation_dir):
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

            # If it's a folder in "Clustered_proof_curation", add a label with the file count
            if source_dir.endswith("Clustered_proof_curation") and os.path.isdir(os.path.join(source_dir, filename)):
                file_count = len(os.listdir(os.path.join(source_dir, filename)))
                file_count_label = Label(frame, text=f"({file_count} files)", bg='white')
                file_count_label.grid(row=i - start, column=3)

            # Don't show "Consensus" and "Extension" button for low copy and clustered proof curation
            if not source_dir.endswith("TE_low_copy") and not source_dir.endswith("Clustered_proof_curation")\
                    and not source_dir.endswith("TE_skipped"):

                # Create "Consensus" button inside button_frame
                copy_button = Button(button_frame, text="Consensus", bg='white', fg='black')
                copy_button.grid(row=0, column=0, padx=5)
                # Bind "Consensus" button with copy_file function with specific source and destination folder
                copy_button.bind('<Button-1>', copy_file(filename, copy_button, consensus_folder,
                                                         source_dir, root))

                # Define "Extension" button
                more_extend_button = Button(button_frame, text="Copy for extension", bg='white', fg='black')
                more_extend_button.grid(row=0, column=1, padx=5)
                # Bind "Extension" button with copy_file function with different destination folder
                more_extend_button.bind('<Button-1>', copy_file(filename, more_extend_button, need_more_extension,
                                                                source_dir, root))
            # Add rescue button for skipped TE
            if os.path.basename(source_dir) == "TE_skipped":
                rescue_skip_button = Button(button_frame, text="Rescue skip", bg='white', fg='black')
                rescue_skip_button.grid(row=0, column=1, padx=5)
                rescue_skip_button.bind('<Button-1>', copy_file(filename, rescue_skip_button, rescue_skip_elements,
                                                                source_dir, root))

            # Add "Save to low copy" button for low copy TE
            if os.path.basename(source_dir) == "TE_low_copy":
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
                open_cluster_folder(filename, source_dir)
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

    def show_help():
        help_window = Toplevel(root)
        help_window.title("Help")
        help_window.geometry('900x700')

        help_text_widget = Text(help_window, bg='white', font=('Arial', 15), wrap='word')
        help_text_widget.insert('1.0', initial_text)
        help_text_widget.config(state='disabled')  # Make the Text widget read-only
        help_text_widget.pack(fill='both', expand=True)

        help_window.update_idletasks()  # Ensure that all widget sizes are calculated

    #####################################################################################################
    # Code block: Define child canvas
    #####################################################################################################

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

    def child_load_files(start, end, frame, canvas, source_dir, current_win, scroll_position=None, button_states=None):
        if scroll_position is not None:
            canvas.yview_moveto(scroll_position)  # Set scrollbar to saved position
        else:
            canvas.yview_moveto(0)  # Reset scrollbar to top if no position is provided

        if button_states is None:
            button_states = {}

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

            # Get button states
            file_button_bg = 'white'
            consensus_button_bg = 'white'
            extension_button_bg = 'white'
            teaid_button_bg = 'white'
            crop_end_by_div_button_bg = 'white'
            crop_end_by_gap_button_bg = 'white'
            remove_gap_column_button_bg = "white"
            others_button_bg = 'white'

            if filename in button_states:
                states = button_states[filename]
                if len(states) >= 1:
                    file_button_bg = states[0]
                if len(states) >= 2:
                    consensus_button_bg = states[1]
                if len(states) >= 3:
                    extension_button_bg = states[2]
                if len(states) >= 4:
                    teaid_button_bg = states[3]
                if len(states) >= 5:
                    crop_end_by_div_button_bg = states[4]
                if len(states) >= 6:
                    crop_end_by_gap_button_bg = states[5]
                if len(states) >= 7:
                    remove_gap_column_button_bg = states[6]
                if len(states) >= 8:
                    others_button_bg = states[7]

            # Add file name button into canvas frame
            file_button = Button(frame, text=filename, anchor='w', bg=file_button_bg)
            file_button.grid(row=i - start, column=1, sticky='ew')
            # Bind with child_open_file function to open file
            file_button.bind('<Double-Button-1>', child_open_file(filename, file_button, source_dir))

            # Build a child button_frame inside frame
            button_frame = Frame(frame, bg='white')
            button_frame.grid(row=i - start, column=2, sticky='e')

            # Create "Consensus" button inside button_frame
            copy_button = Button(button_frame, text="Cons", bg=consensus_button_bg, fg='black')
            copy_button.grid(row=0, column=0, padx=1)
            # Bind "Consensus" button with copy_file function with specific source and destination folder
            copy_button.bind('<Button-1>', copy_file(filename, copy_button, consensus_folder, source_dir, current_win))

            # Define "Extension" button
            more_extend_button = Button(button_frame, text="Extend", bg=extension_button_bg, fg='black')
            more_extend_button.grid(row=0, column=1, padx=1)
            # Bind "Extension" button with copy_file function with different destination folder
            more_extend_button.bind('<Button-1>', extension_function(filename, more_extend_button, source_dir,
                                                                     temp_folder, current_win, chrom_size, frame,
                                                                     canvas, source_dir, current_win))

            # Define "Plotter" button
            plot_button = Button(button_frame, text="TEAid", bg=teaid_button_bg, fg='black')
            plot_button.grid(row=0, column=2, padx=1)
            # Bind "Plotter" button with plotter_function
            plot_button.bind('<Button-1>',
                             plotter_function(filename, plot_button, source_dir, temp_folder, genome_file, current_win))

            # Define "Crop end by divergence" button
            crop_end_by_div_button = Button(button_frame, text="C_div", bg=crop_end_by_div_button_bg, fg='black')
            crop_end_by_div_button.grid(row=0, column=3, padx=1)
            crop_end_by_div_button.bind('<Button-1>',
                                        crop_end_div_gui(filename, crop_end_by_div_button, source_dir, source_dir,
                                                         current_win, frame, canvas, source_dir, current_win))

            # Define "Crop end by gap" button
            crop_end_by_gap_button = Button(button_frame, text="C_gap", bg=crop_end_by_gap_button_bg, fg='black')
            crop_end_by_gap_button.grid(row=0, column=4, padx=1)

            crop_end_by_gap_button.bind('<Button-1>',
                                        crop_end_gap_gui(filename, crop_end_by_gap_button, source_dir, source_dir,
                                                         current_win, frame, canvas, source_dir, current_win))

            # Define "Clean gap column" button
            clean_gap_column_button = Button(button_frame, text="C_col", bg=remove_gap_column_button_bg, fg='black')
            clean_gap_column_button.grid(row=0, column=5, padx=1)

            clean_gap_column_button.bind('<Button-1>', remove_gaps_with_similarity_check_gui(
                filename, clean_gap_column_button, source_dir, source_dir, current_win,
                frame, canvas, source_dir, current_win))

            # Define "Others" button (discard)
            others_button = Button(button_frame, text="Discard", bg=others_button_bg, fg='black')
            others_button.grid(row=0, column=6, padx=1)
            # Bind "Others" button with copy_file function with different destination folder
            others_button.bind('<Button-1>', copy_file(filename, others_button, others_dir, source_dir, current_win))

            button_frame.grid_columnconfigure(0, weight=1)
            button_frame.grid_rowconfigure(0, weight=1)
            frame.grid_columnconfigure(1, weight=1)
            frame.grid_columnconfigure(2, weight=0)

    def open_cluster_folder(folder_n, source_dir):
        # Create a new top-level window
        folder_window = Toplevel()
        folder_window.title(folder_n)
        folder_window.geometry('1100x600')

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

        # Load file, show maximum 1000 files
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

    annotation_folders = ["Clustered_proof_curation", "TE_low_copy", "TE_skipped"]

    # Create sub-menu
    for annotation in annotation_folders:
        annotationMenu = Menu(menubar, tearoff=0)
        menubar.add_cascade(label=annotation, menu=annotationMenu)
        annotation_path = os.path.join(te_trimmer_proof_curation_dir, annotation)

        # Give hits when folder isn't found
        if not os.path.exists(annotation_path):
            annotationMenu.add_command(label="Folder not detected! Use correct input folder path",
                                       command=partial(messagebox.showerror,
                                                       "Error", "Please use the correct input directory."
                                                                " three folder should be contained in your input path, "
                                                                "including 'Clustered_proof_curation', "
                                                                "'TE_skipped', and 'TE_low_copy'"))
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

    # Enable setting cleaning module parameters
    settings_menu = Menu(menubar, tearoff=0)
    settings_menu.add_command(label="Modify Parameters", command=show_settings_dialog)
    menubar.add_cascade(label="Modify parameters", menu=settings_menu)

    # Add Undo button
    undo_menu = Menu(menubar, tearoff=0)
    # Add Undo child menu
    undo_menu.add_command(label="undo last copy", command=undo_last_copy)
    menubar.add_cascade(label="Undo", menu=undo_menu)

    # Add the Help button to the menu
    help_menu = Menu(menubar, tearoff=0)
    help_menu.add_command(label="Show instruction", command=show_help)
    menubar.add_cascade(label="Help", menu=help_menu)

    root.mainloop()


if __name__ == '__main__':
    proof_curation()
