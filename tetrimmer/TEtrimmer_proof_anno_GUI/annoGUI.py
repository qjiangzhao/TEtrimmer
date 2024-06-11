import subprocess
import sys
import os
import re
import shutil

# Insert the bin directory to sys.path
bin_path = os.path.join(os.path.dirname(__file__), 'bin')
if bin_path not in sys.path:
    sys.path.insert(1, bin_path)


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

from tkinter import Tk, Frame, Button, messagebox, Scrollbar, Canvas, Label, Menu, BooleanVar, \
    Toplevel, simpledialog, Text, Entry, ttk
import click
import traceback
from functools import partial
import platform
from Bio import SeqIO
from Plotter import GUI_plotter, con_generater
# Import cleaning module
from crop_end_divergence import crop_end_div
from crop_end_gap import crop_end_gap
from remove_gap import remove_gaps_with_similarity_check
from GUI_functions import separate_sequences, blast, fasta_header_to_bed, extend_bed_regions, \
    extract_fasta_from_bed
from cialign_plot import drawMiniAlignment

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
              help='TEtrimmer proof curation output path.'
                   'Like <TEtrimmer_output_path>/TEtrimmer_for_proof_curation. Three folders should exist in the '
                   'given path including "Clustered_proof_curation", "TE_low_copy", and "TE_skipped". '
                   'If you start the "annoGUI.py" from TETrimmer output directory, you do not need to use this option')
@click.option('--output_dir', '-o', default=None, type=str,
              help='Output directory. Default: input directory')
@click.option('--genome_file', '-g', required=True, type=str,
              help='Genome fasta file path.')
@click.option('--consensus_lib', '-clib', default=None, type=str,
              help='TE consensus library fasta file. You can check and improve other TE consensus '
                   'library e.g. The TE library directly from EDTA2, RepeatModeler2, and other tools. If you want to '
                   'check the same TE library as last time, you do not need to use this option again.')
def proof_curation(te_trimmer_proof_curation_dir, output_dir, genome_file, consensus_lib):
    """
    This tool can help do quick proof curation

    python ./annoGUI.py -i <TEtrimmer_output_folder> -g <genome_file>
    """

    # Make directory for temporary files
    temp_folder = os.path.join(bin_py_path, "temp_folder")
    os.makedirs(temp_folder, exist_ok=True)

    # Define directories for other consensus library checking
    other_cons_lib_folder = os.path.join(bin_py_path, "temp_consensus_lib_check")
    other_cons_lib_single_file_folder = os.path.join(other_cons_lib_folder, "Single_files")

    # Define empty list to store copy history, which enable undo button
    copy_history = []

    # Define cleaning module, blast, and consensus generation global parameters
    global crop_div_thr_g, crop_div_win_g, crop_gap_thr_g, crop_gap_win_g
    global column_gap_thr_g, simi_check_gap_thr_g, similarity_thr_g, min_nucleotide_g, cons_thre_g, blast_e_value_g

    crop_div_thr_g = 0.65
    crop_div_win_g = 40
    crop_gap_thr_g = 0.05
    crop_gap_win_g = 150
    column_gap_thr_g = 0.8
    simi_check_gap_thr_g = 0.4
    similarity_thr_g = 0.7
    min_nucleotide_g = 5
    cons_thre_g = 0.8
    blast_e_value_g = 1e-40

    # Define global variable to determine the current canvas contens
    global current_canvas_content

    current_canvas_content = "tetrimmer_out"

    # If the -i option is None define the default input directory
    if te_trimmer_proof_curation_dir is None:
        te_trimmer_proof_curation_dir = os.path.abspath(os.path.join(bin_py_path, os.pardir))

    # If the -o option is not given, use the parent directory of -i as output directory.
    if output_dir is None:
        output_dir = os.path.join(te_trimmer_proof_curation_dir, "TEtrimmer_proof_anno_results")
    # Define output folders, create them when they are not found
    consensus_folder = os.path.abspath(os.path.join(output_dir, "Proof_curation_consensus_folder"))
    others_dir = os.path.abspath(os.path.join(output_dir, "Proof_curation_discard"))

    other_cons_lib_result_folder = os.path.join(te_trimmer_proof_curation_dir, "Consensus_lib_check_results")

    for dir_path in [consensus_folder, others_dir, other_cons_lib_result_folder]:
        os.makedirs(dir_path, exist_ok=True)

    #####################################################################################################
    # Code block: build TKinter window
    #####################################################################################################

    # Initialize Tk window
    root = Tk()
    root.title("TEtrimmer Proof Curation Tool")
    # width * height
    root.geometry('900x800')

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
    # current_win contains child_canvas
    # child_canvas contains child_frame
    # for each row of child_frame, it contains three elements: label, file name frame (include file name button), and
    # another frame (contains function buttons like "Cons", "Extend", "TEAid".)
    # fresh_save_path is used to define save button path. By defaule it save to other_cons_lib_result_folder
    def fresh_canvas(child_frame, child_canvas, child_source_dir, current_win, file_start=0, file_end=500,
                     update_child_canvas=True, fresh_save_path=None):
        # Record the current scroll position and button states
        scroll_position = child_canvas.yview()[0]

        # button_states dictionary stores all button background and text color for each file name
        button_states = {}

        for row in range(len(child_frame.grid_slaves(column=0))):
            row_widgets = child_frame.grid_slaves(row=row)

            # Before sort looks like
            # [<tkinter.Frame object .!toplevel.!canvas.!frame.!frame>, <tkinter.Button object .!toplevel.!canvas.!frame.!button>,
            # <tkinter.Label object .!toplevel.!canvas.!frame.!label>]
            # After sort looks like
            # [<tkinter.Label object .!toplevel.!canvas.!frame.!label>, <tkinter.Button object .!toplevel.!canvas.!frame.!button>,
            # <tkinter.Frame object .!toplevel.!canvas.!frame.!frame>]
            row_widgets.sort(key=lambda widget: widget.grid_info()["column"])

            # This corresponds to file name button
            filename = row_widgets[1].cget("text")  # Get text content
            file_bg = row_widgets[1].cget("bg")  # Get button background color
            file_fg = row_widgets[1].cget("fg")  # Get text color

            # Update background dictionary
            button_states[filename] = [[file_bg, file_fg]]

            # Get the button frame which is the third widget in the row
            button_frame = row_widgets[2]

            # Iterate over the children of the button frame
            for button in button_frame.winfo_children():
                # Get the button text color and background color
                button_bg = button.cget("bg")
                button_fg = button.cget("fg")
                button_list = [button_bg, button_fg]

                # Update button_states to include other button background
                button_states[filename].append(button_list)

        # Destroy all widgets in child_frame before reload it
        for widget in child_frame.winfo_children():
            widget.destroy()

        if update_child_canvas:
            # Reload the child canvas to show the new file
            child_load_files(file_start, file_end, child_frame, child_canvas, child_source_dir, current_win,
                             scroll_position=scroll_position, button_states=button_states)
        else:
            other_cons_load_files(file_start, file_end, child_frame, child_canvas, child_source_dir, current_win,
                                  save_path=fresh_save_path,
                                  scroll_position=scroll_position, button_states=button_states)

    # Combined extension module
    # source_dir is the folder path contains all loaded files, like the Cluster folder path
    # output_dir is the path to store all intermediate files like temp_folder
    # child_canvas contain child_frame
    # current_win contains child_frame and child_canvas
    def extension_function(input_fasta_n, ext_button, source_dir, output_dir, current_win, chrom_s, child_frame,
                           child_canvas, update_child_canvas=True, file_start=0, file_end=500, save_path=None):
        def _extension_function(event):
            if input_fasta_n.lower().endswith(('.fa', '.fasta')):

                left_ex = simpledialog.askinteger("Input", "Enter left extension length (bp):",
                                                  parent=current_win, minvalue=0, initialvalue=1000)
                if left_ex is None:  # If the user clicks cancel, stop the function
                    return

                right_ex = simpledialog.askinteger("Input", "Enter right extension length (bp):",
                                                   parent=current_win, minvalue=0, initialvalue=1000)
                if right_ex is None:  # If the user clicks cancel, stop the function
                    return

                progress = show_progress_bar(child_canvas)  # Show progress bar

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
                        fresh_canvas(child_frame, child_canvas, source_dir, current_win,
                                     update_child_canvas=update_child_canvas, file_start=file_start, file_end=file_end,
                                     fresh_save_path=save_path)

                    except Exception as e:
                        click.echo(f"An error occurred during extension: \n {traceback.format_exc()}")
                        messagebox.showerror("Error", f"An error occurred during extension: {str(e)}", parent=current_win)

                    finally:
                        progress.stop()
                        progress.pack_forget()  # Hide the progress bar
            else:
                return

        return _extension_function

    #####################################################################################################
    # Code block: set plotter function
    #####################################################################################################
    def plotter_function(input_fasta_n, button, source_dir, output_dir, genome_file, child_canvas, current_win):
        def _plotter_function(event):

            if input_fasta_n.lower().endswith(('.fa', '.fasta')):
                progress = show_progress_bar(child_canvas)  # Show progress bar
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
                    messagebox.showerror("Error", f"Plotting failed, please align sequences first: {str(e)}",
                                         parent=current_win)
                finally:
                    progress.stop()
                    progress.pack_forget()  # Hide the progress bar
            else:
                return
        return _plotter_function

    #####################################################################################################
    # Code block: define cleaning functions
    #####################################################################################################

    # Define cleaning functions using global parameters
    def crop_end_div_gui(input_fasta_n, button, source_dir, output_dir, current_win, child_frame, child_canvas,
                         update_child_canvas=True, file_start=0, file_end=500, save_path=None):
        def _crop_end_div_gui(event):
            if input_fasta_n.lower().endswith(('.fa', '.fasta')):
                progress = show_progress_bar(child_canvas)  # Show progress bar

                try:
                    input_file = os.path.join(source_dir, input_fasta_n)
                    output_file = os.path.join(output_dir, f"{input_fasta_n}_CDiv.fa")

                    crop_end_div(input_file, output_file, threshold=crop_div_thr_g, window_size=crop_div_win_g)

                    if os_type == "Darwin":
                        button.config(fg='red')  # Change button text color under macOS system
                        button.update_idletasks()  # Update UI immediately
                    else:
                        button.config(bg='light green')  # Change button color
                        button.update_idletasks()  # Update UI immediately
                    button.update_idletasks()

                    # Fresh child canvas
                    fresh_canvas(child_frame, child_canvas, source_dir, current_win,
                                 update_child_canvas=update_child_canvas, file_start=file_start, file_end=file_end,
                                 fresh_save_path=save_path)

                except Exception as e:
                    click.echo(f"An error occurred for crop end by divergence: \n {traceback.format_exc()}")
                    messagebox.showerror("Error", f"MSA cleaning crop end by divergence failed, please align sequences first: {str(e)}",
                                         parent=current_win)

                finally:
                    progress.stop()
                    progress.pack_forget()  # Hide the progress bar

        return _crop_end_div_gui

    def crop_end_gap_gui(input_fasta_n, button, source_dir, output_dir, current_win, child_frame, child_canvas,
                         update_child_canvas=True, file_start=0, file_end=500, save_path=None):
        def _crop_end_gap_gui(event):
            if input_fasta_n.lower().endswith(('.fa', '.fasta')):
                progress = show_progress_bar(child_canvas)  # Show progress bar
                try:
                    input_file = os.path.join(source_dir, input_fasta_n)
                    output_file = os.path.join(output_dir, f"{input_fasta_n}_CGap.fa")

                    crop_end_gap(input_file, output_file, gap_threshold=crop_gap_thr_g, window_size=crop_gap_win_g)

                    if os_type == "Darwin":
                        button.config(fg='red')  # Change button text color under macOS system
                        button.update_idletasks()  # Update UI immediately
                    else:
                        button.config(bg='light green')  # Change button color
                        button.update_idletasks()  # Update UI immediately
                    button.update_idletasks()

                    # Fresh child canvas
                    fresh_canvas(child_frame, child_canvas, source_dir, current_win,
                                 update_child_canvas=update_child_canvas, file_start=file_start, file_end=file_end,
                                 fresh_save_path=save_path)

                except Exception as e:
                    click.echo(f"An error occurred for crop end by gap: \n {traceback.format_exc()}")
                    messagebox.showerror("Error", f"MSA cleaning crop end by gap failed, please align sequences first: {str(e)}",
                                         parent=current_win)
                finally:
                    progress.stop()
                    progress.pack_forget()  # Hide the progress bar

        return _crop_end_gap_gui

    def remove_gaps_with_similarity_check_gui(input_fasta_n, button, source_dir, output_dir, current_win,
                                              child_frame, child_canvas, update_child_canvas=True,
                                              file_start=0, file_end=500, save_path=None):

        def _remove_gaps_with_similarity_check_gui(event):
            if input_fasta_n.lower().endswith(('.fa', '.fasta')):
                progress = show_progress_bar(child_canvas)  # Show progress bar
                try:
                    input_file = os.path.join(source_dir, input_fasta_n)
                    output_file = os.path.join(output_dir, f"{input_fasta_n}_CCol.fa")

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
                    fresh_canvas(child_frame, child_canvas, source_dir, current_win,
                                 update_child_canvas=update_child_canvas, file_start=file_start, file_end=file_end,
                                 fresh_save_path=save_path)

                except Exception as e:
                    click.echo(f"An error occurred for cleaning gap columns: \n {traceback.format_exc()}")
                    messagebox.showerror("Error", f"Remove gappy column failed, please align sequences first, : {str(e)}",
                                         parent=current_win)
                finally:
                    progress.stop()
                    progress.pack_forget()  # Hide the progress bar

        return _remove_gaps_with_similarity_check_gui

    #####################################################################################################
    # Code block: functions related with checking TE consensus library
    #####################################################################################################
    def check_other_cons_lib(consensus_lib_file):

        # Check cons_lib_folder exist, if so remove it when new consensus file is given
        if os.path.isdir(other_cons_lib_folder) and consensus_lib_file is not None:
            shutil.rmtree(other_cons_lib_folder)

        if consensus_lib_file is not None:
            # Create other_cons_lib_folder
            os.makedirs(other_cons_lib_folder, exist_ok=True)

            # Create folder to store separated single files
            os.makedirs(other_cons_lib_single_file_folder, exist_ok=True)

            # Separate other consensus library into single files
            separate_sequences(consensus_lib_file, other_cons_lib_single_file_folder)

    check_other_cons_lib(consensus_lib)

    def check_other_cons_lib_blast(input_fasta_n, blast_button, source_dir, output_dir, current_win,
                                   child_frame, child_canvas, genome_f, e_value=1e-40, update_child_canvas=True,
                                   file_start=0, file_end=500, save_path=None):

        def _check_other_cons_lib_blast(event):
            if input_fasta_n.lower().endswith(('.fa', '.fasta')):

                input_fasta_file = os.path.join(source_dir, input_fasta_n)

                try:
                    # Count the number of sequences in the input file
                    sequence_count = sum(1 for _ in SeqIO.parse(input_fasta_file, "fasta"))

                    if sequence_count != 1:
                        messagebox.showerror("Error",
                                             "Blast can only be performed on files containing exactly one sequence.",
                                             parent=current_win)
                        return

                    # Perform the blast operation if the sequence count is exactly one
                    other_cons_bed, other_cons_blast = blast(input_fasta_file, genome_f, output_dir, e_value=e_value,
                                                             bed_file=True)

                    if other_cons_bed and other_cons_blast:

                        # Allow to cover the old blast file
                        destination_other_cons_blast_file = os.path.join(source_dir, os.path.basename(other_cons_blast))
                        if os.path.exists(destination_other_cons_blast_file):
                            os.remove(destination_other_cons_blast_file)

                        # Move other_cons_blast to source_dir
                        shutil.move(other_cons_blast, source_dir)
                        # Define check other consensus library fasta file derived from the bed file
                        other_cons_fasta = os.path.join(source_dir, f"{os.path.basename(other_cons_blast)}.fa")
                        extract_fasta_from_bed(genome_f, other_cons_bed, other_cons_fasta)

                        if os_type == "Darwin":
                            blast_button.config(fg='red')  # Change button text color under macOS system
                        else:
                            blast_button.config(bg='light green')  # Change button color
                        blast_button.update_idletasks()

                        # Fresh child canvas
                        fresh_canvas(child_frame, child_canvas, source_dir, current_win,
                                     update_child_canvas=update_child_canvas, file_start=file_start, file_end=file_end,
                                     fresh_save_path=save_path)

                except Exception as e:
                    click.echo(f"An error occurred during blast: \n {traceback.format_exc()}")
                    messagebox.showerror("Error", f"Blast failed. Refer to terminal for more information: {str(e)}",
                                         parent=current_win)

        return _check_other_cons_lib_blast

    #####################################################################################################
    # Code block: define consensus sequence generation function
    #####################################################################################################
    def generate_cons(input_fasta_n, cons_button, source_dir, current_win, child_frame, child_canvas,
                      cons_thre=0.8, ambiguous="N", update_child_canvas=True, file_start=0, file_end=500,
                      save_path=None):
        def _generate_cons(envet):

            if input_fasta_n.lower().endswith(('.fa', '.fasta')):

                try:
                    input_fasta_file = os.path.join(source_dir, input_fasta_n)
                    # The consensus file will be added _co.fa at the end of the input file name
                    # use source_dir here to redirect the consensus file to the same folder with the input file
                    con_generater(input_fasta_file, source_dir, threshold=cons_thre, ambiguous=ambiguous)

                    if os_type == "Darwin":
                        cons_button.config(fg='red')  # Change button text color under macOS system
                    else:
                        cons_button.config(bg='light green')  # Change button color
                    cons_button.update_idletasks()

                except Exception as e:
                    click.echo(f"An error for consensus sequence genration: \n {traceback.format_exc()}")
                    messagebox.showerror("Error",
                                         f"Consensus sequence generation failed. Refer to terminal for more information: {str(e)}",
                                         parent=current_win)
                    return

                try:
                    # draw CIAlign style plot
                    # Define output plot file
                    cialign_plot_file = os.path.join(source_dir, f"{os.path.basename(input_fasta_file)}_CIAlign.png")
                    drawMiniAlignment(input_fasta_file, cialign_plot_file, dpi=300)

                except Exception as e:
                    click.echo(f"An error for CIAlign plotting: \n {traceback.format_exc()}")
                    messagebox.showerror("Error",
                                         f"CIAlign plotting failed. please align sequences first. "
                                         f"Refer to terminal for more information. : {str(e)}",
                                         parent=current_win)

                # Fresh child canvas
                fresh_canvas(child_frame, child_canvas, source_dir, current_win,
                             update_child_canvas=update_child_canvas, file_start=file_start, file_end=file_end,
                             fresh_save_path=save_path)

        return _generate_cons

    #####################################################################################################
    # Code block: define consensus sequence generation function
    #####################################################################################################
    # Define search function
    def search_files_with_pattern(filename_pattern, dirs, current_canvas_status="tetrimmer_out"):

        # dirs: path contains all files need to search, like te_trimmer_proof_curation_dir
        # and other_cons_lib_single_file_folder
        # current_canvas_status can be "tetrimmer_out" and "cons_lib"
        # "tetrimmer" means to check TEtrimmer output files. "cons_lib" is to search TE consensus library files
        # te_trimmer_proof_curation_dir
        found_paths = []
        print(current_canvas_status)

        if current_canvas_status == "tetrimmer_out":
            # direct dir to
            check_paths = [os.path.join(dirs, "Clustered_proof_curation"),
                           os.path.join(dirs, "TE_low_copy"),
                           os.path.join(dirs, "TE_skipped")]

            # Define found_path to store search result
            for directory in check_paths:
                if not os.path.isdir(directory):
                    continue

                for root, subdirs, files in os.walk(directory):
                    # Check files in the current directory
                    for file in files:
                        if filename_pattern in file:
                            # Only show the last folder name of root and the file name
                            found_paths.append(os.path.join(os.path.basename(root), file))

                    # Check files in each subdirectory
                    # subdir will be a empty directory if no directory is found under root
                    for subdir in subdirs:
                        subdir_path = os.path.join(root, subdir)
                        for subdir_root, _, subdir_files in os.walk(subdir_path):
                            for subdir_file in subdir_files:
                                if filename_pattern in subdir_file:
                                    # Only show the cluster number and the file name
                                    found_paths.append(os.path.join(os.path.basename(subdir_root), subdir_file))
        else:
            # Then current_canvas_status is "cons_lib"
            # The dirs variable should be "other_cons_lib_single_file_folder" to search TE consensus library
            if os.path.isdir(dirs):
                sorted_cons_files = [f for f in sorted(os.listdir(dirs))]
                for index, file in enumerate(sorted_cons_files):
                    if filename_pattern in file:
                        found_paths.append((f"Position number {index + 1}", file))

        return found_paths

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

    log_text = """
        ████████\\ ████████\\ ██\\               ██\\
        \\__██  __|██  _____|██ |              \\__|
           ██ |   ██ |    ██████\\    ██████\\  ██\\ ██████\\████\\  ██████\\████\\   ██████\\   ██████\\
           ██ |   █████\\  \\_██  _|  ██  __██\\ ██ |██  _██  _██\\ ██  _██  _██\\ ██  __██\\ ██  __██\\
           ██ |   ██  __|   ██ |    ██ |  \\__|██ |██ / ██ / ██ |██ / ██ / ██ |████████ |██ |  \\__|
           ██ |   ██ |      ██ |██\\ ██ |      ██ |██ | ██ | ██ |██ | ██ | ██ |██   ____|██ |
           ██ |   ████████\\ \\████  |██ |      ██ |██ | ██ | ██ |██ | ██ | ██ |\\███████\\ ██ |
           \\__|   \\________| \\____/ \\__|      \\__|\\__| \\__| \\__|\\__| \\__| \\__| \\_______|\\__|
        """


    initial_text = "Manual proof curation is highly recommended to improve the quality of TE annotations.\n\n" \
                   "1, Click the <Clustered_proof_curation> button in the menu bar.\n\n" \
                   "   TEs with more than 90% identity are grouped into one cluster.\n\n" \
                   "2, Click each <Cluster> button.\n\n" \
                   "   For each TE, you can find four files: \n" \
                   "     <seq_name.cluster.fa>          multiple sequence alignment (MSA) before clustering\n" \
                   "     <seq_name.raw.fa>               MSA file before cleaning\n" \
                   "     <seq_name.fa>                      MSA file corresponding to TE consensus sequence\n" \
                   "     <seq_name.pdf>                    report plots file for evaluating annotation quality\n\n" \
                   "3, Double click <seq_name.pdf> to evaluate annotation quality.\n\n" \
                   "4, If satisfied, click <Cons> button next to <seq_name.fa>.\n" \
                   "   This MSA file will be copied to <Proof_curation_consensus_folder>.\n\n" \
                   "5, If not satisfied, check <seq_name.fa> or <seq_name.raw.fa>. \n\n" \
                   "6, For more extension, click <Extend> button next to the fasta file.\n\n" \
                   "   Use <CropDiv> <CropGap> and <CleanCol> to clean the MSA. \n" \
                   "   Use <TEAid> to generate interactive report plots. \n\n" \
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

    def show_progress_bar(parent, height=5, speed=1000):
        # Create a custom style
        style = ttk.Style()
        style.configure("custom.Horizontal.TProgressbar", thickness=height)

        progress = ttk.Progressbar(parent, mode='indeterminate', style="custom.Horizontal.TProgressbar")
        progress.pack(side='bottom', fill='x')
        parent.update_idletasks()  # Ensure the layout is updated

        def update_progress():
            if progress.winfo_exists():  # Check if the widget exists before updating
                progress.step(1)
                parent.after(speed, update_progress)  # Adjust the speed as needed

        progress.start()
        update_progress()  # Start updating the progress bar

        return progress

    # Set function to allow to modify MSA cleaning parameters
    def show_settings_dialog():
        settings_window = Toplevel(root)
        settings_window.title("Modify Parameters")
        if os_type == "Darwin":
            settings_window.geometry('550x400')
        else:
            settings_window.geometry('600x450')

        def save_settings():
            global crop_div_thr_g, crop_div_win_g, crop_gap_thr_g, crop_gap_win_g
            global column_gap_thr_g, simi_check_gap_thr_g, similarity_thr_g, min_nucleotide_g
            global cons_thre_g, blast_e_value_g

            try:
                crop_div_thr_g = float(crop_div_thr_entry.get())
                crop_div_win_g = int(crop_div_win_entry.get())
                crop_gap_thr_g = float(crop_gap_thr_entry.get())
                crop_gap_win_g = int(crop_gap_win_entry.get())
                column_gap_thr_g = float(column_gap_thr_entry.get())
                simi_check_gap_thr_g = float(simi_check_gap_thr_entry.get())
                similarity_thr_g = float(similarity_thr_entry.get())
                min_nucleotide_g = int(min_nucleotide_entry.get())
                cons_thre_g = float(cons_thres_entry.get())
                blast_e_value_g = float(e_value_entry.get())

                settings_window.destroy()

            except ValueError as e:
                messagebox.showerror("Invalid input", f"Please enter valid values. Error: {str(e)}")

        # Create a frame to hold the labels and entries
        frame = Frame(settings_window, bg='white')
        frame.pack(fill='both', expand=True, padx=10, pady=10)

        # Add labels and entries with grid for alignment
        Label(frame, text="Cons: Consensus generation Threshold (0-1):", anchor='w', bg='white').grid(row=0, column=0,
                                                                                                      sticky='w',
                                                                                                      pady=2)
        cons_thres_entry = Entry(frame, bg='white')
        cons_thres_entry.insert(0, str(cons_thre_g))
        cons_thres_entry.grid(row=0, column=1, sticky='w', pady=2)

        Label(frame, text="Blast: BLASTN e-value:", anchor='w', bg='white').grid(row=1, column=0, sticky='w', pady=2)
        e_value_entry = Entry(frame, bg='white')
        e_value_entry.insert(0, str(blast_e_value_g))
        e_value_entry.grid(row=1, column=1, sticky='w', pady=2)

        Label(frame, text="CropDiv: Crop End Divergence Threshold (0-1):", anchor='w', bg='white').grid(row=2, column=0,
                                                                                                        sticky='w',
                                                                                                        pady=2)
        crop_div_thr_entry = Entry(frame, bg='white')
        crop_div_thr_entry.insert(0, str(crop_div_thr_g))
        crop_div_thr_entry.grid(row=2, column=1, sticky='w', pady=2)

        Label(frame, text="CropDiv: Crop End Divergence Window Size:", anchor='w', bg='white').grid(row=3, column=0,
                                                                                                    sticky='w', pady=2)
        crop_div_win_entry = Entry(frame, bg='white')
        crop_div_win_entry.insert(0, str(crop_div_win_g))
        crop_div_win_entry.grid(row=3, column=1, sticky='w', pady=2)

        Label(frame, text="CropGap: Crop End Gap Threshold (0-1):", anchor='w', bg='white').grid(row=4, column=0,
                                                                                                 sticky='w', pady=2)
        crop_gap_thr_entry = Entry(frame, bg='white')
        crop_gap_thr_entry.insert(0, str(crop_gap_thr_g))
        crop_gap_thr_entry.grid(row=4, column=1, sticky='w', pady=2)

        Label(frame, text="CropGap: Crop End Gap Window Size:", anchor='w', bg='white').grid(row=5, column=0,
                                                                                             sticky='w', pady=2)
        crop_gap_win_entry = Entry(frame, bg='white')
        crop_gap_win_entry.insert(0, str(crop_gap_win_g))
        crop_gap_win_entry.grid(row=5, column=1, sticky='w', pady=2)

        Label(frame, text="CleanCol: Clean Column Threshold (0-1):", anchor='w', bg='white').grid(row=6, column=0,
                                                                                                  sticky='w', pady=2)
        column_gap_thr_entry = Entry(frame, bg='white')
        column_gap_thr_entry.insert(0, str(column_gap_thr_g))
        column_gap_thr_entry.grid(row=6, column=1, sticky='w', pady=2)

        Label(frame, text="CleanCol: Similarity Check Gap Threshold (0-1):", anchor='w', bg='white').grid(row=7,
                                                                                                          column=0,
                                                                                                          sticky='w',
                                                                                                          pady=2)
        simi_check_gap_thr_entry = Entry(frame, bg='white')
        simi_check_gap_thr_entry.insert(0, str(simi_check_gap_thr_g))
        simi_check_gap_thr_entry.grid(row=7, column=1, sticky='w', pady=2)

        Label(frame, text="CleanCol: Similarity Threshold (0-1):", anchor='w', bg='white').grid(row=8, column=0,
                                                                                                sticky='w', pady=2)
        similarity_thr_entry = Entry(frame, bg='white')
        similarity_thr_entry.insert(0, str(similarity_thr_g))
        similarity_thr_entry.grid(row=8, column=1, sticky='w', pady=2)

        Label(frame, text="CleanCol: Min Nucleotide number:", anchor='w', bg='white').grid(row=9, column=0, sticky='w',
                                                                                           pady=2)
        min_nucleotide_entry = Entry(frame, bg='white')
        min_nucleotide_entry.insert(0, str(min_nucleotide_g))
        min_nucleotide_entry.grid(row=9, column=1, sticky='w', pady=2)

        Button(frame, text="Save", command=save_settings).grid(row=10, columnspan=2, pady=10)

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

            button_frame.grid_columnconfigure(0, weight=1)
            button_frame.grid_rowconfigure(0, weight=1)
            frame.grid_columnconfigure(1, weight=1)
            frame.grid_columnconfigure(2, weight=0)

    def load_files_with_destroy(start, end, frame, canvas, path):

        # update current_canvas_content to help search function
        global current_canvas_content
        current_canvas_content = "tetrimmer_out"

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

    # source_dir is the folder path contains all files need to be loaded
    # canvas is inside current_win
    # frame is inside canvas
    # frame contains child-frames and buttons

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

            """
            # Get button states
        button_defaults = {
            'file_button': ('white', 'black'),
            'save_button': ('white', 'black'),
            'extension_button': ('white', 'black'),
            'teaid_button': ('white', 'black'),
            'cons_button': ('white', 'black'),
            'others_button': ('white', 'black')
        }

        if filename in button_states:
            states = button_states[filename]
            for j, key in enumerate(button_defaults.keys()):
                if j < len(states):
                    button_defaults[key] = states[j]
                    
        Button(button_frame, text="Save", bg=button_defaults['save_button'][0], fg=button_defaults['save_button'][1])
            """

            # Get button states
            file_button_bg = 'white'
            save_button_bg = 'white'
            extension_button_bg = 'white'
            teaid_button_bg = 'white'
            crop_end_by_div_button_bg = 'white'
            crop_end_by_gap_button_bg = 'white'
            remove_gap_column_button_bg = "white"
            cons_button_bg = 'white'
            others_button_bg = 'white'

            file_button_fg = 'black'
            save_button_fg = 'black'
            extension_button_fg = 'black'
            teaid_button_fg = 'black'
            crop_end_by_div_button_fg = 'black'
            crop_end_by_gap_button_fg = 'black'
            remove_gap_column_button_fg = "black"
            cons_button_fg = "black"
            others_button_fg = 'black'

            if filename in button_states:
                states = button_states[filename]
                if len(states) >= 1:
                    file_button_bg = states[0][0]
                    file_button_fg = states[0][1]
                if len(states) >= 2:
                    save_button_bg = states[1][0]
                    save_button_fg = states[1][1]
                if len(states) >= 3:
                    extension_button_bg = states[2][0]
                    extension_button_fg = states[2][1]
                if len(states) >= 4:
                    teaid_button_bg = states[3][0]
                    teaid_button_fg = states[3][1]
                if len(states) >= 5:
                    crop_end_by_div_button_bg = states[4][0]
                    crop_end_by_div_button_fg = states[4][1]
                if len(states) >= 6:
                    crop_end_by_gap_button_bg = states[5][0]
                    crop_end_by_gap_button_fg = states[5][1]
                if len(states) >= 7:
                    remove_gap_column_button_bg = states[6][0]
                    remove_gap_column_button_fg = states[6][1]
                if len(states) >= 8:
                    cons_button_bg = states[7][0]
                    cons_button_fg = states[7][1]
                if len(states) >= 9:
                    others_button_bg = states[8][0]
                    others_button_fg = states[8][1]

            # Add file name button into canvas frame
            file_button = Button(frame, text=filename, anchor='w', bg=file_button_bg, fg=file_button_fg)
            file_button.grid(row=i - start, column=1, sticky='ew')
            # Bind with open_file function to open file
            file_button.bind('<Double-Button-1>', open_file(filename, file_button, source_dir))

            # Build a child button_frame inside frame
            button_frame = Frame(frame, bg='white')
            button_frame.grid(row=i - start, column=2, sticky='e')

            # Create "Consensus" button inside button_frame
            copy_button = Button(button_frame, text="Save", bg=save_button_bg, fg=save_button_fg)
            copy_button.grid(row=0, column=0, padx=1)
            # Bind "Consensus" button with copy_file function with specific source and destination folder
            copy_button.bind('<Button-1>', copy_file(filename, copy_button, consensus_folder, source_dir, current_win))

            # Define "Extension" button
            more_extend_button = Button(button_frame, text="Extend", bg=extension_button_bg, fg=extension_button_fg)
            more_extend_button.grid(row=0, column=1, padx=1)
            # Bind "Extension" button with extension function
            more_extend_button.bind('<Button-1>', extension_function(filename, more_extend_button, source_dir,
                                                                     temp_folder, current_win, chrom_size, frame,
                                                                     canvas))

            # Define "Plotter" button
            plot_button = Button(button_frame, text="TEAid", bg=teaid_button_bg, fg=teaid_button_fg)
            plot_button.grid(row=0, column=2, padx=1)
            # Bind "Plotter" button with plotter_function
            plot_button.bind('<Button-1>',
                             plotter_function(filename, plot_button, source_dir, temp_folder, genome_file,
                                              canvas, current_win))

            # Define "Crop end by divergence" button
            crop_end_by_div_button = Button(button_frame, text="CropDiv", bg=crop_end_by_div_button_bg, fg=crop_end_by_div_button_fg)
            crop_end_by_div_button.grid(row=0, column=3, padx=1)
            crop_end_by_div_button.bind('<Button-1>',
                                        crop_end_div_gui(filename, crop_end_by_div_button, source_dir, source_dir,
                                                         current_win, frame, canvas))

            # Define "Crop end by gap" button
            crop_end_by_gap_button = Button(button_frame, text="CropGap", bg=crop_end_by_gap_button_bg, fg=crop_end_by_gap_button_fg)
            crop_end_by_gap_button.grid(row=0, column=4, padx=1)

            crop_end_by_gap_button.bind('<Button-1>',
                                        crop_end_gap_gui(filename, crop_end_by_gap_button, source_dir, source_dir,
                                                         current_win, frame, canvas))

            # Define "Clean gap column" button
            clean_gap_column_button = Button(button_frame, text="CleanCol", bg=remove_gap_column_button_bg, fg=remove_gap_column_button_fg)
            clean_gap_column_button.grid(row=0, column=5, padx=1)

            clean_gap_column_button.bind('<Button-1>', remove_gaps_with_similarity_check_gui(
                filename, clean_gap_column_button, source_dir, source_dir, current_win,
                frame, canvas))

            # Define "Cons" button
            cons_button = Button(button_frame, text="Cons", bg=cons_button_bg, fg=cons_button_fg)
            cons_button.grid(row=0, column=6, padx=1)

            cons_button.bind('<Button-1>', generate_cons(filename, cons_button, source_dir, current_win, frame,
                                                         canvas, cons_thre=cons_thre_g, ambiguous="N"))

            # Define "Others" button (discard)
            others_button = Button(button_frame, text="Discard", bg=others_button_bg, fg=others_button_fg)
            others_button.grid(row=0, column=7, padx=1)
            # Bind "Others" button with copy_file function with different destination folder
            others_button.bind('<Button-1>', copy_file(filename, others_button, others_dir, source_dir, current_win))

            button_frame.grid_columnconfigure(0, weight=1)
            button_frame.grid_rowconfigure(0, weight=1)
            frame.grid_columnconfigure(1, weight=1)
            frame.grid_columnconfigure(2, weight=0)

    # Build child canvas window
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
    # Code block: Define Low_copy_TE, Skipped_TE and other TE consensus library check canvas
    #####################################################################################################
    # load_files(start, end, frame, canvas, source_dir=te_trimmer_proof_curation_dir)
    # load_files(start, end, frame, canvas, source_dir=path)
    # source_dir is the folder path contains all files need to be loaded
    # canvas is inside current_win
    # frame is inside canvas
    # frame contains child-frames and buttons
    def other_cons_load_files(start, end, frame, canvas, source_dir, current_win, save_path=None,
                              scroll_position=None, button_states=None):

        # Update global current_canvas_content variable to enable search function
        global current_canvas_content
        current_canvas_content = "cons_lib"

        destroy_initial_label()
        clear_frame()
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
            save_button_bg = 'white'
            blast_button_bg = 'white'
            extension_button_bg = 'white'
            teaid_button_bg = 'white'
            crop_end_by_div_button_bg = 'white'
            crop_end_by_gap_button_bg = 'white'
            remove_gap_column_button_bg = "white"
            cons_button_bg = "white"

            file_button_fg = 'black'
            save_button_fg = 'black'
            blast_button_fg = 'black'
            extension_button_fg = 'black'
            teaid_button_fg = 'black'
            crop_end_by_div_button_fg = 'black'
            crop_end_by_gap_button_fg = 'black'
            remove_gap_column_button_fg = "black"
            cons_button_fg = "black"

            if filename in button_states:
                states = button_states[filename]
                if len(states) >= 1:
                    file_button_bg = states[0][0]
                    file_button_fg = states[0][1]
                if len(states) >= 2:
                    save_button_bg = states[1][0]
                    save_button_fg = states[1][1]
                if len(states) >= 3:
                    blast_button_bg = states[2][0]
                    blast_button_fg = states[2][1]
                if len(states) >= 4:
                    extension_button_bg = states[3][0]
                    extension_button_fg = states[3][1]
                if len(states) >= 5:
                    teaid_button_bg = states[4][0]
                    teaid_button_fg = states[4][1]
                if len(states) >= 6:
                    crop_end_by_div_button_bg = states[5][0]
                    crop_end_by_div_button_fg = states[5][1]
                if len(states) >= 7:
                    crop_end_by_gap_button_bg = states[6][0]
                    crop_end_by_gap_button_fg = states[6][1]
                if len(states) >= 8:
                    remove_gap_column_button_bg = states[7][0]
                    remove_gap_column_button_fg = states[7][1]
                if len(states) >= 9:
                    cons_button_bg = states[8][0]
                    cons_button_fg = states[8][1]

            # Add file name button into canvas frame
            file_button = Button(frame, text=filename, anchor='w', bg=file_button_bg, fg=file_button_fg)
            file_button.grid(row=i - start, column=1, sticky='ew')
            # Bind with open_file function to open file
            file_button.bind('<Double-Button-1>', open_file(filename, file_button, source_dir))

            # Build a child button_frame inside frame
            button_frame = Frame(frame, bg='white')
            button_frame.grid(row=i - start, column=2, sticky='e')

            # Create "Save" button inside button_frame
            copy_button = Button(button_frame, text="Save", bg=save_button_bg, fg=save_button_fg)
            copy_button.grid(row=0, column=0, padx=1)
            # Bind "Consensus" button with copy_file function with specific source and destination folder
            copy_button.bind('<Button-1>', copy_file(filename, copy_button, save_path,
                                                     source_dir, current_win))

            # Define "Blast" button
            blast_button = Button(button_frame, text="Blast", bg=blast_button_bg, fg=blast_button_fg)
            blast_button.grid(row=0, column=1, padx=1)
            # Bind blast button with blast functions
            blast_button.bind('<Button-1>', check_other_cons_lib_blast(filename, blast_button, source_dir,
                                                                       other_cons_lib_folder, current_win, frame,
                                                                       canvas, genome_file, e_value=blast_e_value_g,
                                                                       update_child_canvas=False,
                                                                       file_start=start, file_end=end,
                                                                       save_path=save_path))

            # Define "Extension" button
            more_extend_button = Button(button_frame, text="Extend", bg=extension_button_bg, fg=extension_button_fg)
            more_extend_button.grid(row=0, column=2, padx=1)
            # Bind "Extension" button with extension function
            more_extend_button.bind('<Button-1>', extension_function(filename, more_extend_button, source_dir,
                                                                     temp_folder, current_win, chrom_size, frame,
                                                                     canvas, update_child_canvas=False,
                                                                     file_start=start, file_end=end,
                                                                     save_path=save_path))

            # Define "Plotter" button
            plot_button = Button(button_frame, text="TEAid", bg=teaid_button_bg, fg=teaid_button_fg)
            plot_button.grid(row=0, column=3, padx=1)
            # Bind "Plotter" button with plotter_function
            plot_button.bind('<Button-1>',
                             plotter_function(filename, plot_button, source_dir, temp_folder, genome_file,
                                              canvas, current_win))

            # Define "Crop end by divergence" button
            crop_end_by_div_button = Button(button_frame, text="CropDiv", bg=crop_end_by_div_button_bg,
                                            fg=crop_end_by_div_button_fg)
            crop_end_by_div_button.grid(row=0, column=4, padx=1)
            crop_end_by_div_button.bind('<Button-1>',
                                        crop_end_div_gui(filename, crop_end_by_div_button, source_dir, source_dir,
                                                         current_win, frame, canvas, update_child_canvas=False,
                                                         file_start=start, file_end=end, save_path=save_path))

            # Define "Crop end by gap" button
            crop_end_by_gap_button = Button(button_frame, text="CropGap", bg=crop_end_by_gap_button_bg,
                                            fg=crop_end_by_gap_button_fg)
            crop_end_by_gap_button.grid(row=0, column=5, padx=1)
            crop_end_by_gap_button.bind('<Button-1>',
                                        crop_end_gap_gui(filename, crop_end_by_gap_button, source_dir, source_dir,
                                                         current_win, frame, canvas, update_child_canvas=False,
                                                         file_start=start, file_end=end, save_path=save_path))

            # Define "Clean gap column" button
            clean_gap_column_button = Button(button_frame, text="CleanCol", bg=remove_gap_column_button_bg,
                                             fg=remove_gap_column_button_fg)
            clean_gap_column_button.grid(row=0, column=6, padx=1)
            clean_gap_column_button.bind('<Button-1>', remove_gaps_with_similarity_check_gui(
                filename, clean_gap_column_button, source_dir, source_dir, current_win,
                frame, canvas, update_child_canvas=False, file_start=start, file_end=end, save_path=save_path))

            # Define "Cons" button
            cons_button = Button(button_frame, text="Cons", bg=cons_button_bg, fg=cons_button_fg)
            cons_button.grid(row=0, column=7, padx=1)

            cons_button.bind('<Button-1>', generate_cons(filename, cons_button, source_dir, current_win, frame,
                                                         canvas, cons_thre=cons_thre_g, ambiguous="N",
                                                         update_child_canvas=False, file_start=start, file_end=end,
                                                         save_path=save_path))

            button_frame.grid_columnconfigure(0, weight=1)
            button_frame.grid_rowconfigure(0, weight=1)
            frame.grid_columnconfigure(1, weight=1)
            frame.grid_columnconfigure(2, weight=0)

    #####################################################################################################
    # Code block: Define search window
    #####################################################################################################
    def open_search_window():
        search_window = Toplevel(root)
        search_window.title("Search Files")
        search_window.geometry('800x600')

        search_frame = Frame(search_window)
        search_frame.pack(fill='x', padx=10, pady=5)

        search_entry = Entry(search_frame, width=40)
        search_entry.pack(side='left', padx=5)

        search_button = Button(search_frame, text="Search",
                               command=lambda: perform_search_in_window(search_entry.get(), result_frame, result_canvas))
        search_button.pack(side='left', padx=5)

        # Create canvas for displaying search results
        result_canvas = Canvas(search_window, bg='white')
        result_canvas.pack(side='left', fill='both', expand=True)

        result_scrollbar = Scrollbar(search_window, orient='vertical', command=result_canvas.yview)
        result_canvas.configure(yscrollcommand=result_scrollbar.set)
        result_scrollbar.pack(side='right', fill='y')

        search_window.bind("<MouseWheel>", lambda event: scroll(event, result_canvas))
        search_window.bind("<Button-4>", lambda event: scroll(event, result_canvas))
        search_window.bind("<Button-5>", lambda event: scroll(event, result_canvas))

        result_canvas.yview_moveto(0)  # Reset scrollbar to top

        result_frame = Frame(result_canvas, bg='white')
        result_canvas.create_window((0, 0), window=result_frame, anchor='nw')

        result_frame.bind('<Configure>', lambda e: result_canvas.configure(scrollregion=result_canvas.bbox('all')))

        search_window.mainloop()

    def perform_search_in_window(pattern, frame, canvas):
        # te_trimmer_proof_curation_dir, current_canvas_content, other_cons_lib_single_file_folder are global variables
        if pattern:
            if current_canvas_content == "tetrimmer_out":
                search_results = search_files_with_pattern(pattern, te_trimmer_proof_curation_dir, current_canvas_content)
            else:
                search_results = search_files_with_pattern(pattern, other_cons_lib_single_file_folder,
                                                           current_canvas_content)
            display_search_results_in_window(search_results, frame)
            canvas.yview_moveto(0)  # Reset scrollbar to top

    def display_search_results_in_window(search_results, frame):
        for widget in frame.winfo_children():
            widget.destroy()

        if not search_results:
            label = Label(frame, text="No files found matching the pattern.", bg='white')
            label.pack(pady=20)
            return

        for i, result in enumerate(search_results):
            line_number = Label(frame, text=str(i + 1), bg='white')
            line_number.grid(row=i, column=0)

            file_label = Label(frame, text=result, anchor='w', bg='white')
            file_label.grid(row=i, column=1, sticky='ew')

    #####################################################################################################
    # Code block: Add menu bar for mother window
    #####################################################################################################

    # Create mother menu
    menubar = Menu(root)

    # Show menu on the window
    root.configure(menu=menubar)

    annotation_folders = ["Clustered_proof_curation", "TE_low_copy", "TE_skipped", "Check_consensus_lib"]

    # Create sub-menu
    for i, annotation in enumerate(annotation_folders):
        annotationMenu = Menu(menubar, tearoff=0)
        menubar.add_cascade(label=annotation, menu=annotationMenu)

        # Check which menu button is selected
        if i <= 2:  # Means "Clustered_proof_curation", "TE_low_copy", "TE_skipped"
            annotation_path = os.path.join(te_trimmer_proof_curation_dir, annotation)
        else:  # Means "Check_consensus_lib"
            annotation_path = other_cons_lib_single_file_folder

        # Give hits when folder isn't found
        if not os.path.exists(annotation_path):

            if i == 0:
                annotationMenu.add_command(label="Cluster folder not detected! Use correct input folder path",
                                           command=partial(messagebox.showerror,
                                                           "Error", "Please use the correct input directory by '-i' option."
                                                                    " three folder should be contained in your input path, "
                                                                    "including 'Clustered_proof_curation', "
                                                                    "'TE_skipped', and 'TE_low_copy'"))
            elif i == 1:
                annotationMenu.add_command(label="No low copy TEs")
            elif i == 2:
                annotationMenu.add_command(label="No skipped TEs")

            elif annotation == "Check_consensus_lib":
                annotationMenu.add_command(label="TE consensus sequences not found",
                                           command=partial(messagebox.showerror,
                                                           "Error", "Please use option '--consensus_lib' to define the "
                                                                    "TE consensus library path you want to check"))
            continue

        # Sort files inside each folder
        sorted_files_annotation = sorted(os.listdir(annotation_path))

        # Give No files found when the folder is empty
        if not sorted_files_annotation:
            annotationMenu.add_command(label="No file is found.")
        else:
            # Show maximum 150 files on the window
            if i == 0:  # Means "Clustered_proof_curation"
                for j in range(0, len(sorted_files_annotation), 150):
                    end = min(j + 150, len(sorted_files_annotation))

                    annotationMenu.add_command(label=f"{j + 1}-{end}",
                                               command=partial(load_files_with_destroy,
                                                               j, end, frame, canvas, annotation_path))
            elif i == 1 or i == 2:  # Means "TE_low_copy", "TE_skipped"
                for j in range(0, len(sorted_files_annotation), 100):
                    end = min(j + 100, len(sorted_files_annotation))

                    # Add 1000 to the end when checking the last page to enable to show added files by fresh canvas
                    if end == len(sorted_files_annotation):
                        annotationMenu.add_command(label=f"{j + 1}-{end}",
                                                   command=partial(other_cons_load_files,
                                                                   j, end+1000, frame, canvas, annotation_path, canvas,
                                                                   save_path=consensus_folder))
                    else:
                        annotationMenu.add_command(label=f"{j + 1}-{end}",
                                                   command=partial(other_cons_load_files,
                                                                   j, end, frame, canvas, annotation_path,
                                                                   canvas,
                                                                   save_path=consensus_folder))
            elif i == 3:  # Means "Check_consensus_lib"

                for j in range(0, len(sorted_files_annotation), 100):
                    end = min(j + 100, len(sorted_files_annotation))

                    if end == len(sorted_files_annotation):

                        annotationMenu.add_command(label=f"{j + 1}-{end}",
                                                   command=partial(other_cons_load_files,
                                                                   j, end+1000, frame, canvas, annotation_path, canvas,
                                                                   save_path=other_cons_lib_result_folder))
                    else:
                        annotationMenu.add_command(label=f"{j + 1}-{end}",
                                                   command=partial(other_cons_load_files,
                                                                   j, end, frame, canvas, annotation_path,
                                                                   canvas,
                                                                   save_path=other_cons_lib_result_folder))

    # Add confirm menu button
    settings_menu = Menu(menubar, tearoff=0)
    # Set confirm_menu child menu
    # Initialize the BooleanVar to store "Show Confirmation Window" status
    show_confirmation = BooleanVar(value=True)
    settings_menu.add_checkbutton(label="Show Confirmation Window When Click 'Save'", onvalue=True, offvalue=False,
                                  variable=show_confirmation)
    settings_menu.add_command(label="Modify Function Parameters", command=show_settings_dialog)
    menubar.add_cascade(label="Setting", menu=settings_menu)

    # Add Undo button
    undo_menu = Menu(menubar, tearoff=0)
    # Add Undo child menu
    undo_menu.add_command(label="undo last copy", command=undo_last_copy)
    menubar.add_cascade(label="Undo", menu=undo_menu)

    # Add the Search button to the menu
    search_menu = Menu(menubar, tearoff=0)
    search_menu.add_command(label="Search Files", command=lambda: open_search_window())
    menubar.add_cascade(label="Search", menu=search_menu)

    # Add the Help button to the menu
    help_menu = Menu(menubar, tearoff=0)
    help_menu.add_command(label="Show instruction", command=show_help)
    menubar.add_cascade(label="Help", menu=help_menu)

    # Add double confirmation on window close
    def on_closing():
        if messagebox.askokcancel("Quit", "Do you really want to quit?"):

            # Remove all files in the temp_folder
            try:
                click.echo("TEtrimmer is cleaning temporary files")
                shutil.rmtree(temp_folder)
                click.echo("All temporary files are cleaned")
            except Exception as e:
                pass
            root.destroy()

    root.protocol("WM_DELETE_WINDOW", on_closing)

    root.mainloop()


if __name__ == '__main__':
    proof_curation()
