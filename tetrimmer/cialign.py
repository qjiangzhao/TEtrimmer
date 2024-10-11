import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from PIL import Image
import os
import tempfile
from Bio import SeqIO

matplotlib.use('Agg')


#########################################################################
# This code is derived from CIAlign: https://github.com/KatyBrown/CIAlign
#########################################################################


def png_to_pdf(png_path, pdf_path):
    """
    Convert a PNG file to a PDF.

    Parameters:
    - png_path (str): Path for the input PNG file.
    - pdf_path (str): Path for the output PDF file.
    """
    # Open the PNG image
    image = Image.open(png_path)

    # Convert RGBA to RGB
    if image.mode == 'RGBA':
        rgb_image = Image.new('RGB', image.size, (255, 255, 255))  # White background
        rgb_image.paste(image, mask=image.split()[3])  # Paste using alpha channel as mask
        image = rgb_image

        # Release RAM memory
        del rgb_image

    # Convert and save as PDF
    image.save(pdf_path, "PDF", resolution=100.0)

    # Close the image and release memory
    image.close()
    del image


def FastaToArray(infile):
    '''
    Convert an alignment into a numpy array.

    Parameters
    ----------
    infile: string
        path to input alignment file in FASTA format

    Returns
    -------
    arr: np.array
        2D numpy array in the same order as fasta_dict where each row
        represents a single column in the alignment and each column a
        single sequence.
    nams: list
        List of sequence names in the same order as in the input file
    '''

    nams = []
    seqs = []

    valid_chars = {'A', 'G', 'C', 'T', 'N', '-'}

    for record in SeqIO.parse(infile, "fasta"):
        nams.append(record.id)
        seq = [base if base in valid_chars else '-' for base in str(record.seq).upper()]
        seqs.append(seq)

    # Return False if there is only one sequence
    seq_len = len(seqs)
    if seq_len <= 1:
        return None, None, None

    # Check if all sequences are of the same length
    seq_lengths = {len(seq) for seq in seqs}
    if len(seq_lengths) > 1:
        raise ValueError(
            "ERROR: The sequences you provided may not be aligned - not all the sequences are the same length")

    arr = np.array(seqs)
    return arr, nams, seq_len


def arrNumeric(arr):
    '''
    Converts the sequence array into a numerical matrix and a colour map
    which matplotlib can interpret as an image (similar to
                                                https://bit.ly/2CIKOEr)
    The rows in the array are inverted so that the output image has the rows
    in the same order as the input alignment.

    Parameters
    ----------
    arr: np.array
        The alignment stored as a numpy array

    typ: str
        Either 'aa' - amino acid - or 'nt' - nucleotide

    palette: str
        Colour palette, CBS or Bright

    Returns
    -------
    arr2: np.array
        The flipped alignment as an array of integers
    cmap: matplotlib.colors.ListedColormap
        A colour map with the colours corresponding to each base
        or amino acid
    '''

    # turn the array upside down
    arr = np.flip(arr, axis=0)

    color_pattern_tetrimmer = {'A': "#00CC00",
                               'G': "#949494",
                               'T': "#FF6666",
                               'C': "#6161ff",
                               'N': "#c7d1d0",
                               'n': "#c7d1d0",
                               '-': "#FFFFFF",
                               'X': "#c7d1d0"}

    color_pattern_bright = {'A': "#f20707",
                            'G': "#ffd500",
                            'T': "#64bc3c",
                            'C': "#0907f2",
                            'N': "#c7d1d0",
                            'n': "#c7d1d0",
                            '-': "#FFFFFF",
                            'X': "#c7d1d0"}

    # Color-blind people friendly color pattern
    color_pattern = {'A': "#56ae6c",
                     'G': "#c9c433",
                     'T': "#a22c49",
                     'C': "#0038a2",
                     'N': "#6979d3",
                     'n': "#6979d3",
                     '-': "#FFFFFF",
                     'X': "#6979d3"}

    # retrieve the colours for the colour map
    keys = list(color_pattern.keys())
    ali_height, ali_width = np.shape(arr)

    # make a dictionary where each integer corresponds to a base or nt
    i = 0
    nD = dict()
    colours = []
    for key in keys:
        if key in arr:
            nD[key] = i
            colours.append(color_pattern[key])
            i += 1

    arr2 = np.empty([ali_height, ali_width])

    for x in range(ali_width):
        for y in range(ali_height):
            # numeric version of the alignment array
            arr2[y, x] = nD[arr[y, x]]

    cmap = matplotlib.colors.ListedColormap(colours)
    return arr2, cmap


def drawMiniAlignment(input_file, outfile, start_point, end_point,
                      dpi=500, title=None, width=20, height=20, orig_nams=[],
                      keep_numbers=False, force_numbers=False):
    '''
    Draws a "mini alignment" image showing a small representation of the
    whole alignment so that gaps and poorly aligned regions are visible.

    Parameters:
    - arr: np.array
        The alignment stored as a numpy array
    - outfile: str
        Path to the output PDF file (No need to give .pdf extension)
    - ... [rest of your parameters] ...

    Returns:
    - None
    '''

    arr, nams, seq_len = FastaToArray(input_file)

    # arr will be false when only one sequence is found in the fasta file
    if arr is None:
        return False

    if seq_len <= 60:
        height = seq_len * 0.3

    ali_height, ali_width = np.shape(arr)

    fontsize = 14

    # Determine tick interval
    if force_numbers:
        tickint = 1
    elif ali_height <= 11:
        tickint = 1
    elif ali_height <= 500:
        tickint = 10
    else:
        tickint = 100
    lineweight_h = 10 / ali_height
    lineweight_v = 10 / ali_width

    f = plt.figure(figsize=(width, height), dpi=dpi)
    a = f.add_subplot(1, 1, 1)
    a.set_xlim(-0.5, ali_width)
    a.set_ylim(-0.5, ali_height - 0.5)

    arr2, cm = arrNumeric(arr)
    a.imshow(arr2, cmap=cm, aspect='auto', interpolation='nearest')

    # Plot start and end points arrows
    arrow_length = 20
    a.annotate('Start crop Point', xy=(start_point, ali_height - 0.5), xytext=(0, arrow_length),
               textcoords='offset points',
               arrowprops=dict(facecolor='red', edgecolor='red', width=2, headwidth=5, headlength=4),
               ha='center', color='red', va='bottom', fontsize=15)
    a.annotate('End crop Point', xy=(end_point, ali_height - 0.5), xytext=(0, arrow_length),
               textcoords='offset points',
               arrowprops=dict(facecolor='blue', edgecolor='blue', width=2, headwidth=5, headlength=4),
               ha='center', color='blue', va='bottom', fontsize=15)

    f.subplots_adjust(top=0.85, bottom=0.15, left=0.1, right=0.95)
    a.hlines(np.arange(-0.5, ali_height), -0.5, ali_width, lw=lineweight_h, color='white', zorder=100)
    a.vlines(np.arange(-0.5, ali_width), -0.5, ali_height, lw=lineweight_v, color='white', zorder=100)

    a.spines['right'].set_visible(False)
    a.spines['top'].set_visible(False)
    a.spines['left'].set_visible(False)

    if title:
        f.suptitle(title, fontsize=fontsize * 1.5, y=0.92)

    for t in a.get_xticklabels():
        t.set_fontsize(fontsize)

    a.set_yticks(np.arange(ali_height - 1, -1, -tickint))
    x = 1
    if tickint == 1:
        if keep_numbers:
            labs = []
            for nam in orig_nams:
                if nam in nams:
                    labs.append(x)
                x += 1
            a.set_yticklabels(labs, fontsize=fontsize * 0.75)
        else:
            a.set_yticklabels(np.arange(1, ali_height + 1, tickint), fontsize=fontsize * 0.75)
    else:
        a.set_yticklabels(np.arange(0, ali_height, tickint), fontsize=fontsize)

    # Save the plot to a temporary PNG file
    with tempfile.NamedTemporaryFile(suffix='.png', delete=True) as tmp_file:
        temp_png = tmp_file.name
        f.savefig(temp_png, dpi=dpi, bbox_inches='tight')

        # Convert the temporary PNG to the desired PDF file
        png_to_pdf(temp_png, outfile)

    # Explicitly close the plot to release resources
    plt.close()
    del arr, arr2, nams

    if os.path.exists(outfile):
        return outfile
    else:
        return False

