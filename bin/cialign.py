import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
from PIL import Image
import os
import tempfile

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


def getPalette(palette='CBS'):
    '''
    Generates a dictionary that assigns a name to each color using a colorblind-
    friendly palette, generated using
    https://medialab.github.io/iwanthue/
    
    Parameters
    ----------
    palette: str
        The ID of the palette to be used, currently only colorblind-friendly
        (CBS) has been implemented.

    Returns
    -------
    dict
        Dictionary where keys are names of colors and
        values are hexadecimal codes for colors
    '''
    if palette.lower() == 'cbs':
        p = CBSafe()
    if palette.lower() == 'bright':
        p = Bright()
    if palette.lower() == 'te_trimmer':
        p = te_trimmer()
    return p


def getNtColours(palette='CBS'):
    '''
    Generates a dictionary that assigns a color to each nucleotide (including grey
    for "N" and white for "-")
    
    Parameters
    ----------
    pal: str
        A string designating which palette to use, currently only colorblind-friendly
        (CBS) has been implemented.

    Returns
    -------
    dict
        Dictionary where keys are single letter nucleotide codes and
        values are hexadecimal codes for colors
    '''
    pal = getPalette(palette=palette)
    return {'A': pal['green_nt'],
            'G': pal['yellow_nt'],
            'T': pal['red_nt'],
            'C': pal['blue_nt'],
            'N': pal['grey_nt'],
            "-": pal['white'],
            "U": pal['red_nt'],
            "R": pal['grey_nt'],
            "Y": pal['grey_nt'],
            "S": pal['grey_nt'],
            "W": pal['grey_nt'],
            "K": pal['grey_nt'],
            "M": pal['grey_nt'],
            "B": pal['grey_nt'],
            "D": pal['grey_nt'],
            "H": pal['grey_nt'],
            "V": pal['grey_nt'],
            "X": pal['grey_nt']}


def FastaToArray(infile, log=None):
    '''
    Convert an alignment into a NumPy array.

    Parameters
    ----------
    infile: string
        path to input alignment file in FASTA format
    log: logging.Logger
        An open log file object

    Returns
    -------
    arr: np.array
        2D NumPy array in the same order as fasta_dict where each row
        represents a single column in the alignment and each column a
        single sequence.
    nams: list
        List of sequence names in the same order as in the input file
    '''

    formatErrorMessage = "The MSA file needs to be in FASTA format."
    nams = []
    seqs = []
    nam = ""
    seq = ""
    psl = 0
    nseq = 0
    with open(infile) as input:
        for line in input:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == ">":
                sl = len([s.upper() for s in seq])
                if sl != psl and nseq > 1:
                    print(nseq, sl, psl)
                    raise ValueError("""
ERROR: The sequences you provided may not be aligned - not all sequences \
are the same length""")

                psl = sl
                nseq += 1
                seqs.append([s.upper() for s in seq])
                nams.append(nam)
                seq = []
                nam = line.replace(">", "")
            else:
                if len(nams) == 0:
                    if log:
                        log.error(formatErrorMessage)
                    print(formatErrorMessage)
                    exit(1)
                seq += list(line)
    sl = len([s.upper() for s in seq])
    if sl != psl and nseq > 1:
        print(nseq, sl, psl)
        raise ValueError("""
ERROR: The sequences you provided may not be aligned - not all sequences \
are the same length""")
    seqs.append(np.array([s.upper() for s in seq]))
    nams.append(nam)
    arr = np.array(seqs[1:])
    return arr, nams


def arrNumeric(arr, palette='CBS'):
    '''
    Converts the sequence array into a numerical matrix and a color map,
    which matplotlib can interpret as an image (similar to
                                                https://bit.ly/2CIKOEr)
    The rows in the array are inverted so that the output image has the rows
    in the same order as the input alignment.

    Parameters
    ----------
    arr: np.array
        The alignment stored as a NumPy array

    typ: str
        Either 'aa' - amino acid - or 'nt' - nucleotide

    palette: str
        Color palette, CBS or Bright

    Returns
    -------
    arr2: np.array
        The flipped alignment as an array of integers
    cmap: matplotlib.colors.ListedColormap
        A color map with the colors corresponding to each base
        or amino acid
    '''

    # turn the array upside down
    arr = np.flip(arr, axis=0)
    D = getNtColours(palette)

    # retrieve the colors for the color map
    keys = list(D.keys())
    ali_height, ali_width = np.shape(arr)

    # make a dictionary where each integer corresponds to a base or nt
    i = 0
    nD = dict()
    colours = []
    for key in keys:
        if key in arr:
            nD[key] = i
            colours.append(D[key])
            i += 1

    arr2 = np.empty([ali_height, ali_width])

    for x in range(ali_width):
        for y in range(ali_height):
            # numeric version of the alignment array
            arr2[y, x] = nD[arr[y, x]]

    cmap = matplotlib.colors.ListedColormap(colours)
    return arr2, cmap


def drawMiniAlignment(arr, nams, outfile, start_point, end_point,
                      dpi=500, title=None, width=5, height=3, orig_nams=[],
                      keep_numbers=False, force_numbers=False, palette="te_trimmer"):
    '''
    Draws a "mini-alignment" image showing a small representation of the
    whole alignment so that gaps and poorly aligned regions are visible.

    Parameters:
    - arr: np.array
        The alignment stored as a NumPy array
    - outfile: str
        Path to the output PDF file (No need to provide .pdf extension)
    - ... [rest of your parameters] ...

    Returns:
    - None
    '''

    ali_height, ali_width = np.shape(arr)

    fontsize = 2000 / dpi
    om = math.floor(math.log10(ali_height))
    tickint = 1 if om == 0 or force_numbers else 10 if om == 1 else 100
    lineweight_h = 10 / ali_height
    lineweight_v = 10 / ali_width

    f = plt.figure(figsize=(width, height), dpi=dpi)
    a = f.add_subplot(1, 1, 1)
    a.set_xlim(-0.5, ali_width)
    a.set_ylim(-0.5, ali_height - 0.5)

    arr2, cm = arrNumeric(arr, palette)
    a.imshow(arr2, cmap=cm, aspect='auto', interpolation='nearest')

    arrow_length = 13
    a.annotate('Start crop Point', xy=(start_point, ali_height - 0.5), xytext=(0, arrow_length),
               textcoords='offset points',
               arrowprops=dict(facecolor='red', edgecolor='red', width=0.3, headwidth=2, headlength=2),
               ha='center', color='red', va='bottom', fontsize=3)
    a.annotate('End crop Point', xy=(end_point, ali_height - 0.5), xytext=(0, arrow_length),
               textcoords='offset points',
               arrowprops=dict(facecolor='blue', edgecolor='blue', width=0.3, headwidth=2, headlength=2),
               ha='center', color='blue', va='bottom', fontsize=3)

    f.subplots_adjust(top=0.85, bottom=0.1, left=0.1, right=0.95)
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



#########################################################################
# This code is derived from CIAlign: https://github.com/KatyBrown/CIAlign
#########################################################################

def base():
    '''
    Returns the hexadecimal values for black and white

    Parameters
    ----------
    None

    Returns
    -------
    dict
        A dictionary containing the hexadecimal values for black and white
    '''

    return {'black': '#000000',
            'white': '#FFFFFF'}


def CBSafe():
    '''
    Returns the hexadecimal values for a colorblind-friendly color palette

    Parameters
    ----------
    None

    Returns
    -------
    dict
        A dictionary containing the hexadecimal values for the colors used
        in the CIAlign mini-alignments
    '''

    b = base()
    b.update({'yellow_nt': "#c9c433",
              'green_nt': "#56ae6c",
              'red_nt': "#a22c49",
              'blue_nt': "#0038a2",
              'grey_nt': "#6979d3",
              'red_aa': "#a22c49",
              'yellow_aa': "#c9c433",
              'blue_aa':  "#0038a2",
              'orange_aa': "#e57700",
              'midblue_aa': "#589aab",
              'cyan_aa': "#50d3cb",
              'lightgrey_aa': '#eae2ea',
              'green_aa': "#56ae6c",
              'darkgrey_aa': "#888988",
              'purple_aa': '#89236a',
              'paleblue_aa': '#e669ca',
              'peach_aa': "#ffc4a9",
              'tan_aa': "#936e23",
              'remove_insertions': "#9db341",
              'remove_divergent': "#7066bc",
              'crop_ends': '#020545',
              'remove_gap_only': '#f9c1d2',
              'remove_short': "#c85133",
              'crop_divergent': '#ff00d1',
              'user': '#fff4a3'})
    return b


def Bright():
    '''
    Returns the hexadecimal values for a colorblind-friendly color palette

    Parameters
    ----------
    None

    Returns
    -------
    dict
        A dictionary containing the hexadecimal values for the colors used
        in the CIAlign mini-alignments
    '''
    b = base()
    b.update({'yellow_nt': "#ffd500",
              'red_nt': "#64bc3c",
              'green_nt': "#f20707",
              'blue_nt': "#0907f2",
              'grey_nt': "#c7d1d0",
              'red_aa': "#f20707",
              'yellow_aa': "#ffd500",
              'blue_aa':  "#0907f2",
              'orange_aa': "#f4aa03",
              'midblue_aa': "#03b5f4",
              'cyan_aa': "#03f4dd",
              'lightgrey_aa': '#f0f4f4',
              'green_aa': "#64bc3c",
              'darkgrey_aa': "#808080",
              'purple_aa': '#6f0cac',
              'paleblue_aa': '#cde3f8',
              'peach_aa': "#f8e7cd",
              'tan_aa': "#d2a867",
              'remove_insertions': "#9db341",
              'remove_divergent': "#7066bc",
              'crop_ends': '#020545',
              'remove_gaponly': '#f9c1d2',
              'remove_short': "#c85133",
              'crop_divergent': '#ff00d1',
              'user': '#fff4a3'})
    return b

def te_trimmer():
    b = base()
    b.update({'yellow_nt': "#949494",
              'red_nt': "#FF6666",
              'green_nt': "#00CC00",
              'blue_nt': "#6161ff",
              'grey_nt': "#c7d1d0",
              'red_aa': "#f20707",
              'yellow_aa': "#ffd500",
              'blue_aa': "#0907f2",
              'orange_aa': "#f4aa03",
              'midblue_aa': "#03b5f4",
              'cyan_aa': "#03f4dd",
              'lightgrey_aa': '#f0f4f4',
              'green_aa': "#64bc3c",
              'darkgrey_aa': "#808080",
              'purple_aa': '#6f0cac',
              'paleblue_aa': '#cde3f8',
              'peach_aa': "#f8e7cd",
              'tan_aa': "#d2a867",
              'remove_insertions': "#9db341",
              'remove_divergent': "#7066bc",
              'crop_ends': '#020545',
              'remove_gaponly': '#f9c1d2',
              'remove_short': "#c85133",
              'crop_divergent': '#ff00d1',
              'user': '#fff4a3'})
    return b

