#! /usr/bin/env python
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import Cialign_palettes as palettes
import math
matplotlib.use('Agg')


def getPalette(palette='CBS'):
    '''
    Generates a dictionary which assigns a name to each colour using a colour
    blindness safe palette, generated using
    https://medialab.github.io/iwanthue/
    Parameters
    ----------
    palette: str
        The ID of the palette to be used, currently only colour blind safe
        (CBS) is implemented.

    Returns
    -------
    dict
        Dictionary where keys are names of colours and
        values are hexadecimal codes for colours
    '''
    if palette.lower() == 'cbs':
        p = palettes.CBSafe()
    if palette.lower() == 'bright':
        p = palettes.Bright()
    if palette.lower() == 'te_trimmer':
        p = palettes.te_trimmer()
    return p

def getNtColours(palette='CBS'):
    '''
    Generates a dictionary which assigns a colour to each nucleotide (plus grey
    for "N" and white for "-")
    Parameters
    ----------
    pal: str
        A string designating which palette to use, currently only colour blind
        safe (CBS) is implemented.

    Returns
    -------
    dict
        Dictionary where keys are single letter nucleotide codes and
        values are hexadecimal codes for colours
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
    Convert an alignment into a numpy array.

    Parameters
    ----------
    infile: string
        path to input alignment file in FASTA format
    log: logging.Logger
        An open log file object

    Returns
    -------
    arr: np.array
        2D numpy array in the same order as fasta_dict where each row
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
                    print (nseq, sl, psl)
                    raise ValueError ("""
ERROR: The sequences you provided may not be aligned - all the sequences \
are not the same length""")

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
        print (nseq, sl, psl)
        raise ValueError ("""
ERROR: The sequences you provided may not be aligned - all the sequences \
are not the same length""")
    seqs.append(np.array([s.upper() for s in seq]))
    nams.append(nam)
    arr = np.array(seqs[1:])
    return arr, nams[1:]

def arrNumeric(arr, palette='CBS'):
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
    D = getNtColours(palette)

    # retrieve the colours for the colour map
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


def drawMiniAlignment(arr, nams, outfile,
                      dpi=300, title=None, width=5, height=3, ret=False, orig_nams=[],
                      keep_numbers=False, force_numbers=False, palette="te_trimmer"):
    '''
    Draws a "mini alignment" image showing a small representation of the
    whole alignment so that gaps and poorly aligned regions are visible.

    Parameters
    ----------
    arr: np.array
        The alignment stored as a numpy array
    nams: list
        The names of the sequences in the alignment
    log: logging.Logger
        The open log file object
    outfile: str
        Path to the output file
    typ: str
        Either 'aa' - amino acid - or 'nt' - nucleotide
    dpi: int
        DPI for the output image
    title: str
        Title for the output image
    width: int
        Width of the output image
    height: int
        Height of the output image
    markup: bool
        Should the deleted rows and columns be marked on the output image?
    markupdict: dict
        Dictionary where the keys are function names and the values are
        lists of columns, rows or positions which have been removed
    ret: bool
        Return the subplot as a matplotlib object, used to make plots when
        using this function directly rather than the CIAlign workflow
    orig_nams: list
        List of names in the original input plot, used if keep_numbers is
        switched on to keep the original numbering scheme
    keep_numbers: bool
        Number the sequences (rows) based on the original CIAlign input rather
        than renumbering.
    palette: str
        Colour palette, CBS or Bright
    Returns
    -------
    None

    '''
    ali_height, ali_width = np.shape(arr)

    # font size needs to scale with DPI
    fontsize = 1500 / dpi

    # what is the order of magnitude of the number of sequences (rows)
    # 0 = 1 - 10 sequences - label every row
    # 1 = 10 - 100 sequences - label every 10th row
    # 2+ = 100+ sequences - label every 100th row

    om = math.floor(math.log10(ali_height))
    tickint = 1 if om == 0 or force_numbers else 10 if om == 1 else 100

    # use thinner lines for bigger alignments
    lineweight_h = 10 / ali_height
    lineweight_v = 10 / ali_width

    f = plt.figure(figsize=(width, height), dpi=dpi)
    a = f.add_subplot(1, 1, 1)
    a.set_xlim(-0.5, ali_width)
    a.set_ylim(-0.5, ali_height-0.5)

    # generate the numeric version of the array
    arr2, cm = arrNumeric(arr, palette)
    # display it on the axis
    a.imshow(arr2, cmap=cm, aspect='auto', interpolation='nearest')

    # these are white lines between the bases in the alignment - as the
    # image is actually solid
    a.hlines(np.arange(-0.5, ali_height), -0.5,
             ali_width, lw=lineweight_h, color='white',
             zorder=100)
    a.vlines(np.arange(-0.5, ali_width), -0.5,
             ali_height, lw=lineweight_v, color='white',
             zorder=100)

    # aesthetics
    a.spines['right'].set_visible(False)
    a.spines['top'].set_visible(False)
    a.spines['left'].set_visible(False)

    if title:
        f.suptitle(title, fontsize=fontsize*1.5, y=0.92)
    for t in a.get_xticklabels():
        t.set_fontsize(fontsize)
    a.set_yticks(np.arange(ali_height-1, -1, -tickint))
    x = 1
    if tickint == 1:
        if keep_numbers:
            labs = []
            for nam in orig_nams:
                if nam in nams:
                    labs.append(x)
                x += 1
            a.set_yticklabels(labs,
                              fontsize=fontsize*0.75)
        else:
            a.set_yticklabels(np.arange(1, ali_height+1, tickint),
                              fontsize=fontsize*0.75)
    else:
        a.set_yticklabels(np.arange(0, ali_height, tickint), fontsize=fontsize)

    f.tight_layout()
    f.savefig(outfile, dpi=dpi, bbox_inches='tight')
    if ret:
        return f
    plt.close()


input_file = "/Users/panstrugamacbook/Documents/PhD_project_files/TE_Trimmer/Select_gap_block_test/nest_1.fasta.blast.bed.uniq.bed.fil.bed_0_0_cl.fa_maf.fa_gap_rm.fa_cl.fa"
output_file = "/Users/panstrugamacbook/Documents/PhD_project_files/TE_Trimmer/Cluster_TKlinker/test_cialign"
test1, test2 = FastaToArray(input_file)
test3, test4 = arrNumeric(test1)
print(test3)
print(test4)

drawMiniAlignment(test1, test2, output_file)