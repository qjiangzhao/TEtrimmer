import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
import subprocess
import click
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
import warnings
from Bio import BiopythonDeprecationWarning

# Suppress all deprecation warnings
warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)


def con_generater(input_file, output_dir, threshold=0.8, ambiguous="N"):
    # Read input file
    alignment = AlignIO.read(input_file, "fasta")

    # Generate a summary of the alignment
    summary_align = AlignInfo.SummaryInfo(alignment)

    # Calculate consensus sequence with the specified threshold
    consensus = summary_align.dumb_consensus(threshold=threshold, ambiguous=ambiguous)

    # Create SeqRecord for consensus sequence
    consensus_record = SeqRecord(consensus, id=os.path.basename(input_file), description="")

    # Determine the length of the consensus sequence
    consensus_length = len(consensus)

    # Write consensus sequence to a FASTA file
    output_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_co.fa")
    with open(output_file, "w") as file:
        SeqIO.write(consensus_record, file, "fasta")

    return output_file, consensus_length


def check_database(genome_file, output_dir=None):

    # Define database path
    genome_n = os.path.basename(genome_file)
    if output_dir is None:
        database_path = genome_file
    else:
        database_path = os.path.join(output_dir, genome_n)

    database_file = f"{database_path}.nin"

    if not os.path.isfile(database_file):

        try:
            makeblastdb_cmd = f"makeblastdb -in {genome_file} -dbtype nucl -out {database_path} "
            subprocess.run(makeblastdb_cmd, shell=True, check=True, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, text=True)
            return database_path

        except FileNotFoundError:
            click.echo("'makeblastdb' command not found. Please ensure 'makeblastdb' is installed correctly.")
            raise Exception

        except subprocess.CalledProcessError as e:
            click.echo(f"makeblastdb failed with exit code {e.returncode}")
            click.echo(e.stderr)
            return False
    else:
        return database_path


def blast(seq_file, genome_file, output_dir, e_value=1e-40):
    """
    Runs BLAST calling a specified task type and saves the results as a BED file.

    :param seq_file: str, path to input FASTA file
    :param genome_file: str, path to genome FASTA file
    :param output_dir: str, prefix for output files
    """
    input_file = seq_file
    input_file_n = os.path.basename(input_file)
    blast_hits_count = 0

    # define blast outfile
    blast_out_file = os.path.join(output_dir, f"{os.path.basename(input_file)}.b")

    # Modify the blast command to include the specified task
    blast_cmd = (f"blastn -query {input_file} -db {genome_file} "
                 f"-outfmt \"10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore\" "
                 f"-evalue {str(e_value)} | awk 'BEGIN{{print \"V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13\"}};{{print}}' | sed 's/#/-/g' > {blast_out_file}")

    try:
        subprocess.run(blast_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    except FileNotFoundError:
        click.echo("'blastn' command not found. Please ensure 'blastn' is installed correctly.")
        raise Exception

    except subprocess.CalledProcessError as e:
        click.echo(f"\nBLAST failed for {input_file_n} with error code {e.returncode}")
        click.echo(e.stderr)
        raise Exception

    # Check the number of BLAST hits
    with open(blast_out_file) as blast_file:
        for _ in blast_file:
            blast_hits_count += 1

    if blast_hits_count > 1:
        return blast_out_file
    else:
        click.echo("Blast hits number is 0, blast plots won't show")
        return False


#####################################################################################################
# Code block: blast with genome
#####################################################################################################

def blast_with_genome(input_file, output_dir, genome_file):

    # Generate consensus sequence
    con_seq, cons_len = con_generater(input_file, output_dir, threshold=0.6)

    # Check genome blast database availability
    blast_database = check_database(genome_file)

    # Run blast
    if blast_database:
        blast_file = blast(con_seq, blast_database, output_dir)

        if blast_file:
            return blast_file, cons_len
        else:
            return False
    else:
        return False


#####################################################################################################
# Code block: self blast
#####################################################################################################

def self_blast(input_file, output_dir):

    input_n = os.path.basename(input_file)

    output_dir_self_blast = os.path.join(output_dir, f"{input_n}_self_blast")

    os.makedirs(output_dir_self_blast, exist_ok=True)

    # Generate consensus sequence. Use smaller threshold for better dot plot
    con_seq, cons_len = con_generater(input_file, output_dir_self_blast, threshold=0.5)

    # Build blast database
    blast_database = check_database(con_seq, output_dir_self_blast)

    # Run blast. Use bigger e_value for self blast
    if blast_database:
        blast_file = blast(con_seq, blast_database, output_dir_self_blast, e_value=0.05)

        if blast_file:
            return blast_file
        else:
            return False

    else:
        return False

#####################################################################################################
# Code block: genome blast plot
#####################################################################################################


# Function to create hover text
def create_hovertext(row):
    divergence = 100 - row['V3']
    return f"Start: {row['V7']}<br>End: {row['V8']}<br>Divergence: {divergence:.2f}%"


def blast_plot(df, cons_len, full_len_thr=0.8):
    # Calculate full matches based on the threshold and consensus length
    full = df[(abs(df['V7'] - df['V8']) >= full_len_thr * cons_len)]

    # Calculate the max number of divergence for the y-axis range
    y_max = max(100 - df['V3'])

    # Create a figure
    fig = go.Figure()

    # Add segments for all blast hits
    for i in df.index:
        fig.add_trace(go.Scatter(
            x=[df.loc[i, 'V7'], df.loc[i, 'V8']],
            y=[100 - df.loc[i, 'V3'], 100 - df.loc[i, 'V3']],
            mode='lines',
            line=dict(color='black'),
            hoverinfo='text',
            text=create_hovertext(df.loc[i]),
            showlegend=False
        ))

    # Add segments for full matches
    for i in full.index:
        fig.add_trace(go.Scatter(
            x=[full.loc[i, 'V7'], full.loc[i, 'V8']],
            y=[100 - full.loc[i, 'V3'], 100 - full.loc[i, 'V3']],
            mode='lines',
            line=dict(color='red', width=3),
            hoverinfo='text',
            text=create_hovertext(full.loc[i]),
            showlegend=False
        ))

    # Define plot title
    plot_title = f"TE: {df.iloc[0,0]} | consensus size: {cons_len}bp;<br> fragments number: {len(df)}; full length: {len(full)} (>={int(full_len_thr*cons_len)}bp)"

    # Update the layout to add titles and labels with adjusted axis ranges
    fig.update_layout(
        title=plot_title,
        xaxis_title="TE consensus (bp)",
        yaxis_title="Divergence to consensus (%)",
        xaxis=dict(
            title_font=dict(size=16, family="Arial, sans-serif", color="black"),
            showgrid=False,
            range=[0, cons_len + 10],
            linewidth=1,
            linecolor='black',
            mirror=True
        ),
        yaxis=dict(
            showgrid=False,
            title_font=dict(size=16, family="Arial, sans-serif", color="black"),
            range=[-(y_max / 30), y_max + (y_max / 30)],
            linewidth=1,
            linecolor='black',
            mirror=True
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin=dict(l=60, r=60, t=60, b=60),  # Adjust margins to ensure frame visibility
        font=dict(family="Arial, sans-serif", size=12, color="black"),
        width=600,
        height=600
    )
    #fig.show()
    return fig, y_max, plot_title

#####################################################################################################
# Code block: coverage plot
#####################################################################################################


def coverage_plot(df, cons_len):
    # Initialize coverage array
    coverage = np.zeros((len(df), cons_len))

    # Populate coverage array
    for i in range(len(df)):
        start = df.iloc[i]['V7'] - 1  # Convert to zero-indexed
        end = df.iloc[i]['V8']
        coverage[i, start:end] = 1  # Set covered positions to 1

    # Sum coverage over all rows to get total coverage at each position
    total_coverage = np.sum(coverage, axis=0)
    max_coverage = np.max(total_coverage)

    # Create a figure for plotting
    fig = go.Figure()

    # Add the coverage line
    fig.add_trace(go.Scatter(
        x=np.arange(1, cons_len+1),  # X-axis from 1 to cons_len
        y=total_coverage,
        mode='lines',
        line=dict(color='black', width=3),
        showlegend=False
    ))

    # Update layout to add titles and labels
    fig.update_layout(
        title="TE Consensus Genomic Coverage Plot",
        xaxis_title="Position along TE consensus (bp)",
        yaxis_title="Coverage (bp)",
        xaxis=dict(
            title_font=dict(size=16, family="Arial, sans-serif", color="black"),
            showgrid=False,
            range=[0, cons_len + 10],
            linewidth=1,
            linecolor='black',
            mirror=True
        ),
        yaxis=dict(
            showgrid=False,
            title_font=dict(size=16, family="Arial, sans-serif", color="black"),
            range=[-(max_coverage / 30), max_coverage + (max_coverage / 30)]
            ,
            linewidth=1,
            linecolor='black',
            mirror=True
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin=dict(l=60, r=60, t=60, b=60),  # Adjust margins to ensure frame visibility
        font=dict(family="Arial, sans-serif", size=12, color="black"),
        width=600,
        height=600
    )
    #fig.show()
    return fig, max_coverage


#####################################################################################################
# Code block: self dot plot
#####################################################################################################


def dot_plot(df):
    fig = go.Figure()

    # Loop through data and add segments
    for i in range(len(df)):
        segment = df.iloc[i]
        if segment['V10'] > segment['V9']:
            color = "black"
        else:
            color = "#009E73"  # Use different color for TIR
        fig.add_trace(go.Scatter(
            x=[segment['V7'], segment['V8']],
            y=[segment['V9'], segment['V10']],
            mode='lines',
            line=dict(color=color, width=3),
            hoverinfo='text',
            text=f"x left: {segment['V7']}, y left: {segment['V9']}<br>x right: {segment['V8']}, y right: {segment['V10']}",
            showlegend=False
        ))

    fig.update_layout(
        title=f"Self dot plot",
        xaxis_title="TE consensus self dotplot (bp)",
        yaxis_title="TE consensus self dotplot (bp)",
        xaxis=dict(
            title_font=dict(size=16, family="Arial, sans-serif", color="black"),
            showgrid=False,
            range=[0, max(df['V8'])],
            linewidth=1,
            linecolor='black',
            mirror=True
        ),
        yaxis=dict(
            showgrid=False,
            title_font=dict(size=16, family="Arial, sans-serif", color="black"),
            range=[0, max(df['V8'])],
            linewidth=1,
            linecolor='black',
            mirror=True
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin=dict(l=60, r=60, t=60, b=60),  # Adjust margins to ensure frame visibility
        font=dict(family="Arial, sans-serif", size=12, color="black"),
        width=600,
        height=600
    )
    #fig.show()
    return fig


#####################################################################################################
# Code block: TEtrimmer GUI plotter
#####################################################################################################

def GUI_plotter(input_file, output_dir, genome_file):

    fig_blast = None
    fig_coverage = None
    fig_dot = None

    # Do genome blast including database construction and blastn
    genome_blast = blast_with_genome(input_file, output_dir, genome_file)

    if genome_blast:
        genome_blast_out, cons_len = genome_blast
        df_genome_blast_out = pd.read_csv(genome_blast_out)

        # Genome blast plot and blast coverage plot
        fig_blast, divergence_max, blast_plot_title = blast_plot(df_genome_blast_out, cons_len)
        fig_coverage, coverage_max = coverage_plot(df_genome_blast_out, cons_len)

    # Perform self blast
    dotplot_blast = self_blast(input_file, output_dir)

    if dotplot_blast:
        # self blast dot plot
        df_dotplot_blast = pd.read_csv(dotplot_blast)
        fig_dot = dot_plot(df_dotplot_blast)

    # To avoid error caused by single plot problem, check plots availability
    if fig_blast is not None and fig_coverage is not None and fig_dot is not None:
        # Create a 1x3 subplot layout
        fig = make_subplots(rows=1, cols=3, subplot_titles=(blast_plot_title, "BLAST coverage Plot", "Self Dot Plot"))

        # Merge plots into the combined subplot layout
        for i, subplot in enumerate([fig_blast, fig_coverage, fig_dot], start=1):
            for trace in subplot.data:
                fig.add_trace(trace, row=1, col=i)

        # Update the font size of subplot titles
        for annotation in fig['layout']['annotations']:
            annotation['font'] = dict(size=13, family="Arial, sans-serif", color="black")

            # Configuring subplot 1
        fig.update_xaxes(title_text="TE consensus (bp)",
                         showgrid=False,
                         title_font=dict(size=16, family="Arial, sans-serif", color="black"),
                         range=[0, cons_len + 10],
                         linewidth=1,
                         linecolor='black',
                         mirror=True,
                         row=1, col=1)
        fig.update_yaxes(title_text="Divergence to consensus (%)",
                         showgrid=False,
                         title_font=dict(size=16, family="Arial, sans-serif", color="black"),
                         range=[-(divergence_max / 30), divergence_max + (divergence_max / 30)],
                         linewidth=1,
                         linecolor='black',
                         mirror=True,
                         row=1, col=1)

        # Configuring subplot 2
        fig.update_xaxes(title_text="TE consensus (bp)",
                         showgrid=False,
                         title_font=dict(size=16, family="Arial, sans-serif", color="black"),
                         range=[0, cons_len + 10],
                         linewidth=1,
                         linecolor='black',
                         mirror=True,
                         row=1, col=2)
        fig.update_yaxes(title_text="Coverage (bp)",
                         showgrid=False,
                         title_font=dict(size=16, family="Arial, sans-serif", color="black"),
                         range=[-(coverage_max / 30), coverage_max + (coverage_max / 30)],
                         linewidth=1,
                         linecolor='black',
                         mirror=True,
                         row=1, col=2)

        # Configuring subplot 3
        fig.update_xaxes(title_text="TE consensus (bp)",
                         showgrid=False,
                         title_font=dict(size=16, family="Arial, sans-serif", color="black"),
                         range=[0, cons_len + 10],
                         linewidth=1,
                         linecolor='black',
                         mirror=True,
                         row=1, col=3)
        fig.update_yaxes(title_text="TE consensus (bp)",
                         showgrid=False,
                         title_font=dict(size=16, family="Arial, sans-serif", color="black"),
                         range=[0, cons_len + 10],
                         linewidth=1,
                         linecolor='black',
                         mirror=True,
                         row=1, col=3)

        # Update layout if needed
        fig.update_layout(height=600, width=1800, title_text="", plot_bgcolor='white', paper_bgcolor='white')
        fig.show()
    # Show plots separately when one has a problem
    else:
        if fig_blast is not None:
            fig_blast.show()
        if fig_coverage is not None:
            fig_coverage.show()
        if fig_dot is not None:
            fig_dot.show()


"""

input_file = "/Users/panstrugamacbook/Documents/TE_Trimmer/Plotly_figure/clipboard-alignment_3073100289831958694.fasta"
genome = "/Users/panstrugamacbook/Documents/TE_Trimmer/bgh_genome/bgh_dh14_v4.fa"

output_dir = "/Users/panstrugamacbook/Documents/TE_Trimmer/Plotly_figure"

GUI_plotter(input_file, output_dir, genome)
"""

