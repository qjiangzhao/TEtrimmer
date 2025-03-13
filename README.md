[![Bioconda Version](https://img.shields.io/conda/vn/Bioconda/tetrimmer.svg)](https://anaconda.org/bioconda/tetrimmer)
[![Singularity Badge](https://img.shields.io/badge/Singularity-available-brightgreen.svg)](https://quay.io/repository/biocontainers/tetrimmer?tab=tags&tag=1.4.0--hdfd78af_0)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](./LICENSE)
[![Manual Available](https://img.shields.io/badge/Manual-available-brightgreen.svg)](./docs/TEtrimmerv1.4.0Manual.pdf)

# TEtrimmer
Automated and interactive curation of transposable elements.

## Contents

- [Introduction](#introduction)
- [Manual](#manual)
- [Installation](#installation)
  - [1. Conda](#1-conda)
  - [2. Singularity](#2-singularity)
  - [3. Docker](#3-docker)
- [Example Usage](#example-usage)
  - [Inputs](#inputs)
  - [Resume run](#resume-run)
  - [Combining TE libraries](#combining-te-libraries)
  - [Outputs](#outputs)
- [TEtrimmerGUI](#tetrimmergui)
  - [Curate existing TE libraries](#curate-existing-te-libraries)
- [Hardware requirements](#hardware-requirements)
- [Citation](#citation)

## Introduction

Many tools have been developed for the discovery and annotation of transposable elements (TEs).

However, the high-quality TE consensus library construction still requires manual curation of TEs,
which is time-consuming and needs experts with an in-depth understanding of TE biology.


TEtrimmer is a powerful software designed to automate the manual curation of TEs. The input can be a TE library from
*de novo* TE discovery tools, such as EDTA and RepeatModeler2, or a TE library from closely related species.

For each input consensus sequence, TEtrimmer automatically performs BLASTN search,  sequence extraction, extension,
multiple sequence alignment (MSA), MSA clustering, MSA cleaning, TE boundary definition, and TE classification.
TEtrimmer also provides a graphical user interface (GUI) to inspect and improve predicted TEs, which can assist
achieving manual curation-level TE consensus libraries easily.

<a href="https://www.youtube.com/watch?v=PgKw76gmzl8&t=17s&ab_channel=ZhaoJiang">
  <img src="docs/TEtrimmerFlowchart_video.png" alt="TEtrimmerFlowchart" style="width: 800px;">
</a>


## Manual
For detailed instructions, including installation steps, usage options, example outputs, and more,
please refer to [TEtrimmerv1.4.0Manual.pdf](https://github.com/qjiangzhao/TEtrimmer/blob/main/docs/TEtrimmerv1.4.0Manual.pdf)

**Note: Manual is out of date.**

## Installation
TEtrimmer can be installed locally with Conda, or run from a Container with Singularity or Docker.

You will need local copies of the Pfam-A and CDD databases available. See example usage section for download instructions.

### 1. Conda

Create a Conda environment with all requirements for running TEtrimmer.

Install notes:

- You will need to install a conda distribution, such as [miniforge3](https://conda-forge.org/download/).

- Some dependencies are not available for Windows or M-series Macs (arm64 processors)*

- Ignore `ClobberError` or `ClobberWarning` messages during install. BLAST and rmBLAST contain some duplicate tools.

#### 1A. Install from Bioconda

Install TEtrimmer and dependencies directly from Bioconda. (Many thanks to [HangXue](https://github.com/hangxue-wustl))

```bash
# Create a new conda env + bioconda::TEtrimmer
conda create --name TEtrimmer -c bioconda -c conda-forge TEtrimmer

# Activate the env
conda activate TEtrimmer

# Test
TEtrimmer --help
```

#### 1B. Install latest version into local Conda env

Use the [`environment.yml`](https://github.com/qjiangzhao/TE-Trimmer/blob/main/environment.yml) file to create a new env called "TEtrimmer" with dependencies only.

```bash
# Create new conda environment
conda env create -f environment.yml

# Activate TEtrimmer env
conda activate TEtrimmer
```

Now you can install TEtrimmer from this github repo:

```bash
# Clone the latest version of TEtrimmer from Github
git clone https://github.com/qjiangzhao/TEtrimmer.git && cd TEtrimmer

# pip install package
pip install -e '.[tests]'

# Call TEtrimmer
TEtrimmer --help
```

### 2. Singularity

```bash
# Download and generate TEtrimmer "sif" file
singularity pull docker://quay.io/biocontainers/tetrimmer:1.4.0--hdfd78af_0

# Run TEtrimmer based on sif file
# If <your_path_to_store_PFAM_database> doesn't contain PFAM database
# TEtrimmer can automatically download PFAM to <your_path_to_store_PFAM_database>

singularity exec --writable-tmpfs \
--bind <your_path_contain_genome_file>:/genome \
--bind <your_path_contain_input_TE_library_file>:/input \
--bind <your_output_path>:/output \
--bind <your_path_to_store_PFAM_database>:/pfam \
<your_path_contain_sif_file>/tetrimmer_1.4.0--hdfd78af_0.sif \
TEtrimmer \
-i /input/<TE_library_name.fasta> \
-g /genome/<genome_file_name.fasta> \
-o /output \
--pfam_dir /pfam \
-t 20 --classify_all

# The Singularity image includes the TEtrimmer source code, but the version inside may be outdated.
# If you want to run the latest version of TEtrimmer via the singularity image
# Clone the new version of TEtrimmer from Github

git clone https://github.com/qjiangzhao/TEtrimmer.git

singularity exec --writable-tmpfs \
--bind <your_path_to_cloned_TEtrimmer_folder_which_contain_TEtrimmer.py>:/TEtrimmer_cloned \
--bind <your_path_contain_genome_file>:/genome \
--bind <your_path_contain_input_TE_library_file>:/input \
--bind <your_output_path>:/output \
--bind <your_path_to_store_PFAM_database>:/pfam \
<your_path_contain_sif_file>/tetrimmer_1.4.0--hdfd78af_0.sif \
python -m TEtrimmer_cloned_dir.src.TEtrimmer.tetrimmer.TEtrimmer \
-i /input/<TE_library_name.fasta> \
-g /genome/<genome_file_name.fasta> \
-o /output \
--pfam_dir /pfam \
-t 20 --classify_all

```

### 3. Docker

```bash
# Download TEtrimmer docker image
docker pull quay.io/biocontainers/tetrimmer:1.4.0--hdfd78af_0
docker run -it --name TEtrimmer -v <bind_your_path>:/data quay.io/biocontainers/tetrimmer:1.4.0--hdfd78af_0
# Then you can run TEtrimmer inside TEtrimmer container
# Please note: Run TEtrimmer via Docker is relatively slower than Conda and Singularity.
```

## Example Usage

In this example we will clone the development repo and create a fresh Conda environment to run `TEtrimmer`.


```bash
# Clone this fork and navigate into repo dir
git clone git@github.com:qjiangzhao/TEtrimmer.git && cd TEtrimmer

# Create env "TEtrimmer"
conda env create -f environment.yml

# Activate new env
conda activate TEtrimmer

# Install package in editable mode
pip install -e '.[tests]'

# Confirm installation
TEtrimmer --help
```


Fetch Pfam-A database. We will store the database files in the default location. Any database path can be specified with `--pfam_dir`.

``` bash
# Create default db dir
mkdir src/TEtrimmer/pfam_database

# Download files
curl -o src/TEtrimmer/pfam_database/Pfam-A.hmm.gz https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
curl -o src/TEtrimmer/pfam_database/Pfam-A.hmm.dat.gz https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz

# Unzip PFAM files
gunzip src/TEtrimmer/pfam_database/Pfam-A.hmm.gz
gunzip src/TEtrimmer/pfam_database/Pfam-A.hmm.dat.gz
```

Now we can run `TEtrimmer` on the test data files[test_input.fa.gz](https://github.com/qjiangzhao/TEtrimmer/blob/main/tests/data/test_input.fa.gz) and [test_genome.fa.gz](https://github.com/qjiangzhao/TEtrimmer/blob/main/tests/data/test_genome.fa.gz).

``` bash
# Run TEtrimmer
TEtrimmer --input_file tests/data/test_input.fa.gz \
--genome_file tests/data/test_genome.fa.gz \
--output_dir output_directory \
--num_threads 2 \
--classify_all \
--classify_unknown \
--cons_thr 0.8 \
--mini_orf 200 \
--max_msa_lines 100 \
--genome_anno \
--preset divergent \
--pfam_dir src/TEtrimmer/pfam_database
```


After running `TEtrimmer `the processes TE families can be further curated with TEtrimmerGUI.

First we will need to fetch the CDD database of TE domains. Again we will store it in the default location.

```bash
# Create default db path
mkdir src/TEtrimmer/tetrimmerGUI/cdd_database

# Download or copy existing cdd.tar.gz into local dir
curl -o src/TEtrimmer/tetrimmerGUI/cdd_database/cdd.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cdd.tar.gz
```

Finally, we can launch `TEtrimmerGUI`. We provide the path to the CDD database with `--cdd_dir`.

``` bash
# Launch GUI

TEtrimmerGUI -i output_directory/TEtrimmer_for_proof_curation \
-g tests/data/test_genome.fa.gz  \
--output_dir output_dir_GUI \
-cdd_dir src/TEtrimmer/tetrimmerGUI/cdd_database
```


### Inputs

- **Genome file**: The genome sequence in FASTA format (.fa or .fasta). May be GZIP compressed (.fa.gz or .fasta.gz)
- **TE consensus library**: TEtrimmer uses the TE consensus library from *de novo* TE annotation tools, like `RepeatModeler` or `EDTA`, as input.
For this reason, you have to run `RepeatModeler` or other TE annotation software first.

```bash
# RepeatModeler is installed in the TEtrimmer conda environment. Below is an example command for running RepeatModeler.

# Build genome database index files
BuildDatabase -name <genome_file_database_name> <genome_file.fa>

# Run RepeatModeler
RepeatModeler -database <genome_file_database_name> \
              -threads 20 \
              -LTRStruct
# Output: TE_consensus_library.fa
```

Example:

```bash
TEtrimmer --input_file TE_consensus_library.fa \
          --genome_file genome_file.fa \
          --output_dir output_directory \
          --num_threads 20 \
          --classify_all
```

### Resume run

If you want to **continue the analysis based on previous unfinished results in the same directory:**:

```bash
TEtrimmer --input_file TE_consensus_library.fa \
          --genome_file genome_file.fa \
          --output_dir directory_contains_previous_unfinished_results \
          --num_threads 20 \
          --classify_all \
          --continue_analysis
```

### Combining TE libraries

If you want to combine TE predictions from different sources for the input file, we recommend using `--dedup` to remove duplicate sequences before processing.

This step can potentially save overall run time.

**Note: TEtrimmer only accepts a single input file, you have to combine TE libraries in advance.**

```bash
# Combine TE libraries
cat library_1.fa library_2.fa > TE_consensus_library.fa

# Run TEtrimmer with dedup mode
TEtrimmer --input_file TE_consensus_library.fa \
          --genome_file genome_file.fa \
          --output_dir output_directory \
          --num_threads 20 \
          --classify_all
          --dedup
```

### Outputs

- 📁**Classification** - *This folder is used for TE classifications.*
- 📁**Multiple_sequence_alignment** - *All raw files will be stored in this folder if < --debug > is enabled.*
  - 📄**error_file.txt** - *Error file to store all error messages, only visible if errors were found.*
- 📁**Single_fasta_files** - *All sequences in the input file will be separated into single FASTA files and stored here.*
- 📁**TEtrimmer_for_proof_curation** - *This folder contains files used for manual inspection of TEtrimmer annotations.*
  - 📁**Annotation_perfect** - *Four files are associated with each sequence as showed below.*
    - 📄**TE_name.raw.fa** - *Multiple sequence alignment file before TE boundary definition.*
    - 📄**TE_name.fa** - *Multiple sequence alignment file after TE boundary definition, which is used to generate the consensus sequence.*
    - 📄**TE_name.pdf** - *Plot file used to evaluate output.*
    - 📄**TE_name.cluster.fa** - *Multiple sequence alignment file before clustering.*
  - 📁**Annotation_good**
  - 📁**Annotation_check_recommended**
  - 📁**Annotation_check_required**
  - 📁**Clustered_proof_curation** - *This folder contains all the output files from folder "Annotation_perfect", "Annotation_good", "Annotation_check_recommended", and "Annotation_check_required". The difference is TEtrimmer group similar output TEs into one "Cluster", which can make it easier to compare similar outputs.*
  - 📁**TE_low_copy** - *This folder contains low copy TEs.*
  - 📁**TE_skipped** - *Contains TE_Aid plots for all skipped TEs.*
- 📁**HMM** - *This folder is used to store Hidden Markov Model file. Only visible when < --hmm > is enabled.*
- 📄**Sequence_name_mapping.txt** - *This file connects the input sequence names with the modified names from TEtrimmer.*
- 📄**summary.txt** - *Summary file.*
- 📄**TEtrimmer_consensus.fasta** - *TE consensus library file before de-duplication.*
- 📄**TEtrimmer_consensus_merged.fasta** - *TE consensus library file after de-duplication.*

## Example report plots for each output TE consensus sequence

For each TEtrimmer output TE consensus sequence. You will get a report plot file like this:
![Reportplots](docs/TEtrimmer_report_plots_for_each_output.jpg)

## TEtrimmerGUI

You can use the TEtrimmerGUI tool to inspect and improve TEtrimmer generated TE consensus library.
This step is optional! TEtrimmer output can be used for genome-wide TE annotation directly.
But if you want to get a traditional manual-curation level TE consensus library, you should perform this step.

### Manual inspection of TEtrimmer outputs with TEtrimmerGUI

```bash
# Use --help to see all options
TEtrimmerGUI --help

# To start the manual inspection GUI tool
TEtrimmerGUI -i output_directory/TEtrimmer_for_proof_curation -g test_genome_file.fa.gz
```

<a href="https://www.youtube.com/watch?v=52GYZUQyzSE&t=1608s&ab_channel=ZhaoJiang">
  <img src="docs/TEtrimmer_GUI_work_space_video.png" alt="Proof_curation_GUI_work_page" style="width: 800px;">
</a>

### Curate existing TE libraries

The TEtrimmer GUI can also be used to check other TE consensus libraries like the TE library directly from EDTA,
RepeatModeler2, REPET, and other tools.

```bash
# Use --help to see all options
TEtrimmerGUI --help

# Open other consensus lib of TE predictions
TEtrimmerGUI -g <genome_file.fa> -clib <TE_consensus_library.fa>
```


#### Other TEtrimmerGUI Options

```code

Options:
  -i, --te_trimmer_proof_curation_dir TEXT
                                  TEtrimmer proof curation output path.Like <T
                                  Etrimmer_output_path>/TEtrimmer_for_proof_cu
                                  ration. Three folders should exist in the
                                  given path including "TE_clustered",
                                  "TE_low_copy", and "TE_skipped". If you
                                  start the "annoGUI.py" from TETrimmer output
                                  directory, you do not need to use this
                                  option
  -o, --output_dir TEXT           Output directory. Default: input directory
  -g, --genome_file TEXT          Genome fasta file path.
  -clib, --consensus_lib TEXT     TE consensus library FASTA file. You can
                                  check and improve other TE consensus library
                                  e.g. The TE library directly from EDTA2,
                                  RepeatModeler2, and other tools. If you want
                                  to check the same TE library as last time,
                                  you do not need to use this option again.
  -cdd, --cdd_dir TEXT            NCBI cdd database path.
  --max_msa_lines INTEGER         Set the maximum number of sequences to be
                                  included when click "Blast" button. Default:
                                  100
  --top_msa_lines INTEGER         If the sequence number after "Blast" is
                                  greater than <max_msa_lines>, TEtrimmerGUI
                                  will first sort sequences by length and
                                  choose <top_msa_lines> number of sequences.
                                  Then, TEtrimmerGUI will randomly select
                                  sequences from all remaining BLAST hits
                                  until <max_msa_lines>sequences are found.
                                  Default: 100
  --help                          Show this message and exit.

  ```


## All TEtrimmer available options
```commandline
Options:
  -i, --input_file TEXT           Path to TE consensus file (FASTA format). Use the output from
                                  RepeatModeler, EDTA, REPET, et al.  [required]

  -g, --genome_file TEXT          Path to genome FASTA file (FASTA format).  [required]

  -o, --output_dir TEXT           Path to output directory. Default: current working directory.

  -s, --preset [conserved|divergent]
                                  Choose one preset config (conserved or divergent).

  -t, --num_threads INTEGER       Thread number used for TEtrimmer. Default: 10

  --classify_unknown              Use RepeatClassifier to classify the consensus sequence if the input
                                  sequence is not classified or is unknown or the processed sequence
                                  length by TEtrimmer is 2000 bp longer or shorter than the query
                                  sequence.

  --classify_all                  Use RepeatClassifier to classify every consensus sequence. WARNING:
                                  This may take a long time.

  -ca, --continue_analysis        Continue from previous unfinished TEtrimmer run and would use the
                                  same output directory.
  --dedup                         Remove duplicate sequences in input file.

  -ga, --genome_anno              Perform genome TE annotation using RepeatMasker with the TEtrimmer
                                  curated TE libraries.

  --hmm                           Generate HMM files for each processed consensus sequence.

  --debug                         debug mode. This will keep all raw files. WARNING: Many files will be
                                  generated.

  -pd, --pfam_dir TEXT            Pfam database directory. TE Trimmer will download the database
                                  automatically. Only turn on this option if you want to use a local
                                  PFAM database or the automatic download fails.

  --cons_thr FLOAT                The minimum level of agreement required at a given position in the
                                  alignment for a consensus character to be called. Default: 0.8

  --mini_orf INTEGER              Define the minimum ORF length to be predicted by TEtrimmer. Default:
                                  200

  --max_msa_lines INTEGER         Set the maximum number of sequences to be included in a multiple
                                  sequence alignment. Default: 100

  --top_msa_lines INTEGER         If the sequence number of multiple sequence alignment (MSA) is
                                  greater than <max_msa_lines>, TEtrimmer will first sort sequences by
                                  length and choose <top_msa_lines> number of sequences. Then,
                                  TEtrimmer will randomly select sequences from all remaining BLAST
                                  hits until <max_msa_lines>sequences are found for the multiple
                                  sequence alignment. Default: 100

  --min_seq_num INTEGER           The minimum blast hit number required for the input sequence. We do
                                  not recommend decreasing this number. Default: 10

  --min_blast_len INTEGER         The minimum sequence length for blast hits to be included for further
                                  analysis. Default: 150

  --max_cluster_num INTEGER       The maximum number of clusters assigned in each multiple sequence
                                  alignment. Each multiple sequence alignment can be grouped into
                                  different clusters based on alignment patterns WARNING: using a
                                  larger number will potentially result in more accurate consensus
                                  results but will significantly increase the running time. We do not
                                  recommend increasing this value to over 5. Default: 2

  --ext_thr FLOAT                 The threshold to call “N” at a position. For example, if the most
                                  conserved nucleotide in a MSA columnhas proportion smaller than
                                  <ext_thr>, a “N” will be called at this position. Used with
                                  <ext_check_win>. The lower the value of <ext_thr>, the more likely to
                                  get longer the extensions on both ends. You can try reducing
                                  <ext_thr> if TEtrimmer fails to find full-length TEs. Default: 0.7

  --ext_check_win TEXT            the check windows size during defining start and end of the consensus
                                  sequence based on the multiple sequence alignment. Used with
                                  <ext_thr>. If <ext_check_win> bp at the end of multiple sequence
                                  alignment has “N” present (ie. positions have similarity proportion
                                  smaller than <ext_thr>), the extension will stop, which defines the
                                  edge of the consensus sequence. Default: 150

  --ext_step INTEGER              the number of nucleotides to be added to the left and right ends of
                                  the multiple sequence alignment in each extension step. TE_Trimmer
                                  will iteratively add <ext_step> nucleotides until finding the TE
                                  boundary or reaching <max_ext>. Default: 1000

  --max_ext INTEGER               The maximum extension in nucleotides at both ends of the multiple
                                  sequence alignment. Default: 7000

  --gap_thr FLOAT                 If a single column in the multiple sequence alignment has a gap
                                  proportion larger than <gap_thr> and the proportion of the most
                                  common nucleotide in this column is less than <gap_nul_thr>, this
                                  column will be removed from the consensus. Default: 0.4

  --gap_nul_thr FLOAT             The nucleotide proportion threshold for keeping the column of the
                                  multiple sequence alignment. Used with the <gap_thr> option. i.e. if
                                  this column has <40% gap and the portion of T (or any other)
                                  nucleotide is >70% in this particular column, this column will be
                                  kept. Default: 0.7

  --crop_end_div_thr FLOAT        The crop end by divergence function will convert each nucleotide in
                                  the multiple sequence alignment into a proportion value. This
                                  function will iteratively choose a sliding window from each end of
                                  each sequence of the MSA and sum up the proportion numbers in this
                                  window. The cropping will continue until the sum of proportions is
                                  larger than <--crop_end_div_thr>. Cropped nucleotides will be
                                  converted to -. Default: 0.7

  --crop_end_div_win INTEGER      Window size used for the end-cropping process. Used with the
                                  <--crop_end_div_thr> option. Default: 40

  --crop_end_gap_thr FLOAT        The crop end by gap function will iteratively choose a sliding window
                                  from each end of each sequence of the MSA and calculate the gap
                                  proportion in this window. The cropping will continue until the sum
                                  of gap proportions is smaller than <--crop_end_gap_thr>. Cropped
                                  nucleotides will be converted to -. Default: 0.1

  --crop_end_gap_win INTEGER      Define window size used to crop end by gap. Used with the
                                  <--crop_end_gap_thr> option. Default: 250

  --start_patterns TEXT           LTR elements always start with a conserved sequence pattern.
                                  TEtrimmer searches the beginning of the consensus sequence for these
                                  patterns. If the pattern is not found, TEtrimmer will extend the
                                  search of <--start_patterns> to up to 15 nucleotides from the
                                  beginning of the consensus sequence and redefine the start of the
                                  consensus sequence if the pattern is found. Note: The user can
                                  provide multiple LTR start patterns in a comma-separated list, like:
                                  TG,TA,TC (no spaces; the order of patterns determines the priority
                                  for the search). Default: TG

  --end_patterns TEXT             LTR elements always end with a conserved sequence pattern. TEtrimmer
                                  searches the end of the consensus sequence for these patterns. If the
                                  pattern is not found, TEtrimmer will extend the search of
                                  <--end_patterns> to up to 15 nucleotides from the end of the
                                  consensus sequence and redefine the end of the consensus sequence if
                                  the pattern is found. Note: The user can provide multiple LTR end
                                  patterns in a comma-separated list, like: CA,TA,GA (no spaces; the
                                  order of patterns determines the priority for the search). Default:
                                  CA

  --help                          Show this message and exit.
```

## Hardware requirements
System: Linux, macOS, Windows WSL (untested)

RAM:
- For HPC Linux user, enough RAM needs to be assigned. We highly recommend running TEtrimmer on HPC with at least 40 threads and assigning at least 5 GB RAM to each thread.


| Threads | RAM    |
|---------|--------|
| 40      | 200 GB |
| 100     | 600 GB |

- Windows and macOS PC users can use Virtual Memory. Simply assign 20 threads to push the CPU to its limits. We did tests on a Macbook Pro (2020 M1 chip, 16 GB RAM) and compared with HPC, you can find the running time here:

| Query sequence number | Platform       | Threads | RAM                    | Run time |
|-----------------------|----------------|---------|------------------------|--------------|
| 1700                  | Macbook Pro M1 | 20      | 16 GB + Virtual Memory | 50 hours     |
| 1700                  | HPC            | 40      | 150 GB                 | 5 hours      |


## Citation
Qian, J., Xue, H., Ou, S., Storer, J., Fürtauer, L., Wildermuth, M. C., Kusch, S., & Panstruga, R. bioRxiv (2024) https://doi.org/10.1101/2024.06.27.600963
TEtrimmer: A novel tool to automate the manual curation of transposable elements.
