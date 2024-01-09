## Contents
- [Introduction](#Introduction)
- [Installation](#Installation)
- [Usage](#Usage)
  - [Test](#Test)
  - [Hardware requirements](#Hardware-requirements)
  - [Inputs](#Inputs) 
  - [Outputs](#Outputs)
  - [Proof annotation](#Proof-annotation)
- [Flowchart](#Flowchart)
- [All available options](#All-available-options)

## Introduction
Many tools have been developed for *de novo* transposable element (TE) identification. However, manual 
curation is still required for high-quality TE annotation by experts. TETrimmer is designed to replace and assist
manual TE curation. You can find more details about TETrimmer below [flowchart](#Flowchart).

## Installation
You can find the required packages [here](https://github.com/qjiangzhao/TE-Trimmer/blob/main/TE_Trimmer_dependencies) 
and install them by yourself.


**or** create `conda` environment in ***Linux*** based on [TETrimmer_env_for_linux.yml](https://github.com/qjiangzhao/TE-Trimmer/blob/main/TETrimmer_env_for_linux.yml).
```commandline
conda env create -f TETrimmer_env_for_linux.yml
```
**or** install using `mamba` based on [TETrimmer_env_for_linux.yml](https://github.com/qjiangzhao/TE-Trimmer/blob/main/TETrimmer_env_for_linux.yml). 
We highly recommend installation with `mamba`, as it is much faster. 

```commandline
# Install mamba first
conda install -c conda-forge mamba

# Create TETrimmer environment with mamba
mamba env create -f TETrimmer_env_for_linux.yml
```
For ***Windows WSL***, you can follow the same instructions used for Linux. 


For ***macOS***, use the following instructions: [TETrimmer_env_for_macOS.yml](https://github.com/qjiangzhao/TE-Trimmer/blob/main/TETrimmer_env_for_macOS.yml).

**Currently, the macOS yml file is unusable. Please install according to [dependency file](https://github.com/qjiangzhao/TE-Trimmer/blob/main/TETrimmer_dependencies)**

We will develop `conda` and `Docker` packages for `TETrimmer`.

## Usage:
Use --help to access all [options](#All-available-options)

```commandline
python {path to TETrimmer}/TETrimmer.py --help
```
## Hardware requirements
System: Linux, macOS

RAM:
- For HPC Linux users, enough RAM needs to be assigned. We highly recommend to run it on HPC using at least 40 threads and 150 GB RAM.

| Threads | RAM    |
|---------|--------|
| 40      | 150 GB |
| 100     | 400 GB |

- PC macOS users can use Virtual Memory. Simply assign 20 threads to push the CPU to its limits. We did
tests on a Macbook Pro (2020 M1 chip, 16 GB RAM) and compared with HPC, you can find the running time here:

| Query sequence number | Platform       | Threads | RAM                    | Run     time |
|-----------------------|----------------|---------|------------------------|--------------| 
| 1700                  | Macbook Pro M1 | 20      | 16 GB + Virtual Memory | 60 hours     |
| 1700                  | HPC            | 40      | 150 GB                 | 7 hours      | 

- We have not tested it on the WLS of Windows, but it should be feasible to run TETrimmer on it as well given sufficient resources. 

## Test
```commandline
# Perform unit test. 
python {path to TETrimmer}/test.py
```
or
```commandline
# The {output directory} must be empty. Unit test was not be performed. 
python {path to TETrimmer}/TETrimmer.py --input_file {path to TETrimmer}/tests/test_input.fa \
                                        --genome_file {path to TETrimmer}/tests/test_genome.fasta \
                                        --output_dir {output directory} \
                                        --num_threads 10
                                        --species fungi
                                        --classify_unknown
                                          
```
## Inputs
 
- **TE consensus library**: TETrimmer uses the TE consensus library from *de novo* TE annotation tools, like `RepeatModeler` or `EDTA`, as input. 
For this reason, you have to run `RepeatModeler` or other TE annotation software first. 
- **Genome file**: The genome sequence in FASTA format (.fa or .fasta).

Example:

```commandline
# The {output directory} must be empty.
python {path to TETrimmer}/TETrimmer.py --input_file {TE consensus library} \
                                        --genome_file {genome file} \
                                        --output_dir {output directory} \
                                        --num_threads 10
                                          
```
If you want to **continue the analysis based on results from a previous unfinished run**:
```commandline
python {path to TETrimmer}/TETrimmer.py --input_file {TE consensus library} \
                                        --genome_file {genome file} \
                                        --output_dir {directory contains previous unfinished results} \
                                        --num_threads 10 \
                                        --continue_analysis
```
If you want to **remove duplicate sequences** in the input file:
```commandline
python {path to TETrimmer}/TETrimmer.py --input_file {TE consensus library} \
                                        --genome_file {genome file} \
                                        --output_dir {output directory} \
                                        --num_threads 10 \
                                        --dedup    
```
More options are available:
```commandline
  --genome_anno                   Perform genome TE annotation using the TETrimmer curated database. Requires RepeatMasker.
  --hmm                           Generate HMM files for each consensus sequences.
  --pfam_dir TEXT                 PFAM database directory. Omit if you do not have a local PFAM database. TETrimmer will download the database automatically.
  --cons_thr FLOAT                Threshold used for the final consensus sequence generation. Default: 0.8
  --mini_orf INTEGER              Define the minimum ORF length that will be predicted by TETrimmer. Default: 200
  --classify_unknown              Use RepeatClassifier to classify the consensus sequence if the input sequence is not
                                  classified or is unknown.
  --classify_all                  Use RepeatClassifier to classify every consensus sequence.  WARNING: This may take a long time.

```
## Outputs
- 📁**Classification** - *This folder is used for TE classifications.*  
- 📁**Multiple_sequence_alignment** - *All raw files will be stored in this folder if < --debug > is enabled.*
  - 📄**error_file.txt** - *Error file to store all error messages, only visible if errors were found.*
- 📁**Single_fasta_files** - *All sequences in the input file will be separated into single FASTA files and be stored here.*
- 📁**TETrimmer_for_proof_annotation** - *This folder contains files used for manual inspection of TETrimmer annotations.* 
  - 📁**Perfect_annotation** - *Three files are associate with each sequence (anno.fa; fa; pdf).*
    - 📄**TE_name.anno.fa** - *Multiple sequence alignment file before cleaning.*
    - 📄**TE_name.fa** - *Multiple sequence alignment file after cleaning.*
    - 📄**TE_name.pdf** - *Plot file used to evaluate output.*
    - 📄**TE_name.bed** - *BED file used for further sequence elongation.*
  - 📁**Good_annotation** 
  - 📁**Recommend_check_annotation**
  - 📁**Need_check_annotation**
  - 📁**Low_copy_TE** - *This folder contains low copy TEs.*
  - 📁**Skipped_TE** - *Contains TE_Aid plots for all skipped TEs.*
- 📁**HMM** - *This folder is used to store Hidden Markov Model file. Only visible when < --hmm > is enabled.*
- 📄**summary.txt** - *Report file.* 
- 📄**TETrimmer_consensus.fasta** - *TE consensus library file before de-duplication.*
- 📄**TETrimmer_consensus_merged.fasta** - *TE consensus library file after de-duplication.*
- 📄**Sequence_name_mapping.txt** - *This file connects the input sequence names with the modified names from TETrimmer.*


## Proof annotation: Manual inspection of TETrimmer annotations
You can use this graphical user interface tool to assist the manual inspection of TETrimmer-generated annotations. We highly recommend to perform
manual inspection of TE annotations in the "Recommend_check_annotation" and "Need_check_annotation" folders to generate a high-quality TE
consensus library. 
```commandline
# To start the manual inspection GUI tool
python {path to TETrimmer}/annoGUI.py -i {path to TETrimmer_for_proof_annotation folder} \
                                      -o {output directory}
```
You can follow these instructions to perform the inspection. 
![TETrimmer_interface1](https://www.dropbox.com/scl/fi/mynrf8mokblq9egslpsti/Screenshot-2023-10-29-at-12.19.27.png?rlkey=pozzit1llyteux2rhwxnxnn99&raw=1)
The following are files deposited in the "Perfect annotation" folder (Click the "Perfect annotation" button in the menu bar to show this.)
![TETrimmer_interfact2](https://www.dropbox.com/scl/fi/4nh0u7xvirieb68c5knnw/Screenshot-2023-10-29-at-12.20.14.png?rlkey=m2nfsevhriennsp5vf9s766zr&raw=1)

## Benchmarking
TETrimmer is 6-times more accurate to annotate the intact TE than RepeatModeler in case of *Blumeria hordei*. 
![Benchmarking1](https://www.dropbox.com/scl/fi/v1ex6txe0mb9200gmtir3/Benchamrking_joined2.png?rlkey=i742b8ykyht0zw885r3mj9u64&raw=1)

## Acknowledgements

## Flowchart
![image](https://www.dropbox.com/scl/fi/4s0sd2e0ndic62pyt22dt/TE_Trimmer_vertical_flowchart.png?rlkey=ixwbo1p7h05xhz80nh2j47y2o&raw=1)

## All available options 
```commandline
Options:
  -i, --input_file TEXT           Path to TE consensus file (FASTA format). Use the output from RepeatModeler or EDTA
                                  et al.  [required]
                                  
  -g, --genome_file TEXT          Path to genome FASTA file.  [required]
  
  -o, --output_dir TEXT           Output directory. Default: current directory.
  
  -s, --species [fungi|plant|animal|powdery_mildew]
                                  Select the species for which you want to run TETrimmer.
                                  
  -ca, --continue_analysis        Continue to analysis after interruption.
  
  --dedup                         Remove duplicate sequences in input file.
  
  --genome_anno                   Perform genome TE annotation using the TETrimmer curated database. Requires
                                  RepeatMasker.
                                  
  --hmm                           Generate HMM files for each consensus sequences.
  
  --debug                         Open debug mode. This will keep all raw files. WARNING: Many files will be produced.
  
  --fast_mode                     Reduce running time but at the cost of lower accuracy and specificity.
    
  --pfam_dir TEXT                 Pfam database directory. Omit this option if you do not have a local PFAM database.
                                  TETrimmer will download the database automatically in this case.
                                  
  --cons_thr FLOAT                Threshold used for the final consensus sequence generation. Default: 0.8
  
  --mini_orf INTEGER              Define the minimum ORF length that will be predicted by TETrimmer. Default: 200
  
  --max_msa_lines INTEGER         Set the maximum sequences number for multiple sequence alignment. Default: 100
  
  --top_mas_lines INTEGER         When the sequence number of multiple sequence alignment (MSA) is greater than
                                  <--max_msa_lines>, TETrimmer will sort sequences by length and choose
                                  <--top_msa_lines> number of sequences. Then, TETrimmer will randomly select
                                  sequences from all remaining BLAST hits until <--max_msa_lines> sequences are found
                                  for the multiple sequence alignment. Default: 70
                                  
  --min_seq_num INTEGER           The minimum sequence number for each multiple sequence alignment. Note: can not
                                  smaller than 10. Default: 10
                                  
  --min_blast_len INTEGER         The minimum sequence length for blast hits. Default: 150
  
  --max_cluster_num INTEGER       The maximum cluster number for each multiple sequence alignment. Each multiple
                                  sequence alignment can be divided into different clusters. TETrimmer will sort
                                  clusters by sequence number and choose the top <--max_cluster_num> of clusters for
                                  further analysis. WARNING: Big number will dramatically increase running time.
                                  Default: 2
                                  
  --ext_thr FLOAT                 Threshold used for define the extension extent. The lower the value of <--ext_thr>,
                                  the easier the extensions on both ends be longer. Reduce <--ext_thr> if TETrimmer
                                  fails to determine the correct ends of repeat elements. Default: 0.7
                                  
  --ext_check_win TEXT            Define check windows size for extension. Default: 150
  
  --ext_step INTEGER              Number of nucleotides to be added to the left and right ends of the multiple
                                  sequence alignment. TETrimmer will iteratively add <--ext_step> number of
                                  nucleotides until finding the boundary. Default: 1000
                                  
  --max_ext INTEGER               The maximum extension in nucleotides at both ends of the multiple sequence
                                  alignment. Default: 7000
                                  
  --gap_thr FLOAT                 If multiple sequence alignment positions (columns) have a gap proportion larger than
                                  <--gap_thr> and the proportion of the most common nucleotide in this column is less
                                  than <--gap_nul_thr>, this column will be removed from the consensus. Default: 0.4
                                  
  --gap_nul_thr FLOAT             Set nucleotide proportion threshold for keeping the column of the multiple sequence
                                  alignment. Used with the <--gap_thr> option. Default: 0.7
                                  
  --crop_end_div_thr FLOAT        The crop end by divergence function will convert each nucleotide in the multiple
                                  sequence alignment into a proportion value. This function will iteratively choose a
                                  sliding window from each end of each sequence of the MSA and sum up the proportion
                                  numbers in this window. The cropping will continue until the sum of proportions is
                                  larger than <--crop_end_div_thr>. Cropped nucleotides will be converted to -.
                                  Default: 0.8
                                  
  --crop_end_div_win INTEGER      Window size used for the end-cropping process. Used with --crop_end_div_thr option.
                                  Default: 20
                                  
  --crop_end_gap_thr FLOAT        The crop end by gap function will iteratively choose a sliding window from each end
                                  of each sequence of the MSA and calculate the gap proportion in this window. The
                                  cropping will continue until the sum of gap proportions is smaller than
                                  <--crop_end_gap_thr>. Cropped nucleotides will be converted to -. Default: 0.1
                                  
  --crop_end_gap_win INTEGER      Define window size used to crop end by gap, used with <--crop_end_gap_thr> option.
                                  Default: 250
                                  
  --start_patterns TEXT           LTR elements always start with a conserved sequence pattern. TETrimmer searches the
                                  beginning of the consensus sequence for these patterns. If the pattern is not found,
                                  it will extend the search of <--start_patterns> to up to 15 nucleotides from the
                                  beginning of the consensus sequence and redefine the start of the consensus sequence
                                  if the pattern is found. Note: The user can provide multiple LTR start patterns in a
                                  comma-separated list, like: TG,TA,TC (no spaces; the order of patterns determines
                                  the priority for the search). Default: TG
                                  
  --end_patterns TEXT             LTR elements always end with a conserved sequence pattern. TETrimmer searches the
                                  end of the consensus sequence for these patterns. If the pattern is not found, it
                                  will extend the search of <--end_patterns> to up to 15 nucleotides from the end of
                                  the consensus sequence and redefine the end of the consensus sequence if the pattern
                                  is found. Note: The user can provide multiple LTR end patterns in a comma-separated
                                  list, like: CA,TA,GA (no spaces; the order of patterns determines the priority for
                                  the search). Default: CA
                                  
  -t, --num_threads INTEGER       Threads numbers used for TETrimmer. Default: 10
  
  --classify_unknown              Use RepeatClassifier to classify the consensus sequence if the input sequence is not
                                  classified or is unknown.
                                  
  --classify_all                  Use RepeatClassifier to classify every consensus sequence.  WARNING: it will take
                                  longer time.
                                  
  --help                          Show this message and exit.
```
