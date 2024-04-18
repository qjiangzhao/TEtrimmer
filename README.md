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
curation is still required for high-quality TE annotation by experts. TEtrimmer is designed to replace and assist
manual TE curation. You can find more details about TEtrimmer below [flowchart](#Flowchart).

## Manual
For detailed instructions, including installation steps, usage options, example outputs, and more, 
please refer to [TEtrimmerv1.2.0Manual.pdf](https://github.com/qjiangzhao/TEtrimmer/blob/main/docs/TEtrimmerv1.2.0Manual.pdf) 

## Installation
Install TEtrimmer dependencies. 
We highly recommend installation with `mamba`, as it is much faster. 

```commandline
# Create TEtrimmer environment
# Install mamba 
conda install -c conda-forge mamba

# Clone the github repository for TEtrimmer.
git clone https://github.com/qjiangzhao/TEtrimmer.git

# Install TEtrimmer by the "yml" file
mambe env create -f <path to/TEtrimmer_env_for_linux.yml>
```
Here is the provided [TEtrimmer_env_for_linux.yml](https://github.com/qjiangzhao/TE-Trimmer/blob/main/TEtrimmer_env_for_linux.yml)



**or** if you are using macOS, you can install TEtrimmer conda package directly. Note: TEtrimmer installation requires python=3.10 
```commandline
conda create --name tetrimmerenv python=3.10
# Install TETrimmer 
mamba install qianjiangzhao::tetrimmer
# Display options of TETrimmer 
TETrimmer --help
```
TEtrimmer conda package only works for macOS currently. We are developing conda package for Linux. 
We are working on uploading the package to the Bioconda channel and dockerize it. 

**or** See required dependencies [TEtrimmer_dependencies](https://github.com/qjiangzhao/TEtrimmer/blob/main/docs/TEtrimmer_dependencies). 

## Usage:
Use --help to access all [options](#All-available-options)
```commandline
python {path to TEtrimmer}/TEtrimmer.py --help
```
**or** If using TEtrimmer conda package
```commandline
TEtrimmer --help
```

## Hardware requirements
System: Linux, macOS, Windows WSL

RAM:
- For HPC Linux user, enough RAM needs to be assigned. We highly recommend running TEtrimmer on HPC with at least 40 threads and assigning at least 5 GB ram to each thread.


| Threads | RAM    |
|---------|--------|
| 40      | 200 GB |
| 100     | 600 GB |

- PC macOS users can use Virtual Memory. Simply assign 20 threads to push the CPU to its limits. We did tests on a Macbook Pro (2020 M1 chip, 16 GB RAM) and compared with HPC, you can find the running time here:

| Query sequence number | Platform       | Threads | RAM                    | Run time |
|-----------------------|----------------|---------|------------------------|--------------| 
| 1700                  | Macbook Pro M1 | 20      | 16 GB + Virtual Memory | 50 hours     |
| 1700                  | HPC            | 40      | 150 GB                 | 5 hours      | 

- We have not tested it on the WLS of Windows, but it should be feasible to run TEtrimmer on it as well given sufficient resources. 

## Test

- Download the test files [test_input.fa](https://github.com/qjiangzhao/TEtrimmer/blob/main/tests/test_input.fa) and [test_genome.fasta](https://github.com/qjiangzhao/TEtrimmer/blob/main/tests/test_genome.fasta).

```commandline
TEtrimmer --input_file {path to test_input.fa} \
          --genome_file {path to test_genome.fasta} \
          --output_dir {output directory} \
          --num_threads 10
          --classify_unknown                                          
```
## Inputs
- **Genome file**: The genome sequence in FASTA format (.fa or .fasta).
- **TE consensus library**: TEtrimmer uses the TE consensus library from *de novo* TE annotation tools, like `RepeatModeler` or `EDTA`, as input. 
For this reason, you have to run `RepeatModeler` or other TE annotation software first. 

Example:

```commandline
TEtrimmer --input_file {TE consensus library} \
          --genome_file {genome file} \
          --output_dir {output directory} \
          --num_threads 10
                                          
```
If you want to **continue the analysis based on previous unfinished results in the same directory:**:
```commandline
TEtrimmer --input_file {TE consensus library} \
          --genome_file {genome file} \
          --output_dir {directory contains previous unfinished results} \
          --num_threads 10 \
          --continue_analysis
```
If you want to **combine files from different sources for the input file, we recommend removing duplicate sequences before processing. This step can potentially save overall run time in the input file**:
```commandline
TEtrimmer --input_file {TE consensus library} \
          --genome_file {genome file} \
          --output_dir {output directory} \
          --num_threads 10 \
          --dedup    
```
More options are available:
```commandline
  -s, --preset [conserved|divergent]
                                  Choose one preset config (conserved or divergent). Default: conserved.
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

```
## Outputs
- 📁**Classification** - *This folder is used for TE classifications.*  
- 📁**Multiple_sequence_alignment** - *All raw files will be stored in this folder if < --debug > is enabled.*
  - 📄**error_file.txt** - *Error file to store all error messages, only visible if errors were found.*
- 📁**Single_fasta_files** - *All sequences in the input file will be separated into single FASTA files and stored here.*
- 📁**TEtrimmer_for_proof_annotation** - *This folder contains files used for manual inspection of TEtrimmer annotations.* 
  - 📁**Annotation_perfect** - *Four files are associated with each sequence (anno.fa; fa; pdf).*
    - 📄**TE_name.anno.fa** - *Multiple sequence alignment file before cleaning.*
    - 📄**TE_name.fa** - *Multiple sequence alignment file after cleaning.*
    - 📄**TE_name.pdf** - *Plot file used to evaluate output.*
    - 📄**TE_name.bed** - *BED file used for further sequence elongation.*
  - 📁**Annotation_good** 
  - 📁**Annotation_check_recommended**
  - 📁**Annotation_check_required**
  - 📁**Clustered_proof_annotation** - *The folder group processed TEs to different clusters based on TE consensus sequence similarity.*
  - 📁**TE_low_copy** - *This folder contains low copy TEs.*
  - 📁**TE_skipped** - *Contains TE_Aid plots for all skipped TEs.*
  - 📁**TEtrimmer_proof_anno_GUI** - *The folder contains graphical user interface tools for manual proof annotation.*
    - 📄**annoGUI.py** - *Use python ./annoGUI.py to start manual proof annotation GUI.*
- 📁**HMM** - *This folder is used to store Hidden Markov Model file. Only visible when < --hmm > is enabled.*
- 📄**Sequence_name_mapping.txt** - *This file connects the input sequence names with the modified names from TEtrimmer.*
- 📄**summary.txt** - *Summary file.* 
- 📄**TEtrimmer_consensus.fasta** - *TE consensus library file before de-duplication.*
- 📄**TEtrimmer_consensus_merged.fasta** - *TE consensus library file after de-duplication.*


## Proof annotation: Manual inspection of TEtrimmer annotations
You can use this graphical user interface tool to assist the manual inspection of TEtrimmer-generated annotations. We highly recommend to perform
manual inspection of TE annotations in the "Recommend_check_annotation" and "Need_check_annotation" folders to generate a high-quality TE
```commandline
# To start the manual inspection GUI tool
# Open your Linux, macOS, or Windows terminal and type
python <path to your output_directory>/TEtrimmer_for_proof_annotation/TEtrimmer_proof_anno_GUI/annoGUI.py
```
The following are clusters and files (Click the "Clustered_proof_annotation" button in the menu bar to show this.)
![Proof annotation GUI](https://github.com/qjiangzhao/TEtrimmer/blob/main/docs/Proof_annotation.pdf)
                        
## Benchmarking
TEtrimmer is 6-times more accurate to annotate the intact TE than RepeatModeler in case of *Blumeria hordei*. 
![Benchmarking1](https://www.dropbox.com/scl/fi/v1ex6txe0mb9200gmtir3/Benchamrking_joined2.png?rlkey=i742b8ykyht0zw885r3mj9u64&raw=1)

## Acknowledgements

## Flowchart
![image](https://github.com/qjiangzhao/TEtrimmer/blob/main/docs/TEtrimmerFlowchart.pdf)

## All available options 
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

  --fast_mode                     Reduce running time at the cost of lower accuracy and specificity.

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