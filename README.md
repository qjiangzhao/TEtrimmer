## Contents
- [Introduction](#Introduction)
- [Installation](#Installation)
- [Usage](#Usage)
  - [Test](#Test)
  - [Hardware requirement](#Hardware-requirement)
  - [Inputs](#Inputs) 
  - [Outputs](#Outputs)
  - [Proof annotation](#Proof-annotation)
- [Flowchart](#Flowchart)

## Introduction
Many tools have been developed for *de novo* transposable element (TE) identification. But manual 
curation is still required for high quality TE annotation by experts. TE Trimmer is designed to replace and assistant
TE manual curation. You can find more details for TE Trimmer [flowchart](#Flowchart).


## Installation
You can find required dependencies from [here](https://github.com/qjiangzhao/TE-Trimmer/blob/Jiangzhao_TE_trimmer/TE_Trimmer_dependencies). 
We will develop Conda and Docker packages for TE Trimmer.

## Usage:
Use --help to access all options

```commandline
python {path to TE Trimmer}/TE_Trimmer.py --help
```
### Hardware requirement
System: Linux, macOS

RAM:
- For HPC Linux user, enough RAM has to be assigned. We highly recommend to run it on HPC with at least 40 threads.

| Threads | RAM    |
|---------|--------|
| 40      | 250 GB |
| 100     | 500 GB |

- For PC macOS user, because Virtual Memory can be used. You can simply use 20 threads to push the CPU to its limits. We did
test on Macbook Pro (2020 M1 chip 16 GB) and compared with HPC, you can find the running time here:

| Query sequence number | Platform       | Threads | RAM                    | Running time |
|-----------------------|----------------|---------|------------------------|--------------| 
| 1700                  | Macbook Pro M1 | 20      | 16 GB + Virtual Memory | 60 hours     |
| 1700                  | HPC            | 40      | 250 GB                 | 7 hours      | 

- We haven't tested it on WLS of Windows, it should be feasible to run TE Trimmer on it too. 

### Test
Test file is still not available.

### Inputs
 
- **TE consensus library**: TE Trimmer use the TE consensus library from *de novo* TE annotation tools like RepeatModeler or EDTA as input. 
For this reason, you have to run RepeatModeler or other TE annotation software first. 
- **Genome file**

Example:

```commandline
# The {output directory} must be empty.
python {path to TE Trimmer}/TE_Trimmer.py --input_file {TE consensus library} \
                                          --genome_file {genome file} \
                                          --output_dir {output directory} \
                                          --num_threads 10
                                          
```
If you want to continue the analysis based on previous unfinished result:
```commandline
python {path to TE Trimmer}/TE_Trimmer.py --input_file {TE consensus library} \
                                          --genome_file {genome file} \
                                          --output_dir {directory contains previous unfinished result} \
                                          --num_threads 10 \
                                          --continue_analysis
```
If you want to remove duplicate sequences in the input file:
```commandline
python {path to TE Trimmer}/TE_Trimmer.py --input_file {TE consensus library} \
                                          --genome_file {genome file} \
                                          --output_dir {output directory} \
                                          --num_threads 10 \
                                          --dedup    
```
More options are available:
```commandline
  --genome_anno                   Perform genome TE annotation using the TE Trimmer curated database. Requires RepeatMasker.
  --hmm                           Generate HMM files for each consensus sequences.
  --fast_mode                     Reduce running time but at the cost of lower accuracy and specificity.
  --plot_skip                     Perform TE_Aid plot for skipped elements
  --pfam_dir TEXT                 Pfam database directory. Omit this if you do not have a local PFAM database. TE Trimmer will download the database automatically.
  --cons_thr FLOAT                Threshold used for the final consensus sequence generation. Default: 0.8
  --mini_orf INTEGER              Define the minimum ORF length that will be predicted by TE Trimmer. Default: 200
  --classify_unknown              Use RepeatClassifier to classify the consensus sequence if the input sequence is not
                                  classified or is unknown.
  --classify_all                  Use RepeatClassifier to classify every consensus sequence.  WARNING: it will take
                                  longer time.
```
### Outputs
- 📁**Classification** - Folder used for TEs classification  
- 📁**Multiple_sequence_alignment** - This folder contain all the raw files. < --debug > will keep all raw files here.
- 📁**Single_fasta_files** - All sequences in the input file will be separated to single fasta file and be stored here.
- 📁**TE_Trimmer_for_proof_annotation** - Folder contains files used for proof annotation. 
  - 📁**Perfect_annotation** - For each sequence, three files are corresponded with it (anno.fa; fa; pdf)
  - - 📄TE_name.anno.fa - Multiple sequence alignment file before cleaning
  - - 📄TE_name.fa - Multiple sequence alignment file after cleaning.
  - - 📄TE_name.pdf - Plot file used to evaluate output
  - 📁**Good_annotation** 
  - 📁**Recommend_check_annotation**
  - 📁**Need_check_annotation**
  - 📁**Low_copy_TE** - Contains low copy TEs
  - 📁**Skipped_TE** - Contains TE_Aid plots for all skipped TEs. Only available when < --plot_skip > in enabled.
- 📁**HMM** - Folder used to store Hidden Markov Model file. Only visible when < --hmm > is enabled.
- 📄**Finished_sequence_recording.txt** - Report file 
- 📄**TE_Trimmer_consensus.fasta** - Consensus library file before de-duplication
- 📄**TE_Trimmer_consensus_merged.fasta** - Consensus library file after de-duplication
- 📄**error_file.txt** - Error file to store all error messages. 

### Proof annotation

## Benchmark

## Acknowledgements

## Flowchart
![image](https://private-user-images.githubusercontent.com/75024559/278702015-5f01f336-f9bd-4827-89e4-7b1b055b7e3f.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTEiLCJleHAiOjE2OTg0MjQxMjcsIm5iZiI6MTY5ODQyMzgyNywicGF0aCI6Ii83NTAyNDU1OS8yNzg3MDIwMTUtNWYwMWYzMzYtZjliZC00ODI3LTg5ZTQtN2IxYjA1NWI3ZTNmLnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFJV05KWUFYNENTVkVINTNBJTJGMjAyMzEwMjclMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjMxMDI3VDE2MjM0N1omWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTliZWY1ODkzOWU0Yzc4NDA1ZjY1MDY0NDYzMGJmZDMzN2IwYmZiMWQwYWMwNWFlZTQ0NGNhNmFmMGJhMWIxMWEmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0JmFjdG9yX2lkPTAma2V5X2lkPTAmcmVwb19pZD0wIn0.pdVoWAxEymwxSfrFfZwlIw1mDvD7GQ5sgjL2Ya-xGWs)



