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
You can find required packages from [here](https://github.com/qjiangzhao/TE-Trimmer/blob/Jiangzhao_TE_trimmer/package-list.txt). 
To create the conda enviroment, do
```commandline
conda env create -f {path to TE Trimmer}/tetrimmerenv.yml
```
or 
```commandline
conda create --name tetrimmerenv --file {path to TE Trimmer}/package-list.txt
```
or 
```commandline
conda create -n tetrimmerenv -c conda-forge mamba
conda activate tetrimmerenv
mamba install -c conda-forge -c bioconda python perl r-base biopython matplotlib multiprocess pandas pypdf2 numpy scikit-learn urllib3 regex tk click requests bedtools blast cd-hit emboss hmmer mafft pfam_scan 'repeatmasker<4.1.5' 'rmblast<2.11' samtools dataclasses 'repeatmodeler=2.0.1' 'muscle=3.8'
```
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
| 100     | 600 GB |

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
If you want to **continue the analysis based on previous unfinished result**:
```commandline
python {path to TE Trimmer}/TE_Trimmer.py --input_file {TE consensus library} \
                                          --genome_file {genome file} \
                                          --output_dir {directory contains previous unfinished result} \
                                          --num_threads 10 \
                                          --continue_analysis
```
If you want to **remove duplicate sequences** in the input file:
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
- 📁**Classification** - This folder is used for TEs classification.  
- 📁**Multiple_sequence_alignment** - All raw files will be kept in this folder when < --debug > is enables.
- 📁**Single_fasta_files** - All sequences in the input file will be separated to single fasta files and be stored here.
- 📁**TE_Trimmer_for_proof_annotation** - This folder contains files used for proof annotation. 
  - 📁**Perfect_annotation** - For each sequence, three files are associate with it (anno.fa; fa; pdf)
    - 📄TE_name.anno.fa - Multiple sequence alignment file before cleaning.
    - 📄TE_name.fa - Multiple sequence alignment file after cleaning.
    - 📄TE_name.pdf - Plot file used to evaluate output.
  - 📁**Good_annotation** 
  - 📁**Recommend_check_annotation**
  - 📁**Need_check_annotation**
  - 📁**Low_copy_TE** - This folder contains low copy TEs.
  - 📁**Skipped_TE** - Contains TE_Aid plots for all skipped TEs. Only visible when < --plot_skip > is enabled.
- 📁**HMM** - This folder is used to store Hidden Markov Model file. Only visible when < --hmm > is enabled.
- 📄**Finished_sequence_recording.txt** - Report file. 
- 📄**TE_Trimmer_consensus.fasta** - TE consensus library file before de-duplication.
- 📄**TE_Trimmer_consensus_merged.fasta** - TE consensus library file after de-duplication.
- 📄**error_file.txt** - Error file to store all error messages. 

### Proof annotation
You can use this graphical user interface tool to assistant your proof annotation. We highly recommend to do proof 
annotation for the TEs in "Recommend_check_annotation" and "Need_check_annotation" folder to achieve high quality TE
consensus library. 
```commandline
# To start the proof annotation GUI tool
python {path to TE Trimmer}/Proof_annotation_GUI.py -i {path to TE_Trimmer_for_proof_annotation folder} \
                                                    -o {output directory}
```
You can follow the instruction to perform the proof annotation. 
![TE_Trimmer_interface1](https://www.dropbox.com/scl/fi/mynrf8mokblq9egslpsti/Screenshot-2023-10-29-at-12.19.27.png?rlkey=pozzit1llyteux2rhwxnxnn99&raw=1)
Those are all files in "Perfect annotation" folder (Click "Perfect annotation" button in the menu bar to show this.)
![TE_Trimmer_interfact2](https://www.dropbox.com/scl/fi/4nh0u7xvirieb68c5knnw/Screenshot-2023-10-29-at-12.20.14.png?rlkey=m2nfsevhriennsp5vf9s766zr&raw=1)

## Benchmarking
TE Trimmer is 6 time more accurate to annotate the intact TE than RepeatModeler for *B.hordei*. 
![Benchmarking1](https://www.dropbox.com/scl/fi/v1ex6txe0mb9200gmtir3/Benchamrking_joined2.png?rlkey=i742b8ykyht0zw885r3mj9u64&raw=1)

## Acknowledgements

## Flowchart
![image](https://www.dropbox.com/scl/fi/4s0sd2e0ndic62pyt22dt/TE_Trimmer_vertical_flowchart.png?rlkey=ixwbo1p7h05xhz80nh2j47y2o&raw=1)
