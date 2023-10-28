## Contents
- [Introduction](#Introduction)
- [Installation](#Installation)
- [Usage](#Usage)
- - [Test](#Test)
- - [Inputs](#Inputs)
- - [Outputs](#Outputs)
- - [Proof annotation](#Proof-annotation)
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

### Test
Test file is still not available.

### Inputs
 
- TE consensus library: TE Trimmer use the TE consensus library from *de novo* TE annotation tools like RepeatModeler or EDTA as input. 
For this reason, you have to run RepeatModeler or other TE annotation software first. 
- Genome file

Example:

```commandline
python {path to TE Trimmer}/TE_Trimmer.py --input_file {TE consensus library} \
                                          --genome_file {genome file} \
                                          --output_dir {output directory}
```
If you want to continue analysis based on previous unfinished result:
```commandline
python {path to TE Trimmer}/TE_Trimmer.py --input_file {TE consensus library} \
                                          --genome_file {genome file} \
                                          --output_dir {directory contains previous unfinished result} \
                                          -ca
```
If you want to remove duplicate sequences in input file (cd-hit-est is used):
```commandline
python {path to TE Trimmer}/TE_Trimmer.py --input_file {TE consensus library} \
                                          --genome_file {genome file} \
                                          --output_dir {output directory} \
                                          --dedup
```
More options are available:
```commandline
  --genome_anno                   Perform genome TE annotation using the TE Trimmer curated database. Requires RepeatMasker.
  --hmm                           Generate HMM files for each consensus sequences.
  --debug                         Open debug mode. This will keep all raw files. WARNING: Many files will be produced.
  --fast_mode                     Reduce running time but at the cost of lower accuracy and specificity.
  --plot_skip                     Perform TE_Aid plot for skipped elements
  --pfam_dir TEXT                 Pfam database directory. Omit this if you do not have a local PFAM database. TE Trimmer will download the database automatically.
  --cons_thr FLOAT                Threshold used for the final consensus sequence generation. Default: 0.8
  --mini_orf INTEGER              Define the minimum ORF length that will be predicted by TE Trimmer. Default: 200
  --classify_unknown              Use RepeatClassifier to classify the consensus sequence if the input sequence is not
                                  classified or is unknown.
  --classify_all                  Use RepeatClassifier to classify every consensus sequence.  WARNING: it will take
                                  longer time.
  -t, --num_threads INTEGER       Threads numbers used for TE Trimmer. Default: 10
```
### Outputs




### Proof annotation

## Benchmark

## Acknowledgements

## Flowchart
![image](https://private-user-images.githubusercontent.com/75024559/278702015-5f01f336-f9bd-4827-89e4-7b1b055b7e3f.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTEiLCJleHAiOjE2OTg0MjQxMjcsIm5iZiI6MTY5ODQyMzgyNywicGF0aCI6Ii83NTAyNDU1OS8yNzg3MDIwMTUtNWYwMWYzMzYtZjliZC00ODI3LTg5ZTQtN2IxYjA1NWI3ZTNmLnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFJV05KWUFYNENTVkVINTNBJTJGMjAyMzEwMjclMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjMxMDI3VDE2MjM0N1omWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTliZWY1ODkzOWU0Yzc4NDA1ZjY1MDY0NDYzMGJmZDMzN2IwYmZiMWQwYWMwNWFlZTQ0NGNhNmFmMGJhMWIxMWEmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0JmFjdG9yX2lkPTAma2V5X2lkPTAmcmVwb19pZD0wIn0.pdVoWAxEymwxSfrFfZwlIw1mDvD7GQ5sgjL2Ya-xGWs)



