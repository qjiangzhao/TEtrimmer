# TE-Trimmer
TE Trimmer: a tool to replace transposable element manual curation.

Usage:

# See options and explanation
python ./path_to_TE_Trimmer_folder/bin/main.py --help 

# Example of running TE Trimmer
python ./path_to_TE_Trimmer_folder/bin/main.py --input_file [your_TE_consensus_file_path] \
                                               --genome_file [your_genome_file_path] \
                                               --output_dir [output_directory]
                                               
# Example of graphical user interface-based proof annotation
python ./path_to_TE_Trimmer_folder/bin/Class_TKinter_proof_annotation.py --te_trimmer_output_dir [your_TE_Trimmer_output_directory/Multiple_sequence_alignment_less] 

Options:
  -i, --input_file TEXT       Input fasta file. TE consensus sequences  [required]
  -o, --output_dir TEXT       Output directory. Default: current directory
  -g, --genome_file TEXT      Provide the genome file path  [required]
  --continue_analysis         If --continue_analysis is provided on the command line, TE Trimmer
                              will continue the analysis based on the existing data. Otherwise, it
                              will overlap the existing files
  --pfam_dir TEXT             Pfam database directory. Leave this option when you don't have Pfam
                              database, TE Trimmer will download automatically
  --max_msa_lines INTEGER     Set the maximum sequences number for multiple sequence alignment.
                              Default: 100
  --min_blast_len INTEGER     Blast sequence length lower than --min_blast_len will be ignored
  --top_mas_lines INTEGER     When the sequence number of multiple sequence alignment (MSA) is
                              greater than "max_mas_lines". It will order sequences by length and
                              choose "top_msa_lines" number of sequences. Then randomly choose the
                              rest number of sequences Default: 100
  --min_seq_num INTEGER       The minimum sequence number for each multiple sequence alignment.
                              Default: 10
  --max_cluster_num INTEGER   The maximum cluster number for each multiple sequence alignment.
                              Each multiple sequence alignment can be divided into different
                              groups. The group is called cluster
  --min_el INTEGER            If the query sequence length is smaller than the given number,
                              1500bp will be added to the left and right sides of this sequence.
                              Meanwhile, redundant elongation will be cleaned
  --min_el_dna INTEGER        Same like --min_el option. But set different given numbers for DNA
                              elements. Default: 500
  --min_el_sine INTEGER       Same like --min_el option. But set different given number for SINE
                              element. Default: 200
  --cons_thr FLOAT            Threshold used for consensus sequence generation for multiple
                              sequence alignment
  --ex_step INTEGER           Number of nucleotides will be added to the left or right side of
                              multiple sequence alignment. TE_Trimmer will iteratively add
                              --ex_step number of nucleotide until finding the boundary
  --max_extension INTEGER     The maximum extension number for the right and left side. For
                              example, if --ex_step is 1000, it can only add seven times to the
                              MSA left or right side
  --gap_thr FLOAT             If columns have a larger gap proportion than --gap_thr and the most
                              common nucleotide proportion in this column is less than
                              --gap_nul_thr, this column will be removed
  --gap_nul_thr FLOAT         Set nucleotide proportion to decide if remove this column. Coupled
                              with --gap_nul_thr option
  --crop_end_thr INTEGER      Crop end function will convert each nucleotide in MSA into
                              proportion number. This function will check from the beginning and
                              end of each sequence from MSA by iteratively choosing a slide window
                              and sum up the proportion numbers. It will stop until the summary of
                              proportion is larger than --crop_end_thr. The nucleotide do not
                              match --crop_end_thr will be converted to -
  --crop_end_win INTEGER      Window size used for crop end process. Coupled with --crop_end_thr
                              option
  --crop_end_gap_thr FLOAT    Crop end by gap function will check from the beginning and end of
                              each sequence from MSA by iteratively choosing a slide window. I
                              will stop until the gap proportion in this slide window is smaller
                              than --crop_end_gap_thr. The nucleotide do not match the requirement
                              will be converted to -
  --crop_end_gap_win INTEGER  Define window size used to crop end by gap, coupled with
                              --crop_end_gap_thr option
  --start_patterns TEXT       LTR elements will always start with fixed pattern. TE Trimmer will
                              check if it starts with those patterns. If not, it will seek around
                              the start point, if the pattern is found, the start point will be
                              converted to there. Note: if you want to give multiple start
                              patterns, separate them by comma. Like: TG,TA,TC (No space between
                              them). The order of the given patterns is important
  --end_patterns TEXT         LTR elements will always end with fixed pattern. TE Trimmer will
                              check if it end with those patterns. If not, it will seek around the
                              end point, if the pattern is found, the start point will be
                              converted to there. Note: if you want to give multiple end patterns,
                              separate them by comma. Like: CA,TA,GA (No space between them). The
                              order of the given patterns is important
  --mini_orf INTEGER          Set the minimum ORF length that will be predicted by TE Trimmer
  -t, --num_threads INTEGER   Threads numbers used for TE Trimmer. Default: 10
  --help                      Show this message and exit.
