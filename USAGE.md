##########################################################################################
TE Trimmer v1.1 (09/OCT/2023)
https://github.com/qjiangzhao/TE-Trimmer

Developers of TE Trimmer:
Jiangzhao Qian. Email: jqian@bio1.rwth-aachen.de
Hang Xue. Email: hang_xue@berkeley.edu

Funding source:
Panstruga's Lab. Website: https://www.bio1.rwth-aachen.de/PlantMolCellBiology/index.html
RWTH Aachen University
Many thanks to Dr. Stefan Kusch
##########################################################################################
python ./path_to_TE_Trimmer_bin/main.py -i <TE_consensus_file> -o <genome_file>

TE Trimmer is designed to replace transposable element (TE) manual curation. Two mandatory arguments are required including <genome file> and <TE consensus file> from TE annotation software like RepeatModeler, EDTA, and REPET et al. TE Trimmer can do blast, extension, multiple sequence alignment, and defining TE boundaries.

Options:
 -i, --input_file TEXT      TE consensus fasta file. Use the output of RepeatModeler, EDTA, or REPET et al.
                 [required]
 -g, --genome_file TEXT     Genome file path. [required]
 -o, --output_dir TEXT      Output directory. Default: current directory.
 -s, --species [fungi|plant|animal|powdery_mildew]
                 Select the species for which you want to run TE Trimmer. [required]
 --continue_analysis       Continue to analysis based on interrupted results.
 --merge             Merge input file to remove duplicate sequences.
 --genome_anno          Perform genome TE annotation based on TE Trimmer curated database at the end.
 --hmm              Generate HMM files for each consensus sequences.
 --keep_intermediate       Keep all raw files. WARNING: Many files will be produced.
 --fast_mode           Use less running time but lower accuracy and specificity
 --pfam_dir TEXT         Pfam database directory. Leave this option when you don't have Pfam database, TE
                 Trimmer will download automatically
 --cons_thr FLOAT        Threshold used for the final consensus sequence generation. Default: 0.8
 --max_msa_lines INTEGER     Set the maximum sequences number for multiple sequence alignment. Default: 100
 --top_mas_lines INTEGER     When the sequence number of multiple sequence alignment (MSA) is greater than
                 "max_mas_lines". It will order sequences by length and choose "top_msa_lines" number
                 of sequences. Then randomly choose the rest number of sequences Default: 100
 --min_seq_num INTEGER      The minimum sequence number for each multiple sequence alignment. Default: 10
 --min_blast_len INTEGER     The minimum hit sequence length for blast. Default: 150
 --max_cluster_num INTEGER    The maximum cluster number for each multiple sequence alignment. Each multiple
                 sequence alignment can be divided into different clusters. TE Trimmer will sort
                 cluster by sequence number and choosethe top --max_cluster_num of clusters for the
                 further analysis. Default: 2
 --ext_thr FLOAT         threshold used for define the extension extent. The smaller number means it become
                 easier to have a final longer extension for each side of the sequence. Default: 0.7
 --ex_step INTEGER        Number of nucleotides will be added to the left or right side of multiple sequence
                 alignment. TE_Trimmer will iteratively add --ex_step number of nucleotide until
                 finding the boundary. Default: 1000
 --max_extension INTEGER     The maximum extension number for the right and left side. For example, if --ex_step
                 is 1000, it can only add seven times to the MSA left or right side. Default: 7000
 --gap_thr FLOAT         If columns have a larger gap proportion than --gap_thr and the most common
                 nucleotide proportion in this column is less than --gap_nul_thr, this column will be
                 removed. Default: 0.4
 --gap_nul_thr FLOAT       Set nucleotide proportion to decide if remove this column. Coupled with
                 --gap_nul_thr option. Default: 0.7
 --crop_end_win INTEGER     Window size used for crop end process. Coupled with --crop_end_thr option. Default:
                 20
 --crop_end_thr INTEGER     Crop end function will convert each nucleotide in MSA into proportion number. This
                 function will check from the beginning and end of each sequence from MSA by
                 iteratively choosing a slide window and sum up the proportion numbers. It will stop
                 until the summary of proportion is larger than --crop_end_thr. The nucleotide do not
                 match --crop_end_thr will be converted to -. The recommended number is 0.8 *
                 --crop_end_win. Default: 16
 --crop_end_gap_win INTEGER   Define window size used to crop end by gap, coupled with --crop_end_gap_thr option.
                 Default: 150
 --crop_end_gap_thr FLOAT    Crop end by gap function will check from the beginning and end of each sequence from
                 MSA by iteratively choosing a slide window. It will stop until the gap proportion in
                 this slide window is smaller than --crop_end_gap_thr. The nucleotide do not match
                 the requirement will be converted to - Default: 0.1
 --start_patterns TEXT      LTR elements will always start with fixed pattern. TE Trimmer will check if it
                 starts with those patterns. If not, it will seek around the start point, if the
                 pattern is found, the start point will be converted to there. Note: if you want to
                 give multiple start patterns, separate them by comma. Like: TG,TA,TC (No space
                 between them). The order of the given patterns is matter. Default: TG
 --end_patterns TEXT       LTR elements will always end with fixed pattern. TE Trimmer will check if it end
                 with those patterns. If not, it will seek around the end point, if the pattern is
                 found, the start point will be converted to there. Note: if you want to give
                 multiple end patterns, separate them by comma. Like: CA,TA,GA (No space between
                 them). The order of the given patterns is matter. Default: AC
 --mini_orf INTEGER       Set the minimum ORF length that will be predicted by TE Trimmer. Default: 200
 --check_extension_win TEXT   Define check windows size for extension. Deafault: 150
 -t, --num_threads INTEGER    Threads numbers used for TE Trimmer. Default: 10
 --classify_unknown       Use RepeatClassfier to classify the consensus sequence if the input sequence is not
                 classfied or is unknown. Default: False
 --classify_all         Use RepeatClassfier to classify every consensus sequence. WARNING: it will take
                 longer. Default: False
 --help             Show this message and exit.
