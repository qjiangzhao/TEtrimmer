##########################################################################################
TE Trimmer v1.1 (09/OCT/2023)
https://github.com/qjiangzhao/TE-Trimmer

Developers of TE Trimmer:
Jiangzhao Qian. Email: jqian@bio1.rwth-aachen.de
Hang Xue. Email: hang_xue@berkeley.edu

Funding source:
Panstruga Lab. Website: https://www.bio1.rwth-aachen.de/PlantMolCellBiology/index.html
RWTH Aachen University

Many thanks to Dr. Stefan Kusch
##########################################################################################

General usage:
```ShellSession
python ./path_to_TE_Trimmer_bin/TE_Trimmer.py -i <TE_consensus_file> -o <genome_file>
```

TE Trimmer is designed to replace transposable element (TE) manual curation. Two mandatory arguments are required, which are <genome file> and <TE consensus file> from TE annotation software like RepeatModeler, EDTA, and REPET et al. TE Trimmer can perform BLAST searches, extend matching sequence regions, multiple sequence alignment, and define TE boundaries.

Options:
```
 -i, --input_file STRING      Path to TE consensus file (FASTA format). Use the output from RepeatModeler, EDTA, or REPET et al.
                 [required]
 -g, --genome_file STRING     Path to genome FASTA file. [required]
 -o, --output_dir STRING      Output directory. Default: current working directory.
 -s, --species [fungi|plant|animal|powdery_mildew]
                 Select the type of organism for which you want to run TE Trimmer. [required]
 --continue_analysis       Continue analysis after interruption.
 --dedup             Remove duplicate sequences in input file.
 --genome_anno          Perform genome TE annotation using the TE Trimmer curated database. Requires RepeatMasker.
 --hmm              Generate HMM files for each consensus sequence.
 --keep_intermediate       Keep all raw files. WARNING: Many files will be produced.
 --fast_mode           Reduce running time but at the cost of lower accuracy and specificity.
 --pfam_dir STRING         PFAM database directory. Omit this option if you do not have a local PFAM database - TE
                 Trimmer will download the database automatically in this case.
 --cons_thr FLOAT        Threshold used for the final consensus sequence generation. Default: 0.8
 --max_msa_lines INTEGER     Set the maximum sequence number for multiple sequence alignment. Default: 100
 --top_msa_lines INTEGER     When the sequence number of multiple sequence alignment (MSA) is greater than
                 <max_msa_lines>, TE Trimmer will sort sequences by length and choose <top_msa_lines> number
                 of sequences. Then, TE Trimmer will randomly select sequences from all remaining BLAST hits until
                 <max_msa_lines> sequences are found for the multiple sequence alignment. Default: 100
 --min_seq_num INTEGER      The minimum sequence number for each multiple sequence alignment. Default: 10
 --min_blast_len INTEGER     The minimum sequence length for BLAST hits. Default: 150
 --max_cluster_num INTEGER    The maximum cluster number for each multiple sequence alignment. Each multiple
                 sequence alignment can be divided into different clusters. TE Trimmer will sort
                 clusters by sequence number and choose the top <max_cluster_num> of clusters for
                 further analysis. Default: 2
 --ext_thr FLOAT            Sequence similarity threshold used for defining start and end of the consensus sequence,
                 based on the multiple sequence alignment. Nucleotides in each position with a similarity proportion
                 smaller than <ext_thr> will be assigned the value N. If no N values are found, the algorithm will
                 extend the multiple sequence alignment to determine the limits of the consensus sequence. The lower
                 the value of <ext_thr>, the longer the extensions on both ends. Reduce <ext_thr> if TE Trimmer fails
                 to determine the correct ends of repeat elements. Default: 0.7
 --ext_step INTEGER        Number of nucleotides to be added to the left and right ends of the multiple sequence
                 alignment. TE_Trimmer will iteratively add <ex_step number> of nucleotides until
                 finding the boundary. Default: 1000
 --max_extension INTEGER     The maximum extension in nucleotides at both ends of the multiple sequence alignment. 
                 Default: 7000
 --gap_thr FLOAT         If multiple sequence alignment positions (columns) have a gap proportion larger than <gap_thr>
                 and the proportion of the most common nucleotide in this column is less than <gap_nul_thr>, this column
                 will be removed from the consensus. Default: 0.4
 --gap_nul_thr FLOAT       Set nucleotide proportion threshold for keeping the column of the multiple sequence alignment. 
                 Used with the --gap_nul_thr option. Default: 0.7
 --crop_end_win INTEGER     Window size used for the end-cropping process. Used with the --crop_end_thr option. Default:
                 20
 --crop_end_thr FLOAT     The crop end function will convert each nucleotide in the multiple sequence alignment into a
                 proportion value. This function will iteratively choose a slide window and sum up the proportion numbers 
                 at each end of each sequence of the MSA. The algorithm will continue until the sum of proportions at the
                 respective end of the MSA is larger than <crop_end_thr>. Nucleotides in the consensus sequence below
                 <crop_end_thr> will be converted to -. The recommended value is 0.8 * <crop_end_win>. Default: 16
 --crop_end_gap_win INTEGER   Define window size used to crop end by gap. Used with the --crop_end_gap_thr option.
                 Default: 150
 --crop_end_gap_thr FLOAT    The crop end by gap function iteratively determines the gap proportion in sliding windows at
                 both ends of the multiple sequence alignment. The function continues until the gap proportion in
                 the sliding window is below <crop_end_gap_thr>. Nucleotide positions in the consensus sequence below
                 <crop_end_gap_thr> will be converted to -. Used with the --crop_end_gap_win option. Default: 0.1
 --start_patterns STRING      LTR elements always start with a conserved sequence pattern. TE Trimmer searches the 
                 beginning of the consensus sequence for these patterns. If the pattern is not found, it will extend the
                 search of <start_patterns> to up to 15 nucleotides from the beginning of the consensus sequence, and
                 redefine the start of the consensus element if the pattern is found. Note: The user can provide multiple
                 LTR start patterns in a comma-separated list, like: TG,TA,TC (No spaces; the order of patterns determines
                 the priority for the search). Default: TG
 --end_patterns STRING       LTR elements always end with a conserved sequence pattern. TE Trimmer searches the 
                 end of the consensus sequence for these patterns. If the pattern is not found, it will extend the
                 search of <start_patterns> to up to 15 nucleotides from the end of the consensus sequence, and
                 redefine the end of the consensus element if the pattern is found. Note: The user can provide multiple
                 LTR end patterns in a comma-separated list, like: CA,TA,GA (No spaces; the order of patterns determines
                 the priority for the search). Default: CA
 --mini_orf INTEGER       Define the minimum ORF length that will be predicted by TE Trimmer. Default: 200
 --check_extension_win INTEGER   Define window size for extension. Deafault: 150
 -t, --num_threads INTEGER    Thread number used for TE Trimmer. Default: 10
 --classify_unknown       Use RepeatClassfier to classify the consensus sequence if the input sequence is not
                 classfied or is unknown. Default: False
 --classify_all         Use RepeatClassfier to classify every consensus sequence. WARNING: This may take a long
                 time. Default: False
 --help             Show this message and exit.
```
