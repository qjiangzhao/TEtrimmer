# TE-Trimmer
TE Trimmer: a tool to replace transposable element manual curation.

Usage:

# Use --help for a detailed list of options
python ./path_to_TE_Trimmer_folder/bin/main.py --help 

# Example of running TE Trimmer
python ./path_to_TE_Trimmer_folder/bin/main.py --input_file [your_TE_consensus_file_path] \
                                               --genome_file [your_genome_file_path] \
                                               --output_dir [output_directory]
                                               
# Example of graphical user interface-based proof annotation
python ./path_to_TE_Trimmer_folder/bin/Class_TKinter_proof_annotation.py --te_trimmer_output_dir [your_TE_Trimmer_output_directory/Multiple_sequence_alignment_less] 
