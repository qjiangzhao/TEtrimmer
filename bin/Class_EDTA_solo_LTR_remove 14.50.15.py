#!/bin/python

import re


def read_gff_data(file_path):
    with open(file_path, 'r') as gff_file:
        return gff_file.read().split('\n')


def write_gff_data(file_path, data):
    with open(file_path, 'w') as gff_file:
        gff_file.write('\n'.join(data))


def process_gff_data(gff_data):
    remove_id = []

    for idx, line in enumerate(gff_data):
        fields = line.split()
        if not fields:
            continue
        chrom = fields[0]
        seq_type = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        LTR_ID = fields[8].split(";")[1]
        if "LTR_retrotransposon" in seq_type:
            for inner_line in gff_data:
                inner_fields = inner_line.split()
                if not inner_fields:
                    continue
                inner_chrom = inner_fields[0]
                inner_seq_type = inner_fields[2]
                inner_start = int(inner_fields[3])
                inner_end = int(inner_fields[4])

                if "LTR_retrotransposon" in inner_seq_type and inner_start == start and inner_end < end and inner_chrom == chrom:
                    remove_id.append(LTR_ID.replace("Parent=repeat_", ""))
                    break
                elif "LTR_retrotransposon" in inner_seq_type and inner_start > start and inner_end == end and inner_chrom == chrom:
                    remove_id.append(LTR_ID.replace("Parent=repeat_", ""))
                    break
    return remove_id


input_file_path = '/work/ur376715/TE_analysis/Second_time_EDTA_result_based_on_bgh_gene_model_v4.2/bgh_dh14_v4.fa.mod.EDTA.intact.gff3'
output_file_path = '/work/ur376715/TE_analysis/Second_time_EDTA_result_based_on_bgh_gene_model_v4.2/bgh_dh14_v4.fa.mod.EDTA.intact_modified.gff3'

gff_data = read_gff_data(input_file_path)
new_gff_data = process_gff_data(gff_data)
write_gff_data(output_file_path, new_gff_data)
