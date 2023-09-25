import click
import os
from Bio import AlignIO


def calc_conservation(col):
    """
    Calculate the conservation of a column as the fraction of the most common nucleotide.
    :param col: one single column content for multiple sequence alignment
    :return: most abundant nucleotide proportion for this column (gap isn't included)
    """
    # Convert all nucleotides to lowercase
    col = [nucleotide.lower() for nucleotide in col]
    nucleotide_counts = {'a': 0, 'c': 0, 'g': 0, 't': 0}

    for nucleotide in col:
        if nucleotide in nucleotide_counts:
            nucleotide_counts[nucleotide] += 1

    total_nucleotides = sum(nucleotide_counts.values())
    max_count = max(nucleotide_counts.values())

    return max_count / total_nucleotides


@click.command()
@click.option("--input_file", "-i", required="True", type=str,
              help="Multiple sequence alignment fasta file path")
@click.option("--output_file", "-o", required="True", type=str,
              help="Output file")
@click.option("--gap_threshold", default=0.8, type=float,
              help="Columns with gap percentage greater (not include equal) than threshold will be removed directly")
@click.option("--simi_check_gap_thr", default=0.4, type=float,
              help="Columns with gap percentage between --simi_check_gap_thr and --gap_threshold will be deleted "
                   "when the most abundant nucleotide proportion is smaller than --similarity_thr")
@click.option("--similarity_thr", default=0.7, type=float,
              help="If the most abundant nucleotide proportion is smaller than --similarity_thr and gap proportion is"
                   "greater than --simi_check_gap_thr, this colum will be removed")
@click.option("--min_nucleotide", default=5, type=int,
              help="Columns with less than --min_nucleotide nucleotides will be removed")
def remove_gaps_with_similarity_check(input_file, output_file, gap_threshold,
                                      simi_check_gap_thr, similarity_thr, min_nucleotide):
    """
    Remove columns directly when gap percentage is bigger than --gap_threshold. Remove columns when nucleotide number
    is less than --min_nucleotide.
    When gap percentage is equal or bigger than --simi_check_gap_thr and smaller than --gap_threshold,
    it will calculate most abundant nucleotide proportion for this column. If it is less than --similarity_thr
    this column will be removed

    """
    if os.path.isfile(input_file):
        keep_list = []
        MSA_mafft = AlignIO.read(input_file, "fasta")

        column_mapping = {}  # Stores the mapping of column indices from filtered MSA to original MSA

        if len(MSA_mafft) < 5:
            raise ValueError("Number of sequences is less than 5. Cannot remove gaps.")

        for col_idx in range(MSA_mafft.get_alignment_length()):
            col = MSA_mafft[:, col_idx]
            gap_count = col.count('-')
            gap_fraction = gap_count / len(col)

            # Check if the number of nucleotides in a column is less than 5
            nt_count = len(col) - gap_count
            if nt_count < min_nucleotide:
                continue

            # If the gap fraction is less than the gap_threshold, consider the column for further analysis
            if gap_fraction <= gap_threshold:
                # If the gap fraction is between 0.4 and 0.8, check for the similarity of the remaining nucleotides
                if simi_check_gap_thr > gap_fraction:
                    keep_list.append(col_idx)
                    column_mapping[len(keep_list) - 1] = col_idx
                elif simi_check_gap_thr <= gap_fraction:
                    nt_fraction = calc_conservation(col)

                    # If the nucleotides are less than similarity_thr, skip this column
                    if nt_fraction < similarity_thr:
                        continue
                    else:
                        # If the column passes all checks, add it to the keep_list
                        keep_list.append(col_idx)
                        # Store the mapping of original MSA index to filtered MSA index
                        column_mapping[len(keep_list) - 1] = col_idx

        # The index in python won't include the second value. That it is to say the returned end_posit from
        # DefineBoundary() will one more index than the final len(keep_list) -1] for this reason, you have to
        # add one more value
        column_mapping[len(keep_list)] = column_mapping[len(keep_list) - 1] + 1

        # Keep the columns
        MSA_mafft_filtered = MSA_mafft[:, keep_list[0]:keep_list[0] + 1]
        for i in keep_list[1:]:
            MSA_mafft_filtered += MSA_mafft[:, i:i + 1]

        # Write the filtered MSA to the output file
        with open(output_file, 'w') as f:
            AlignIO.write(MSA_mafft_filtered, f, 'fasta')
        return output_file

    else:
        raise FileNotFoundError(f"The file '{input_file}' does not exist.")


if __name__ == '__main__':
    remove_gaps_with_similarity_check()
