from collections import defaultdict


class RepeatMaskerAnalyzer:
    def __init__(self, file_path):
        # Read and parse the RepeatMasker output file upon instantiation
        self.rm_data = self.read_repeatmasker_output(file_path)

    def read_repeatmasker_output(self, file_path):
        # This function reads a RepeatMasker output file and stores the data as a list of dictionaries
        data = []
        with open(file_path, 'r') as file:
            for line in file:
                # Ignore lines starting with "SW", "score", and empty lines
                if line.startswith("SW") or line.startswith("score") or line.startswith("  SW") or line.strip() == "":
                    continue

                # Split each line into columns
                cols = line.split()

                # Create a dictionary for each row with the appropriate column names and values
                row = {
                    'SW_score': float(cols[0]),
                    'perc_div': float(cols[1]),
                    'perc_del': float(cols[2]),
                    'perc_ins': float(cols[3]),
                    'seq_name': cols[4],
                    'seq_start': int(cols[5]),
                    'seq_end': int(cols[6]),
                    'seq_remaining': cols[7],
                    'orientation': cols[8],
                    'repeat_name': cols[9],
                    'repeat_class_family': cols[10],
                    'repeat_start': int(cols[11]),
                    'repeat_end': int(cols[12])
                }
                data.append(row)
        return data

    def weighted_average_divergences(self):
        # This function calculates the weighted average divergence for all available TE consensus families

        # Initialize dictionaries to store the weighted sum of divergence values and total alignment lengths for each TE family
        weighted_sums = defaultdict(float)
        total_lengths = defaultdict(float)

        # Loop through each row in the RepeatMasker data
        for row in self.rm_data:
            # Get the TE family name, divergence, and alignment length
            te_family = row['repeat_name']
            divergence = row['perc_div']
            alignment_length = row['seq_end'] - row['seq_start'] + 1

            # Accumulate the weighted sum of divergence values and total alignment lengths for each TE family
            weighted_sums[te_family] += divergence * alignment_length
            total_lengths[te_family] += alignment_length

        # Calculate the weighted average divergence for each TE family
        weighted_avg_divergences = {}
        for te_family in weighted_sums.keys():
            if total_lengths[te_family] == 0:
                weighted_avg_divergences[te_family] = None
            else:
                weighted_avg_divergences[te_family] = weighted_sums[te_family] / total_lengths[te_family]

        return weighted_avg_divergences
