import os.path
from Bio.Align import AlignInfo
from Bio import AlignIO


# Check if file name contains the string LTR
def is_LTR(input_file):
    input_file_name = os.path.basename(input_file)
    return 'LTR' in input_file_name


# Generate consensus sequence
def generate_consensus_sequence(input_file, threshold, ambiguous):
    alignment = AlignIO.read(input_file, "fasta")
    summary = AlignInfo.SummaryInfo(alignment)
    consensus_seq = list(summary.dumb_consensus(threshold=threshold, ambiguous=ambiguous).upper())
    return consensus_seq


# Check if start and end are matchable with the given pattern
def check_start_end(consensus_seq, start, end, start_patterns, end_patterns):
    start_matched = end_matched = True

    if start_patterns:
        for start_pattern in start_patterns:
            if consensus_seq[start: start + len(start_pattern)] == list(start_pattern):
                break
        else:
            start_matched = False

    if end_patterns:
        for end_pattern in end_patterns:
            if consensus_seq[end - len(end_pattern): end] == list(end_pattern):
                break
        else:
            end_matched = False

    return start_matched, end_matched


def check_and_update(consensus_seq, start, end, start_patterns, end_patterns):
    """
    Check if start and end of the MSA are equal to the given patterns
    """
    # Ensure all patterns are uppercase
    start_patterns = [pattern.upper() for pattern in start_patterns.split(',')] if start_patterns else None  
    end_patterns = [pattern.upper() for pattern in end_patterns.split(',')] if end_patterns else None

    start_matched, end_matched = check_start_end(consensus_seq, start, end, start_patterns, end_patterns)

    if not start_matched or not end_matched:
        # If exact position matching fails, try sliding window approach
        if not start_matched and start_patterns:
            for start_pattern in start_patterns:
                start_window = consensus_seq[max(0, start - 15): start + 15]
                for i in range(len(start_window) - len(start_pattern) + 1):
                    if start_window[i: i + len(start_pattern)] == list(start_pattern):
                        start = max(0, start - 15) + i
                        break
                break

        if not end_matched and end_patterns:
            for end_pattern in end_patterns:
                end_window = consensus_seq[max(0, end - 15): end + 15]
                for i in reversed(range(len(end_window) - len(end_pattern) + 1)):
                    if end_window[i: i + len(end_pattern)] == list(end_pattern):
                        end = max(0, end - 15) + i + len(end_pattern)
                        break
                break

    # Check again after updating positions
    start_matched, end_matched = check_start_end(consensus_seq, start, end, start_patterns, end_patterns)
    return start_matched, end_matched, start, end
