import collections
from Bio import AlignIO
from collections import defaultdict
from itertools import combinations
from sklearn.cluster import DBSCAN
import numpy as np


def get_nt_radius_indices(sequence, center_idx, nt_radius, direction="both"):
    """Finds MSA column indices containing nt_radius nucleotides, skipping gaps."""
    seq_len = len(sequence)
    start_col = center_idx
    end_col = center_idx

    if direction in ["left", "both"]:
        count = 0
        while start_col > 0 and count < nt_radius:
            start_col -= 1
            if sequence[start_col] != '-':
                count += 1
    if direction in ["right", "both"]:
        count = 0
        while end_col < seq_len and count < nt_radius:
            if sequence[end_col] != '-':
                count += 1
            end_col += 1
    return start_col, end_col


# --- K-MER SIMILARITY HELPERS ---
def get_sub_kmers(sequence, k=2):
    """Breaks a TSD into small overlapping k-mers (default 2-mers)."""
    return set(sequence[i:i + k] for i in range(len(sequence) - k + 1))


def calculate_jaccard_similarity(seq1, seq2):
    """Calculates similarity based on shared sub-k-mers."""
    set1 = get_sub_kmers(seq1.upper())
    set2 = get_sub_kmers(seq2.upper())
    if not set1 or not set2: return 0
    return len(set1.intersection(set2)) / len(set1.union(set2))


def identify_tsds_with_boundary_logic(file_path, left_boundary, right_boundary,
                                      nt_radius=100, k_min=4, k_max=7):
    alignment = AlignIO.read(file_path, "fasta")
    num_seqs = len(alignment)
    location_hits = defaultdict(list)

    for record in alignment:
        seq_str = str(record.seq)

        # Map: Nucleotide_Index <-> MSA_Column for the current sequence
        all_nt_pos = [i for i, char in enumerate(seq_str) if char != '-']
        nt_to_msa = {i: col for i, col in enumerate(all_nt_pos)}
        msa_to_nt = {col: i for i, col in enumerate(all_nt_pos)}

        # 1. Define Windows
        l_col_start, l_col_end = get_nt_radius_indices(seq_str, left_boundary, nt_radius, direction="left")
        r_col_start, r_col_end = get_nt_radius_indices(seq_str, right_boundary, nt_radius, direction="right")

        # 2. Map nucleotides in windows to indices within the all_nt_pos list
        l_nt_pos_map = [c for c in range(l_col_start, l_col_end) if seq_str[c] != '-']
        r_nt_pos_map = [c for c in range(r_col_start, r_col_end) if seq_str[c] != '-']

        # --- 3. Exact K-mer Search & Extension (No Mismatch) ---
        # Logic: Start from k_min to find all seeds, then extend to maximal match.
        k = k_min
        left_kmers_seeds = defaultdict(list)

        # Build Seed Map for Left Window
        for i in range(len(l_nt_pos_map) - k + 1):
            kmer_cols = l_nt_pos_map[i: i + k]
            kmer = "".join(seq_str[c] for c in kmer_cols)

            # This will record the position for each set of k-mer
            # i is the position number of the k-mer first letter
            left_kmers_seeds[kmer].append(i)

        # Track unique hits for this record to avoid sub-kmer redundancy
        record_unique_hits = set()

        # store temporary hits first
        record_hits = []

        # Search & Extend in Right Window
        for j in range(len(r_nt_pos_map) - k + 1):
            kmer_cols_r = r_nt_pos_map[j: j + k]
            kmer = "".join(seq_str[c] for c in kmer_cols_r)

            # When find the identical k-mer from the right sequence to the left sequence
            if kmer in left_kmers_seeds:

                # Go through each position of this k-mer in the left sequence
                for seed_idx_l in left_kmers_seeds[kmer]:
                    # 1. Extend Leftwards
                    l_shift = 0
                    while (seed_idx_l - l_shift - 1 >= 0 and j - l_shift - 1 >= 0):
                        char_l = seq_str[l_nt_pos_map[seed_idx_l - l_shift - 1]]
                        char_r = seq_str[r_nt_pos_map[j - l_shift - 1]]
                        if char_l == char_r:
                            l_shift += 1
                        else:
                            break

                    # 2. Extend Rightwards
                    r_shift = 0
                    while (seed_idx_l + k + r_shift < len(l_nt_pos_map) and
                           j + k + r_shift < len(r_nt_pos_map)):
                        char_l = seq_str[l_nt_pos_map[seed_idx_l + k + r_shift]]
                        char_r = seq_str[r_nt_pos_map[j + k + r_shift]]
                        if char_l == char_r:
                            r_shift += 1
                        else:
                            break

                    # 3. Get Final Maximal Boundaries
                    l_indices = l_nt_pos_map[seed_idx_l - l_shift: seed_idx_l + k + r_shift]
                    r_indices = r_nt_pos_map[j - l_shift: j + k + r_shift]

                    # Deduplicate using coordinates (l_start, l_end, r_start, r_end)
                    hit_coords = (l_indices[0], l_indices[-1], r_indices[0], r_indices[-1])

                    if hit_coords not in record_unique_hits:
                        record_unique_hits.add(hit_coords)

                        # Apply TSD Boundary Logic
                        last_nt_l_col = l_indices[-1]
                        nt_idx_l = msa_to_nt[last_nt_l_col]

                        if nt_idx_l + 1 in nt_to_msa:
                            te_start_col = nt_to_msa[nt_idx_l + 1]

                            first_nt_r_col = r_indices[0]
                            nt_idx_r = msa_to_nt[first_nt_r_col]

                            if nt_idx_r - 1 in nt_to_msa:
                                te_end_col = nt_to_msa[nt_idx_r - 1]

                                final_kmer_seq = "".join(seq_str[c] for c in l_indices)
                                # if the k-mer is longer than >=5 but all with the same nucleotide delete it
                                if len(final_kmer_seq) >= 5 and len(set(final_kmer_seq.lower())) == 1:
                                    continue

                                record_hits.append(
                                    ((te_start_col, te_end_col),
                                     record.id,
                                     final_kmer_seq)
                                )

        # Keep longest hit first
        record_hits.sort(key=lambda x: len(x[2]), reverse=True)

        filtered_hits = []

        for hit in record_hits:
            coords, rec_id, seq = hit
            keep = True

            for kept_hit in filtered_hits:
                kept_coords, kept_id, kept_seq = kept_hit

                # same record + same TE end + shorter sequence covered by longer one
                if (rec_id == kept_id and
                    coords[1] == kept_coords[1] and
                    seq in kept_seq):
                    keep = False
                    break

            if keep:
                filtered_hits.append(hit)

        # save final non-redundant hits
        for coords, rec_id, seq in filtered_hits:
            location_hits[coords].append((rec_id, seq))

    return location_hits

# --- 4. ENTIRE FILTERING PIPELINE (K-mer Profile + Nucleotide Cleaning) ---
def process_tsd_hits_old(hits, max_sim_threshold=0.35, nt_threshold=0.9):
    """
    1. Diversity Filtering: Ensures TSDs at a location are diverse.
    2. Nucleotide Cleaning: Trims homopolymer-like boundaries (poly-A/poly-T).
    """
    final_processed_hits = defaultdict(list)

    for (te_start, te_end), matches in hits.items():
        if not matches: continue

        # --- PHASE 1: Diversity Filtering ---
        kmers = [m[1].upper() for m in matches]

        hit_n = len(kmers)

        if hit_n >= 2:
            sim_scores = [calculate_jaccard_similarity(s1, s2) for s1, s2 in combinations(kmers, 2)]
            avg_similarity = sum(sim_scores) / len(sim_scores)

            # If similarity is too high, skip this location entirely
            if avg_similarity > max_sim_threshold:
                continue

        # --- PHASE 2: Nucleotide Cleaning (Poly-N Detection) ---
        first_letters = [k[0] for k in kmers if len(k) > 0]
        last_letters = [k[-1] for k in kmers if len(k) > 0]

        shift_start = 0
        shift_end = 0

        if hit_n >= 5:
            # Check if start column is > 90% same nucleotide
            if first_letters:
                counts_f = collections.Counter(first_letters)
                _, freq_f = counts_f.most_common(1)[0]
                if (freq_f / hit_n) >= nt_threshold:
                    shift_end = 1  # Adjust coordinate

            # Check if end column is > 85% same nucleotide
            if last_letters:
                counts_l = collections.Counter(last_letters)
                _, freq_l = counts_l.most_common(1)[0]
                if (freq_l / hit_n) >= nt_threshold:
                    shift_start = -1  # Adjust coordinate

            # Recalculate coordinates and trim the TSD strings
            new_loc = (te_start + shift_start, te_end + shift_end)
            cleaned_matches = []
            for rec_id, seq in matches:
                t_seq = seq
                if shift_start != 0: t_seq = t_seq[1:]
                if shift_end != 0: t_seq = t_seq[:-1]
                cleaned_matches.append((rec_id, t_seq))

            final_processed_hits[new_loc].extend(cleaned_matches)

        else:
            final_processed_hits[te_start, te_end].extend(matches)

    return final_processed_hits


def process_tsd_hits(hits, max_sim_threshold=0.35, nt_threshold=0.9, complexity_threshold=0.8):
    """
    1. Diversity Filtering: Ensures TSDs at a location are diverse.
    2. Nucleotide Cleaning: Trims homopolymer-like boundaries (poly-A/poly-T).
    3. Complexity Filtering: Removes sequences dominated by a single nucleotide.
    """
    final_processed_hits = defaultdict(list)

    for (te_start, te_end), matches in hits.items():
        if not matches: continue

        # --- PHASE 1: Diversity Filtering ---
        kmers = [m[1].upper() for m in matches]
        hit_n = len(kmers)

        if hit_n >= 2:
            sim_scores = [calculate_jaccard_similarity(s1, s2) for s1, s2 in combinations(kmers, 2)]
            avg_similarity = sum(sim_scores) / len(sim_scores)

            if avg_similarity > max_sim_threshold:
                continue

        # --- PHASE 2: Nucleotide Cleaning & Complexity Check ---
        cleaned_matches = []

        # We check each sequence in the match list
        for rec_id, seq in matches:
            seq_upper = seq.upper()
            seq_len = len(seq_upper)

            if seq_len == 0: continue

            # 1. Complexity check: count the most frequent nucleotide in this specific k-mer
            counts = collections.Counter(seq_upper)
            _, most_freq_count = counts.most_common(1)[0]

            # If the sequence is > 80% (complexity_threshold) one nucleotide, delete/skip it
            if (most_freq_count / seq_len) > complexity_threshold:
                continue

            cleaned_matches.append((rec_id, seq))

        # Update hit count after complexity filtering
        hit_n = len(cleaned_matches)
        if hit_n == 0: continue

        # --- PHASE 3: Poly-N Detection (Coordinate Shifting) ---
        first_letters = [m[1][0].upper() for m in cleaned_matches if len(m[1]) > 0]
        last_letters = [m[1][-1].upper() for m in cleaned_matches if len(m[1]) > 0]

        shift_start = 0
        shift_end = 0

        if hit_n >= 5:
            # Check if start column is > 90% same nucleotide
            if first_letters:
                counts_f = collections.Counter(first_letters)
                _, freq_f = counts_f.most_common(1)[0]
                if (freq_f / hit_n) >= nt_threshold:
                    shift_end = 1

                    # Check if end column is > 90% same nucleotide
            if last_letters:
                counts_l = collections.Counter(last_letters)
                _, freq_l = counts_l.most_common(1)[0]
                if (freq_l / hit_n) >= nt_threshold:
                    shift_start = -1

                    # Apply coordinate shifts and final trimming
            new_loc = (te_start + shift_start, te_end + shift_end)
            final_matches = []
            for rec_id, seq in cleaned_matches:
                t_seq = seq
                if shift_start != 0 and len(t_seq) > 1: t_seq = t_seq[1:]
                if shift_end != 0 and len(t_seq) > 1: t_seq = t_seq[:-1]
                final_matches.append((rec_id, t_seq))

            final_processed_hits[new_loc].extend(final_matches)

        else:
            # If cluster is small, keep cleaned matches without coordinate shifting
            final_processed_hits[te_start, te_end].extend(cleaned_matches)

    return final_processed_hits


def cluster_tsds_ml(filtered_hits, eps=1.1):
    """
    DBSCAN with Euclidean distance.
    - eps=1.1: Groups hits if coordinates differ by 1 in one direction.
    - eps=1.5: Groups hits if coordinates differ by 1 in both directions.
    """
    if not filtered_hits:
        return []

    locs = list(filtered_hits.keys())
    data = np.array(locs)

    # Use Euclidean to allow for fractional distance tuning
    db = DBSCAN(eps=eps, min_samples=1, metric='manhattan').fit(data)

    clusters = defaultdict(list)
    for idx, label in enumerate(db.labels_):
        clusters[label].append(locs[idx])

    results = []
    for label, member_locs in clusters.items():
        all_matches_with_coords = []
        for l in member_locs:
            for match in filtered_hits[l]:
                all_matches_with_coords.append({
                    "id": match[0],
                    "seq": match[1],
                    "coords": l
                })

        results.append({
            "cluster_id": label,
            "all_locations": member_locs,
            "hits": all_matches_with_coords,
            "unique_count": len(set(m["id"] for m in all_matches_with_coords))
        })

    results.sort(key=lambda x: x['unique_count'], reverse=True)
    return results


def TSD_identifier(left_n, right_n, extend_n, input_file):
    location_hits = identify_tsds_with_boundary_logic(input_file, left_n, right_n, nt_radius=extend_n)

    filtered_location_hits = process_tsd_hits(location_hits, max_sim_threshold=0.3, nt_threshold=0.95)


    final_clusters = cluster_tsds_ml(filtered_location_hits, eps=1)

# desired output format

# if there are TSD near the defined boundary
# if possible to use TSD for the boundary definition?
# Give the TSD position