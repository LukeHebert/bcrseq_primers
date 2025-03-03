'''
This script is meant to take an input FASTA of C gene sequences (CH1 for heavy 
chain and coding region for light chain) as input. It performs pairwise 
comparisons across all inputs to determine if any 20 bp-long windows overlap
between two or more sequences. The goal is to outline bad areas to avoid when
designin OE-RT-PCR or Nested primers targetting the C genes. If you prime those
exact 20 bp subsequences, you will not be able to avoid an off-target primed
sequence from another C gene.
'''

#!/usr/bin/env python3
import argparse
import itertools
import os
from Bio import SeqIO

def extend_match(seq1, seq2, pos1, pos2, initial_match):
    """
    Extend the initial exact match (of length k) to the right as far as possible.
    """
    match = initial_match
    while (pos1 + len(match) < len(seq1)) and (pos2 + len(match) < len(seq2)) \
          and seq1[pos1 + len(match)] == seq2[pos2 + len(match)]:
        match += seq1[pos1 + len(match)]
    return match

def main():
    parser = argparse.ArgumentParser(
        description="Identify shared subsequences (>= min_match bases) among immunoglobulin C gene sequences."
    )
    parser.add_argument('-i', '--input', required=True,
                        help="Input FASTA file containing immunoglobulin C gene sequences.")
    parser.add_argument('-n', '--num_bases', type=int, default=45,
                        help="Number of bases to consider from the start of each sequence (default 45).")
    parser.add_argument('--min_match', type=int, default=20,
                        help="Minimum length of shared subsequence (default 20).")
    args = parser.parse_args()

    # Determine output file name based on input file name and argument abbreviations.
    input_path = os.path.abspath(args.input)
    input_dir = os.path.dirname(input_path)
    input_base = os.path.splitext(os.path.basename(input_path))[0]
    # Create an abbreviated string from the arguments (e.g. "i45_n20")
    arg_str = f"i{args.num_bases}_n{args.min_match}"
    output_file = os.path.join(input_dir, f"{input_base}_analysis_{arg_str}.txt")

    # Read and trim sequences
    records = list(SeqIO.parse(args.input, "fasta"))
    trimmed_sequences = []
    for record in records:
        # Use at most the first N bases of the sequence
        trimmed_seq = str(record.seq)[:args.num_bases]
        trimmed_sequences.append((record.description, trimmed_seq))

    k = args.min_match
    num_bases = args.num_bases
    num_seqs = len(trimmed_sequences)
    num_windows = max(0, num_bases - k + 1)  # number of windows per sequence

    # Calculate the theoretical number of comparisons.
    # For each pair of sequences (num_seqs choose 2), every window in one is compared to every window in the other.
    theoretical_comparisons = (num_seqs * (num_seqs - 1) // 2) * (num_windows ** 2)

    # Build dictionary mapping each k-mer (of length k) to its occurrences.
    # Each occurrence is stored as a tuple: (sequence index, start position in trimmed sequence)
    kmer_dict = {}
    for idx, (desc, seq) in enumerate(trimmed_sequences):
        # Only consider positions where a full k-mer can be extracted
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            kmer_dict.setdefault(kmer, []).append((idx, i))

    # To avoid duplicate reports, store reported matches.
    reported = set()

    problematic_matches = 0

    # List to store output lines.
    output_lines = []

    # Include run parameters in the output file.
    output_lines.append("Run parameters:")
    output_lines.append(f"Input file: {args.input}")
    output_lines.append(f"Number of bases considered (num_bases): {args.num_bases}")
    output_lines.append(f"Minimum match length (min_match): {args.min_match}")
    output_lines.append("")

    # For each k-mer that occurs in more than one sequence, check for extended matches.
    for kmer, occurrences in kmer_dict.items():
        if len(occurrences) < 2:
            continue  # k-mer unique to a single sequence; nothing to compare.
        # Consider all unique pairs among the occurrences.
        for (idx1, pos1), (idx2, pos2) in itertools.combinations(occurrences, 2):
            if idx1 == idx2:
                continue  # Skip comparisons within the same sequence.
            # Canonically order the pair (by sequence index) for consistent reporting.
            if idx1 > idx2:
                idx1, pos1, idx2, pos2 = idx2, pos2, idx1, pos1

            desc1, seq1 = trimmed_sequences[idx1]
            desc2, seq2 = trimmed_sequences[idx2]

            # Extend the match to the right from the k-mer start.
            match = extend_match(seq1, seq2, pos1, pos2, kmer)
            if len(match) >= k:
                # Create a unique key for this reported match to avoid duplicates.
                key = (idx1, pos1, idx2, pos2, len(match))
                if key not in reported:
                    reported.add(key)
                    problematic_matches += 1
                    # Report positions as 1-indexed.
                    start1 = pos1 + 1
                    end1 = pos1 + len(match)
                    start2 = pos2 + 1
                    end2 = pos2 + len(match)
                    # Append the match information to the output.
                    output_lines.append(f"{desc1} (positions {start1} to {end1}) and {desc2} (positions {start2} to {end2}) share a subsequence: {match}")

    # Append summary of comparisons.
    output_lines.append("\nSummary:")
    output_lines.append(f"Total theoretical window comparisons: {theoretical_comparisons}")
    output_lines.append(f"Total problematic subsequences found: {problematic_matches}")

    # Write results to the output file.
    with open(output_file, "w") as out_f:
        out_f.write("\n".join(output_lines))

    print(f"Analysis complete. Results written to {output_file}")

if __name__ == "__main__":
    main()
