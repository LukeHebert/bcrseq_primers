'''
Reads a single FASTA file of target sequences.

Parses a user parameters for primer design such as
    -primer length range (for example “20-23”), 
    -a maximum number of allowed degenerate positions, and 
    -whether to target the “start” or “end” of a sequence with optional target 
    offset/wiggle room

For each allowed offset (from 0 to the target offset) and for each allowed 
primer length the script tries seeding a group with one unassigned target 
and then greedily adds other primable sequence targets if the consensus 
(computed column by column and “degeneratized” via IUPAC codes when needed) does
not exceed the allowed number of degenerate positions. (Recall that a 
“degenerate” column is one where the target nucleotides are not unanimous.)

It then picks the candidate primers that cover the most primable sequences (and,
when tied, prefers the longer primer) and removes those targets from further 
consideration.

Finally, it writes 
    -one FASTA file per group (with the group’s targets)
    -one FASTA file of the primer sequences (with a header that notes the unique 
    group/primer id and the list of target ids covered), and 
    -a TSV file that lists each target, its group assignment, and the final 
    primer sequence.

(Note that primer design is a minimal set–cover problem, which is NP-hard; 
this script uses one of multiple possible heuristics to solve the problem.)
'''

#!/usr/bin/env python3
import argparse
import os
from collections import OrderedDict

def parse_args():
    parser = argparse.ArgumentParser(description="Design degenerate primers covering target sequences with minimal primer groups.")
    parser.add_argument("input_fasta", help="Input FASTA file with target sequences")
    parser.add_argument("--primer_length_range", required=True, help="Allowed primer lengths, e.g. '20-23'")
    parser.add_argument("--max_degenerate", type=int, required=True, help="Maximum number of degenerate positions allowed in any primer")
    parser.add_argument("--target", choices=["start", "end"], required=True,
                        help="Whether the primer should target the 'start' or 'end' of the sequence")
    parser.add_argument("--target_offset", type=int, default=0,
                        help="Optional shift (in nucleotides) allowed from the designated position (default: 0)")
    parser.add_argument("--output_prefix", default="output", help="Prefix for output files")
    return parser.parse_args()

def parse_fasta(filename):
    sequences = OrderedDict()
    with open(filename, 'r') as f:
        seq_id = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id is not None:
                    sequences[seq_id] = ''.join(seq_lines)
                # take first word after ">" as id
                seq_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if seq_id is not None:
            sequences[seq_id] = ''.join(seq_lines)
    return sequences

def get_candidate_region(seq, target, offset, length):
    """Extract a candidate region given a sequence, target position, offset and primer length.
       For 'start': region = seq[offset:offset+length]
       For 'end': if offset==0, region = seq[-length:], else region = seq[-(length+offset):-offset]
    """
    if target == "start":
        if offset + length <= len(seq):
            return seq[offset:offset+length]
        else:
            return None
    elif target == "end":
        if length + offset <= len(seq):
            if offset == 0:
                return seq[-length:]
            else:
                return seq[-(length+offset):-offset]
        else:
            return None
    else:
        return None

def degenerate_code(nucs):
    """Return the IUPAC degenerate nucleotide code for a set of nucleotides."""
    mapping = {
        frozenset(['A']): 'A',
        frozenset(['C']): 'C',
        frozenset(['G']): 'G',
        frozenset(['T']): 'T',
        frozenset(['A', 'G']): 'R',
        frozenset(['C', 'T']): 'Y',
        frozenset(['G', 'T']): 'K',
        frozenset(['A', 'C']): 'M',
        frozenset(['C', 'G']): 'S',
        frozenset(['A', 'T']): 'W',
        frozenset(['A', 'C', 'G']): 'V',
        frozenset(['A', 'C', 'T']): 'H',
        frozenset(['A', 'G', 'T']): 'D',
        frozenset(['C', 'G', 'T']): 'B',
        frozenset(['A', 'C', 'G', 'T']): 'N'
    }
    return mapping.get(frozenset(nucs), 'N')

def compute_consensus(seqs):
    """Given a list of sequences (all same length), compute the consensus,
       count how many positions are degenerate, and return (consensus, degenerate_count).
    """
    if not seqs:
        return "", 0
    length = len(seqs[0])
    consensus = []
    degenerate_count = 0
    for i in range(length):
        column = {seq[i] for seq in seqs}
        if len(column) == 1:
            consensus.append(next(iter(column)))
        else:
            degenerate_count += 1
            consensus.append(degenerate_code(column))
    return ''.join(consensus), degenerate_count

def design_groups(sequences, target, allowed_offset, min_length, max_length, max_degenerate):
    """
    Greedily form groups of targets that can share a primer.
    For each allowed offset (0..allowed_offset) and each allowed length,
    try seeding a group from each unassigned sequence and then add others if the
    consensus (computed over the candidate regions) has <= max_degenerate degenerate positions.
    Returns a list of groups. Each group is a dict with keys:
      - 'primer': the consensus primer sequence (with degenerate bases)
      - 'offset': the chosen offset
      - 'length': the chosen primer length
      - 'members': list of target IDs in the group
      - 'degenerate_count': number of degenerate positions in the consensus.
    """
    ungrouped = set(sequences.keys())
    groups = []
    
    while ungrouped:
        best_group = None
        best_consensus = None
        best_off = None
        best_len = None
        best_deg = None
        
        # Try all allowed offset and length combinations
        for off in range(allowed_offset + 1):
            for length in range(min_length, max_length + 1):
                # For each seed in ungrouped, build a group
                for seed in list(ungrouped):
                    candidate_seed = get_candidate_region(sequences[seed], target, off, length)
                    if candidate_seed is None:
                        continue
                    current_group = {seed: candidate_seed}
                    # Try adding every other ungrouped target (order is arbitrary)
                    for other in list(ungrouped):
                        if other in current_group:
                            continue
                        candidate_other = get_candidate_region(sequences[other], target, off, length)
                        if candidate_other is None:
                            continue
                        # Build temporary list of candidate regions for current group plus this candidate
                        temp_candidates = list(current_group.values()) + [candidate_other]
                        consensus, deg_count = compute_consensus(temp_candidates)
                        if deg_count <= max_degenerate:
                            current_group[other] = candidate_other
                    group_size = len(current_group)
                    # Choose best group: first by size, then by longer primer length if tied.
                    if best_group is None or group_size > len(best_group) or (group_size == len(best_group) and length > best_len):
                        best_group = current_group.copy()
                        best_consensus, best_deg = compute_consensus(list(best_group.values()))
                        best_off = off
                        best_len = length
        if best_group is None:
            # Fallback: assign one target if nothing worked (should not occur in practice)
            seed = next(iter(ungrouped))
            candidate = get_candidate_region(sequences[seed], target, 0, max_length)
            if candidate is None:
                candidate = sequences[seed]
            best_group = {seed: candidate}
            best_consensus, best_deg = compute_consensus([candidate])
            best_off = 0
            best_len = max_length
        
        # Record the best group found
        groups.append({
            'primer': best_consensus,
            'offset': best_off,
            'length': best_len,
            'members': list(best_group.keys()),
            'degenerate_count': best_deg
        })
        # Remove grouped sequences from further consideration
        for s in best_group.keys():
            ungrouped.discard(s)
    return groups

def write_outputs(groups, sequences, output_prefix):
    """
    Writes three outputs:
      1. One FASTA file per group of target sequences (named group_<i>.fasta).
      2. A single FASTA file of primer sequences (named primers.fasta). Each header includes a unique group id and the list of target ids covered.
      3. A TSV file (primer_assignments.tsv) enumerating for each target: target id, group id, primer consensus.
    """
    # Create output directory if needed.
    outdir = output_prefix + "_groups"
    os.makedirs(outdir, exist_ok=True)
    
    # Prepare primer FASTA and TSV data
    primer_fasta_lines = []
    tsv_lines = ["Target\tGroup\tPrimer"]
    
    for i, group in enumerate(groups, start=1):
        group_id = f"primer{i}"
        # Write group FASTA file with the target sequences
        group_filename = os.path.join(outdir, f"{group_id}.fasta")
        with open(group_filename, "w") as gf:
            for tid in group['members']:
                gf.write(f">{tid}\n{sequences[tid]}\n")
                tsv_lines.append(f"{tid}\t{group_id}\t{group['primer']}")
        # Prepare primer FASTA record: include list of target IDs in header.
        header = f">{group_id} | covers {','.join(group['members'])}"
        primer_fasta_lines.append(header)
        primer_fasta_lines.append(group['primer'])
    
    # Write primer FASTA file
    with open(output_prefix + "_primers.fasta", "w") as pf:
        pf.write("\n".join(primer_fasta_lines) + "\n")
    # Write TSV file
    with open(output_prefix + "_primer_assignments.tsv", "w") as tf:
        tf.write("\n".join(tsv_lines) + "\n")

def main():
    args = parse_args()
    # Parse primer length range (e.g., "20-23")
    try:
        min_length, max_length = map(int, args.primer_length_range.split("-"))
    except Exception as e:
        print("Error parsing primer_length_range. Please provide in the format min-max (e.g. 20-23).")
        return

    # Read input sequences
    sequences = parse_fasta(args.input_fasta)
    if not sequences:
        print("No sequences found in the input FASTA.")
        return

    # Design groups using the greedy heuristic.
    groups = design_groups(sequences, args.target, args.target_offset, min_length, max_length, args.max_degenerate)
    
    # Write outputs: one FASTA per group, one primer FASTA, and a TSV mapping.
    write_outputs(groups, sequences, args.output_prefix)
    print(f"Designed {len(groups)} primer group(s).")
    print("Outputs written: one primer FASTA file, one TSV file, and one FASTA per primer group (in folder '{}_groups').".format(args.output_prefix))

if __name__ == "__main__":
    main()
