#!/usr/bin/env python3
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

import argparse
import os
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_args():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Design degenerate primers covering target sequences with minimal primer groups."
    )
    parser.add_argument("input_fasta", help="Input FASTA file with target sequences")
    parser.add_argument("--primer_length_range", required=True, help="Allowed primer lengths, e.g. '20-23'")
    parser.add_argument("--max_degenerate", type=int, required=True, help="Maximum number of degenerate positions allowed in any primer")
    parser.add_argument("--target", choices=["start", "end"], required=True,
                        help="Whether the primer should target the 'start' or 'end' of the sequence")
    parser.add_argument("--target_offset", type=int, default=0,
                        help="Optional shift (in nucleotides) allowed from the designated position (default: 0)")
    parser.add_argument("--orientation", choices=["fwd", "rev"], default="fwd",
                        help="Primer orientation: 'fwd' (default) for forward or 'rev' for reverse complement output")
    parser.add_argument("--output_prefix", default="output", help="Prefix for output files (placed in the input FASTA directory)")
    return parser.parse_args()

def parse_fasta(filename):
    """
    Read sequences from a FASTA file.
    """
    sequences = OrderedDict()
    for record in SeqIO.parse(filename, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def get_candidate_region(seq, target, offset, length):
    """
    Extract a candidate region given a sequence, target position, offset, and primer length.
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
    """
    Return the IUPAC degenerate nucleotide code for a set of nucleotides.
    """
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
    """
    Given a list of sequences (all of the same length), compute the consensus sequence.
    In each column, if nucleotides are unanimous, that nucleotide is used; otherwise, the corresponding 
    degenerate nucleotide code is used. Also counts the number of degenerate positions.
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
    
    For each allowed offset (0..allowed_offset) and each allowed primer length,
    a seed target is chosen from the unassigned set. Then other targets are added if the
    consensus of their candidate regions (computed column-by-column) remains within the
    maximum allowed degenerate positions. The best group is the one that covers the most targets,
    preferring longer primer lengths when tied.
    """
    ungrouped = set(sequences.keys())
    groups = []
    
    while ungrouped:
        best_group = None
        best_off = None
        best_len = None
        best_deg = None
        
        # Try all allowed offset and primer length combinations.
        for off in range(allowed_offset + 1):
            for length in range(min_length, max_length + 1):
                # For each seed in ungrouped targets, try building a group.
                for seed in list(ungrouped):
                    candidate_seed = get_candidate_region(sequences[seed], target, off, length)
                    if candidate_seed is None:
                        continue
                    current_group = {seed: candidate_seed}
                    # Try adding every other ungrouped target.
                    for other in list(ungrouped):
                        if other in current_group:
                            continue
                        candidate_other = get_candidate_region(sequences[other], target, off, length)
                        if candidate_other is None:
                            continue
                        # Build temporary list of candidate regions for current group plus this candidate.
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
        
        groups.append({
            'primer': best_consensus,
            'offset': best_off,
            'length': best_len,
            'members': list(best_group.keys()),
            'degenerate_count': best_deg
        })
        # Remove grouped sequences from further consideration.
        for s in best_group.keys():
            ungrouped.discard(s)
    return groups

def write_outputs(groups, sequences, args):
    """
    Writes output files:
      1. One FASTA file per group of target sequences (named primer<i>.fasta in a groups directory).
      2. A single FASTA file of primer sequences (named <output_prefix>_primers.fasta). 
         Each primer's record header includes a unique group id, its length, number of degenerate positions,
         number of targets it covers, and the list of target IDs.
      3. A TSV file (<output_prefix>_primer_assignments.tsv) enumerating each target, the group it belongs to,
         and the final primer sequence.
    """
    # Determine the directory of the input FASTA and prepare output paths.
    input_dir = os.path.dirname(os.path.abspath(args.input_fasta))
    output_prefix = os.path.join(input_dir, args.output_prefix)
    groups_dir = output_prefix + "_groups"
    os.makedirs(groups_dir, exist_ok=True)
    
    # If orientation is "rev", update primer sequences to be reverse complements.
    if args.orientation == "rev":
        for group in groups:
            group['primer'] = str(Seq(group['primer']).reverse_complement())
    
    # Prepare primer FASTA records and TSV data.
    primer_records = []
    tsv_lines = ["Target\tGroup\tPrimer"]
    
    for i, group in enumerate(groups, start=1):
        group_id = f"primer{i}"
        num_targets = len(group['members'])
        length = group['length']
        degenerate_count = group['degenerate_count']
        target_list = ",".join(group['members'])
        
        # Create a SeqRecord for the primer with extended header information.
        rec = SeqRecord(Seq(group['primer']),
                        id=group_id,
                        description=f"len:{length} degenerate:{degenerate_count} targets:{num_targets} covers:{target_list}")
        primer_records.append(rec)
        
        # Write group FASTA file with the target sequences.
        group_filename = os.path.join(groups_dir, f"{group_id}.fasta")
        group_records = []
        for tid in group['members']:
            target_record = SeqRecord(Seq(sequences[tid]),
                                      id=tid,
                                      description="")
            group_records.append(target_record)
            tsv_lines.append(f"{tid}\t{group_id}\t{group['primer']}")
        SeqIO.write(group_records, group_filename, "fasta")
    
    # Write the primer FASTA file.
    primer_fasta_file = output_prefix + "_primers.fasta"
    SeqIO.write(primer_records, primer_fasta_file, "fasta")
    
    # Write the TSV file.
    tsv_file = output_prefix + "_primer_assignments.tsv"
    with open(tsv_file, "w") as tf:
        tf.write("\n".join(tsv_lines) + "\n")

def main():
    """
    Main function:
      - Parses arguments.
      - Reads input sequences from the FASTA file.
      - Designs primer groups using a greedy heuristic.
      - Writes output files (group FASTA files, primer FASTA, and a TSV assignment file).
    """
    args = parse_args()
    # Parse primer length range (e.g. "20-23")
    try:
        min_length, max_length = map(int, args.primer_length_range.split("-"))
    except Exception:
        print("Error parsing primer_length_range. Please provide in the format min-max (e.g. 20-23).")
        return

    sequences = parse_fasta(args.input_fasta)
    if not sequences:
        print("No sequences found in the input FASTA.")
        return

    groups = design_groups(sequences, args.target, args.target_offset, min_length, max_length, args.max_degenerate)
    write_outputs(groups, sequences, args)
    
    print(f"Designed {len(groups)} primer group(s).")
    print("Outputs written to the same directory as the input FASTA:")
    print(f" - Primer FASTA file: {os.path.join(os.path.dirname(os.path.abspath(args.input_fasta)), args.output_prefix + '_primers.fasta')}")
    print(f" - TSV assignments: {os.path.join(os.path.dirname(os.path.abspath(args.input_fasta)), args.output_prefix + '_primer_assignments.tsv')}")
    print(f" - Group FASTA files are in: {os.path.join(os.path.dirname(os.path.abspath(args.input_fasta)), args.output_prefix + '_groups')}")

if __name__ == "__main__":
    main()
