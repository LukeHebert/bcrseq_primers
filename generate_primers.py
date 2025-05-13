#!/usr/bin/env python3
'''
Reads a single FASTA file of same-length DNA sequences. For these sequences, 
attempts to generate the minimum number of required primers that will target
all input sequences in a multiplexed PCR.

User specifices primer design parameters such as primer position (and position 
tolerances aka "wiggle"), length(s), and maximum allowed degenerate nuclotides.

For each allowed wiggle (0-primer_wiggle) and each allowed primer length the
script tries seeding a group with one unassigned target sequence & then greedily adds
other primable sequence targets if the consensus (computed column by column and
“degeneratized” via IUPAC codes when needed) does not exceed the allowed number
of degenerate positions and those degenerate positions are not near the 3' end.

It then picks the candidate primers that cover the most primable sequences
(and, when tied, prefers the longer primer) and removes those target seqs from
further consideration.

Finally, it writes
    – one FASTA file of the primer sequences (with a header that notes the
     the list of input sequences targeted), 
    – one FASTA file per input target sequence group primed by a shared primer  
    – a TSV file that lists each target, its group assignment, and its group's 
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
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Design degenerate primers covering target sequences with minimal primer groups."
    )
    parser.add_argument("input_fasta", help="Input FASTA file with target sequences")
    parser.add_argument("--primer_length_range", required=True,
                        help="Allowed primer lengths, e.g. '20-23'")
    parser.add_argument("--max_degenerate", type=int, required=True,
                        help="Maximum number of degenerate positions allowed in any primer")
    parser.add_argument("--prime_side", choices=["start", "end"], required=True,
                        help="Whether the primer should anneal to the 'start' or 'end' of the target sequences")
    parser.add_argument("--prime_offset", type=int, default=0,
                        help="Number of bases to skip from the chosen end before primer design "
                             "(1-indexed to the user; e.g. --prime_side start --prime_offset 2 "
                             "ignores positions 1 and 2).")
    parser.add_argument("--prime_wiggle", type=int, default=0,
                        help="How far the primer may be shifted inward from the offset anchor.")
    parser.add_argument("--orientation", choices=["fwd", "rev"], default="fwd",
                        help="Primer orientation: 'fwd' (default) for forward or 'rev' for reverse complement output")
    parser.add_argument("--nondegenerate_tail", type=int, default=0,
                        help="Require the last X (3′) positions to be non-degenerate (A/C/G/T only).")
    parser.add_argument("--output_prefix", default="output",
                        help="Prefix for output files (placed in the input FASTA directory)")
    return parser.parse_args()

def parse_fasta(filename):
    """
    Read sequences from a FASTA file.
    """
    sequences = OrderedDict()
    for record in SeqIO.parse(filename, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def get_candidate_region(seq, side, offset_skip, wiggle, length):
    """
    Extract a candidate region given a sequence, side of target to prime, an initial
    offset to skip, an additional inward wiggle, and primer length.

    offset_skip – bases excluded from the sequence end  
    wiggle      – additional shift (0-prime_wiggle) applied by the algorithm
    """
    if side == "start":
        pos = offset_skip + wiggle          # 0-based index of primer start
        if pos + length <= len(seq):
            return seq[pos:pos + length]
        return None
    elif side == "end":
        end_idx   = len(seq) - offset_skip - wiggle      # index after last primer base
        start_idx = end_idx - length
        if start_idx >= 0:
            return seq[start_idx:end_idx]
        return None
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
    Given a list of sequences (all the same length), compute the consensus
    sequence and count how many positions are degenerate.
    """
    if not seqs:
        return "", 0
    length = len(seqs[0])
    consensus = []
    degenerate_count = 0
    for i in range(length):
        col = {s[i] for s in seqs}
        if len(col) == 1:
            consensus.append(next(iter(col)))
        else:
            degenerate_count += 1
            consensus.append(degenerate_code(col))
    return ''.join(consensus), degenerate_count

def tail_is_clean(primer, tail_len):
    """
    Check whether the last `tail_len` bases of the primer are non-degenerate.
    """
    if tail_len == 0:
        return True
    tail_len = min(tail_len, len(primer))
    return all(b in "ACGT" for b in primer[-tail_len:])

def design_groups(sequences, side, offset_skip, allowed_wiggle,
                  min_length, max_length, max_degenerate, tail_len):
    """
    Greedily form groups of target sequences that can share a primer.

    offset_skip     – user-specified --prime_offset  
    allowed_wiggle  – user-specified --prime_wiggle
    """
    ungrouped = set(sequences.keys())
    groups = []

    while ungrouped:
        best_group = None
        best_off   = None
        best_len   = None
        best_deg   = None

        # Try all allowed shifts (from 0 to allowed_wiggle) and primer lengths.
        for off in range(allowed_wiggle + 1):
            for length in range(min_length, max_length + 1):
                for seed in list(ungrouped):
                    cand_seed = get_candidate_region(sequences[seed], side,
                                                    offset_skip, off, length)
                    if cand_seed is None:
                        continue
                    current = {seed: cand_seed}
                    for other in list(ungrouped):
                        if other in current:
                            continue
                        cand_other = get_candidate_region(sequences[other], side,
                                                          offset_skip, off, length)
                        if cand_other is None:
                            continue
                        tmp = list(current.values()) + [cand_other]
                        consensus, deg = compute_consensus(tmp)
                        if deg <= max_degenerate and tail_is_clean(consensus, tail_len):
                            current[other] = cand_other
                    size = len(current)
                    if (best_group is None or size > len(best_group) or
                       (size == len(best_group) and length > best_len)):
                        best_group = current.copy()
                        best_cons, best_deg = compute_consensus(list(best_group.values()))
                        best_off  = off
                        best_len  = length

        if best_group is None:  # Fallback (should rarely happen)
            seed = next(iter(ungrouped))
            cand = get_candidate_region(sequences[seed], side,
                                        offset_skip, 0, max_length)
            if cand is None:
                cand = sequences[seed]
            best_group = {seed: cand}
            best_cons, best_deg = compute_consensus([cand])
            best_off  = 0
            best_len  = max_length

        groups.append({
            'primer'          : best_cons,
            'offset'          : best_off,      # wiggle actually used
            'length'          : best_len,
            'members'         : list(best_group.keys()),
            'degenerate_count': best_deg
        })
        for s in best_group:
            ungrouped.discard(s)
    return groups

def write_outputs(groups, sequences, args):
    """
    Write one FASTA per target sequence group, a primer FASTA, and a TSV map.
    """
    input_dir = os.path.dirname(os.path.abspath(args.input_fasta))
    prefix = os.path.join(input_dir, args.output_prefix)
    grp_dir = prefix + "_groups"
    os.makedirs(grp_dir, exist_ok=True)

    if args.orientation == "rev":
        for g in groups:
            g['primer'] = str(Seq(g['primer']).reverse_complement())

    primer_records = []
    tsv_lines = ["Target\tGroup\tPrimer"]

    for idx, g in enumerate(groups, start=1):
        gid = f"primer{idx}"
        rec = SeqRecord(Seq(g['primer']),
                        id=gid,
                        description=f"len:{g['length']} degenerate:{g['degenerate_count']} "
                                    f"targets:{len(g['members'])} covers:{','.join(g['members'])}")
        primer_records.append(rec)

        grp_file = os.path.join(grp_dir, f"{gid}.fasta")
        grp_records = [SeqRecord(Seq(sequences[tid]), id=tid, description="")
                       for tid in g['members']]
        SeqIO.write(grp_records, grp_file, "fasta")
        for tid in g['members']:
            tsv_lines.append(f"{tid}\t{gid}\t{g['primer']}")

    SeqIO.write(primer_records, prefix + "_primers.fasta", "fasta")
    with open(prefix + "_primer_assignments.tsv", "w") as tf:
        tf.write("\n".join(tsv_lines) + "\n")

def main():
    """
    Main function:
      – Parse arguments  
      – Read input FASTA  
      – Design primer groups  
      – Write outputs
    """
    args = parse_args()
    try:
        min_len, max_len = map(int, args.primer_length_range.split("-"))
    except Exception:
        print("Error parsing primer_length_range. Use min-max (e.g. 20-23).")
        return

    sequences = parse_fasta(args.input_fasta)
    if not sequences:
        print("No sequences found in the input FASTA.")
        return

    groups = design_groups(sequences, args.prime_side,
                           args.prime_offset, args.prime_wiggle,
                           min_len, max_len,
                           args.max_degenerate, args.nondegenerate_tail)
    write_outputs(groups, sequences, args)

    base = os.path.dirname(os.path.abspath(args.input_fasta))
    print(f"Designed {len(groups)} primer group(s).")
    print("Outputs written to the same directory as the input FASTA:")
    print(f" - Primer FASTA : {os.path.join(base, args.output_prefix + '_primers.fasta')}")
    print(f" - TSV mapping  : {os.path.join(base, args.output_prefix + '_primer_assignments.tsv')}")
    print(f" - Group FASTAs : {os.path.join(base, args.output_prefix + '_groups')}")

if __name__ == "__main__":
    main()
