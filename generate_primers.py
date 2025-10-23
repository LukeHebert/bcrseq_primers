#!/usr/bin/env python3
"""
Design degenerate primers that cover a set of equal-length DNA targets using as
few primer groups as possible (greedy set-cover heuristic).
"""

import argparse
import os
from collections import OrderedDict
from typing import Dict, List, Tuple, Iterable, Optional
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
    parser.add_argument("input_fasta", help="Input FASTA file with target sequences (all same length).")
    parser.add_argument("--primer_length_range", required=True,
                        help="Allowed primer lengths, e.g. '20-23'")
    parser.add_argument("--max_degenerate", type=int, required=True,
                        help="Maximum number of degenerate positions allowed in any primer")
    parser.add_argument("--prime_side", choices=["start", "end"], required=True,
                        help="Which end of the target the primer anneals to ('start' = 5′ end, 'end' = 3′ end of the target)")
    parser.add_argument("--prime_offset", type=int, default=0,
                        help="Number of bases to skip from the chosen end before primer design "
                             "(1-indexed to the user; e.g. --prime_side start --prime_offset 2 "
                             "ignores positions 1 and 2).")
    parser.add_argument("--prime_wiggle", type=int, default=0,
                        help="How far the primer may be shifted inward from the offset anchor.")
    parser.add_argument("--orientation", choices=["fwd", "rev"], default="fwd",
                        help="Primer orientation to output: 'fwd' (default) or reverse complement 'rev'")
    parser.add_argument("--nondegenerate_tail", type=int, default=0,
                        help="Require the last X (3′) positions of the **primer** to be non-degenerate (A/C/G/T only).")
    parser.add_argument("--output_prefix", default="output",
                        help="Prefix for output files (placed in the input FASTA directory)")
    return parser.parse_args()


def parse_fasta(filename: str) -> Dict[str, str]:
    """
    Read sequences from a FASTA file into an ordered dict of {id: sequence}.
    """
    sequences = OrderedDict()
    for record in SeqIO.parse(filename, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences


def get_candidate_region(seq: str, side: str, offset_skip: int, wiggle: int, length: int) -> Optional[str]:
    """
    Extract a candidate region to prime against given:
      - side: which end of the target (start or end),
      - offset_skip: how many bases to skip from that end,
      - wiggle: extra inward shift from the offset,
      - length: primer length.
    Returns None if the window would fall outside the sequence.
    """
    if side == "start":
        # 0-based index where the candidate window begins
        start_idx = offset_skip + wiggle
        end_idx = start_idx + length
        if end_idx <= len(seq):
            return seq[start_idx:end_idx]
        return None
    else:  # side == "end"
        # end_idx is one past the last base included
        end_idx = len(seq) - offset_skip - wiggle
        start_idx = end_idx - length
        if start_idx >= 0:
            return seq[start_idx:end_idx]
        return None


# IUPAC code map for degenerate nucleotides
_IUPAC_MAP = {
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


def degenerate_code(nucs: Iterable[str]) -> str:
    """Return the IUPAC degenerate code for a set of nucleotides."""
    return _IUPAC_MAP.get(frozenset(nucs), 'N')


def compute_consensus(seqs: List[str]) -> Tuple[str, int]:
    """
    Given equal-length sequences, compute a column-wise IUPAC consensus string
    and return (consensus, degenerate_count).
    """
    if not seqs:
        return "", 0
    length = len(seqs[0])
    consensus_chars: List[str] = []
    degenerate_count = 0
    for i in range(length):
        col = {s[i] for s in seqs}
        if len(col) == 1:
            consensus_chars.append(next(iter(col)))
        else:
            degenerate_count += 1
            consensus_chars.append(degenerate_code(col))
    return ''.join(consensus_chars), degenerate_count


def has_clean_3prime_tail(consensus_forward: str, tail_len: int, orientation: str) -> bool:
    """
    Enforce a non-degenerate 3′ tail of length `tail_len` on the ACTUAL primer
    that will be used:
      - For forward primers, check the 3′ end of the consensus as-is.
      - For reverse primers, check the 3′ end of the reverse-complement.
    """
    if tail_len <= 0:
        return True
    # The primer sequence used depends on orientation
    primer_seq = (consensus_forward if orientation == "fwd"
                  else str(Seq(consensus_forward).reverse_complement()))
    tail_len = min(tail_len, len(primer_seq))
    return all(b in "ACGT" for b in primer_seq[-tail_len:])


def iter_design_space(allowed_wiggle: int,
                      min_len: int,
                      max_len: int,
                      seeds: Iterable[str]) -> Iterable[Tuple[int, int, str]]:
    """
    Yield all (wiggle, length, seed_id) combinations explored by the algorithm.
    This replaces three nested loops with a single, readable iterator.
    """
    for wiggle in range(allowed_wiggle + 1):
        for length in range(min_len, max_len + 1):
            for seed in seeds:
                yield (wiggle, length, seed)


def design_groups(sequences: Dict[str, str],
                  side: str,
                  offset_skip: int,
                  allowed_wiggle: int,
                  min_length: int,
                  max_length: int,
                  max_degenerate: int,
                  tail_len: int,
                  orientation: str) -> List[Dict]:
    """
    Greedily form groups of target sequences that can share a primer.

    We repeatedly:
      1) explore candidate design settings (wiggle, length, seed),
      2) grow a seed-based group by adding sequences whose combined consensus
         remains within `max_degenerate` and passes the non-degenerate 3′ tail
         rule in the ACTUAL primer orientation,
      3) select the candidate that covers the most sequences (breaking ties
         by preferring longer primers), then remove those sequences and repeat.
    """
    ungrouped = set(sequences.keys())
    groups: List[Dict] = []

    while ungrouped:
        best_group_map: Optional[Dict[str, str]] = None  # {seq_id: candidate_window}
        best_wiggle: Optional[int] = None
        best_len: Optional[int] = None
        best_deg: Optional[int] = None
        best_consensus: Optional[str] = None

        # Single loop over a generator instead of three nested loops
        for wiggle, length, seed in iter_design_space(
                allowed_wiggle, min_length, max_length, list(ungrouped)):

            # Candidate window for the seed sequence
            seed_window = get_candidate_region(sequences[seed], side, offset_skip, wiggle, length)
            if seed_window is None:
                continue

            # Start a tentative group with the seed
            candidate_map: Dict[str, str] = {seed: seed_window}

            # Try to greedily add other ungrouped sequences that keep consensus valid
            for other in list(ungrouped):
                if other in candidate_map:
                    continue
                other_window = get_candidate_region(sequences[other], side, offset_skip, wiggle, length)
                if other_window is None:
                    continue

                tmp = list(candidate_map.values()) + [other_window]
                consensus, deg = compute_consensus(tmp)

                # Tail rule must apply to the ACTUAL primer orientation.
                if deg <= max_degenerate and has_clean_3prime_tail(consensus, tail_len, orientation):
                    candidate_map[other] = other_window

            # Evaluate this candidate
            size = len(candidate_map)
            if best_group_map is None or size > len(best_group_map) or (size == len(best_group_map) and length > (best_len or 0)):
                best_group_map = candidate_map.copy()
                best_consensus, best_deg = compute_consensus(list(best_group_map.values()))
                best_wiggle = wiggle
                best_len = length

        # Fallback: if nothing worked (should be rare), make a singleton with a reasonable window
        if best_group_map is None:
            seed = next(iter(ungrouped))
            cand = get_candidate_region(sequences[seed], side, offset_skip, 0, max_length)
            if cand is None:
                cand = sequences[seed]
            best_group_map = {seed: cand}
            best_consensus, best_deg = compute_consensus([cand])
            best_wiggle = 0
            best_len = max_length

        groups.append({
            "primer": best_consensus,            # consensus in forward orientation
            "offset": best_wiggle,               # wiggle actually used
            "length": best_len,
            "members": list(best_group_map.keys()),
            "degenerate_count": best_deg
        })

        # Remove assigned targets
        for sid in best_group_map:
            ungrouped.discard(sid)

    return groups


def write_outputs(groups: List[Dict], sequences: Dict[str, str], args) -> None:
    """
    Write:
      - one FASTA per target sequence group,
      - a primer FASTA,
      - a TSV map of target -> group -> primer.
    """
    input_dir = os.path.dirname(os.path.abspath(args.input_fasta))
    prefix = os.path.join(input_dir, args.output_prefix)
    grp_dir = prefix + "_groups"
    os.makedirs(grp_dir, exist_ok=True)

    # Convert stored forward-consensus primers to the requested orientation on output
    if args.orientation == "rev":
        for g in groups:
            g["primer"] = str(Seq(g["primer"]).reverse_complement())

    primer_records: List[SeqRecord] = []
    tsv_lines = ["Target\tGroup\tPrimer"]

    for idx, g in enumerate(groups, start=1):
        gid = f"primer{idx}"
        rec = SeqRecord(
            Seq(g["primer"]),
            id=gid,
            description=f"len:{g['length']} degenerate:{g['degenerate_count']} "
                        f"targets:{len(g['members'])} covers:{','.join(g['members'])}"
        )
        primer_records.append(rec)

        grp_file = os.path.join(grp_dir, f"{gid}.fasta")
        grp_records = [SeqRecord(Seq(sequences[tid]), id=tid, description="") for tid in g["members"]]
        SeqIO.write(grp_records, grp_file, "fasta")

        for tid in g["members"]:
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

    groups = design_groups(
        sequences=sequences,
        side=args.prime_side,
        offset_skip=args.prime_offset,
        allowed_wiggle=args.prime_wiggle,
        min_length=min_len,
        max_length=max_len,
        max_degenerate=args.max_degenerate,
        tail_len=args.nondegenerate_tail,
        orientation=args.orientation
    )

    write_outputs(groups, sequences, args)

    base = os.path.dirname(os.path.abspath(args.input_fasta))
    print(f"Designed {len(groups)} primer group(s).")
    print("Outputs written to the same directory as the input FASTA:")
    print(f" - Primer FASTA : {os.path.join(base, args.output_prefix + '_primers.fasta')}")
    print(f" - TSV mapping  : {os.path.join(base, args.output_prefix + '_primer_assignments.tsv')}")
    print(f" - Group FASTAs : {os.path.join(base, args.output_prefix + '_groups')}")


if __name__ == "__main__":
    main()
