#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script groups DNA sequences from an input FASTA file by their primer 
region, extracted from either the beginning or the end of each sequence. The 
user specifies an allowable primer size range (e.g. 20-23) and a target region 
("start" or "end"), and the script selects the primer length within that range 
that minimizes the number of distinct primer groups (thereby reducing the total
number of primers required). The output consists of separate FASTA files—named 
partly after the input file—for each group of sequences sharing an identical 
primer region.

Use group_primer_targets.py --help for instructions
"""

import argparse
import os
from Bio import SeqIO

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Group sequences by primer regions based on user-specified parameters.")
    parser.add_argument("fasta_file", help="Input FASTA file containing DNA sequences.")
    parser.add_argument("--region", choices=["start", "end"], required=True,
                        help="Region of the sequence to target for primers ('start' or 'end').")
    parser.add_argument("--primer_range", required=True,
                        help="Allowable primer sizes as a range (e.g. 20-23).")
    return parser.parse_args()

def get_primer_groups(records, region, primer_length):
    """Return a dictionary mapping primer sequences to lists of records based on the specified region and primer length."""
    groups = {}
    for record in records:
        if len(record.seq) < primer_length:
            # Skip sequences that are shorter than the required primer length.
            continue
        primer = record.seq[:primer_length] if region == "start" else record.seq[-primer_length:]
        primer = str(primer).upper()
        groups.setdefault(primer, []).append(record)
    return groups

def choose_optimal_primer_length(records, region, min_length, max_length):
    """Return the primer length within the given range that minimizes the number of groups, along with its grouping."""
    best_length = None
    best_groups = None
    best_group_count = None
    for L in range(min_length, max_length + 1):
        groups = get_primer_groups(records, region, L)
        group_count = len(groups)
        if best_group_count is None or group_count < best_group_count:
            best_length = L
            best_groups = groups
            best_group_count = group_count
    return best_length, best_groups

def write_fasta_files(groups, input_basename, primer_length, region):
    """Write each group of sequences to a separate FASTA file named based on the input file and primer region."""
    for primer, records in groups.items():
        out_filename = f"{input_basename}_{region}_{primer_length}bp_{primer}.fasta"
        with open(out_filename, "w") as out_handle:
            SeqIO.write(records, out_handle, "fasta")
        print(f"Wrote {len(records)} records to {out_filename}")

def main():
    """Main function to group sequences by primer regions and write output FASTA files."""
    args = parse_arguments()
    try:
        min_length, max_length = map(int, args.primer_range.split("-"))
    except ValueError:
        print("Error: primer_range must be in the format 'min-max' (e.g. 20-23).")
        return

    records = list(SeqIO.parse(args.fasta_file, "fasta"))
    if not records:
        print("No records found in the input FASTA file.")
        return

    optimal_length, groups = choose_optimal_primer_length(records, args.region, min_length, max_length)
    print(f"Optimal primer length chosen: {optimal_length} (resulting in {len(groups)} groups)")
    
    input_basename = os.path.splitext(os.path.basename(args.fasta_file))[0]
    write_fasta_files(groups, input_basename, optimal_length, args.region)

if __name__ == "__main__":
    main()
