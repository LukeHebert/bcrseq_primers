#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script groups DNA sequences from an input FASTA file by their primer region 
(either at the start or the end) and writes separate FASTA files for each target 
group, placing them in the same directory as the input file. The user specifies 
an allowable primer size range (e.g. 20-23) and a region ("start" or "end"), and 
the script determines an optimal primer length that minimizes the number of primer 
groups while outputting progress messages. In addition, all primers (the shared 
subsequences) are collected into a single FASTA file, and any groups that remain 
unchanged when extended to the maximum allowed length are “upgraded” to that 
longer primer.
"""

import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Group sequences by primer regions and output target groups and a combined primer FASTA file based on user parameters.")
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
        print(f"Trying primers of length {L}...")
        groups = get_primer_groups(records, region, L)
        group_count = len(groups)
        print(f"Found {group_count} target sequence groups.")
        if best_group_count is None or group_count < best_group_count:
            best_length = L
            best_groups = groups
            best_group_count = group_count
    return best_length, best_groups

def upgrade_primer_groups(groups, region, max_length):
    """
    Upgrade groups that remain unchanged when extended to max_length so that their primer is the max allowed length.
    Returns the upgraded groups and the count of groups that were upgraded.
    """
    upgraded_groups = {}
    upgrade_count = 0
    for primer, rec_list in groups.items():
        # Compute the set of max_length primers for all records in the group.
        extended = { str(rec.seq[:max_length] if region == "start" else rec.seq[-max_length:]).upper() for rec in rec_list }
        if len(extended) == 1:
            new_primer = extended.pop()
            # Count as an upgrade if the original primer is shorter than max_length and differs from the max_length primer.
            if len(primer) < max_length and new_primer != primer:
                upgrade_count += 1
            upgraded_groups[new_primer] = rec_list
        else:
            upgraded_groups[primer] = rec_list
    return upgraded_groups, upgrade_count

def write_output_files(groups, input_basename, region, output_dir):
    """
    For each primer group, write a FASTA file containing the target sequences and accumulate primer records.
    All target group files are written to the input file's directory and all primer records are collected into a single FASTA file.
    The target FASTA filenames now start with the number of contained sequences followed by an underscore.
    """
    primer_records = []
    group_counter = 1
    for primer, records in groups.items():
        effective_length = len(primer)
        group_count = len(records)
        # Construct the target group filename in the input file's directory with the count as prefix.
        target_filename = os.path.join(output_dir, f"{group_count}_{input_basename}_{region}_{effective_length}bp_{primer}.fasta")
        with open(target_filename, "w") as out_handle:
            SeqIO.write(records, out_handle, "fasta")
        print(f"Wrote {group_count} records to {target_filename}")
        
        # Create a primer SeqRecord with a unique ID and description mapping to the target file.
        primer_id = f"primer_{group_counter}"
        target_filename_base = os.path.basename(target_filename)
        primer_record = SeqRecord(Seq(primer),
                                  id=primer_id,
                                  description=f"targets file: {target_filename_base}")
        primer_records.append(primer_record)
        group_counter += 1

    # Write all primer records to a single FASTA file in the input file's directory.
    primers_filename = os.path.join(output_dir, f"{input_basename}_{region}_primers.fasta")
    with open(primers_filename, "w") as p_handle:
        SeqIO.write(primer_records, p_handle, "fasta")
    print(f"Wrote {len(primer_records)} primer records to {primers_filename}")

def main():
    """Main function to group sequences by primer regions, upgrade primer groups when possible, and write output FASTA files."""
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

    # Choose the primer length that minimizes the number of groups.
    optimal_length, groups_opt = choose_optimal_primer_length(records, args.region, min_length, max_length)
    print(f"Optimal primer length chosen: {optimal_length} (resulting in {len(groups_opt)} groups)")
    
    # Upgrade groups that remain unchanged when using the maximum allowed primer length.
    groups_final, upgrade_count = upgrade_primer_groups(groups_opt, args.region, max_length)
    if optimal_length < max_length:
        if upgrade_count:
            print(f"{upgrade_count} groups remain unchanged when testing primer length {optimal_length} as compared to primer length {max_length}. These groups are now upgraded to have primer lengths of {max_length} bp.")
        else:
            print(f"No groups remained unchanged when extending from primer length {optimal_length} to {max_length}.")
    else:
        print("Optimal primer length is already the maximum allowed; no upgrades performed.")

    # Determine the directory of the input file and the basename.
    output_dir = os.path.dirname(os.path.abspath(args.fasta_file))
    input_basename = os.path.splitext(os.path.basename(args.fasta_file))[0]
    write_output_files(groups_final, input_basename, args.region, output_dir)

if __name__ == "__main__":
    main()
