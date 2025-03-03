#!/usr/bin/env python3
import sys
import os

def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence.
    """
    # Create a translation table for nucleotide complement
    trans = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq[::-1].translate(trans)

def process_fasta(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        header = None
        seq_lines = []
        for line in infile:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # If there's an existing record, process and write it out
                if header is not None:
                    sequence = "".join(seq_lines)
                    revcomp_seq = reverse_complement(sequence)
                    # Append 'revcomp' to the description line
                    new_header = header + " revcomp"
                    outfile.write(new_header + "\n")
                    outfile.write(revcomp_seq + "\n")
                # Start a new record
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        # Process the final record
        if header is not None:
            sequence = "".join(seq_lines)
            revcomp_seq = reverse_complement(sequence)
            new_header = header + " revcomp"
            outfile.write(new_header + "\n")
            outfile.write(revcomp_seq + "\n")

def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: python revcomp.py <input_fasta_file>")
    input_file = sys.argv[1]
    if not os.path.isfile(input_file):
        sys.exit(f"Error: File '{input_file}' does not exist.")
    # Build output filename by inserting _revcomp before the extension.
    base, ext = os.path.splitext(input_file)
    output_file = base + "_revcomp.fasta"
    process_fasta(input_file, output_file)
    print(f"Reverse complement FASTA written to: {output_file}")

if __name__ == "__main__":
    main()
