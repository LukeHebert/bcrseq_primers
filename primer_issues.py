'''
This script evaluates an input FASTA file of candidate primer sequences for 
issues or ineffeciencies. Specifically it looks at:
- GC clamps
- Melting temperatures
- Self & inter-primer complementarity
- Nucleotide repeat runs

Based largely on ThermoFisher primer guidelines:
https://www.thermofisher.com/blog/behindthebench/pcr-primer-design-tips/

Use python primer_issues.py --help for instructions
'''

#!/usr/bin/env python3
import argparse
import os
import re
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

def resolve_iupac(seq):
    """
    Replace ambiguous IUPAC degenerate symbols in the sequence with a specific base.
    The mapping is defined arbitrarily:
      R -> A, Y -> T, S -> G, W -> A, K -> G, M -> A,
      B -> C, D -> A, H -> A, V -> A, N -> A.
    Returns the resolved sequence.
    """
    iupac_dict = {
        'R': 'A',
        'Y': 'T',
        'S': 'G',
        'W': 'A',
        'K': 'G',
        'M': 'A',
        'B': 'C',
        'D': 'A',
        'H': 'A',
        'V': 'A',
        'N': 'A'
    }
    resolved = []
    for base in seq:
        if base in iupac_dict:
            resolved.append(iupac_dict[base])
        else:
            resolved.append(base)
    return "".join(resolved)

def check_gc_clamp(seq):
    """
    Returns True if the primer has a sufficient GC clamp:
      - Ends with a G or C, OR
      - Has at least two G/C bases in its last six bases.
    """
    last_six = seq[-6:] if len(seq) >= 6 else seq
    return (seq[-1] in "GC") or (sum(1 for base in last_six if base in "GC") >= 2)

def compute_tm(seq):
    """Compute the melting temperature (Tm) using Tm_GC."""
    return mt.Tm_GC(seq)

def count_tm_neighbors(tm, all_tms):
    """
    Count how many primers have a melting temperature within ±5°C of tm.
    Returns a tuple: (count, total, percentage).
    """
    total = len(all_tms)
    count = sum(1 for t in all_tms if abs(t - tm) <= 5)
    percentage = count / total * 100
    return count, total, percentage

def find_self_homology(seq):
    """
    Check for self-homology by scanning all 3-base windows.
    If the reverse complement of a triplet is found at another position,
    record the triplet and its positions (1-indexed).
    Returns a list of tuples: (triplet, pos1, pos2).
    """
    matches = []
    for i in range(len(seq) - 2):
        triplet = seq[i:i+3]
        rev_comp_triplet = str(Seq(triplet).reverse_complement())
        for j in range(len(seq) - 2):
            if j == i:
                continue
            if seq[j:j+3] == rev_comp_triplet:
                match = (triplet, i + 1, j + 1)  # positions are 1-indexed
                if match not in matches:
                    matches.append(match)
    return matches

def check_interprimer_homology(seq, other_seqs):
    """
    Check if the given primer's sequence is the same as the reverse complement
    of any other primer. Returns a list of primer IDs (from other_seqs) that match.
    """
    homologous_with = []
    for other_id, other_seq in other_seqs:
        other_rev_comp = str(Seq(other_seq).reverse_complement())
        if seq == other_rev_comp:
            homologous_with.append(other_id)
    return homologous_with

def find_repeats(seq):
    """
    Look for:
      - Mononucleotide repeats: four or more of the same base.
      - Dinucleotide repeats: a two-base pattern repeated four or more times.
    Returns a list of tuples: (repeat_type, repeat_sequence, start, end)
    where start and end are 1-indexed positions.
    """
    repeats = []
    # Mononucleotide repeat: e.g., AAAA, CCCC, etc.
    mono_pattern = re.compile(r'([ACGT])\1{3,}')
    for match in mono_pattern.finditer(seq):
        start = match.start() + 1
        end = match.end()
        repeats.append(("mono", match.group(), start, end))
    
    # Dinucleotide repeat: e.g., ATATATAT (AT repeated 4 times)
    di_pattern = re.compile(r'(([ACGT]{2}))(\2){3,}')
    for match in di_pattern.finditer(seq):
        start = match.start() + 1
        end = match.end()
        repeats.append(("di", match.group(), start, end))
    
    return repeats

def plot_histogram(data, xlabel, ylabel, title, output_file):
    """Plot and save a histogram from the data."""
    plt.figure()
    plt.hist(data, bins=20)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.savefig(output_file, dpi=600)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Multiplex PCR Primer Evaluation")
    parser.add_argument("fasta", help="Input FASTA file with candidate primers")
    args = parser.parse_args()

    fasta_file = args.fasta
    parent_dir = os.path.dirname(os.path.abspath(fasta_file))
    base_name = os.path.splitext(os.path.basename(fasta_file))[0]
    # Create an output subdirectory named after the input FASTA file.
    output_dir = os.path.join(parent_dir, base_name)
    os.makedirs(output_dir, exist_ok=True)

    # Define output file names.
    output_txt = os.path.join(output_dir, f"{base_name}_analysis.txt")
    tm_hist_file = os.path.join(output_dir, f"{base_name}_tm_histogram.png")

    # Read primers from FASTA.
    primers = list(SeqIO.parse(fasta_file, "fasta"))
    if not primers:
        print("No primers found in the FASTA file.")
        return

    # Preprocess primers: create a list of tuples (id, original_sequence, resolved_sequence)
    primer_info = []
    for record in primers:
        original_seq = str(record.seq).upper()
        resolved_seq = resolve_iupac(original_seq)
        primer_info.append((record.id, original_seq, resolved_seq))
    
    # For analysis, use the resolved sequences.
    resolved_primer_info = [(pid, res) for pid, orig, res in primer_info]

    # Pre-calculate Tm for all primers using the resolved sequences.
    tm_values = {}
    for pid, orig, res in primer_info:
        tm_values[pid] = compute_tm(res)

    # Summary counters.
    no_gc_clamp_count = 0
    tm_outside_count = 0
    self_homology_count = 0
    interprimer_homology_count = 0
    repeats_count = 0

    # Prepare a list to hold analysis results.
    analysis_results = []

    for pid, original_seq, resolved_seq in primer_info:
        # GC Clamp.
        gc_clamp = check_gc_clamp(resolved_seq)
        if not gc_clamp:
            no_gc_clamp_count += 1

        # Melting Temperature.
        tm_val = tm_values[pid]
        tm_within_range = (65 <= tm_val <= 75)
        if not tm_within_range:
            tm_outside_count += 1

        similar_count, total, percentage = count_tm_neighbors(tm_val, list(tm_values.values()))

        # Self-homology.
        self_hom_matches = find_self_homology(resolved_seq)
        if self_hom_matches:
            self_homology_count += 1

        # Interprimer homology.
        other_info = [item for item in resolved_primer_info if item[0] != pid]
        interprimer_matches = check_interprimer_homology(resolved_seq, other_info)
        if interprimer_matches:
            interprimer_homology_count += 1

        # Repeats.
        repeats_found = find_repeats(resolved_seq)
        if repeats_found:
            repeats_count += 1

        analysis_results.append({
            "id": pid,
            "original_sequence": original_seq,
            "resolved_sequence": resolved_seq,
            "gc_clamp": gc_clamp,
            "tm": tm_val,
            "tm_similar_count": similar_count,
            "tm_total": total,
            "tm_similar_percentage": percentage,
            "self_homology": self_hom_matches,
            "interprimer_homology": interprimer_matches,
            "repeats": repeats_found
        })

    # Write summary and detailed analysis to the output text file.
    with open(output_txt, "w") as f:
        # Include user arguments at the top.
        summary_header = "User Arguments:\n"
        for arg, value in vars(args).items():
            summary_header += f"  {arg}: {value}\n"
        summary_header += "\n"
        
        summary = (
            f"Summary:\n"
            f"Total primers: {len(primers)}\n"
            f"Primers with no GC clamp: {no_gc_clamp_count}\n"
            f"Primers with melting temp outside 65-75°C: {tm_outside_count}\n"
            f"Primers with (3-mer) self-homology: {self_homology_count}\n"
            f"Primers with interprimer homology: {interprimer_homology_count}\n"
            f"Primers with one or more (4-mer or compound 4-mer) repeats: {repeats_count}\n\n"
        )
        f.write(summary_header)
        f.write(summary)
        
        for res in analysis_results:
            f.write(f"Primer ID: {res['id']}\n")
            f.write(f"Original Sequence: {res['original_sequence']}\n")
            f.write(f"IUPAC-resolved Sequence: {res['resolved_sequence']}\n")
            f.write(f"GC Clamp: {'Sufficient' if res['gc_clamp'] else 'Insufficient'}\n")
            f.write(f"Melting Temperature (Tm): {res['tm']:.2f} °C\n")
            f.write(f"Primers with Tm within ±5°C: {res['tm_similar_count']} / {res['tm_total']} ({res['tm_similar_percentage']:.2f}%)\n")
            
            if res['self_homology']:
                f.write("Self-homology (3-base complementary matches):\n")
                for match in res['self_homology']:
                    f.write(f"    Triplet '{match[0]}' found at positions {match[1]} and {match[2]}\n")
            else:
                f.write("No self-homology detected.\n")
            
            if res['interprimer_homology']:
                f.write("Interprimer homology with primer(s): " + ", ".join(res['interprimer_homology']) + "\n")
            else:
                f.write("No interprimer homology detected.\n")
            
            if res['repeats']:
                f.write("Repeats detected:\n")
                for rep in res['repeats']:
                    rep_type, rep_seq, start, end = rep
                    f.write(f"    {rep_type}-repeat: '{rep_seq}' from position {start} to {end}\n")
            else:
                f.write("No repeats detected.\n")
            f.write("\n")

    # Plot histogram for melting temperatures.
    tm_vals = list(tm_values.values())
    plot_histogram(tm_vals,
                   xlabel="Melting Temperature (°C)",
                   ylabel="Frequency",
                   title="Distribution of Melting Temperatures",
                   output_file=tm_hist_file)

    # Print file locations and summary to the terminal.
    print("Analysis complete.")
    print(f"Output text file: {output_txt}")
    print(f"Melting temperature histogram: {tm_hist_file}")
    print("\n" + summary)

if __name__ == "__main__":
    main()
