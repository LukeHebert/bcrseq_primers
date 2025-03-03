## Overview
This repository holds script(s) used in the generation of OE-RT-PCR primers for any organism for which the user has heavy and light chain V and C gene sequences. Here, ferret sequences are used as example input and outputs, contained in the subdirectories.

Important note: these scripts cannot be used alone for primer generation. They serve as a first step. The "overlap" subsequences of primers required for overlap-extension-RT-PCR must be added, then melting temperatures and primer compatibility must be simulated with e.g. [IDT tools](https://www.idtdna.com/pages/tools/primerquest?utm_source=google&utm_medium=cpc&utm_campaign=00588_1i_03&utm_content=search&gad_source=1&gclid=Cj0KCQiA2oW-BhC2ARIsADSIAWqLLeU8xxYsdb8FC8H3mHMXEwgz3zivyRlAaNtvaHjE65R4VGVlCs8aAobeEALw_wcB). Finally, primers must be tested to the user's satisfaction *in vitro* to ensure amplicon length and depth of gene capture.

## Background

"BCRseq" or B cell receptor sequencing traditionally involves the isolation of bulk mRNA from peripheral blood mononuclear cells (which are ~10% B cells in humans) followed by the reverse transcription & ultimate high-throughput sequencing of only BCR mRNA. Authors of the two publications linked below pioneered a method to not only acquire heavy chain or light chain BCR transcript sequences, but capture the natively paired transcripts of both chains in one sequence. 
- [McDaniel JR, DeKosky BJ, Tanno H, Ellington AD, Georgiou G. Ultra-high-throughput sequencing of the immune receptor repertoire from millions of lymphocytes. Nat Protoc. 2016 Mar;11(3):429-42. doi: 10.1038/nprot.2016.024. Epub 2016 Feb 4. PMID: 26844430.](https://pubmed.ncbi.nlm.nih.gov/26844430/)
- [Tanno H, McDaniel JR, Stevens CA, Voss WN, Li J, Durrett R, Lee J, Gollihar J, Tanno Y, Delidakis G, Pothukuchy A, Ellefson JW, Goronzy JJ, Maynard JA, Ellington AD, Ippolito GC, Georgiou G. A facile technology for the high-throughput sequencing of the paired VH:VL and TCRβ:TCRα repertoires. Sci Adv. 2020 Apr 22;6(17):eaay9093. doi: 10.1126/sciadv.aay9093. PMID: 32426460; PMCID: PMC7176429.](https://pmc.ncbi.nlm.nih.gov/articles/PMC7176429/)

These papers focused on human BCR genes. For myriad reasons including the application of BCRseq for biomedical model organisms such as the ferret (the example data used here), this technique would be useful in non-human organisms. This repository's scripts aid in the initial design of overlap-extension reverse-transcription PCR (OE-RT-PCR) as described in the papers above.


## Contents
- `generate_primers.py`: This script takes a FASTA file of some target V gene sequences (e.g. heavy chain V gene framework 1 sequences) and a handful of user provided parameters such as range of allowable primer length, number of allowable degenerate nucleotide positions, etc. It outputs suggested primers and their corresponding target sequences.

- `bad_c_areas.py`: This script takes a FASTA file of all target C gene sequences and searches for areas to avoid priming due to perfect shared homology between two or more C genes.

- Four-letter directories such as `IGHV` include both the necessary input files for designing ferret primers (e.g. `ighv_fam1.fasta`) as well as the output from `generate_primers.py`

- `all_Cs`: Includes the relevant input ferret C gene sequences tested for undesirable priming areas using `bad_c_areas.py` as well as the output and resulting manually chosen C gene primer sequence choies.

- `revcomp.py`: Small script that simply takes an input FASTA and creates a new FASTA with reverse compliment versions of the input's sequences. Can be helpful when manually designing C gene primers.

- `environment.yml`: Conda environment i.e. dependencies file for ease of setup

## Helpful resources

- **Getting target gene sequences for your organism**: Assuming your organism's immunoglobulin genes have been sufficiently annotated. [IMGT's gene lookup tool](https://www.imgt.org/genedb/) is highly useful for acquiring the necessary V and C gene sequences/regions used as input for this repository's script(s). 

- **Double checking primer targets**: The [Clustal Omega multiple sequence alignment tool](https://www.ebi.ac.uk/jdispatcher/msa/clustalo) is useful for a quick double-check that the primers suggested by `generate_primers.py` do target the user's desired location, and can provide information that informed any degenerate primer positions. 

- **Simulating primer behavior/issues**: [IDT](https://www.idtdna.com/pages/tools/primerquest?utm_source=google&utm_medium=cpc&utm_campaign=00588_1i_03&utm_content=search&gad_source=1&gclid=Cj0KCQiA2oW-BhC2ARIsADSIAWqJSldnoqmAcMCrlltLo1xX62u5B97-xaWZbliHxkra0yiwiJ2StboaAhzLEALw_wcB) provides a number of tools for analyzing primer melting temperatures, undesirable secondary structures, etc. for multiplex PCR primer sets.