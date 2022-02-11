# BIOT670_COVIDAnalysisProject

import tempfile
import itertools
import numpy as np
import biotite.sequence as seq
import biotite.sequence.io.genbank as gb
import biotite.sequence.io.fasta as fasta
import biotite.database.entrez as entrez


# Get the E. coli wuhan genome as annotated sequence
gb_file = gb.GenBankFile.read(
    entrez.fetch("MN908947.3", tempfile.gettempdir(), "gb", "nuccore", "gb")
)
# We are only interested in CDS features
wuhan_genome = gb.get_annotated_sequence(gb_file, include_only=["CDS"])


# This dictionary will count how often each codon occurs in the genome
# For increased performance the dictionary uses symbol codes ([0 3 2])
# instead of symbols (['A' 'T' 'G']) as keys
codon_counter = {
    codon: 0 for codon
    in itertools.product( *([range(len(wuhan_genome.sequence.alphabet))] * 3) )
}
# For demonstration purposes print the 64 codons in symbol code form
print(list(codon_counter.keys()))

########################################################################
# As expected the dictionary encodes each codon as tuple of 3 numbers,
# where ``0`` represents ``A``, ``1`` ``C``, ``2`` ``G`` and ``3`` ``T``.
# These mappings are defined by the alphabet of the genomic sequence.
#
# In the next step the occurrences of codons in the coding sequences
# are counted, the relative frequencies are calculated and the codon
# table is printed.

# Iterate over all CDS features
for cds in wuhan_genome.annotation:
    # Obtain the sequence of the CDS
    cds_seq = wuhan_genome[cds]
    if len(cds_seq) % 3 != 0:
        # A CDS' length should be a multiple of 3,
        # otherwise the CDS is malformed
        continue
    # Iterate over the sequence in non-overlapping frames of 3
    # and count the occurence of each codon
    for i in range(0, len(cds_seq), 3):
        codon_code = tuple(cds_seq.code[i:i+3])
        codon_counter[codon_code] += 1

# Convert the total frequencies into relative frequencies
# for each amino acid
# The NCBI codon table with ID 11 is the bacterial codon table
table = seq.CodonTable.load(11)
# As the script uses symbol codes, each amino acid is represented by a
# number between 0 and 19, instead of the single letter code
for amino_acid_code in range(20):
    # Find all codons coding for the amino acid
    # The codons are also in symbol code format, e.g. ATG -> (0, 3, 2)
    codon_codes_for_aa = table[amino_acid_code]
    # Get the total amount of codon occurrences for the amino acid
    total = 0
    for codon_code in codon_codes_for_aa:
        total += codon_counter[codon_code]
    # Replace the total frequencies with relative frequencies
    # and print it
    for codon_code in codon_codes_for_aa:
        # Convert total frequencies into relative frequencies
        codon_counter[codon_code] /= total
        # The rest of this code block prints the codon usage table
        # in human readable format
        amino_acid = seq.ProteinSequence.alphabet.decode(amino_acid_code)
        # Convert the symbol codes for each codon into symbols...
        # ([0,3,2] -> ['A' 'T' 'G'])
        codon = wuhan_genome.sequence.alphabet.decode_multiple(codon_code)
        # ...and represent as string
        # (['A' 'T' 'G'] -> "ATG")
        codon = "".join(codon)
        freq = codon_counter[codon_code]
        print(f"{amino_acid}   {codon}   {freq:.2f}")
    print()

