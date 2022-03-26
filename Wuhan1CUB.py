# =============================================================================
# Group 2 COVID Analysis Project
# BIOT 670, Spring 2022
# =============================================================================

# Please refer for calculations:
# Sharp, P. M., & Li, W.-H. (1987). Nucleic Acids Research, 15(3), 1281–1295. 
# https://doi.org/10.1093/nar/15.3.1281

# =============================================================================
# Imports
# =============================================================================
from Bio import SeqIO
from Bio.Seq import Seq
from math import exp, log
import pandas as pd
# =============================================================================


# Define amino acids for each codon
gencode = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

# Define synonymous_codons
synonymous_codons = {
"C": ["TGT", "TGC"],
"D": ["GAT", "GAC"],
"S": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"],
"Q": ["CAA", "CAG"],
"M": ["ATG"],
"N": ["AAC", "AAT"],
"P": ["CCT", "CCG", "CCA", "CCC"],
"K": ["AAG", "AAA"],
"_": ["TAG", "TGA", "TAA"],
"T": ["ACC", "ACA", "ACG", "ACT"],
"F": ["TTT", "TTC"],
"A": ["GCA", "GCC", "GCG", "GCT"],
"G": ["GGT", "GGG", "GGA", "GGC"],
"I": ["ATC", "ATA", "ATT"],
"L": ["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"],
"H": ["CAT", "CAC"],
"R": ["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"],
"W": ["TGG"],
"V": ["GTA", "GTC", "GTG", "GTT"],
"E": ["GAG", "GAA"],
"Y": ["TAT", "TAC"]}

def gb_parser(file):
# Return biopython genbank record from a genbank formatted file

    gb_record = SeqIO.read(open(file,"r"), "genbank")
    return gb_record


def in_locs(gene_locs, test):
# Determine if a locus (test) falls within a set of gene loci

    found = False
    for i in gene_locs:
        if test in i:
            found = True
    return found


def genes(gb_record):
# From a genbank record, return a list of all CDS names, loci, and sequences

    names = []
    locs = []
    seqs = []
    for feature in gb_record.features:
        # Extract gene names, loci, and sequences from gb record
        if feature.type == "CDS":
            try:
                name = feature.qualifiers.get("gene")[0]
            except TypeError:
                # If gb record does not include gene names, use product name
                name = feature.qualifiers.get("product")[0]
            if name in names:
                # Avoid duplicate entries
                continue
            if (in_locs(locs, feature.location.start) == True and in_locs(locs, 
                feature.location.end) == True):
                    # Avoid duplicate sequences (i.e. ORF1a & ORF1ab)
                    continue
            names.append(name)
            locs.append(feature.location)
            seqs.append(feature.location.extract(gb_record).seq)

    # Assemble list of lists containing all names, loci, and sequences
    genes_all = [names, locs, seqs]
    return genes_all


def cds_codons(cds_list):
# Accept a list of CDS sequences and return as a list of codons

    cds_codons = []
    for cds in cds_list:
        for i in range(0,len(cds),3):
            cds_codons.append(cds[i:i+3])
    return cds_codons


def cds_divide(cds_codons, x):
# Accept a list of codons, return as a list of sequence fragments of x length

    cds_fragments = []
    i=0

    # Split input list of codons into lists of x codons each
    while i < len(cds_codons):
        cds_fragments.append(cds_codons[i:i + x])
        i += x

    # Join codons into sequences
    for f in range(0,(len(cds_fragments))):
        cds_fragments[f] = Seq("").join(cds_fragments[f])
            
    return cds_fragments    


def codon_count(seq):
# Accept a sequence and return codon counts as a dictionary

    codon_dict = {
    "TTT": 0, "TTC": 0, "TTA": 0, "TTG": 0,
    "CTT": 0, "CTC": 0, "CTA": 0, "CTG": 0,
    "ATT": 0, "ATC": 0, "ATA": 0, "ATG": 0,
    "GTT": 0, "GTC": 0, "GTA": 0, "GTG": 0,
    "TAT": 0, "TAC": 0, "TAA": 0, "TAG": 0,
    "CAT": 0, "CAC": 0, "CAA": 0, "CAG": 0,
    "AAT": 0, "AAC": 0, "AAA": 0, "AAG": 0,
    "GAT": 0, "GAC": 0, "GAA": 0, "GAG": 0,
    "TCT": 0, "TCC": 0, "TCA": 0, "TCG": 0,
    "CCT": 0, "CCC": 0, "CCA": 0, "CCG": 0,
    "ACT": 0, "ACC": 0, "ACA": 0, "ACG": 0,
    "GCT": 0, "GCC": 0, "GCA": 0, "GCG": 0,
    "TGT": 0, "TGC": 0, "TGA": 0, "TGG": 0,
    "CGT": 0, "CGC": 0, "CGA": 0, "CGG": 0,
    "AGT": 0, "AGC": 0, "AGA": 0, "AGG": 0,
    "GGT": 0, "GGC": 0, "GGA": 0, "GGG": 0}
    
    for i in range(0,len(seq),3):
        codon_dict[seq[i:i+3]] += 1            
    return codon_dict


def codon_aafreq(codon_counts):
# Calculate codon frequencies per amino acid from dictionary of codon counts

    codon_aafreq = {
    "TTT": 0, "TTC": 0, "TTA": 0, "TTG": 0,
    "CTT": 0, "CTC": 0, "CTA": 0, "CTG": 0,
    "ATT": 0, "ATC": 0, "ATA": 0, "ATG": 0,
    "GTT": 0, "GTC": 0, "GTA": 0, "GTG": 0,
    "TAT": 0, "TAC": 0, "TAA": 0, "TAG": 0,
    "CAT": 0, "CAC": 0, "CAA": 0, "CAG": 0,
    "AAT": 0, "AAC": 0, "AAA": 0, "AAG": 0,
    "GAT": 0, "GAC": 0, "GAA": 0, "GAG": 0,
    "TCT": 0, "TCC": 0, "TCA": 0, "TCG": 0,
    "CCT": 0, "CCC": 0, "CCA": 0, "CCG": 0,
    "ACT": 0, "ACC": 0, "ACA": 0, "ACG": 0,
    "GCT": 0, "GCC": 0, "GCA": 0, "GCG": 0,
    "TGT": 0, "TGC": 0, "TGA": 0, "TGG": 0,
    "CGT": 0, "CGC": 0, "CGA": 0, "CGG": 0,
    "AGT": 0, "AGC": 0, "AGA": 0, "AGG": 0,
    "GGT": 0, "GGC": 0, "GGA": 0, "GGG": 0}
    
    for codon in codon_counts:
        # Find total counts for all synonymous codons
        sumni = 0
        for aa in synonymous_codons:
            if codon in synonymous_codons[aa]:
                for i in synonymous_codons[aa]:
                    sumni += codon_counts[i]
        # Calculate frequency for each synonymous codon
        if sumni > 0:
            codon_aafreq[codon] = codon_counts[codon] / sumni
        else:
            # Avoid divide by zero error
            codon_aafreq[codon] = 0

    return codon_aafreq


def get_RSCU(codon_freq):
# Calculate RSCU values from a dictionary of codon frequencies

    rscu = {}
    for val in synonymous_codons.values():
        # Calculate expected frequency if equal usage among synonymous codons
        exp_freq = 1/len(val)

        # Calculate rscu (observed/expected)
        for i in val:
            rscu[i] = codon_freq[i]/exp_freq

    return rscu


def codon_to_num(codon_dict):
# Accept a dictionary with codons as keys and convert the codons to numbers

    # Define numerical values for codons
    codon_num = {
    'ATA':1, 'ATC':2, 'ATT':3, 'ATG':4,
    'ACA':5, 'ACC':6, 'ACG':7, 'ACT':8,
    'AAC':9, 'AAT':10, 'AAA':11, 'AAG':12,
    'AGC':13, 'AGT':14, 'AGA':15, 'AGG':16,
    'CTA':17, 'CTC':18, 'CTG':19, 'CTT':20,
    'CCA':21, 'CCC':22, 'CCG':23, 'CCT':24,
    'CAC':25, 'CAT':26, 'CAA':27, 'CAG':28,
    'CGA':29, 'CGC':30, 'CGG':31, 'CGT':32,
    'GTA':33, 'GTC':34, 'GTG':35, 'GTT':36,
    'GCA':37, 'GCC':38, 'GCG':39, 'GCT':40,
    'GAC':41, 'GAT':42, 'GAA':43, 'GAG':44,
    'GGA':45, 'GGC':46, 'GGG':47, 'GGT':48,
    'TCA':49, 'TCC':50, 'TCG':51, 'TCT':52,
    'TTC':53, 'TTT':54, 'TTA':55, 'TTG':56,
    'TAC':57, 'TAT':58, 'TAA':59, 'TAG':60,
    'TGC':61, 'TGT':62, 'TGA':63, 'TGG':64}

    num_dict = {}
    for key in codon_dict.keys():
        num_dict[codon_num[key]] = codon_dict[key]

    return num_dict


def file_writer(any_dict, name):
# Write dictionary items into a file

    with open(name + ".csv","w") as f:
        for key, val in any_dict.items():
            f.write(str(key) + "," + str(round(val,3)) + "\n")
    print(name + ".csv is written!")
    
    
def file_writer2(sample_frags, name):
# Write multiple dictionary items into the same file

    key_list=[]
    freq_list=[]
    for key, val in sample_frags.items():
        key_list.append(key) #alternate 1
        freq_list.append(val)
    
    df = pd.DataFrame(freq_list, index=key_list)
    df.to_csv(name + ".csv", sep=",")
    print (name + ".csv is written!")


# =============================================================================

def main():
# From a gb file, calculate RSCU for whole genome and for each gene, output to
# csv. Option to also use fragments of arbitrary size (commented out).
    
    # Get genome sequence from gb file
    file_name = input("Please enter your GB filename: ")


    # Parse genbank file
    sample_gb_record = gb_parser(file_name)
    sample_genes = genes(sample_gb_record)
    sample_codons = cds_codons(sample_genes[2])


    # Analyze whole genome      
    counts_genome = codon_count(Seq("").join(sample_genes[2]))
    freq_genome = codon_aafreq(counts_genome)
    rscu_genome = get_RSCU(freq_genome)


    # Analyze by gene
    freq_genes = {}
    rscu_genes = {}
    for i in range(0,len(sample_genes[0])):
        name = sample_genes[0][i]
        count = codon_count(sample_genes[2][i])
        freq = codon_aafreq(count)
        rscu = get_RSCU(freq)
        freq_genes[name] = freq
        rscu_genes[name] = rscu

  
##    # Analyze by fragments of any number of codons
##    frag_number = int(input("To divide genome into regions of arbitrary size, "
##                            "enter the desired number of codons per region: "))
##    # Divide genome into fragments and calculate RSCU
##    sample_fragments = cds_divide(sample_codons, frag_number)
##    freq_frags = {}
##    rscu_frags = {}
##    pos = 1
##    for i in range(0,len(sample_fragments)):
##        count = codon_count(sample_fragments[i])
##        freq = codon_aafreq(count)
##        rscu = get_RSCU(freq)
##        freq_frags[str(pos)+"-"+str(pos+frag_number-1)] = freq
##        rscu_frags[str(pos)+"-"+str(pos+frag_number-1)] = rscu
##        pos += frag_number


    # Option to convert codons to numbers
    convert = input("Convert codons to numeric values? Enter Y/N: ")
    if convert.lower() == "y":
        rscu_genome = codon_to_num(rscu_genome)
        for key in rscu_genes.keys():
            rscu_genes[key] = codon_to_num(rscu_genes[key])
        for key in rscu_frags.keys():
            rscu_frags[key] = codon_to_num(rscu_frags[key])


    # Write genome rscu table to csv file
    file_writer(rscu_genome, file_name + "_rscu_genome")

    # Write gene rscu table to csv file
    file_writer2(rscu_genes, file_name + "_rscu_genes")

##    # Write fragment rscu table to csv file
##    file_writer2(rscu_frags, file_name + "_rscu_frags")


# =============================================================================

if __name__ == "__main__":
    main()
