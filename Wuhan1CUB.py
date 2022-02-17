# ------------ PROJECT DEFINITION ------------------
# 1.	Identify the codon set most used in the Wuhan-1 strain
# 2.	Do a full count of all codons by AA
# 3.	Create a usage by AA table (all actual codons used by frequency for that AA)
#     a.	Look at the distribution of usage throughout the genome, e.g. by plotting the usage frequencies over the genome
#     b.	Look for the non-highest-frequency codons within the genome
#     c.	See if these are clustered in regions
# 4.	Compare to known non-human-pathogenic strain 
# ---------------------------------------------------

from Bio import SeqIO

# This function gets a genbank file and returns biopython genbank record
def gb_parser(file):
    gb_record = SeqIO.read(open(file,"r"), "genbank")
    return gb_record

# This fucntion gets a genbank record and returns a raw sequence
def gb_to_seq(gb_record):
    return gb_record.seq

# This function extracts CDSs from genbank record. It returns a CDS sequence list
def gb_to_cds(gb_record):
    cds_seqs =[]
    for feature in gb_record.features:
        if feature.type == "CDS":
            cds_seqs.append(feature.location.extract(gb_record).seq)
    
    return(cds_seqs)

# This fuction calculates codon frequencies and returns a codon frequency dictionary
def codon_freq(cds_seqs):
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
    
    count=0
    for seq in cds_seqs:
        for i in range(0,len(seq),3):
            codon_dict[seq[i:i+3]] += 1
            count += 1             

    for key in codon_dict.keys():
        codon_dict[key] /= count
    
    return codon_dict

# This function calculates RSCU from codon frequencies        
def get_RSCU(codon_freq):
    synonymous_codons = {
    "C": ["TGT", "TGC"],
    "D": ["GAT", "GAC"],
    "S": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"],
    "Q": ["CAA", "CAG"],
    "M": ["ATG"],
    "N": ["AAC", "AAT"],
    "P": ["CCT", "CCG", "CCA", "CCC"],
    "K": ["AAG", "AAA"],
    "*": ["TAG", "TGA", "TAA"],
    "T": ["ACC", "ACA", "ACG", "ACT"],
    "F": ["TTT", "TTC"],
    "A": ["GCA", "GCC", "GCG", "GCT"],
    "Z": ["GGT", "GGG", "GGA", "GGC"],
    "I": ["ATC", "ATA", "ATT"],
    "L": ["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"],
    "H": ["CAT", "CAC"],
    "R": ["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"],
    "W": ["TGG"],
    "V": ["GTA", "GTC", "GTG", "GTT"],
    "E": ["GAG", "GAA"],
    "Y": ["TAT", "TAC"]}
    
    rscu = {}
    for val in synonymous_codons.values():
        # calculate expected frequence (assume equal usage)
        exp_freq = 0
        for codon in val:
            exp_freq += codon_freq[codon]
        exp_freq /= len(val)

        #calculate rscu (observed/expected)
        for i in val:
            rscu[i] = codon_freq[i]/exp_freq
    
    return rscu

# This function creates CAI from RSCU
def get_CAI(rscu):
    cai={}
    max_rscu = max(rscu.values())
    for key, val in rscu.items():        
        cai[key] = val/max_rscu
        
    return cai
    
# define a fuction for plotting usage frequencies over the genome, codon_freq_plot()

# get wuhan-1 strain genome sequence (from file or fetch from NCBI)
w1_file = "wuhan-hu-1_sequence.gb"

# parse genbank file
w1_record = gb_parser(w1_file)
w1_seq = gb_to_seq(w1_record)
w1_cds = gb_to_cds(w1_record)

# -- test --       
w1_codon_freq = codon_freq(w1_cds)
w1_rscu = get_RSCU(w1_codon_freq)
w1_cai = get_CAI(w1_rscu)

# -- print --
for key,val in w1_cai.items():
    print(key, round(val,2))

# define a fuction for plotting usage frequencies over the genome, codon_freq_plot()

# get known non-human-pathogenic strain sequence 

# count codons on wuhan-1

# create codon usage by AA table

# plot codon usage

# find non-highest-frequency codons

# comparison with other strains
