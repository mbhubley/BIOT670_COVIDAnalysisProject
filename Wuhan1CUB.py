# ------------ PROJECT DEFINITION ------------------
# 1.	Identify the codon set most used in the Wuhan-1 strain
# 2.	Do a full count of all codons by AA
# 3.	Create a usage by AA table (all actual codons used by frequency for that AA)
#     a.	Look at the distribution of usage throughout the genome, e.g. by plotting the usage frequencies over the genome
#     b.	Look for the non-highest-frequency codons within the genome
#     c.	See if these are clustered in regions
# 4.	Compare to known non-human-pathogenic strain 
# ---------------------------------------------------

# Please refer for calculations:
# Sharp, P. M., & Li, W.-H. (1987). Nucleic Acids Research, 15(3), 1281–1295. 
# https://doi.org/10.1093/nar/15.3.1281

from Bio import SeqIO
from math import exp, log

# define aamino acid codon table
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

# define synonymous_codons
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

# This function calculates relative adaptiveness(W) of codons
def get_w(rscu):
    w={}
    for key, val in rscu.items():
        aa = gencode[key]
        rscu_list=[]
        for codon in synonymous_codons[aa]:
            rscu_list.append(rscu[codon])
        w[key] = val/max(rscu_list)
        
    return w

# This function camculates CAI (geometric mean of the RSCU values)
def calc_cai(sequence, w):
    L = len(sequence)/3 # number of codon
    cai = 0
    for i in range(0,len(sequence),3):
        cai += log(w[sequence[i:i+3]])

    return round(exp((1/L)*cai),3)
    

# This function writes dictionary items into a file
def file_writer(any_dict, name):
    with open(name + ".csv","w") as f:
        for key, val in any_dict.items():
            f.write(key + "," + str(round(val,3)) + "\n")
    print(name + ".csv is written!")
    
# define a fuction for plotting usage frequencies over the genome, codon_freq_plot()

# get wuhan-1 strain genome sequence (from file or fetch from NCBI)
# file_name = "wuhan-hu-1_sequence.gb"
file_name = input("Please enter your GB filename: ")

# parse genbank file
sample_gb_record = gb_parser(file_name)
sample_seq = gb_to_seq(sample_gb_record)
sample_cds = gb_to_cds(sample_gb_record)

# -- test --       
freq = codon_freq(sample_cds)
rscu = get_RSCU(freq)
w = get_w(rscu)

# -- print --
#print(rscu)
#print(w)

# write frequency table
file_writer(freq, file_name + "_freq" )
# write rscu table  
file_writer(rscu, file_name + "_rscu")

# example cai calculation
# seq = "AAATTT"
seq = str(sample_cds[0])

cai = calc_cai(seq,w)
print("CAI: " + str(cai))

# define a fuction for plotting usage frequencies over the genome, codon_freq_plot()

# get known non-human-pathogenic strain sequence 

# count codons on wuhan-1

# create codon usage by AA table

# plot codon usage

# find non-highest-frequency codons

# comparison with other strains
