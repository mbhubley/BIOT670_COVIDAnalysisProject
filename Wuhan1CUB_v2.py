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
import pandas as pd

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

# define a empty codon dictionary
empty_codon_dict ={
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


# This function calculates codon counts and returns a codon count dictionary
def codon_count(cds_seqs):
    # get a copy of codon dictionary
    codon_dict = empty_codon_dict.copy()
    
    count=0
    for seq in cds_seqs:
        for i in range(0,len(seq),3):
            codon_dict[seq[i:i+3]] += 1
            count += 1             
    
    return codon_dict


# This function calculates codon frequencies per AA
def codon_aafreq(codon_counts):
    # get a copy of codon dictionary
    codon_aafreq = empty_codon_dict.copy()
    
    for codon in codon_counts:
        #Find total counts for all synonymous codon
        sumni = 0
        for aa in synonymous_codons:
            if codon in synonymous_codons[aa]:
                for i in synonymous_codons[aa]:
                    sumni += codon_counts[i]
        if sumni > 0:
            codon_aafreq[codon] = codon_counts[codon] / sumni
        else:
            codon_aafreq[codon] = 0

    return codon_aafreq


# This function calculates RSCU from codon frequencies per AA        
def get_RSCU(codon_freq): 
    rscu = {}
    for val in synonymous_codons.values():
        # calculate expected frequence (assume equal usage)
        exp_freq = 1/len(val)

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
    

# This function returns a list of fragments 
def genes(gb_record):
    genes = {}
    for feature in gb_record.features:
        if feature.type == "gene":
            start = feature.location.start
            end = feature.location.end
            key = feature.qualifiers.get("gene")[0]
            genes[key] = [start, end]

    return genes


# this function returns gene start positions
def start_list(gb_record):
    slist = []
    for feature in gb_record.features:
        if feature.type == "gene":
            slist.append(feature.location.start)
            
    return slist


# This function breaks sequence into given size fragments and returns a dictinary
def seq_split(sequence, size):  
    frag={}
    # force triplets
    size -= size%3
    for i in range(0,len(sequence),size):
        locus = str(i)+"-"+str(i+size)
        frag[locus]=(sequence[i:i+size])
        
    return frag


# This function tells related genes with given positions 
def detect_location(begin, end, gene_dict):
    start=""
    stop=""
    for key, val in gene_dict.items():
        if begin >= val[0] and begin <= val[1]:
            start = key
        if end >= val[0] and end <= val[1]:
            stop = key
    return start + "," + stop
        

# this function calculates codon frequencies of the fragments
def frag_codon_freq(sequence):
    seq_len = len(sequence) - len(sequence)%3
    codon_dict = empty_codon_dict.copy()
    
    count=0
    for i in range(0,seq_len,3):
        codon_dict[sequence[i:i+3]] += 1
        count += 1             
    
    for key, val in codon_dict.items():
        codon_dict[key]=  round(val/count,3)
   
    return codon_dict


# generates artifical genome from cds
def cds_to_genome(sample_cds):
    # genome = ""
    # for feature in gb_record.features:
    #     if feature.type == "CDS":
    #         genome += feature.location.extract(gb_record).seq
    # return genome
    genome = ""
    for cds in sample_cds:
        genome += cds
        
    return genome


# finds gene name in artifical genome
def find_gene(gb_record):
    name=[]
    start=[]
    end=[]
    genes={}  
    for feature in gb_record.features:
        if feature.type == "CDS":
            name.append(feature.qualifiers.get("gene")[0])
            start.append(feature.location.start)
            end.append(feature.location.end)

    n =  1
    for i in range(len(name)):
        len_cds= end[i]-start[i]
        if name[i] not in genes.keys():
            genes[name[i]] =[n , n+len_cds-1]
            n = n+len_cds
        else:
            genes[name[i]+"."] =[n , n+len_cds-1]
            n = n+len_cds
    return genes


# This function writes dictionary items into a file
def file_writer(any_dict, name):
    with open(name + ".csv","w") as f:
        for key, val in any_dict.items():
            f.write(key + "," + str(round(val,3)) + "\n")
    print(name + ".csv is written!")


# file witer for multiple fragment frequency dictionary
def file_writer2(sample_frags, name):
    key_list=[]
    freq_list=[]
    for key, val in sample_frags.items():
        pos = key.rsplit("-") # alternate 2
        gene = detect_location(int(pos[0]), int(pos[1]), art_genes) # alternate 2
        key_list.append(gene + "\t" + key) #alternate 2
        #key_list.append(key) #alternate 1
        freq_list.append(frag_codon_freq(val))
    
    df = pd.DataFrame(freq_list, index=key_list)
    df.to_csv(name + ".tsv", sep="\t")
    print (name + ".tsv is written!")
    


# get user inputs
file_name = "wuhan-hu-1_sequence.gb"
frag_size = 1000
#file_name = input("Please enter your GB filename: ")
#frag_size = eval(input("Please enter the size of fragments: "))

# parse genbank file
sample_gb_record = gb_parser(file_name)
sample_seq = gb_to_seq(sample_gb_record)
sample_cds = gb_to_cds(sample_gb_record)
sample_genes = genes(sample_gb_record)


art_genome = cds_to_genome(sample_cds)
sample_frags = seq_split(art_genome, frag_size)
art_genes = find_gene(sample_gb_record)
  
counts = codon_count(sample_cds)
freq = codon_aafreq(counts)
rscu = get_RSCU(freq)
#w = get_w(rscu)


# write frequency table
file_writer(freq, file_name + "_freq" )
# write rscu table  
file_writer(rscu, file_name + "_rscu")
# write fragment codon frequecies
file_writer2(sample_frags, file_name + "_frags_fregs")


# example cai calculation
# seq = "AAATTT"
# seq = str(sample_cds[0])

# cai = calc_cai(seq,w)
# print("CAI: " + str(cai))

  


    
    

# define a fuction for plotting usage frequencies over the genome, codon_freq_plot()

# get known non-human-pathogenic strain sequence 

# count codons on wuhan-1

# create codon usage by AA table

# plot codon usage

# find non-highest-frequency codons

# comparison with other strains
