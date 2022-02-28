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
from Bio.Seq import Seq
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

# This function accepts a genbank file and returns biopython genbank record
def gb_parser(file):
    gb_record = SeqIO.read(open(file,"r"), "genbank")
    return gb_record


# This function determines if a location is inside a set of gene locations
def in_locs(gene_locs, test):
    found = False
    for i in gene_locs:
        if test in i:
            found = True
    return found


# This function creates a list of lists containing all CDS with name, coordinates, and sequence
def genes(gb_record):
    names = []
    locs = []
    seqs = []
    for feature in gb_record.features:
        if feature.type == "CDS":
            try:
                name = feature.qualifiers.get("gene")[0]
            except TypeError:
                name = feature.qualifiers.get("product")[0]
            if name in names:
                continue
            if in_locs(locs, feature.location.start) == True and in_locs(locs, feature.location.end) == True:
                    continue
            names.append(name)
            locs.append(feature.location)
            seqs.append(feature.location.extract(gb_record).seq)
    genes_all = [names, locs, seqs]
    return genes_all


# This function accepts a genbank record and returns a list containing the CDS sequences
def gb_to_cds(gb_record):
    cds_seqs =[]
    for feature in gb_record.features:
        if feature.type == "CDS":
            cds_seqs.append(feature.location.extract(gb_record).seq)
    
    return(cds_seqs)


# This function accepts a list of CDS data and returns the sequence as a list of codons
def cds_codons(cds_list):
    cds_codons = []
    for cds in cds_list:
        for i in range(0,len(cds),3):
            cds_codons.append(cds[i:i+3])
    return cds_codons


# This function accepts a sequence as a list of codons and returns sequence fragments of x codons in length
def cds_divide(cds_codons, x):
    cds_fragments = []
    i=0

    # create list of specified numbers of codons
    while i < len(cds_codons):
        cds_fragments.append(cds_codons[i:i + x])
        i += x

    # join codons into combined sequence
    for f in range(0,(len(cds_fragments))):
        cds_fragments[f] = Seq("").join(cds_fragments[f])

    with open("fragments.txt","w") as f:
        for frag in cds_fragments:
            f.write(str(frag) + "\n")
            
    return cds_fragments    


# This function accepts a sequence and returns codon counts as a dictionary
def codon_count(seq):
    
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


# This function accepts a codon count dictionary and returns frequencies per AA as a dictionary
def codon_aafreq(codon_counts):
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
        #Find total counts for all synonymous codons
        sumni = 0
        for aa in synonymous_codons:
            if codon in synonymous_codons[aa]:
                for i in synonymous_codons[aa]:
                    sumni += codon_counts[i]
        #Calculate frequency for each synonymous codon
        if sumni > 0:
            codon_aafreq[codon] = codon_counts[codon] / sumni
        else:
            #To avoid divide by zero error, if total codon counts for aa = 0, just return 0
            codon_aafreq[codon] = 0

    return codon_aafreq


# This function accepts a dictionary of codon frequencies and returns a dictionary of RSCU values
def get_RSCU(codon_freq): 
    rscu = {}
    for val in synonymous_codons.values():
        # calculate expected frequency (assume equal usage)
        exp_freq = 1/len(val)

        #calculate rscu (observed/expected)
        for i in val:
            rscu[i] = codon_freq[i]/exp_freq

    return rscu


# This function accepts a dictionary with codons as keys and converts codons to numbers
def codon_to_num(codon_dict):
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


# This function calculates CAI (geometric mean of the RSCU values)
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
            f.write(str(key) + "," + str(round(val,3)) + "\n")
    print(name + ".csv is written!")
    
    
# file witer for multiple fragment frequency dictionary
def file_writer2(sample_frags, name):
    key_list=[]
    freq_list=[]
    for key, val in sample_frags.items():
        #pos = key.rsplit("-") # alternate 2
        #gene = detect_location(int(pos[0]), int(pos[1]), sample_genes) # alternate 2
        #key_list.append(gene + " " + key) #alternate 2
        key_list.append(key) #alternate 1
        freq_list.append(val)
    
    df = pd.DataFrame(freq_list, index=key_list)
    df.to_csv(name + ".csv", sep=",")
    print (name + ".csv is written!")


def main():
    # get genome sequence from file
    file_name = input("Please enter your GB filename: ")

    # parse genbank file
    sample_gb_record = gb_parser(file_name)
    sample_cds = gb_to_cds(sample_gb_record)
    sample_genes = genes(sample_gb_record)
    sample_codons = cds_codons(sample_genes[2])


    # analyze whole genome      
    counts_genome = codon_count(Seq("").join(sample_genes[2]))
    freq_genome = codon_aafreq(counts_genome)
    rscu_genome = get_RSCU(freq_genome)

    # analyze by gene
    freq_genes = {}
    rscu_genes = {}
    for i in range(0,len(sample_genes[0])):
        name = sample_genes[0][i]
        count = codon_count(sample_genes[2][i])
        freq = codon_aafreq(count)
        rscu = get_RSCU(freq)
        freq_genes[name] = freq
        rscu_genes[name] = rscu
    
    # analyze by fragments of any number of codons
    # get number of codons per fragment to divide the genome into from user input
    frag_number = int(input("To divide genome into regions of arbitrary size, enter the desired number of codons per region: "))

    # divide genome into fragments
    sample_fragments = cds_divide(sample_codons, frag_number)

    freq_frags = {}
    rscu_frags = {}
    pos = 1
    for i in range(0,len(sample_fragments)):
        count = codon_count(sample_fragments[i])
        freq = codon_aafreq(count)
        rscu = get_RSCU(freq)
        freq_frags[str(pos)+"-"+str(pos+frag_number-1)] = freq
        rscu_frags[str(pos)+"-"+str(pos+frag_number-1)] = rscu
        pos += frag_number

    # ask to convert codons to numbers
    #convert = input("Convert codons to numeric values? Enter Y/N: ")
    #if convert.lower() == "y":
        #rscu_genome = codon_to_num(rscu_genome)
        #for i in range(0,len(rscu_genes)):
            #rscu_genes[i] = codon_to_num(rscu_genes[i])
        #for i in range(0,len(rscu_genes)):
            #rscu_fragments[i] = codon_to_num(rscu_fragments[i])

    # write genome frequency table
    # file_writer(freq_genome, file_name + "_freq_genome" )
    # write genome rscu table  
    file_writer(rscu_genome, file_name + "_rscu_genome")

    # write gene rscu tables and save in genes folder
    file_writer2(rscu_genes, file_name + "_rscu_genes")

    # write fragment rscu tables and save in fragment folder
    file_writer2(rscu_frags, file_name + "_rscu_frags")
    

    # comparison with other strains
if __name__ == "__main__":
    main()
