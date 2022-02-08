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

# This function extracts CDSs from genbank record. It returns a CDS dictionary
def gb_to_cds(gb_record):
    cds_dict={}

    for feature in gb_record.features:    
        start = feature.location.start.position
        end = feature.location.end.position
        if feature.type == "5'UTR":        
            cds_dict["5'UTR"] = [start, end]
        elif feature.type == "CDS":
            locus_tag = feature.qualifiers['protein_id'][0]
            cds_dict[locus_tag] = [start, end]
            
    return(cds_dict)

# define a fuction for codon count, codon_count()
    # define codon:AA dictionary

# define a function for creating an codon usage table, codon_use_table()
    
# define a fuction for plotting usage frequencies over the genome, codon_freq_plot()

# get wuhan-1 strain genome sequence (from file or fetch from NCBI)
w1_file = "wuhan-hu-1_sequence.gb.txt"

# parse genbank file
w1_record = gb_parser(w1_file)
w1_seq = gb_to_seq(w1_record)
w1_cds = gb_to_cds(w1_record)

# -- test --            
print(w1_seq)
print(w1_cds)
      
# get known non-human-pathogenic strain sequence 

# count codons on wuhan-1

# create codon usage by AA table

# plot codon usage

# find non-highest-frequency codons

# comparison with other strains
