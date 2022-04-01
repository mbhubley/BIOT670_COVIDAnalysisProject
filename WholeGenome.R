# R code for visualization of codon frequencies of different SARS-CoV-2 viral strains
install.packages("ggplot2")
library(ggplot2)
codon_usage <- read.csv("Wuhan_codon_list.csv")
ggplot(codon_usage, aes(x=Codon, y=Frequency, fill = AA)) +
facet_wrap(~ AA, scale="free_x") + theme_bw()+
geom_bar(stat = "identity")+
ylab("Codon Usage Frequency")+xlab("Codon")+
labs(title = "Wuhan-Hu-1 Codon Frequency Per Amino Acid")+
theme(axis.text.x = element_text(size = 7))+
theme(legend.position="none")
