# R code for visualization of RSCU values of specific genes and gene fragments 
# for the different SARS-CoV-2 viral strains
install.packages("ggplot2")
library(ggplot2)
install.packages("tidyr")
library(tidyr)
codon_usage <- read.csv("test1.csv")
ggplot(codon_usage, aes(x=Codon, y=RSCU_Values, fill = Measure))+
geom_bar(stat='identity')+
facet_wrap(~Measure, ncol=1, strip.position = "left")+
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
theme(legend.position="none")+
labs(title = "input_graph_title_here")
