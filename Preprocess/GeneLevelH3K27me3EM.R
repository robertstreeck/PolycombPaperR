library(tidyverse)
library(rtracklayer)
source("EM/EMcalc.R")

GFF = readGFF("/Users/streeck/Genomes/DmelBDGP6.91/Drosophila_melanogaster.BDGP6.91.gtf", filter = list(type = "gene"))

resting_H3K27me3 = read.delim("data/fig1/H3K27me3ChIP3772.count", sep = "\t", skip = 1)
colnames(resting_H3K27me3)[7:10] = c("H3K27me3_A", "input_A", "H3K27me3_B", "input_B")
resting_H3K27me3 = resting_H3K27me3 %>% select(Geneid, H3K27me3_A:input_B)

S = as.matrix(resting_H3K27me3[,c("H3K27me3_A","H3K27me3_B")])
N = as.matrix(resting_H3K27me3[,c("H3K27me3_A","H3K27me3_B")] + 
                resting_H3K27me3[,c("input_A","input_B")])
GeneLevelH3K27me3Cluster = BinomEMwrapperParallel(N, S, k=3, ncores = 10)

save(GeneLevelH3K27me3Cluster, file = "Data/fig1/GeneLevelH3K27me3Cluster.Rdata")

resting_H3K27me3 = resting_H3K27me3 %>%
  left_join(GFF %>%
              select(gene_id, gene_name, gene_biotype),
            by = c(Geneid = "gene_id")) %>%
  mutate(class = case_when(GeneLevelH3K27me3Cluster$Group == 1 ~ "non-Pc",
                           GeneLevelH3K27me3Cluster$Group == 2 ~ "Pc-I",
                           GeneLevelH3K27me3Cluster$Group == 3 ~ "Pc-H"))

resting_H3K27me3$class = factor(resting_H3K27me3$class, 
                                levels = c("non-Pc", "Pc-I", "Pc-H"))

save(GeneLevelH3K27me3Cluster, file = "data/fig1/GeneLevelH3K27me3Cluster.Rdata")
save(resting_H3K27me3, file = "data/fig1/resting_H3K27me3.Rdata")
