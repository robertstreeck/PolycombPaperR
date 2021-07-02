library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggthemes)


chromosomes = c("X", "Y", "2L", "2R", "3L", "3R")

broasPeak_cols = c("chr", "start", "end", "name", "signal", "strand", "FC", "log_pval", "log_qval")
H3K27me3_resting_macs2 = read.delim("data/fig1/MACS2Test/H3K27me3_resting_peaks.broadPeak", sep = "\t", header = F,
                                    col.names = broasPeak_cols)

H3K27me3_resting_macs2 = H3K27me3_resting_macs2 %>% 
  filter(chr %in% chromosomes)

GFF_dmel6.91 = readGFF("/Users/streeck/Genomes/DmelBDGP6.91/Drosophila_melanogaster.BDGP6.91.gtf", filter = list(type = c("exon")))
GFF_dmel6.91 = GFF_dmel6.91 %>%
  filter(seqid %in% chromosomes) %>%
  droplevels()

GR_GFF_dmel6.91 = makeGRangesFromDataFrame(GFF_dmel6.91, keep.extra.columns = TRUE)

H3K27me3_resting_macs2_sub = H3K27me3_resting_macs2 %>% 
  filter(log_qval > 30) 

GR_H3K27me3_resting_macs2_high = makeGRangesFromDataFrame(H3K27me3_resting_macs2_sub, keep.extra.columns = TRUE)

overlapping_genes_high = subsetByOverlaps(GR_GFF_dmel6.91, GR_H3K27me3_resting_macs2_high)
overlapping_genes_high = unique(overlapping_genes_high$gene_id)

GR_H3K27me3_resting_macs2_any = makeGRangesFromDataFrame(H3K27me3_resting_macs2, keep.extra.columns = TRUE)

overlapping_genes_any = subsetByOverlaps(GR_GFF_dmel6.91, GR_H3K27me3_resting_macs2_any)
overlapping_genes_any = unique(overlapping_genes_any$gene_id)


load("data/fig1/resting_H3K27me3.Rdata")


log_ratio = function(chip, input){
  chip = (chip)/sum(chip)
  input = (input)/sum(input)
  return(log2(chip) - log2(input))
}

hist_split = function(x, n = 100){
  x_min = min(x[!is.nan(x) & !is.infinite(x)])
  x_max = max(x[!is.nan(x) & !is.infinite(x)])
  breaks = seq(x_min, x_max, length.out = n+1)
  assignments = as.numeric(cut(x, breaks = breaks))
  centers = breaks[1:n] - mean(breaks[1:n] - breaks[2:(n+1)])
  return(centers[assignments])
}

resting_H3K27me3 = resting_H3K27me3 %>%
  mutate(lr = log_ratio(H3K27me3_A + H3K27me3_B, input_A + input_B)) %>%
  mutate(bin = hist_split(lr))

resting_summary = resting_H3K27me3 %>%
  group_by(bin) %>%
  summarise(n = n(),
            non_pc = sum(class == "non-Pc")/n(),
            pc_h = sum(class == "Pc-H")/n(),
            pc_m = sum(class == "Pc-I")/n(),
            macs_high = sum(Geneid %in% overlapping_genes_high)/n(),
            macs_any = sum(Geneid %in% overlapping_genes_any)/n()) 

group_labels = c("Pc-H",
                 "Pc-M",
                 "non-Pc",
                 "MACS2 All",
                 "MACS2 High")

labs = c("eve", "en", "abd-A", "Ubx", "Antp", "Abd-B", "eya")

a = resting_summary %>%
  mutate(frac = n/sum(n)) %>%
  ggplot(aes(bin, frac/2, height = frac)) + geom_tile() +
  scale_y_continuous(breaks = c(0, .025, .05), labels = c("0", "2.5%", "5%"),
                     limits = c(0,.05))+
  theme_bw() + ylab("Percent of Genes") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        plot.margin = unit(c(5,5,5,36.7), "points")) +
  geom_label_repel(aes(bin, .001, label = gene_name),
                   data = resting_H3K27me3 %>%
                     filter(gene_name %in% labs),
                   inherit.aes = F, max.overlaps = 10, nudge_y = .01)

b = resting_summary %>%
  pivot_longer(!bin:n) %>%
  mutate(name = factor(name, levels = c("pc_h", "pc_m", "non_pc", 
                                        "macs_any", "macs_high"))) %>%
  ggplot(aes(bin, name, alpha = value, fill = name)) + geom_tile() + 
  geom_point(aes(color = value), alpha = 0) + 
  scale_color_gradient(high = "black", low = "white",
                       name = "Fraction of Genes",
                       guide = guide_colorbar(direction = "horizontal",
                                              title.position = "top",
                                              barwidth = 5, barheight = .5)) +
  theme_bw() + xlab("H3K27me3 Log2 Ratio") + ylab(" ") +
  scale_alpha_continuous(range = c(0,1)) + 
  guides(fill = FALSE, alpha = F) + 
  scale_y_discrete(labels = group_labels) +
  scale_fill_manual(values = rev(c("#cc6677", "#ddcc77", "#6388b4", "#ffae34", "#ef6f6a"))) +
  theme(legend.position = c(.12, .2), legend.title = element_text(size = 8),
        legend.text = element_text(size = 8))

gridExtra::grid.arrange(a, b, ncol = 1, heights = c(.6, .40))
