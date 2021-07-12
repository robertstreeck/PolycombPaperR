library(rtracklayer)
library(tidyverse)
library(viridis)
library(ComplexHeatmap)
library(ggplot2)
library(ggthemes)



TSSGFF1kb = import.gff("/Users/streeck/Genomes/DmelBDGP6.91/Drosophila_melanogaster.BDGP6.91.gtf")
TSSGFF1kbup = TSSGFF1kb[TSSGFF1kb$type == "transcript" & TSSGFF1kb$gene_biotype == "protein_coding"]
end(TSSGFF1kbup[strand(TSSGFF1kbup)=="+",]) = start(TSSGFF1kbup[strand(TSSGFF1kbup)=="+",]) + 1000
start(TSSGFF1kbup[strand(TSSGFF1kbup)=="+",]) = start(TSSGFF1kbup[strand(TSSGFF1kbup)=="+",]) -3000
start(TSSGFF1kbup[strand(TSSGFF1kbup)=="-",]) = end(TSSGFF1kbup[strand(TSSGFF1kbup)=="-",]) - 1000
end(TSSGFF1kbup[strand(TSSGFF1kbup)=="-",]) = end(TSSGFF1kbup[strand(TSSGFF1kbup)=="-",]) + 3000


load(file = "data/fig1/SevenClassGenomeModel.Rdata")
load("data/fig1/resting_H3K27me3.Rdata")
GFF = readGFF("/Users/streeck/Genomes/DmelBDGP6.91/Drosophila_melanogaster.BDGP6.91.gtf", filter = list(type = "gene"))

TSSGFF1kbup$H3K27me3State = resting_H3K27me3$class[match(TSSGFF1kbup$gene_id, resting_H3K27me3$Geneid)]

state_genome = multi_chip_fit$genome

TSSCluster = matrix(nrow = 21, ncol = length(TSSGFF1kbup))

for (i in 1:length(TSSGFF1kbup)) {
  if(as.logical(TSSGFF1kbup@strand[i] == "+")){
    t = subsetByOverlaps(state_genome, TSSGFF1kbup[i], minoverlap = 1)$name
  }
  if(as.logical(TSSGFF1kbup@strand[i] == "-")){
    t = rev(subsetByOverlaps(state_genome, TSSGFF1kbup[i], minoverlap = 1)$name)
  }
  if(length(t) == 21){
    TSSCluster[,i] = t
  }
}


TSSCluster = t(TSSCluster)
save(TSSCluster, file = "data/fig1/2021_06_09_PerGenePromotorState.Rdata")


dist_matrices = list()
for (i in c("non-Pc", "Pc-I", "Pc-H")) {
  dist_matrices[[i]] = StatMatch::gower.dist(TSSCluster[TSSGFF1kbup$H3K27me3State == i & rowSums(is.na(TSSCluster)) == 0,])
}

save(dist_matrices, file = "data/fig1/2021_06_09_PerGenePromotorStateDistData.Rdata")

df.list = list()
for (i in c("non-Pc", "Pc-I", "Pc-H")) {
  in_clust_order = fastcluster::hclust(as.dist(dist_matrices[[i]]))$order
  ts_ids = TSSGFF1kbup$transcript_id[TSSGFF1kbup$H3K27me3State == i & rowSums(is.na(TSSCluster)) == 0]
  df.list[[i]] = cbind(
    data.frame(
      gene_id = TSSGFF1kbup$gene_id[TSSGFF1kbup$H3K27me3State == i & rowSums(is.na(TSSCluster)) == 0], 
      transcript_id = factor(ts_ids, levels = ts_ids[in_clust_order]),
      gene_state = i
    ),
    as.data.frame(
      TSSCluster[TSSGFF1kbup$H3K27me3State == i & rowSums(is.na(TSSCluster)) == 0,]
    )
  )
}


class_coding = c(`non-Pc` = "non-Pc", `Pc-H` = "Pc-H", 
                   `Pc-I` = "Pc-M")

state_coding = c(`1` = "EnhW", `2` = "Pc-I", `3` = "TEl", `4` = "EnhS", 
                 `5` = "Pc-H", `6` = "TSS", `7` = "Het", `None` = "None")

TSS_heatmap_plot = data.table::rbindlist(df.list) %>%
  filter(!is.na(transcript_id)) %>%
  dplyr::mutate(gene_state = recode(gene_state, !!!class_coding)) %>%
  dplyr::mutate(gene_state = factor(gene_state, 
                                    levels = c("non-Pc", "Pc-M", "Pc-H"))) %>%
  dplyr::mutate(across(V1:V21, ~recode(., !!!state_coding)))
  

named_cols = c(ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value[1:8][c(6,2,1,5,3,4,7)], "#FFFFFF")
names(named_cols) = state_coding

### Fig7StateCoverageHeatmap.pdf

Heatmap(TSS_heatmap_plot[,4:24], 
        column_order = 1:21,
        row_order = order(TSS_heatmap_plot$transcript_id),
        row_split = TSS_heatmap_plot$gene_state,
        show_row_names = F,
        show_column_names = F,
        col = named_cols,
        gap = unit(2, "mm"))

class_dict = c("EnhW", "Pc-I", "TEl", "EnhS", "Pc-H", "TSS", "Het")

TSS_state_plot = cbind(gene_id = TSSGFF1kbup$gene_id, 
                       transcript_id = TSSGFF1kbup$transcript_id,
                       gene_state = TSSGFF1kbup$H3K27me3State,
                       TSSCluster) %>%
  as.data.frame() %>%
  pivot_longer(!gene_id:gene_state, values_to = "chrom_state", names_to = "position") %>%
  dplyr::mutate(position = as.numeric(substr(position, 2, 3)) - 3) %>%
  filter(!is.na(transcript_id))

TSS_profile = TSS_state_plot %>%
  group_by(gene_state, position, chrom_state) %>%
  dplyr::summarise(count = n()) %>%
  group_by(gene_state, position) %>%
  dplyr::mutate(freq = count/sum(count)) %>%
  filter(!is.na(chrom_state)) %>%
  filter(chrom_state != "None") %>%
  dplyr::mutate(chrom_state = class_dict[as.numeric(chrom_state)])

### Fig7StateCoverageProfile.pdf

TSS_profile %>%
  dplyr::mutate(gene_state = factor( c("non-Pc", "Pc-M", "Pc-H")[as.numeric(gene_state)], 
                                     levels = c("non-Pc", "Pc-M", "Pc-H"))) %>%
  ggplot(aes(position, freq, color = gene_state)) + 
  geom_line() + theme_bw() + 
  scale_color_tableau(palette = "Superfishel Stone", name = "Gene State") +
  ylab("Fraction of bins in Gene state") + theme(axis.title.x = element_blank()) +
  scale_x_continuous(breaks = c(1, 15, 21), labels = c("-3 kb", "TSS", "+1 kb")) + 
  facet_wrap(vars(chrom_state), ncol = 1, scales = "free_y", strip.position = "right") +
  scale_y_continuous(limits = c(0,1))




#### Trying proximity


merge_intervals = function(granges.object){
  GenomicRanges::reduce(granges.object + 2) - 2
}

merged_state_genome = merge_intervals(state_genome[state_genome$name == "1"])
merged_state_genome$name = "1"
for (i in 2:7) {
  a = merge_intervals(state_genome[state_genome$name == as.character(i)])
  a$name = as.character(i)
  merged_state_genome = c(merged_state_genome,a)
}

merged_state_genome = sort(merged_state_genome)

TSSquery = TSSGFF1kb[TSSGFF1kb$type == "transcript"]
end(TSSquery[strand(TSSquery)=="+",]) = start(TSSquery[strand(TSSquery)=="+",]) + 1
start(TSSquery[strand(TSSquery)=="-",]) = end(TSSquery[strand(TSSquery)=="-",]) - 1


max.dist = 10000

EnhW_nearest_TSS = distanceToNearest(merged_state_genome[merged_state_genome$name == "1"], TSSquery)
EnhW_nearest_TSS = EnhW_nearest_TSS[EnhW_nearest_TSS@elementMetadata$distance < max.dist]
EnhW_nearest_TSS = data.frame(gene_id = TSSquery$gene_id[EnhW_nearest_TSS@to]) %>%
  group_by(gene_id) %>%
  dplyr::summarise(count = n())


EnhS_nearest_TSS = distanceToNearest(merged_state_genome[merged_state_genome$name == "4"], TSSquery)
EnhS_nearest_TSS = EnhS_nearest_TSS[EnhS_nearest_TSS@elementMetadata$distance < max.dist]
EnhS_nearest_TSS = data.frame(gene_id = TSSquery$gene_id[EnhS_nearest_TSS@to]) %>%
  group_by(gene_id) %>%
  dplyr::summarise(count = n())


TSS_nearest_TSS = distanceToNearest(merged_state_genome[merged_state_genome$name == "6"], TSSquery)
TSS_nearest_TSS = TSS_nearest_TSS[TSS_nearest_TSS@elementMetadata$distance < max.dist]
TSS_nearest_TSS = data.frame(gene_id = TSSquery$gene_id[TSS_nearest_TSS@to]) %>%
  group_by(gene_id) %>%
  dplyr::summarise(count = n())



Proximal_state = resting_H3K27me3 %>%
  dplyr::select(Geneid, gene_name, gene_biotype, class) %>%
  filter(gene_biotype == "protein_coding") %>%
  left_join(EnhS_nearest_TSS,
            by = c("Geneid" = "gene_id")) %>%
  left_join(EnhW_nearest_TSS,
            by = c("Geneid" = "gene_id"),
            suffix = c(".EnhS", ".EnhW")) %>%
  left_join(TSS_nearest_TSS,
            by = c("Geneid" = "gene_id")) %>%
  dplyr::rename(count.TSS = count) %>%
  replace_na(list(count.TSS = 0, count.EnhS = 0, count.EnhW = 0)) %>%
  pivot_longer(count.EnhS:count.TSS) %>%
  dplyr::mutate(bin_group = cut(value, breaks = c(-.5, .5, 1.5, 5.5, Inf), labels = c("0", "1", "2-5", ">5")))

### Fig7StateCoverageProximalEnhancers.pdf

Proximal_state %>%
  dplyr::mutate(name = factor(substr(name, 7,12), levels = c("EnhW", "EnhS", "TSS"))) %>%
  ggplot(aes(x = class, fill = bin_group)) +
  geom_bar(position = "fill") +
  scale_fill_carto_d(palette = "Emrld", name = "number of\nproximal elements") +
  facet_wrap(vars(name), nrow = 3, strip.position = "right") +
  xlab("Gene State") + ylab("Fraction of genes in state") +
  theme_bw() + coord_flip() +
  theme(legend.position = "bottom")


