library(knitr)
library(kableExtra)
library(tibble)
library(dplyr)
library(circlize)
library(reshape2)
library(tidyverse)
library(rcartocolor)
library(data.table)
library(viridis)
library(ggrepel)
library(edgeR)
library(openxlsx)
library(rtracklayer)
library(DESeq2)
library(ggthemes)
library(ComplexHeatmap)

GFF = readGFF("/Users/streeck/Genomes/DmelBDGP6.91/Drosophila_melanogaster.BDGP6.91.gtf", filter = list(type = "gene"))


CountTable = read.delim("data/fig2/unstranded.count", skip = 1)
CountTable = CountTable[,c(1,7:dim(CountTable)[2])]

colnames(CountTable)[2:dim(CountTable)[2]] = sapply(strsplit(colnames(CountTable)[2:dim(CountTable)[2]], ".", fixed = T), function(x) paste0(x[length(x)-2], "_", x[length(x)-1]))

CountTable[,dim(CountTable)[2]+1] = read.delim("data/fig2/unstrandedfeatureCounts.count", skip = 1)[7]
colnames(CountTable)[dim(CountTable)[2]] = "NG_26757_24"


SampleDesctiptors = read.xlsx("data/fig2/2021-04-11MappedSamples.xlsx")
SampleDesctiptors$sample_num = sapply(strsplit(SampleDesctiptors$sample_id, "_", fixed = T), function(x) x[2])
SampleDesctiptors$sample_id = gsub("-", "_", SampleDesctiptors$sample_id, fixed = T)

## exclude samples 1 and 8 (look bad on PCA)
CountTable = CountTable[, !colnames(CountTable) %in% c("NG_26757_1", "NG_26757_8")]
SampleDesctiptors = SampleDesctiptors %>%
  filter(!sample_num %in% c("1", "8"))

groups = SampleDesctiptors$condition
designMatrix = model.matrix(~ 0 + groups)
colnames(designMatrix) = substr(colnames(designMatrix), 7, 100)
colnames(designMatrix) = gsub(" ", "_", colnames(designMatrix), fixed = T)

DGEobject = DGEList(CountTable[,2:dim(CountTable)[2]], group = groups, genes = CountTable[,1])
DGEobject = DGEobject[rowSums(cpm(DGEobject)>1) >= max(table(groups)), , keep.lib.sizes=FALSE]
DGEobject$genes = GFF[match(DGEobject$genes$genes, GFF$gene_id), c(9,10,12)]
DGEobject = calcNormFactors(DGEobject)
DGEobject = estimateDisp(DGEobject, design = designMatrix, robust = T)

## lift data into DESeq and rlog transform

dds = DESeqDataSetFromMatrix(countData = as.matrix(DGEobject$counts), colData = data.frame(row.names=row.names(DGEobject$samples), gr = DGEobject$samples$group), design = ~ gr)
dds = dds[rowSums(counts(dds)) > 1,]
rld = vst(dds, blind = FALSE)


# adjust contrasts
## contrasts and fits
contrasts = makeContrasts(eGFP_induction = EGFP_Ind - EGFP_control,
                          hop_induction = hop_Ind - hop_control,
                          PGRP_induction = PGRP_Ind - PGRP_control,
                          hop_fourway = (hop_Ind - hop_control) - (EGFP_Ind - EGFP_control),
                          PGRP_fourway = (PGRP_Ind - PGRP_control) - (EGFP_Ind - EGFP_control),
                          levels = designMatrix)

fit = glmQLFit(DGEobject, designMatrix, robust = T)
QFresults = list()
for (i in 1:dim(contrasts)[2]) {
  QFresults[[colnames(contrasts)[i]]] = glmQLFTest(fit, contrast = contrasts[,i])
}



selected_contrasts = c("eGFP_induction", "hop_induction", "PGRP_induction", "hop_fourway", "PGRP_fourway")

load("data/fig1/resting_H3K27me3.Rdata")

hypMod = function(qf){
  df = data.frame(gene_id = qf$genes$gene_id,
                  de = as.numeric(decideTests.DGELRT(qf, lfc = log2(1.5))))
  df %>%
    left_join(resting_H3K27me3 %>%
                replace(is.na(.), " ") %>%
                dplyr::select(c(Geneid, class)),
              by = c(gene_id = "Geneid")) %>%
    group_by(de, class) %>%
    summarise(q = n()) %>%
    group_by(de) %>%
    mutate(k = sum(q)) %>%
    group_by(class) %>%
    mutate(m = sum(q)) %>%
    mutate(n = dim(DGEobject$genes)[1] - m) %>%
    tidyr::drop_na() %>%
    mutate(p_val = phyper(q-1, m, n, k, lower.tail = F, log.p = T)) %>%
    mutate(effect_size = (q/k)/(m/(m+n))) %>%
    mutate(p = exp(p_val)) %>%
    filter(de != 0) %>%
    filter(class == "Pc-I") %>%
    arrange(p)
}

hyped = sapply(QFresults[selected_contrasts], simplify = F, hypMod, USE.NAMES = T)

y_coding = c(hop_fourway = "hop-fourway", PGRP_fourway = "PGRP-fourway", 
             eGFP_induction = "eGFP", PGRP_induction = "PGRP-LC", hop_induction = "hop")



# FigSwitchEnrichment.pdf -------------------------------------------------

a = rbindlist(hyped, idcol = "contrast") %>%
  mutate(p.adj = -log10(p.adjust(p))) %>%
  mutate(y = recode(contrast, !!!y_coding)) %>%
  mutate(y = factor(y, levels = y_coding)) %>%
  mutate(greys = as.numeric(de))
a %>%
  mutate(cols = case_when(greys == 1 ~ "Up-regulated", T ~ "Down-regulated")) %>%
  ggplot(aes(y, k)) + coord_flip() +
  geom_bar(aes(y = greys), color = "black", fill = "#AAAAAA", stat = "identity") + 
  geom_bar(aes(y = de*q/k, fill = p.adj), color = "black", stat = "identity") + 
  theme_bw() + ylab("Fraction of Pc-M") +
  scale_fill_viridis(option = "viridis", limits = c(0, max(a$p.adj)), name = "-log10\nadjusted\np-value") + 
  xlab("Contrast") +
  geom_text(aes(y = de*.85, label = paste0("n = ", k))) +
  geom_hline(aes(yintercept = de*m/(m+n)), linetype = 5) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  facet_grid(rows = vars(y), cols = vars(cols), scales = "free")


rbindlist(hyped, idcol = "contrast") %>%
  mutate(q = -log10(p.adjust(p))) %>%
  mutate(y = recode(contrast, !!!y_coding)) %>%
  mutate(y = factor(y, levels = y_coding)) %>%
  dplyr::select(contrast, class, de, effect_size, p, q) %>% kable() %>% 
  kableExtra::kable_styling(bootstrap_options = c("striped", "condensed"))


# volcanos ----------------------------------------------------------------


log_ratio = function(chip, input){
  chip = (chip)/sum(chip)
  input = (input)/sum(input)
  return(log2(chip) - log2(input))
}

resting_H3K27me3 = resting_H3K27me3 %>%
  mutate(lr_A = log_ratio(H3K27me3_A, input_A),
         lr_B = log_ratio(H3K27me3_B, input_B)) %>%
  mutate(h3k27me3_mean_lr = (lr_A + lr_B)/2)

### FigSwtichPGRPLCVolcano.pdf
qf = QFresults$PGRP_induction
colored_volcano = cbind(qf$genes, qf$table)
colored_volcano$de = as.numeric(decideTests.DGELRT(qf, lfc = log2(1.5)))
colored_volcano = colored_volcano %>% 
  left_join(resting_H3K27me3 %>% 
              dplyr::select(Geneid, h3k27me3_mean_lr, class),
            by = c("gene_id" = "Geneid")) %>%
  mutate(col = case_when(de != 0 & class == "non-Pc" ~ "non-Pc",
                         de != 0 & class != "non-Pc" ~ "Pc-M",
                         T ~ "not up-regulated"))

colored_volcano$col = factor(colored_volcano$col,
                             levels = c("non-Pc", "Pc-M", "not up-regulated"))

volcano_labels = colored_volcano %>% 
  mutate(e = logFC^2 + log10(PValue)^2) %>%
  arrange(desc(e)) %>%
  filter(de == 1) %>%
  slice_head(n = 20)

a = ggplot(colored_volcano, aes(logFC, -log10(PValue), color = col)) + 
  geom_point(alpha = .8, shape = 16) +
  geom_label_repel(data = volcano_labels, 
                   aes(logFC, -log10(PValue), label = gene_name, color = col), 
                   size = 2.5, fill = rgb(1,1,1,0.5), inherit.aes = F, show.legend = F,
                   max.overlaps = 20) +
  scale_color_manual(name = "",
    values = ggthemes_data$tableau$`color-palettes`$regular$`Superfishel Stone`$value[c(1,2,10)]) + 
  theme_classic() + xlab("log2 fold change") + ylab("p-value") +
  scale_y_continuous(breaks = c(0, 4,8,12,16,20,24), 
                     labels = c(1, expression(10^-4), expression(10^-8), expression(10^-12), expression(10^-16), expression(10^-20), expression(10^-24))) +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 2))) +
  scale_x_continuous(limits = c(-max(abs(colored_volcano$logFC)), max(abs(colored_volcano$logFC)))) +
  theme(legend.position = c(0.13, 0.2),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  ggtitle("PGRP-LC Induction")


boxplot_data = colored_volcano %>% 
  filter(class %in% c("non-Pc", "Pc-I")) %>%
  mutate(direct = case_when(de == 1 ~ "up",
                            de == -1 ~ "down",
                            T ~ "not de")) %>%
  dplyr::select(logFC, direct, class)

boxplot_data$direct = factor(boxplot_data$direct,
                             levels = c("down", "not de", "up"))


b = ggplot(boxplot_data, aes(logFC, direct, fill = class)) + 
  geom_boxplot(outlier.alpha = 0) +
  scale_fill_manual(name = "",
                     values = ggthemes_data$tableau$`color-palettes`$regular$`Superfishel Stone`$value[c(1,2)]) + 
  theme_classic() + xlab("log2 fold change") + 
  scale_y_discrete(name = "") + 
  scale_x_continuous(limits = c(-max(abs(colored_volcano$logFC)), max(abs(colored_volcano$logFC)))) +
  theme(legend.position = "none")

gridExtra::grid.arrange(a, b, ncol = 1, heights = c(.8, .20))


### FigSwtichhopVolcano.pdf
qf = QFresults$hop_induction
colored_volcano = cbind(qf$genes, qf$table)
colored_volcano$de = as.numeric(decideTests.DGELRT(qf, lfc = log2(1.5)))
colored_volcano = colored_volcano %>% 
  left_join(resting_H3K27me3 %>% 
              dplyr::select(Geneid, h3k27me3_mean_lr, class),
            by = c("gene_id" = "Geneid")) %>%
  mutate(col = case_when(de != 0 & class == "non-Pc" ~ "non-Pc",
                         de != 0 & class != "non-Pc" ~ "Pc-M",
                         T ~ "not up-regulated"))

colored_volcano$col = factor(colored_volcano$col,
                             levels = c("non-Pc", "Pc-M", "not up-regulated"))

volcano_labels = colored_volcano %>% 
  mutate(e = logFC^2 + log10(PValue)^2) %>%
  arrange(desc(e)) %>%
  filter(de == 1) %>%
  slice_head(n = 20)

a = ggplot(colored_volcano, aes(logFC, -log10(PValue), color = col)) + 
  geom_point(alpha = .8, shape = 16) +
  geom_label_repel(data = volcano_labels, 
                   aes(logFC, -log10(PValue), label = gene_name, color = col), 
                   size = 2.5, fill = rgb(1,1,1,0.5), inherit.aes = F, show.legend = F,
                   max.overlaps = 20) +
  scale_color_manual(name = "",
                     values = ggthemes_data$tableau$`color-palettes`$regular$`Superfishel Stone`$value[c(1,2,10)]) + 
  theme_classic() + xlab("log2 fold change") + ylab("p-value") +
  scale_y_continuous(breaks = c(0, 4,8,12,16,20,24), 
                     labels = c(1, expression(10^-4), expression(10^-8), expression(10^-12), expression(10^-16), expression(10^-20), expression(10^-24))) +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 2))) +
  scale_x_continuous(limits = c(-max(abs(colored_volcano$logFC)), max(abs(colored_volcano$logFC)))) +
  theme(legend.position = c(0.13, 0.2),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  ggtitle("hop Induction")


boxplot_data = colored_volcano %>% 
  filter(class %in% c("non-Pc", "Pc-I")) %>%
  mutate(direct = case_when(de == 1 ~ "up",
                            de == -1 ~ "down",
                            T ~ "not de")) %>%
  dplyr::select(logFC, direct, class)

boxplot_data$direct = factor(boxplot_data$direct,
                             levels = c("down", "not de", "up"))


b = ggplot(boxplot_data, aes(logFC, direct, fill = class)) + 
  geom_boxplot(outlier.alpha = 0) +
  scale_fill_manual(name = "",
                    values = ggthemes_data$tableau$`color-palettes`$regular$`Superfishel Stone`$value[c(1,2)]) + 
  theme_classic() + xlab("log2 fold change") + 
  scale_y_discrete(name = "") + 
  scale_x_continuous(limits = c(-max(abs(colored_volcano$logFC)), max(abs(colored_volcano$logFC)))) +
  theme(legend.position = "none")

gridExtra::grid.arrange(a, b, ncol = 1, heights = c(.8, .20))


### FigSwticheGFPVolcano.pdf
qf = QFresults$eGFP_induction
colored_volcano = cbind(qf$genes, qf$table)
colored_volcano$de = as.numeric(decideTests.DGELRT(qf, lfc = log2(1.5)))
colored_volcano = colored_volcano %>% 
  left_join(resting_H3K27me3 %>% 
              dplyr::select(Geneid, h3k27me3_mean_lr, class),
            by = c("gene_id" = "Geneid")) %>%
  mutate(col = case_when(de != 0 & class == "non-Pc" ~ "non-Pc",
                         de != 0 & class != "non-Pc" ~ "Pc-M",
                         T ~ "not up-regulated"))

colored_volcano$col = factor(colored_volcano$col,
                             levels = c("non-Pc", "Pc-M", "not up-regulated"))

volcano_labels = colored_volcano %>% 
  mutate(e = logFC^2 + log10(PValue)^2) %>%
  arrange(desc(e)) %>%
  filter(de == 1) %>%
  slice_head(n = 20)

a = ggplot(colored_volcano, aes(logFC, -log10(PValue), color = col)) + 
  geom_point(alpha = .8, shape = 16) +
  geom_label_repel(data = volcano_labels, 
                   aes(logFC, -log10(PValue), label = gene_name, color = col), 
                   size = 2.5, fill = rgb(1,1,1,0.5), inherit.aes = F, show.legend = F,
                   max.overlaps = 20) +
  scale_color_manual(name = "",
                     values = ggthemes_data$tableau$`color-palettes`$regular$`Superfishel Stone`$value[c(1,2,10)]) + 
  theme_classic() + xlab("log2 fold change") + ylab("p-value") +
  scale_y_continuous(breaks = c(0, 4,8,12,16,20,24), 
                     labels = c(1, expression(10^-4), expression(10^-8), expression(10^-12), expression(10^-16), expression(10^-20), expression(10^-24))) +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 2))) +
  scale_x_continuous(limits = c(-max(abs(colored_volcano$logFC)), max(abs(colored_volcano$logFC)))) +
  theme(legend.position = c(0.13, 0.2),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  ggtitle("eGFP Induction")


boxplot_data = colored_volcano %>% 
  filter(class %in% c("non-Pc", "Pc-I")) %>%
  mutate(direct = case_when(de == 1 ~ "up",
                            de == -1 ~ "down",
                            T ~ "not de")) %>%
  dplyr::select(logFC, direct, class)

boxplot_data$direct = factor(boxplot_data$direct,
                             levels = c("down", "not de", "up"))


b = ggplot(boxplot_data, aes(logFC, direct, fill = class)) + 
  geom_boxplot(outlier.alpha = 0) +
  scale_fill_manual(name = "",
                    values = ggthemes_data$tableau$`color-palettes`$regular$`Superfishel Stone`$value[c(1,2)]) + 
  theme_classic() + xlab("log2 fold change") + 
  scale_y_discrete(name = "") + 
  scale_x_continuous(limits = c(-max(abs(colored_volcano$logFC)), max(abs(colored_volcano$logFC)))) +
  theme(legend.position = "none")

gridExtra::grid.arrange(a, b, ncol = 1, heights = c(.8, .20))

### FigSubSwitchUasInduction.pdf

induction_cpm = cbind(DGEobject$genes, cpm(DGEobject)) %>%
  filter(gene_name %in% c("PGRP-LC", "hop")) %>%
  pivot_longer(!gene_id:gene_biotype, names_to = "sample_id", values_to = "cpm") %>%
  mutate(sample_id = substr(sample_id, 4, 20)) %>%
  mutate(logcpm = log10(cpm)) %>%
  group_by(gene_name) %>%
  mutate(z = (logcpm - mean(logcpm))/sd(logcpm)) %>%
  left_join(SampleDesctiptors %>% 
              dplyr::select(sample_id, condition) %>%
              mutate(sample_id = substr(sample_id, 4, 20)), by = "sample_id") %>%
  mutate(construct = sapply(condition, function(x){strsplit(x, "_", fixed = T)[[1]][1]}, USE.NAMES = F)) %>%
  mutate(induction = sapply(condition, function(x){strsplit(x, "_", fixed = T)[[1]][2]}, USE.NAMES = F)) %>%
  filter(construct != "tkv") %>%
  mutate(construct = case_when(construct == "PGRP" ~ "PGRP-LC",
                               T ~ construct))


ggplot(induction_cpm, aes(construct, cpm, fill = condition)) +
  geom_boxplot() + scale_fill_brewer(palette = "Paired", name = "sample") +
  ylab("cpm") + xlab("genotype") + theme_bw() +
  scale_y_log10() + facet_grid(cols = vars(gene_name), scales = "free")


# heatmap -----------------------------------------------------------------


specific_sets = GFF %>%
  dplyr::select(gene_id) %>%
  mutate(class = case_when(gene_id %in% GFF$gene_id[GeneLevelH3K27me3Cluster$Group == 2] ~ "Pc-M",
                           T ~ "a"))

hyper_gs = DGEobject$genes %>%
  dplyr::select(gene_id) %>%
  left_join(resting_H3K27me3 %>%
              select(Geneid, class), 
            by = c(gene_id = "Geneid")) %>%
  group_by(class) %>%
  dplyr::summarize(m = n()) %>%
  dplyr::mutate(n = dim(DGEobject$genes)[1] - m)

rld_assay = SummarizedExperiment::assay(rld) %>%
  as.data.frame 

rld_rows = DGEobject$genes %>%
  dplyr::select(gene_id) %>%
  left_join(resting_H3K27me3 %>%
              select(Geneid, class), 
            by = c(gene_id = "Geneid")) 

cond_levels = paste0(rep(c("EGFP", "hop", "PGRP"), each = 2), rep(c("_control", "_Ind"), 3))

nk = 6

row_select = ((decideTests.DGELRT(QFresults$eGFP_induction, lfc = 1, p.value = 0.00001) != 0) | 
                (decideTests.DGELRT(QFresults$hop_induction, lfc = 1, p.value = 0.00001) != 0) | 
                (decideTests.DGELRT(QFresults$PGRP_induction, lfc = 1, p.value = 0.00001) != 0))

rld_selected_counts = rld_assay[row_select,]
scaled_rld = rld_selected_counts %>% 
  t %>%
  scale %>%
  t

rld_cluster = DGEobject$genes %>%
  .[row_select,]

rld_cluster$cl = kmeans(scaled_rld, centers = nk, algorithm = 'MacQueen', iter.max = 200, nstart = 40)$cluster

heatmap_cols = SampleDesctiptors %>%
  dplyr::select(condition, sample_id) %>%
  filter(condition %in% cond_levels)%>%
  mutate(condition = factor(condition, levels = cond_levels)) %>%
  mutate(sample_id = gsub("-", "_", sample_id, fixed = T)) %>%
  arrange(condition)

scaled_rld = scaled_rld %>%
  .[,as.character(heatmap_cols$sample_id)]
scaled_rld[scaled_rld > quantile(scaled_rld, 0.999)] = quantile(scaled_rld, 0.999)
scaled_rld[scaled_rld < quantile(scaled_rld, 0.001)] = quantile(scaled_rld, 0.001)

hm_rows = rld_rows[row_select,] %>%
  mutate(Pc.M = case_when(class == "Pc-I" ~ "Pc-M",
                          T ~ "a"),
         Pc.H = case_when(class == "Pc-H" ~ "Pc-H",
                          T ~ "a"))

hm_rows2 = rld_cluster %>%
  left_join(rld_cluster %>%
              left_join(rld_rows, by = "gene_id") %>%
              dplyr::select(cl, class) %>%
              group_by(cl, class) %>%
              summarise(q = n()) %>%
              group_by(cl) %>%
              mutate(k = sum(q)) %>%
              left_join(hyper_gs, by = "class") %>%
              mutate(p_val = phyper(q-1, m, n, k, lower.tail = F, log.p = F)) %>%
              mutate(p.adj = -log10(p.adjust(p_val))) %>%
              mutate(effect_size = (q/k)/(m/(m+n))) %>%
              mutate(class = recode(class, !!!class_recode)) %>%
              pivot_wider(cl, names_from = class, values_from = effect_size, values_fill = 1) %>%
              select(!non.Pc),
            by = "cl")

### FigSwitchGal4FullHeatmap.pdf

ht_list = Heatmap(matrix = scaled_rld,
                  column_order = 1 : dim(scaled_rld)[2],
                  column_split = gsub("_", "\n", heatmap_cols$condition, fixed = T),
                  show_column_names = FALSE,
                  col = colorRampPalette(carto_pal(n = 7, name = "ArmyRose"))(10),
                  show_row_names = FALSE, name = "row\nz score")

ht_list = ht_list + Heatmap(hm_rows$Pc.M, 
                            name = "Pc-M", 
                            col = c("white", "black"),
                            width = unit(.3, "cm"))

ht_list = ht_list + Heatmap(hm_rows2$Pc.M, 
                            name = "effect\nsizenPc-M",
                            col = colorRamp2(c(1, 4), c("#FFFFFF", "#ffae34")),
                            width = unit(.2, "cm"),
                            show_column_names = F)

ht_list = ht_list + Heatmap(hm_rows$Pc.H, 
                            name = "Pc-H", 
                            col = c("white", "black"),
                            width = unit(.3, "cm"),
                            show_heatmap_legend = F)

ht_list = ht_list + Heatmap(hm_rows2$Pc.H, 
                            name = "effect\nsize\nPc-H",
                            col = colorRamp2(c(1, 15), c("#FFFFFF", "#ef6f6a")),
                            width = unit(.2, "cm"),
                            show_column_names = F)


t = draw(ht_list, row_split = rld_cluster$cl,
         show_row_dend = F)

### FigSwitchGal4FullHeatmapTable.png

rld_cluster %>%
  left_join(rld_rows, by = "gene_id") %>%
  dplyr::select(cl, class) %>%
  group_by(cl, class) %>%
  summarise(q = n()) %>%
  group_by(cl) %>%
  mutate(k = sum(q)) %>%
  left_join(hyper_gs, by = "class") %>%
  mutate(p_val = phyper(q-1, m, n, k, lower.tail = F, log.p = F)) %>%
  mutate(p.adj = -log10(p.adjust(p_val))) %>%
  mutate(effect_size = (q/k)/(m/(m+n))) %>%
  mutate(class = recode(class, !!!class_recode)) %>%
  filter(p.adj > -log10(0.05)) %>%
  arrange(-p.adj) %>%
  kable(caption = "p_vals", format.args = list(big.mark = ",")) %>% 
  kable_classic()

### FigSwitchGal4FullHeatmapBoxplot.pdf

as.data.frame(scaled_rld) %>%
  bind_cols(rld_cluster %>% dplyr::select(gene_id, cl)) %>%
  pivot_longer(!gene_id:cl, names_to = "sample_id") %>%
  left_join(heatmap_cols, by = "sample_id") %>%
  group_by(condition, gene_id, cl) %>%
  summarise(mean_z = mean(value)) %>%
  mutate(cl = factor(cl, levels = rev(as.numeric(t@ht_list[[1]]@row_title)))) %>%
  ggplot(aes(as.factor(cl), mean_z, fill = condition)) + 
  geom_boxplot(outlier.alpha = 0) +
  scale_fill_brewer(palette = "Paired", name = "Condition") +
  theme_bw() + ylab("mean z-score") + xlab("Cluster") +
  coord_flip()


### FigSwitchGal4TruncHeatmapV1.pdf

ht_list2 = Heatmap(matrix = scaled_rld[rld_cluster$cl %in% c(3,6,4),],
                  column_order = 1 : dim(scaled_rld)[2],
                  column_split = gsub("_", "\n", heatmap_cols$condition, fixed = T),
                  show_column_names = FALSE,
                  col = colorRampPalette(carto_pal(n = 7, name = "ArmyRose"))(10),
                  show_row_names = FALSE, name = "row\nz score")

ht_list2 = ht_list2 + Heatmap(hm_rows$Pc.M[rld_cluster$cl %in% c(3,6,4)], 
                            name = "Pc-M", 
                            col = c("white", "black"),
                            width = unit(.3, "cm"))

ht_list2 = ht_list2 + Heatmap(hm_rows2$Pc.M[rld_cluster$cl %in% c(3,6,4)], 
                            name = "effect\nsizenPc-M",
                            col = colorRamp2(c(1, 4), c("#FFFFFF", "#ffae34")),
                            width = unit(.2, "cm"),
                            show_column_names = F)

ht_list2 = ht_list2 + Heatmap(hm_rows$Pc.H[rld_cluster$cl %in% c(3,6,4)], 
                            name = "Pc-H", 
                            col = c("white", "black"),
                            width = unit(.3, "cm"),
                            show_heatmap_legend = F)

ht_list2 = ht_list2 + Heatmap(hm_rows2$Pc.H[rld_cluster$cl %in% c(3,6,4)], 
                            name = "effect\nsize\nPc-H",
                            col = colorRamp2(c(1, 15), c("#FFFFFF", "#ef6f6a")),
                            width = unit(.2, "cm"),
                            show_column_names = F)

draw(ht_list2, row_split = rld_cluster$cl[rld_cluster$cl %in% c(3,6,4)],
     show_row_dend = F)

### FigSwitchGal4TruncHeatmapV2.pdf

ht_list2 = Heatmap(matrix = scaled_rld[rld_cluster$cl %in% c(3,6),],
                   column_order = 1 : dim(scaled_rld)[2],
                   column_split = gsub("_", "\n", heatmap_cols$condition, fixed = T),
                   show_column_names = FALSE,
                   col = colorRampPalette(carto_pal(n = 7, name = "ArmyRose"))(10),
                   show_row_names = FALSE, name = "row\nz score")

ht_list2 = ht_list2 + Heatmap(hm_rows$Pc.M[rld_cluster$cl %in% c(3,6)], 
                              name = "Pc-M", 
                              col = c("white", "black"),
                              width = unit(.3, "cm"))

ht_list2 = ht_list2 + Heatmap(hm_rows2$Pc.M[rld_cluster$cl %in% c(3,6)], 
                              name = "effect\nsizenPc-M",
                              col = colorRamp2(c(1, 4), c("#FFFFFF", "#ffae34")),
                              width = unit(.2, "cm"),
                              show_column_names = F)

ht_list2 = ht_list2 + Heatmap(hm_rows$Pc.H[rld_cluster$cl %in% c(3,6)], 
                              name = "Pc-H", 
                              col = c("white", "black"),
                              width = unit(.3, "cm"),
                              show_heatmap_legend = F)

ht_list2 = ht_list2 + Heatmap(hm_rows2$Pc.H[rld_cluster$cl %in% c(3,6)], 
                              name = "effect\nsize\nPc-H",
                              col = colorRamp2(c(1, 15), c("#FFFFFF", "#ef6f6a")),
                              width = unit(.2, "cm"),
                              show_column_names = F)

draw(ht_list2, row_split = rld_cluster$cl[rld_cluster$cl %in% c(3,6)],
     show_row_dend = F)
