library(ComplexHeatmap)
library(rcartocolor)
library(kableExtra)
library(Rfast)
library(edgeR)
library(openxlsx)
library(rtracklayer)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggthemes)

load("data/fig1/resting_H3K27me3.Rdata")

# import RNAi data --------------------------------------------------------

GFF = readGFF("/Users/streeck/Genomes/DmelBDGP6.91/Drosophila_melanogaster.BDGP6.91.gtf", filter = list(type = "gene"))

CountTable = read.delim("data/fig3/2021-03-20featureCounts.count", skip = 1)
CountTable = CountTable[,c(1,7:dim(CountTable)[2])] 

colnames(CountTable)[2:dim(CountTable)[2]] = sapply(strsplit(colnames(CountTable)[2:dim(CountTable)[2]], ".", fixed = T), function(x) x[length(x)-1])

CountTable = CountTable %>%
  dplyr::select(!V_4905)


SampleDesctiptors = read.xlsx("data/fig3/2021-03-21MappedSamples.xlsx") %>%
  filter(sample_id != "V_4905")

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

## select 1000 genes with the highes variance
rv <- rowVars(assay(rld))
select <- order(rv, decreasing = TRUE)[seq_len(3000)]

## apply PCA analysis
pca <- prcomp(t(assay(rld[select, ])))
pcaPlot = as.data.frame(pca$x)
pcaPlot = cbind(groups, pcaPlot)
pcaPlot = pcaPlot[,1:11]
pcaPlot$ID = sapply(strsplit(row.names(pcaPlot), "_", fixed = T), function(x) x[1])
var.explained = round(100*pca$sdev^2/sum(pca$sdev^2), 2)


## contrasts and fits
contrasts = makeContrasts(luc = luc_septic - luc_untreated,
                          nej = nej_septic - nej_untreated,
                          Set1 = Set1_septic - Set1_untreated,
                          trr = trr_septic - trr_untreated,
                          trx = trx_septic - trx_untreated,
                          Utx = Utx_septic - Utx_untreated,
                          nej_baseline = nej_untreated - luc_untreated,
                          Set1_baseline =  Set1_untreated - luc_untreated,
                          trr_baseline = trr_untreated - luc_untreated,
                          trx_baseline = trx_untreated - luc_untreated,
                          Utx_baseline = Utx_untreated - luc_untreated,
                          nej_septic = nej_septic - luc_septic,
                          Set1_septic =  Set1_septic - luc_septic,
                          trr_septic = trr_septic - luc_septic,
                          trx_septic = trx_septic - luc_septic,
                          Utx_septic = Utx_septic - luc_septic,
                          levels = designMatrix)

fit = glmQLFit(DGEobject, designMatrix, robust = T)
QFresults = list()
for (i in 1:dim(contrasts)[2]) {
  QFresults[[colnames(contrasts)[i]]] = glmQLFTest(fit, contrast = contrasts[,i])
}


# Heatmap1 ----------------------------------------------------------------


hyper_gs = DGEobject$genes %>%
  dplyr::select(gene_id) %>%
  left_join(resting_H3K27me3 %>%
              dplyr::select(Geneid, class), 
            by = c(gene_id = "Geneid")) %>%
  group_by(class) %>%
  dplyr::summarize(m = n()) %>%
  dplyr::mutate(n = dim(DGEobject$genes)[1] - m)

cond_levels = c("luc", "trx", "Set1", "trr", "Utx", "nej")

heatmap_cols = SampleDesctiptors %>%
  mutate(treat = sapply(strsplit(condition, split = "_"), function(x){x[2]})) %>%
  filter(treat == "untreated") %>%
  dplyr::select(condition, sample_id) %>%
  mutate(condition = factor(
    sapply(strsplit(condition, split = "_"), function(x){x[1]}), 
    levels = cond_levels)) %>%
  filter(condition != "nej") %>%
  arrange(condition)

rld_assay = SummarizedExperiment::assay(rld) %>%
  as.data.frame() %>%
  dplyr::select(heatmap_cols$sample_id)


rld_rows = DGEobject$genes %>%
  dplyr::select(gene_id) %>%
  left_join(resting_H3K27me3 %>%
              dplyr::select(Geneid, class), 
            by = c(gene_id = "Geneid")) 




l = list()
for (i in 8:11) {
  l[[i]] = log10(QFresults[[i]]$table$PValue)*as.numeric(decideTests.DGELRT(QFresults[[i]]))
}

metric = as.matrix(
  do.call(cbind, l)
)

metric = rowMaxs(metric, value = T)

n_top = 1000
nk = 6


row_select <- order(metric, decreasing = TRUE)[1:n_top]


rld_selected_counts = rld_assay[row_select,]
scaled_rld = rld_selected_counts %>% 
  t %>%
  scale %>%
  t

rld_cluster = DGEobject$genes %>%
  .[row_select,]

rld_cluster$cl = kmeans(scaled_rld, centers = nk, algorithm = 'MacQueen', iter.max = 50, nstart = 10)$cluster

scaled_rld = scaled_rld %>%
  .[,as.character(heatmap_cols$sample_id)]
scaled_rld[scaled_rld > quantile(scaled_rld, 0.999)] = quantile(scaled_rld, 0.999)
scaled_rld[scaled_rld < quantile(scaled_rld, 0.001)] = quantile(scaled_rld, 0.001)

hm_rows = rld_rows[row_select,] %>%
  mutate(Pc.M = case_when(class == "Pc-I" ~ "Pc-M",
                          T ~ "a"),
         Pc.H = case_when(class == "Pc-H" ~ "Pc-H",
                          T ~ "a"))

class_recode = c("non-Pc" = "non.Pc", "Pc-I" = "Pc.M", "Pc-H" = "Pc.H")

hm_rows2 = rld_cluster %>%
  left_join(rld_cluster %>%
              left_join(rld_rows, by = "gene_id") %>%
              dplyr::select(cl, class) %>%
              group_by(cl, class) %>%
              dplyr::summarise(q = n()) %>%
              group_by(cl) %>%
              dplyr::mutate(k = sum(q)) %>%
              left_join(hyper_gs, by = "class") %>%
              dplyr::mutate(p_val = phyper(q-1, m, n, k, lower.tail = F, log.p = F)) %>%
              dplyr::mutate(p.adj = -log10(p.adjust(p_val))) %>%
              dplyr::mutate(effect_size = (q/k)/(m/(m+n))) %>%
              dplyr::mutate(class = recode(class, !!!class_recode)) %>%
              pivot_wider(cl, names_from = class, values_from = effect_size, values_fill = 1) %>%
              dplyr::select(!non.Pc),
            by = "cl")

### FigRNAiBaselineHeatmap.pdf

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

### FigRNAiBaselineHeatmapTable.png

rld_cluster %>%
  left_join(rld_rows, by = "gene_id") %>%
  dplyr::select(cl, class) %>%
  group_by(cl, class) %>%
  dplyr::summarise(q = n()) %>%
  group_by(cl) %>%
  dplyr::mutate(k = sum(q)) %>%
  left_join(hyper_gs, by = "class") %>%
  dplyr::mutate(p_val = phyper(q-1, m, n, k, lower.tail = F, log.p = F)) %>%
  dplyr::mutate(p.adj = -log10(p.adjust(p_val))) %>%
  dplyr::mutate(effect_size = (q/k)/(m/(m+n))) %>%
  dplyr::mutate(class = recode(class, !!!class_recode)) %>%
  filter(p.adj > -log10(0.05)) %>%
  arrange(-p.adj) %>%
  kable(caption = "p_vals", format.args = list(big.mark = ",")) %>% 
  kable_classic()

### FigRNAiBaselineHeatmapBoxplot.pdf

as.data.frame(scaled_rld) %>%
  bind_cols(rld_cluster %>% dplyr::select(gene_id, cl)) %>%
  pivot_longer(!gene_id:cl, names_to = "sample_id") %>%
  left_join(heatmap_cols, by = "sample_id") %>%
  group_by(condition, gene_id, cl) %>%
  dplyr::summarise(mean_z = mean(value)) %>%
  dplyr::mutate(cl = factor(cl, levels = rev(as.numeric(t@ht_list[[1]]@row_title)))) %>%
  ggplot(aes(as.factor(cl), mean_z, fill = condition)) + 
  geom_boxplot(outlier.alpha = 0) +
  scale_fill_tableau(name = "Condition") +
  theme_bw() + ylab("mean z-score") + xlab("Cluster") +
  coord_flip()


# Heatmap2 ----------------------------------------------------------------

hyper_gs = DGEobject$genes %>%
  dplyr::select(gene_id) %>%
  left_join(resting_H3K27me3 %>%
              dplyr::select(Geneid, class), 
            by = c(gene_id = "Geneid")) %>%
  group_by(class) %>%
  dplyr::summarize(m = n()) %>%
  dplyr::mutate(n = dim(DGEobject$genes)[1] - m)


cond_levels = c("luc", "trx", "Set1", "trr", "Utx", "nej")

heatmap_cols = SampleDesctiptors %>%
  mutate(treat = sapply(strsplit(condition, split = "_"), function(x){x[2]})) %>%
  mutate(treat = factor(treat, levels = c("untreated", "septic"))) %>%
  dplyr::select(condition, sample_id, treat) %>%
  mutate(condition = factor(
    sapply(strsplit(condition, split = "_"), function(x){x[1]}), 
    levels = cond_levels)) %>%
  filter(condition != "nej") %>%
  arrange(condition, treat) %>%
  mutate(split = factor(paste0(condition, "\n", treat),
                        levels = paste0(rep(c("luc", "Set1", "trx", "trr", "Utx"), each = 2), "\n",
                                        rep(c("untreated", "septic"), 5))))

rld_assay = SummarizedExperiment::assay(rld) %>%
  as.data.frame() %>%
  dplyr::select(heatmap_cols$sample_id)


rld_rows = DGEobject$genes %>%
  dplyr::select(gene_id) %>%
  left_join(resting_H3K27me3 %>%
              dplyr::select(Geneid, class), 
            by = c(gene_id = "Geneid")) 

l1 = list()
l2 = list()
l3 = rep(0, length(QFresults$luc$genes$gene_id))
for (i in 1:6) {
  l1[[i]] = QFresults[[i]]$table$PValue
  l2[[i]] = QFresults[[i]]$table$logFC
  l3 = l3 + decideTests.DGELRT(QFresults[[i]], lfc = 1)
}

Infection_ptable = cbind(QFresults$luc$genes, 
                         med_p = rowMins(as.matrix(do.call(cbind, l1)), value = T),
                         med_lfc = rowmeans(as.matrix(do.call(cbind, l2))))

metric = -log(Infection_ptable$med_p) * sign(Infection_ptable$med_lfc)

n_top = 1200
nk = 5


row_select <- order(metric, decreasing = TRUE)[1:n_top]


rld_selected_counts = rld_assay[row_select,]
scaled_rld = rld_selected_counts %>% 
  t %>%
  scale %>%
  t

rld_cluster = DGEobject$genes %>%
  .[row_select,]

rld_cluster$cl = kmeans(scaled_rld, centers = nk, algorithm = 'MacQueen', iter.max = 50, nstart = 10)$cluster

scaled_rld = scaled_rld %>%
  .[,as.character(heatmap_cols$sample_id)]
scaled_rld[scaled_rld > quantile(scaled_rld, 0.999)] = quantile(scaled_rld, 0.999)
scaled_rld[scaled_rld < quantile(scaled_rld, 0.001)] = quantile(scaled_rld, 0.001)


hm_rows = rld_rows[row_select,] %>%
  mutate(Pc.M = case_when(class == "Pc-I" ~ "Pc-M",
                          T ~ "a"),
         Pc.H = case_when(class == "Pc-H" ~ "Pc-H",
                          T ~ "a"))

class_recode = c("non-Pc" = "non.Pc", "Pc-I" = "Pc.M", "Pc-H" = "Pc.H")

hm_rows2 = rld_cluster %>%
  left_join(rld_cluster %>%
              left_join(rld_rows, by = "gene_id") %>%
              dplyr::select(cl, class) %>%
              group_by(cl, class) %>%
              dplyr::summarise(q = n()) %>%
              group_by(cl) %>%
              dplyr::mutate(k = sum(q)) %>%
              left_join(hyper_gs, by = "class") %>%
              dplyr::mutate(p_val = phyper(q-1, m, n, k, lower.tail = F, log.p = F)) %>%
              dplyr::mutate(p.adj = -log10(p.adjust(p_val))) %>%
              dplyr::mutate(effect_size = (q/k)/(m/(m+n))) %>%
              dplyr::mutate(class = recode(class, !!!class_recode)) %>%
              pivot_wider(cl, names_from = class, values_from = effect_size, values_fill = 1) %>%
              dplyr::select(!non.Pc),
            by = "cl")

### FigRNAiFullHeatmap.pdf

ht_list = Heatmap(matrix = scaled_rld,
                  column_order = 1 : dim(scaled_rld)[2],
                  column_split = heatmap_cols$split,
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

### FigRNAiFullHeatmapTable.png

rld_cluster %>%
  left_join(rld_rows, by = "gene_id") %>%
  dplyr::select(cl, class) %>%
  group_by(cl, class) %>%
  dplyr::summarise(q = n()) %>%
  group_by(cl) %>%
  dplyr::mutate(k = sum(q)) %>%
  left_join(hyper_gs, by = "class") %>%
  dplyr::mutate(p_val = phyper(q-1, m, n, k, lower.tail = F, log.p = F)) %>%
  dplyr::mutate(p.adj = -log10(p.adjust(p_val))) %>%
  dplyr::mutate(effect_size = (q/k)/(m/(m+n))) %>%
  dplyr::mutate(class = recode(class, !!!class_recode)) %>%
  filter(p.adj > -log10(0.05)) %>%
  arrange(-p.adj) %>%
  kable(caption = "p_vals", format.args = list(big.mark = ",")) %>% 
  kable_classic()

### FigRNAiFullHeatmapBoxplot.pdf

as.data.frame(scaled_rld) %>%
  bind_cols(rld_cluster %>% dplyr::select(gene_id, cl)) %>%
  pivot_longer(!gene_id:cl, names_to = "sample_id") %>%
  left_join(heatmap_cols, by = "sample_id") %>%
  group_by(split, gene_id, cl) %>%
  dplyr::summarise(mean_z = mean(value)) %>%
  dplyr::mutate(cl = factor(cl, levels = rev(as.numeric(t@ht_list[[1]]@row_title)))) %>%
  ggplot(aes(as.factor(cl), mean_z, fill = split)) + 
  geom_boxplot(outlier.alpha = 0) +
  scale_fill_tableau(name = "Condition", palette = "Tableau 20") +
  theme_bw() + ylab("mean z-score") + xlab("Cluster") +
  coord_flip()



# Enrichment plot ---------------------------------------------------------


hypMod = function(qf){
  df = data.frame(gene_id = qf$genes$gene_id,
                  de = as.numeric(decideTests.DGELRT(qf, lfc = log2(1.5))))
  df %>%
    left_join(resting_H3K27me3 %>%
                replace(is.na(.), " ") %>%
                dplyr::select(c(Geneid, class)),
              by = c(gene_id = "Geneid")) %>%
    group_by(de, class) %>%
    dplyr::summarise(q = n()) %>%
    group_by(de) %>%
    dplyr::mutate(k = sum(q)) %>%
    group_by(class) %>%
    dplyr::mutate(m = sum(q)) %>%
    dplyr::mutate(n = dim(DGEobject$genes)[1] - m) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(p_val = phyper(q-1, m, n, k, lower.tail = F, log.p = T)) %>%
    dplyr::mutate(effect_size = (q/k)/(m/(m+n))) %>%
    dplyr::mutate(p = exp(p_val)) %>%
    filter(de != 0) %>%
    filter(class == "Pc-I") %>%
    arrange(p)
}


hyped = sapply(QFresults[8:11], simplify = F, hypMod, USE.NAMES = T)


y_coding = c(Set1_baseline = "Set1", trr_baseline = "trr", 
             trx_baseline = "trx", Utx_baseline = "Utx")


# FigRNAiEnrichment.pdf

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



# KD efficiency -----------------------------------------------------------


GFF_exon = readGFF("/Users/streeck/Genomes/DmelBDGP6.91/Drosophila_melanogaster.BDGP6.91.gtf", filter = list(type = "exon"))

RNAi_list = c("trx", "trr", "Set1", "Utx", "nej")

GFF_exon = GFF_exon %>%
  filter(gene_name %in% RNAi_list)

rtracklayer::export.gff2(GFF_exon, "RNAi_gene_selection.gtf")

t = GFF_exon %>%
  dplyr::select(gene_name, exon_id)

CountTableRNAi = read.delim("data/fig3/RNAi_level_featureCounts.count", skip = 1) %>%
  dplyr::select(!Chr:Strand) %>%
  pivot_longer(!Geneid:Length) %>%
  dplyr::mutate(name = sapply(strsplit(name, ".", fixed = T), function(x){x[6]})) %>%
  left_join(GFF_exon %>%
              dplyr::select(gene_name, exon_id) %>%
              distinct(),
            by = c("Geneid" = "exon_id")) %>%
  left_join(DGEobject$samples %>%
              dplyr::mutate(name = row.names(DGEobject$samples)) %>%
              dplyr::mutate(RNAi = sapply(strsplit(as.character(group), "_", fixed = T), function(x){x[1]})),
            by = "name") %>%
  filter(!is.na(RNAi)) %>%
  filter(RNAi != "nej") %>%
  filter(gene_name != "nej") %>%
  dplyr::mutate(eff_cpm = 1e6*value/(Length*lib.size*norm.factors)) %>%
  filter(value != 0) %>%
  group_by(name, gene_name, RNAi) %>%
  filter(between(eff_cpm, quantile(eff_cpm, 0.251), quantile(eff_cpm, 0.75))) %>%
  dplyr::summarise(mean_tpm = exp(mean(log(eff_cpm))))

kd_efficiency = CountTableRNAi %>%
  filter(RNAi == "luc" | RNAi == gene_name) %>%
  dplyr::mutate(log_tpm = log(mean_tpm)) %>%
  group_by(gene_name, RNAi) %>%
  dplyr::summarise(mean_log_tpm = mean(log_tpm)) %>%
  dplyr::mutate(treat = case_when(RNAi == "luc" ~ "control",
                           T ~ "knock_down")) %>%
  pivot_wider(!RNAi, names_from = treat, values_from = mean_log_tpm) %>%
  dplyr::mutate(kd_eff = exp(control - knock_down)) %>%
  dplyr::mutate(lab = paste0(round(kd_eff, 1), "-fold k.d."))

### FigRNAiKDefficiency.pdf

ggplot(CountTableRNAi, aes(gene_name, mean_tpm, fill = RNAi)) +
  geom_boxplot() + theme_bw() + 
  ggthemes::scale_fill_tableau(palette = "Superfishel Stone") +
  scale_y_log10() + xlab("Gene") + ylab("TPM") + 
  geom_text(data = kd_efficiency, inherit.aes = F, 
            aes(x = gene_name, y = 0.28, label = lab))



# Violins ----------------------------------------------------------------

Infection_ptable$any.de = as.numeric(l3 != 0) * sign(l3)
inf_up_set = Infection_ptable$gene_id[Infection_ptable$any.de == 1]


table_for_violins = as.data.frame(sapply(QFresults[7:11], function(x){x$table$logFC}))
table_for_violins$gene_id = QFresults$luc$genes$gene_id

table_for_violins = table_for_violins %>%
  filter(gene_id %in% inf_up_set) %>%
  pivot_longer(!gene_id) %>%
  filter(name != "nej_baseline") %>%
  left_join(resting_H3K27me3 %>%
              dplyr::select(Geneid, class),
            by = c("gene_id" = "Geneid")) %>%
  dplyr::mutate(chrom_state = case_when(class == "Pc-I" ~ "Pc-M",
                                 T ~ "non-Pc"))


wcx = function(x){wilcox.test(x, alternative = "less")$p.value}

stat_test = table_for_violins %>%
  filter(chrom_state == "Pc-M")%>%
  group_by(name) %>%
  dplyr::summarise(wilcox = wcx(value)) %>%
  filter(wilcox <= 0.05) %>%
  dplyr::mutate(lab = paste0("p = ", sprintf(wilcox, fmt = "%#.2e"))) %>%
  dplyr::mutate(x = sapply(strsplit(name, "_", fixed = T), function(x){x[1]}))

### FigRNAiViolinsResting.pdf

table_for_violins %>%
  filter(chrom_state == "Pc-M") %>%
  dplyr::mutate(x = sapply(strsplit(name, "_", fixed = T), function(x){x[1]})) %>%
  ggplot(aes(x, value, fill = x)) + 
  geom_violin(scale = "width") + 
  geom_boxplot(width = .2, outlier.alpha = 0, notch = T) + theme_bw() + ylab("log2 Fold Change") + 
  theme(axis.title.x = element_blank()) +
  ggthemes::scale_fill_tableau(palette = "Superfishel Stone", name = "RNAi") +
  scale_y_continuous(limits = c(-8, 7)) +
  geom_text(data = stat_test, inherit.aes = F, mapping = aes(x = x, y = 7, label = lab),
            nudge_x = 0) +
  ggtitle("Mean Fold Change in Resting Plasmatocytes")

###

table_for_violins = as.data.frame(sapply(QFresults[12:16], function(x){x$table$logFC}))
table_for_violins$gene_id = QFresults$luc$genes$gene_id

table_for_violins = table_for_violins %>%
  filter(gene_id %in% inf_up_set) %>%
  pivot_longer(!gene_id) %>%
  filter(name != "nej_septic") %>%
  left_join(resting_H3K27me3 %>%
              dplyr::select(Geneid, class),
            by = c("gene_id" = "Geneid")) %>%
  dplyr::mutate(chrom_state = case_when(class == "Pc-I" ~ "Pc-M",
                                        T ~ "non-Pc"))


wcx = function(x){wilcox.test(x, alternative = "less")$p.value}

stat_test = table_for_violins %>%
  filter(chrom_state == "Pc-M")%>%
  group_by(name) %>%
  dplyr::summarise(wilcox = wcx(value)) %>%
  filter(wilcox <= 0.05) %>%
  dplyr::mutate(lab = paste0("p = ", sprintf(wilcox, fmt = "%#.2e"))) %>%
  dplyr::mutate(x = sapply(strsplit(name, "_", fixed = T), function(x){x[1]}))

## FigRNAiViolinsSeptic.pdf

table_for_violins %>%
  filter(chrom_state == "Pc-M") %>%
  dplyr::mutate(x = sapply(strsplit(name, "_", fixed = T), function(x){x[1]})) %>%
  ggplot(aes(x, value, fill = x)) + 
  geom_violin(scale = "width") + 
  geom_boxplot(width = .2, outlier.alpha = 0, notch = T) + theme_bw() + ylab("log2 Fold Change") + 
  theme(axis.title.x = element_blank()) +
  ggthemes::scale_fill_tableau(palette = "Superfishel Stone", name = "RNAi") +
  scale_y_continuous(limits = c(-8, 7)) +
  geom_text(data = stat_test, inherit.aes = F, mapping = aes(x = x, y = 7, label = lab),
            nudge_x = 0) +
  ggtitle("Mean Fold Change after Septic Injury")


# PCA ---------------------------------------------------------------------

library(RColorBrewer)
library(plot3D)

pcaPlot$groups = factor(pcaPlot$groups)

t20 = ggthemes_data$tableau$`color-palettes`$regular$`Tableau 20`$value

pcaPlot = pcaPlot %>%
  filter(!(groups %in% c("nej_untreated", "nej_septic"))) %>%
  mutate(groups = droplevels(groups))


## FigRNAiPCA.pdf

scatter3D(pcaPlot$PC1, pcaPlot$PC2, pcaPlot$PC3, 
          bgvar = as.numeric(as.factor(pcaPlot$groups)), 
          col = "black", bg = t20[as.numeric(pcaPlot$groups)],
          theta = 150, phi = 10, xlab = paste0("PC1: ", var.explained[1], "% var"),
          ylab = paste0("PC2: ", var.explained[2], "% var"), zlab = paste0("PC3: ", var.explained[3], "% var"), 
          bty ="b2", lwd = 1.5, cex = 2.5, colkey = F, alpha = .8, pch = 21)
legend("right",levels(pcaPlot$groups), 
       pt.bg = t20, pt.lwd = 1.2, cex = 1, pt.cex = 2.5, inset = c(0.07,0.07), pch = 21)


# Volcanos ----------------------------------------------------------------

##

log_ratio = function(chip, input){
  chip = (chip)/sum(chip)
  input = (input)/sum(input)
  return(log2(chip) - log2(input))
}

resting_H3K27me3 = resting_H3K27me3 %>%
  mutate(lr_A = log_ratio(H3K27me3_A, input_A),
         lr_B = log_ratio(H3K27me3_B, input_B)) %>%
  mutate(h3k27me3_mean_lr = (lr_A + lr_B)/2)

qf.volcano = function(qf, name){

  
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
    filter(de == -1) %>%
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
    ggtitle(name)
  
  
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
}

### FigRNAiSet1baselineVolcano.pdf
qf.volcano(QFresults$Set1_baseline, "Set1 baseline")

### FigRNAitrrbaselineVolcano.pdf
qf.volcano(QFresults$trr_baseline, "trr baseline")

### FigRNAitrxbaselineVolcano.pdf
qf.volcano(QFresults$trx_baseline, "trx baseline")

### FigRNAiUtxbaselineVolcano.pdf
qf.volcano(QFresults$Utx_baseline, "Utx baseline")



