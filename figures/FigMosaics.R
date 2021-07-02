library(ggplot2)
library(plotly)
library(edgeR)
library(openxlsx)
library(rtracklayer)
library(biomaRt)
library(DESeq2)
library(tidyverse)
library(knitr)
library(kableExtra)
library(ggrepel)
library(ggthemes)
library(ComplexHeatmap)
library(rcartocolor)
library(circlize)
library(tibble)
library(dplyr)
library(reshape2)
library(data.table)
library(viridis)
library(Rfast)
library(ggnewscale)
library(eulerr)

GFF = readGFF("/Users/streeck/Genomes/DmelBDGP6.91/Drosophila_melanogaster.BDGP6.91.gtf", filter = list(type = "gene"))

CountTable2537 = read.delim("data/fig3/2537.count", skip = 1)
CountTable2537 = CountTable2537[,c(1,7:dim(CountTable2537)[2])]
CountTable2537 = CountTable2537[,c(T, grepl("NoDup", colnames(CountTable2537)[2:dim(CountTable2537)[2]]))]
colnames(CountTable2537)[2:dim(CountTable2537)[2]] = sapply(
  strsplit(colnames(CountTable2537)[2:dim(CountTable2537)[2]], "_", fixed = T), 
  function(x){paste0("Lib_2537_", x[2])})

CountTable3506 = read.delim("data/fig3/3506.count", skip = 1)
CountTable3506 = CountTable3506[,c(1,7:dim(CountTable3506)[2])]
colnames(CountTable3506)[2:dim(CountTable3506)[2]] = sapply(
  strsplit(colnames(CountTable3506)[2:dim(CountTable3506)[2]], "_", fixed = T), 
  function(x){paste0("Lib_3506_", x[2])})

CountTable = CountTable2537 %>%
  left_join(CountTable3506, by = "Geneid")

SampleDesctiptors = read.xlsx("data/fig3/MosaicSamples.xlsx") %>%
  mutate(bam = paste0("Lib_", Project, "_", Lib_id))
SampleDesctiptors$bam = factor(SampleDesctiptors$bam, levels = colnames(CountTable)[2:dim(CountTable)[2]])
SampleDesctiptors = SampleDesctiptors[order(SampleDesctiptors$bam),]

groups = SampleDesctiptors$genotype
designMatrix = model.matrix(~ 0 + groups)
colnames(designMatrix) = substr(colnames(designMatrix), 7, 100)
colnames(designMatrix) = gsub(" ", "_", colnames(designMatrix), fixed = T)

DGEobject = DGEList(CountTable[,2:dim(CountTable)[2]], group = groups, genes = CountTable[,1])
DGEobject = DGEobject[rowSums(cpm(DGEobject)>1) >= max(table(groups)), , keep.lib.sizes=FALSE]
DGEobject$genes = GFF[match(DGEobject$genes$genes, GFF$gene_id), c(9,10,12)]
DGEobject = calcNormFactors(DGEobject)
DGEobject = estimateDisp(DGEobject, design = designMatrix, robust = T)

dds = DESeqDataSetFromMatrix(countData = as.matrix(DGEobject$counts), colData = data.frame(row.names=row.names(DGEobject$samples), gr = DGEobject$samples$group), design = ~ gr)
dds = dds[rowSums(counts(dds)) > 1,]
rld = vst(dds, blind = FALSE)



## contrasts and fits
contrasts = makeContrasts(H3.K27R = H3.K27R - H3.wt,
                          Suz12 = Suz12.4 - Suz12.wt,
                          Ez = Ez.731 - Ez.wt,
                          levels = designMatrix)

fit = glmQLFit(DGEobject, designMatrix, robust = T)
QFresults = list()
for (i in 1:dim(contrasts)[2]) {
  QFresults[[colnames(contrasts)[i]]] = glmQLFTest(fit, contrast = contrasts[,i])
}

load("data/fig1/resting_H3K27me3.Rdata")



# volcanos ----------------------------------------------------------------

make_volcano = function(qf){
  colored_volcano = cbind(qf$genes, qf$table)
  colored_volcano$de = as.factor(decideTests.DGELRT(qf))
  colored_volcano = colored_volcano %>% 
    left_join(resting_H3K27me3 %>%
                dplyr::select(Geneid, class),
              by = c("gene_id" = "Geneid")) %>%
    mutate(col = case_when(de != 0 & class == "non-Pc" ~ "non-Pc",
                           de != 0 & class == "Pc-I" ~ "Pc-M",
                           de != 0 & class == "Pc-H" ~ "Pc-H",
                           T ~ "not up-regulated"))
  
  colored_volcano$col = factor(colored_volcano$col,
                               levels = c("non-Pc", "Pc-M", "Pc-H", "not up-regulated"))
  
  volcano_labels = colored_volcano %>% 
    mutate(e = logFC^2 + log10(PValue)^2) %>%
    arrange(desc(e)) %>%
    filter(de == 1) %>%
    slice_head(n = 20)
  
  a = ggplot(colored_volcano, aes(logFC, -log10(PValue), color = col)) + 
    geom_point(alpha = .8, shape = 16) +
    geom_label_repel(data = volcano_labels, 
                     aes(logFC, -log10(PValue), label = gene_name, color = col), 
                     size = 2.5, fill = rgb(1,1,1,0.5), inherit.aes = F, show.legend = F) +
    scale_color_tableau(palette = "Superfishel Stone", name = "") + 
    theme_classic() + xlab("log2 fold change") + ylab("p-value") +
    scale_y_continuous(breaks = c(0, 4, 8, 12, 16), 
                       labels = c(1, expression(10^-4), expression(10^-8), expression(10^-12), expression(10^-16))) +
    guides(color = guide_legend(override.aes = list(shape = 16, size = 2))) +
    scale_x_continuous(limits = c(-max(abs(colored_volcano$logFC)), max(abs(colored_volcano$logFC)))) +
    theme(legend.position = c(0.15, 0.2),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_blank(), axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          legend.background = element_blank())
  
  
  boxplot_data = colored_volcano %>%
    mutate(direct = case_when(de == 1 ~ "up",
                              T ~ "not up")) %>%
    dplyr::select(logFC, direct, class) %>%
    mutate(class = recode(class, "Pc-I" = "Pc-M")) %>%
    mutate(class = factor(class, levels =  c("non-Pc", "Pc-M", "Pc-H")))
  
  boxplot_data$direct = factor(boxplot_data$direct,
                               levels = c("not up", "up"))
  
  
  b = ggplot(boxplot_data, aes(logFC, direct, fill = class)) + 
    geom_boxplot(outlier.alpha = 0, outlier.shape = 16) +
    scale_fill_tableau(palette = "Superfishel Stone", name = "") + 
    theme_classic() + xlab("log2 fold change") + 
    scale_y_discrete(name = "") + 
    scale_x_continuous(limits = c(-max(abs(colored_volcano$logFC)), max(abs(colored_volcano$logFC)))) +
    theme(legend.position = "none")
  
  gridExtra::grid.arrange(a, b, ncol = 1, heights = c(.8, .20))
}


### FigMosaicVolcanoH3K27R.pdf

make_volcano(QFresults$H3.K27R)

### FigMosaicVolcanoEz.pdf

make_volcano(QFresults$Ez)

### FigMosaicVolcanoSuz12.pdf

make_volcano(QFresults$Suz12)



# heatmap -----------------------------------------------------------------

hyper_gs = DGEobject$genes %>%
  dplyr::select(gene_id) %>%
  left_join(resting_H3K27me3 %>%
              dplyr::select(Geneid, class), 
            by = c("gene_id" = "Geneid")) %>%
  group_by(class) %>%
  dplyr::summarize(m = n()) %>%
  dplyr::mutate(n = dim(DGEobject$genes)[1] - m)

rld_assay = SummarizedExperiment::assay(rld) %>%
  as.data.frame 

rld_rows = DGEobject$genes %>%
  dplyr::select(gene_id) %>%
  left_join(resting_H3K27me3 %>%
              dplyr::select(Geneid, class), 
            by = c("gene_id" = "Geneid"))


cond_levels = c("H3.wt", "H3.K27R", "Ez.wt", "Ez.731", "Suz12.wt", "Suz12.4")



n_top = 500
nk = 4

l = list()
for (i in 1:3) {
  l[[i]] = -log10(QFresults[[i]]$table$PValue)
}

metric = as.matrix(
  do.call(cbind, l)
)

metric = rowMedians(metric)

row_select <- order(metric, decreasing = TRUE)[1:n_top]


heatmap_cols = SampleDesctiptors %>%
  dplyr::select(genotype, bam) %>%
  filter(genotype %in% cond_levels)%>%
  mutate(genotype = factor(genotype, levels = rev(cond_levels))) %>%
  arrange(genotype)

scaled_rld = rld_assay[row_select,] %>%
  .[,as.character(heatmap_cols$bam)]

### per condition centering
center_condition = function(count_matrix, grps){
  for(i in unique(grps)){
    selector = grps == i
    count_matrix[selector] = t(scale(t(
      count_matrix[selector]), center = T, scale = F))
  }
  return(as.matrix(count_matrix))
}


scaled_rld = center_condition(
  scaled_rld, substr(heatmap_cols$genotype, 1, 2)
)

rld_cluster = DGEobject$genes %>%
  .[row_select,]

rld_cluster$cl = kmeans(scaled_rld, centers = nk, algorithm = 'MacQueen', iter.max = 200, nstart = 40)$cluster

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




### FigMosaicHeatmap.pdf

ht_list = Heatmap(matrix = scaled_rld,
                  column_order = 1 : dim(scaled_rld)[2],
                  column_split = gsub("_", "\n", heatmap_cols$genotype, fixed = T),
                  show_column_names = FALSE,
                  col = colorRampPalette(carto_pal(n = 7, name = "ArmyRose"))(10),
                  show_row_names = FALSE, name = "mean centerd\nrow signal")


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
                            col = colorRamp2(c(1, 20), c("#FFFFFF", "#ef6f6a")),
                            width = unit(.2, "cm"),
                            show_column_names = F)


t = draw(ht_list, row_split = rld_cluster$cl, show_row_den = F)



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


### FigMosaicHeatmapBars.pdf

as.data.frame(scaled_rld) %>%
  bind_cols(rld_cluster %>% dplyr::select(gene_id, cl)) %>%
  pivot_longer(!gene_id:cl, names_to = "sample_id") %>%
  left_join(heatmap_cols, by = c("sample_id" = "bam")) %>%
  group_by(genotype, gene_id, cl) %>%
  dplyr::summarise(mean_z = mean(value)) %>%
  dplyr::mutate(cl = factor(cl, levels = rev(as.numeric(t@ht_list[[1]]@row_title)))) %>%
  ggplot(aes(as.factor(cl), mean_z, fill = genotype)) + 
  geom_boxplot(outlier.alpha = 0) +
  scale_fill_brewer(palette = "Paired", name = "Condition") +
  theme_bw() + ylab("mean z-score") + xlab("Cluster") +
  coord_flip()


# Enrichment Plot ---------------------------------------------------------


hypMod = function(qf){
  df = data.frame(gene_id = qf$genes$gene_id,
                  de = as.numeric(decideTests.DGELRT(qf, lfc = log2(2))))
  df %>%
    left_join(resting_H3K27me3 %>%
                dplyr::select(c(Geneid, class)),
              by = c("gene_id" = "Geneid")) %>%
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
    arrange(p)
}


hyped = sapply(c(QFresults), simplify = F, hypMod, USE.NAMES = T)

y_coding = c(H3.K27R = "H3 K27R", Suz12 = "Su(z)12", Ez = "E()z")

plot.df = rbindlist(hyped, idcol = "contrast") %>%
  mutate(p.adj = -log10(p.adjust(p))) %>%
  mutate(y = recode(contrast, !!!y_coding)) %>%
  mutate(y = factor(y, levels = y_coding)) %>%
  filter(class != "non-Pc" & de == 1)


### FigMosaicEnrichment.pdf

ggplot(plot.df , aes(y, k)) +
  geom_col(aes(y = 1), color = "black", fill = "#AAAAAA") + 
  geom_col(aes(y = q/k, fill = p.adj), color = "black") + 
  coord_flip() +
  theme_bw() + ylab("Fraction of Pc-M among up-regulated genes") +
  scale_fill_viridis(option = "viridis", limits = c(0, max(plot.df$p.adj)), name = "-log10\nadjusted\np-value") + 
  xlab("Contrast") +
  geom_text(aes(y = .85, label = paste0("n = ", k)))  +
  geom_hline(aes(yintercept = m/(m+n)), linetype = 5) +
  facet_wrap(facets = vars(class), ncol = 1, scales = "free_y", strip.position = "right")


rbindlist(hyped, idcol = "contrast") %>%
  mutate(log.p.adj = -log10(p.adjust(p))) %>%
  mutate(y = recode(contrast, !!!y_coding)) %>%
  mutate(y = factor(y, levels = y_coding)) %>%
  mutate(fraction = q/k) %>%
  mutate(background_frac = m/(n+m)) %>%
  filter(log.p.adj >= -log10(0.05)) %>%
  arrange(-log.p.adj) %>%
  dplyr::select(contrast, de, class, fraction, effect_size, log.p.adj) %>% kable() %>% 
  kableExtra::kable_styling(bootstrap_options = c("striped", "condensed"))



# vs trxG RNAi ------------------------------------------------------------

# RNAi import 

CountTableRNAi = read.delim("data/fig3/2021-03-20featureCounts.count", skip = 1)
CountTableRNAi = CountTableRNAi[,c(1,7:dim(CountTableRNAi)[2])]

colnames(CountTableRNAi)[2:dim(CountTableRNAi)[2]] = sapply(strsplit(colnames(CountTableRNAi)[2:dim(CountTableRNAi)[2]], ".", fixed = T), function(x) x[length(x)-1])


SampleDesctiptorsRNAi = read.xlsx("data/fig3/2021-03-21MappedSamples.xlsx")

groupsRNAi = SampleDesctiptorsRNAi$condition
designMatrixRNAi = model.matrix(~ 0 + groupsRNAi)
colnames(designMatrixRNAi) = substr(colnames(designMatrixRNAi), 11, 100)
colnames(designMatrixRNAi) = gsub(" ", "_", colnames(designMatrixRNAi), fixed = T)

DGEobjectRNAi = DGEList(CountTableRNAi[,2:dim(CountTableRNAi)[2]], group = groupsRNAi, genes = CountTableRNAi[,1])
DGEobjectRNAi = DGEobjectRNAi[rowSums(cpm(DGEobjectRNAi)>1) >= max(table(groupsRNAi)), , keep.lib.sizes=FALSE]
DGEobjectRNAi$genes = GFF[match(DGEobjectRNAi$genes$genes, GFF$gene_id), c(9,10,12)]
DGEobjectRNAi = calcNormFactors(DGEobjectRNAi)
DGEobjectRNAi = estimateDisp(DGEobjectRNAi, design = designMatrixRNAi, robust = T)

## contrasts and fits
contrastsRNAi = makeContrasts(luc = luc_septic - luc_untreated,
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
                              nej_full = nej_septic - nej_untreated - (luc_septic - luc_untreated),
                              Set1_full = Set1_septic - Set1_untreated- (luc_septic - luc_untreated),
                              trr_full = trr_septic - trr_untreated- (luc_septic - luc_untreated),
                              trx_full = trx_septic - trx_untreated- (luc_septic - luc_untreated),
                              Utx_full = Utx_septic - Utx_untreated- (luc_septic - luc_untreated),
                              levels = designMatrixRNAi)

fitRNAi = glmQLFit(DGEobjectRNAi, designMatrixRNAi, robust = T)
QFresultsRNAi = list()
for (i in 1:dim(contrastsRNAi)[2]) {
  QFresultsRNAi[[colnames(contrastsRNAi)[i]]] = glmQLFTest(fitRNAi, contrast = contrastsRNAi[,i])
}


## venns

Mosaic.de = data.frame(
  gene_id = QFresults$H3.K27R$genes$gene_id,
  H3K27R.de = as.numeric(decideTests.DGELRT(QFresults$H3.K27R, lfc = 1)),
  Suz12.de = as.numeric(decideTests.DGELRT(QFresults$Suz12, lfc = 1)),
  Ez.de = as.numeric(decideTests.DGELRT(QFresults$Ez, lfc = 1))
)

RNAi.de = data.frame(
  gene_id = QFresultsRNAi$luc$genes$gene_id
)


for(i in 7:16){
  RNAi.de = cbind(
    RNAi.de,
    as.numeric(decideTests.DGELRT(QFresultsRNAi[[i]], lfc = 1))
  )
  colnames(RNAi.de)[dim(RNAi.de)[2]] = paste0(names(QFresultsRNAi)[i], ".de")
}

Full.de = Mosaic.de %>%
  inner_join(RNAi.de, by = "gene_id") %>%
  mutate(mosaic.any.up = rowAny(as.matrix(.[,2:4] == 1))) %>%
  left_join(resting_H3K27me3 %>%
              dplyr::select(c(Geneid, class)),
            by = c("gene_id" = "Geneid"))


### FigMosaicVennVsRNAiPcMonly.pdf

fit = euler(Full.de %>%
              filter(class == "Pc-I") %>%
              mutate(trr_septic.de = trr_septic.de == -1,
                     trx_septic.de = trx_septic.de == -1,
                     Utx_septic.de = Utx_septic.de == -1) %>%
              dplyr::select(trr_septic.de, trx_septic.de, Utx_septic.de, mosaic.any.up))
plot(fit,
     quantities = list(type = c("counts"), cex = 1.5),
     fills = list(fill = ggthemes_data$tableau$`color-palettes`$regular$`Superfishel Stone`$value[1:4]),
     legend = list(side = "right", vgap = 2, cex = 1.2), main = "Pc-M")

### FigMosaicVennVsRNAinonPconly.pdf

fit = euler(Full.de %>%
              filter(class == "non-Pc") %>%
              mutate(trr_septic.de = trr_septic.de == -1,
                     trx_septic.de = trx_septic.de == -1,
                     Utx_septic.de = Utx_septic.de == -1) %>%
              dplyr::select(trr_septic.de, trx_septic.de, Utx_septic.de, mosaic.any.up))
plot(fit,
     quantities = list(type = c("counts"), cex = 1.5),
     fills = list(fill = ggthemes_data$tableau$`color-palettes`$regular$`Superfishel Stone`$value[1:4]),
     legend = list(side = "right", vgap = 2, cex = 1.2), main = "non-Pc")



# vs Inducible genes ------------------------------------------------------

## switch

CountTable_switch = read.delim("data/fig2/unstranded.count", skip = 1)
CountTable_switch = CountTable_switch[,c(1,7:dim(CountTable_switch)[2])]

colnames(CountTable_switch)[2:dim(CountTable_switch)[2]] = sapply(
  strsplit(colnames(CountTable_switch)[2:dim(CountTable_switch)[2]], ".", fixed = T), 
  function(x) paste0(x[length(x)-2], "_", x[length(x)-1]))

CountTable_switch[,dim(CountTable_switch)[2]+1] = read.delim("data/fig2/unstrandedfeatureCounts.count", skip = 1)[7]
colnames(CountTable_switch)[dim(CountTable_switch)[2]] = "NG_26757_24"


SampleDesctiptors_switch = read.xlsx("data/fig2/2021-04-11MappedSamples.xlsx")
SampleDesctiptors_switch$sample_num = sapply(strsplit(SampleDesctiptors_switch$sample_id, "_", fixed = T), function(x) x[2])
SampleDesctiptors_switch$sample_id = gsub("-", "_", SampleDesctiptors_switch$sample_id, fixed = T)

## exclude samples 1 and 8 (look bad on PCA)
CountTable_switch = CountTable_switch[, !colnames(CountTable_switch) %in% c("NG_26757_1", "NG_26757_8")]
SampleDesctiptors_switch = SampleDesctiptors_switch %>%
  filter(!sample_num %in% c("1", "8"))

groups_switch = SampleDesctiptors_switch$condition
designMatrix_switch = model.matrix(~ 0 + groups_switch)
colnames(designMatrix_switch) = substr(colnames(designMatrix_switch), 14, 100)
colnames(designMatrix_switch) = gsub(" ", "_", colnames(designMatrix_switch), fixed = T)

DGEobject_switch = DGEList(CountTable_switch[,2:dim(CountTable_switch)[2]], group = groups_switch, genes = CountTable_switch[,1])
DGEobject_switch = DGEobject_switch[rowSums(cpm(DGEobject_switch)>1) >= max(table(groups_switch)), , keep.lib.sizes=FALSE]
DGEobject_switch$genes = GFF[match(DGEobject_switch$genes$genes, GFF$gene_id), c(9,10,12)]
DGEobject_switch = calcNormFactors(DGEobject_switch)
DGEobject_switch = estimateDisp(DGEobject_switch, design = designMatrix_switch, robust = T)

contrasts_switch = makeContrasts(eGFP_induction = EGFP_Ind - EGFP_control,
                                 hop_induction = hop_Ind - hop_control,
                                 PGRP_induction = PGRP_Ind - PGRP_control,
                                 levels = designMatrix_switch)

fit_switch = glmQLFit(DGEobject_switch, designMatrix_switch, robust = T)
QFresults_switch = list()
for (i in 1:dim(contrasts)[2]) {
  QFresults_switch[[colnames(contrasts_switch)[i]]] = glmQLFTest(fit_switch, contrast = contrasts_switch[,i])
}


## infection


CountTableInfection = read.delim("data/fig2/OrRInfection.count", skip = 1)
CountTableInfection = CountTableInfection[,c(1,7:dim(CountTableInfection)[2])]

SampleDesctiptorsInfection = read.xlsx("data/fig2/OrRInfectionLibraries.xlsx")
colnames(CountTableInfection)[2:21] = sapply(strsplit(colnames(CountTableInfection)[2:21], "_", fixed = T), function(x){paste0("Lib_", x[2])})
SampleDesctiptorsInfection$bam = factor(SampleDesctiptorsInfection$bam, levels = colnames(CountTableInfection)[2:21])
SampleDesctiptorsInfection = SampleDesctiptorsInfection[order(SampleDesctiptorsInfection$bam),]

groupsInfection = SampleDesctiptorsInfection$condition
designMatrixInfection = model.matrix(~ 0 + groupsInfection)
colnames(designMatrixInfection) = substr(colnames(designMatrixInfection), 16, 100)
colnames(designMatrixInfection) = gsub(" ", "_", colnames(designMatrixInfection), fixed = T)

DGEobjectInfection = DGEList(CountTableInfection[,2:dim(CountTableInfection)[2]], group = groupsInfection, genes = CountTableInfection[,1])
DGEobjectInfection = DGEobjectInfection[rowSums(cpm(DGEobjectInfection)>1) >= max(table(groupsInfection)), , keep.lib.sizes=FALSE]
DGEobjectInfection$genes = GFF[match(DGEobjectInfection$genes$genes, GFF$gene_id), c(9,10,12)]
DGEobjectInfection = calcNormFactors(DGEobjectInfection)
DGEobjectInfection = estimateDisp(DGEobjectInfection, design = designMatrixInfection, robust = T)



## contrasts and fits
contrastsInfection = makeContrasts(de_3h = septic_3h - control,
                                   de_6h = septic_6h - control,
                                   de_18h = septic_18h - control,
                                   levels = designMatrixInfection)

fitInfection = glmQLFit(DGEobjectInfection, designMatrixInfection, robust = T)
QFresultsInfection = list()
for (i in 1:dim(contrastsInfection)[2]) {
  QFresultsInfection[[colnames(contrastsInfection)[i]]] = glmQLFTest(fitInfection, contrast = contrastsInfection[,i])
}


## mosaic vs switch

Mosaic.de = data.frame(
  gene_id = QFresults$H3.K27R$genes$gene_id,
  H3K27R.de = as.numeric(decideTests.DGELRT(QFresults$H3.K27R, lfc = 1)),
  Suz12.de = as.numeric(decideTests.DGELRT(QFresults$Suz12, lfc = 1)),
  Ez.de = as.numeric(decideTests.DGELRT(QFresults$Ez, lfc = 1))
)

InductionSwitch.de = data.frame(
  gene_id = QFresults_switch$eGFP_induction$genes$gene_id,
  eGFP.de = as.numeric(decideTests.DGELRT(QFresults_switch$eGFP_induction, lfc = 1)),
  hop.de = as.numeric(decideTests.DGELRT(QFresults_switch$hop_induction, lfc = 1)),
  PGRP.de =  as.numeric(decideTests.DGELRT(QFresults_switch$PGRP_induction, lfc = 1))
)



Full.de = Mosaic.de %>%
  inner_join(InductionSwitch.de, by = "gene_id") %>%
  mutate(mosaic.any.up = rowAny(as.matrix(.[,2:4] == 1))) %>%
  left_join(resting_H3K27me3 %>%
              dplyr::select(c(Geneid, class)),
            by = c("gene_id" = "Geneid"))

### FigMosaicVennVsSwitchGal4PcMonly.pdf

fit = euler(Full.de %>%
              filter(class == "Pc-I") %>%
              mutate(eGFP.de = eGFP.de == 1,
                     hop.de = hop.de == 1,
                     PGRP.de = PGRP.de == 1) %>%
              dplyr::select(eGFP.de, hop.de, PGRP.de, mosaic.any.up))
plot(fit,
     quantities = list(type = c("counts"), cex = 1.5),
     fills = list(fill = ggthemes_data$tableau$`color-palettes`$regular$`Superfishel Stone`$value[1:4]),
     legend = list(side = "right", vgap = 2, cex = 1.2), main = "Pc-M")

### FigMosaicVennVsSwitchGal4nonPconly.pdf

fit = euler(Full.de %>%
              filter(class == "non-Pc") %>%
              mutate(eGFP.de = eGFP.de == 1,
                     hop.de = hop.de == 1,
                     PGRP.de = PGRP.de == 1) %>%
              dplyr::select(eGFP.de, hop.de, PGRP.de, mosaic.any.up))

plot(fit,
     quantities = list(type = c("counts"), cex = 1.5),
     fills = list(fill = ggthemes_data$tableau$`color-palettes`$regular$`Superfishel Stone`$value[1:4]),
     legend = list(side = "right", vgap = 2, cex = 1.2), main = "non-Pc")

## vs Infections

Mosaic.de = data.frame(
  gene_id = QFresults$H3.K27R$genes$gene_id,
  H3K27R.de = as.numeric(decideTests.DGELRT(QFresults$H3.K27R, lfc = 1)),
  Suz12.de = as.numeric(decideTests.DGELRT(QFresults$Suz12, lfc = 1)),
  Ez.de = as.numeric(decideTests.DGELRT(QFresults$Ez, lfc = 1))
)

InductionInfection.de = data.frame(
  gene_id = QFresultsInfection$de_3h$genes$gene_id,
  si3h.de = as.numeric(decideTests.DGELRT(QFresultsInfection$de_3h, lfc = 1)),
  si6h.de = as.numeric(decideTests.DGELRT(QFresultsInfection$de_6h, lfc = 1))
)



Full.de = Mosaic.de %>%
  inner_join(InductionInfection.de, by = "gene_id") %>%
  mutate(mosaic.any.up = rowAny(as.matrix(.[,2:4] == 1))) %>%
  left_join(resting_H3K27me3 %>%
              dplyr::select(c(Geneid, class)),
            by = c("gene_id" = "Geneid"))

### FigMosaicVennVsInfectionGal4PcMonly.pdf

fit = euler(Full.de %>%
              filter(class == "Pc-I") %>%
              mutate(si3h.de = si3h.de == 1,
                     si6h.de = si6h.de == 1) %>%
              dplyr::select(si3h.de, si6h.de, mosaic.any.up))


plot(fit,
     quantities = list(type = c("counts"), cex = 1.5),
     fills = list(fill = ggthemes_data$tableau$`color-palettes`$regular$`Superfishel Stone`$value[1:3]),
     legend = list(side = "right", vgap = 2, cex = 1.2), main = "Pc-M")


### FigMosaicVennVsInfectionGal4nonPconly.pdf

fit = euler(Full.de %>%
              filter(class == "non-Pc") %>%
              mutate(si3h.de = si3h.de == 1,
                     si6h.de = si6h.de == 1) %>%
              dplyr::select(si3h.de, si6h.de, mosaic.any.up))

plot(fit,
     quantities = list(type = c("counts"), cex = 1.5),
     fills = list(fill = ggthemes_data$tableau$`color-palettes`$regular$`Superfishel Stone`$value[1:4]),
     legend = list(side = "right", vgap = 2, cex = 1.2), main = "non-Pc")

## Enrichment


Mosaic.de = data.frame(
  gene_id = QFresults$H3.K27R$genes$gene_id,
  H3K27R.de = as.numeric(decideTests.DGELRT(QFresults$H3.K27R, lfc = 1)),
  Suz12.de = as.numeric(decideTests.DGELRT(QFresults$Suz12, lfc = 1)),
  Ez.de = as.numeric(decideTests.DGELRT(QFresults$Ez, lfc = 1)))  %>%
  mutate(mosaic.any.up = case_when(rowAny(as.matrix(.[,2:4] == 1)) ~ "Mosaic up",
                                   T ~ "Mosaic not up")) %>%
  left_join(resting_H3K27me3 %>%
              dplyr::select(c(Geneid, class)),
            by = c("gene_id" = "Geneid"))


hypMod = function(qf){
  df = data.frame(gene_id = qf$genes$gene_id,
                  de = as.numeric(decideTests.DGELRT(qf, lfc = log2(1.5))))
  df %>%
    inner_join(Mosaic.de %>%
                 dplyr::select(c(gene_id, mosaic.any.up)) ,
               by = "gene_id") %>%
    group_by(de, mosaic.any.up) %>%
    dplyr::summarise(q = n()) %>%
    group_by(de) %>%
    dplyr::mutate(k = sum(q)) %>%
    group_by(mosaic.any.up) %>%
    dplyr::mutate(m = sum(q)) %>%
    dplyr::mutate(n = dim(DGEobject$genes)[1] - m) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(p_val = phyper(q-1, m, n, k, lower.tail = F, log.p = T)) %>%
    dplyr::mutate(effect_size = (q/k)/(m/(m+n))) %>%
    dplyr::mutate(p = exp(p_val)) %>%
    arrange(p)
}


hyped = sapply(c(QFresultsInfection, QFresults_switch), simplify = F, hypMod, USE.NAMES = T)

y_coding = c(de_18h = "18h p.i.", de_6h = "6h p.i.", de_3h = "3h p.i.", 
             eGFP_induction = "eGFP", PGRP_induction = "PGRP", hop_induction = "hop")

enricht.plot.df = rbindlist(hyped, idcol = "contrast") %>%
  filter(mosaic.any.up == "Mosaic up", de == 1) %>%
  dplyr::mutate(p.adj = -log10(p.adjust(p))) %>%
  dplyr::mutate(y = recode(contrast, !!!y_coding)) %>%
  dplyr::mutate(y = factor(y, levels = y_coding)) %>%
  dplyr::mutate(facet = case_when(contrast %in% c("de_3h", "de_6h", "de_18h") ~ "Septic Injury",
                           T ~ "Hormone induction"))

### FigMosaicEnrichmentOdInducible.pdf

ggplot(enricht.plot.df, aes(y, k)) +
  geom_col(aes(y = 1), color = "black", fill = "#AAAAAA") + 
  geom_col(aes(y = q/k, fill = p.adj), color = "black") + 
  coord_flip() +
  theme_bw() + ylab("Fraction of Mosaic Induced genes among up-regulated genes") +
  scale_fill_viridis(option = "viridis", limits = c(0, max(enricht.plot.df$p.adj)), name = "-log10\nadjusted\np-value") + 
  xlab("Contrast") +
  geom_text(aes(y = .85, label = paste0("n = ", k)))  +
  geom_hline(aes(yintercept = m/(m+n)), linetype = 5) +
  facet_wrap(facets = vars(facet), ncol = 1, scales = "free_y", strip.position = "right")

enricht.plot.df %>% kable() %>% 
  kableExtra::kable_styling(bootstrap_options = c("striped", "condensed"))
