library(ggplot2)
library(gridExtra)
library(rcartocolor)
library(GenomicRanges)
library(data.table)
library(circlize)
library(ggnewscale)
library(gridBase)
library(ComplexHeatmap)
library(ggridges)
library(ggrepel)
library(edgeR)
library(openxlsx)
library(rtracklayer)
library(DESeq2)
library(tidyverse)
library(Rfast)
library(ggthemes)


load("data/fig1/resting_H3K27me3.Rdata")

# data import -------------------------------------------------------------


GFF = readGFF("/Users/streeck/Genomes/DmelBDGP6.91/Drosophila_melanogaster.BDGP6.91.gtf", filter = list(type = "gene"))

CountTable = read.delim("data/fig2/OrRInfection.count", skip = 1)
CountTable = CountTable[,c(1,7:dim(CountTable)[2])]

SampleDesctiptors = read.xlsx("data/fig2/OrRInfectionLibraries.xlsx")
colnames(CountTable)[2:21] = sapply(strsplit(colnames(CountTable)[2:21], "_", fixed = T), function(x){paste0("Lib_", x[2])})
SampleDesctiptors$bam = factor(SampleDesctiptors$bam, levels = colnames(CountTable)[2:21])
SampleDesctiptors = SampleDesctiptors[order(SampleDesctiptors$bam),]

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
rld = rlog(dds, blind = FALSE)

## select 1000 genes with the highes variance
rv <- rowVars(assay(rld))
row_select <- order(rv, decreasing = TRUE)[seq_len(3000)]

## apply PCA analysis
pca <- prcomp(t(assay(rld[row_select, ])))
pcaPlot = as.data.frame(pca$x)
pcaPlot = cbind(groups, pcaPlot)
pcaPlot$ID = substr(row.names(pcaPlot), 1, nchar(row.names(pcaPlot))-4)
var.explained = round(100*pca$sdev^2/sum(pca$sdev^2), 2)


## contrasts and fits
contrasts = makeContrasts(de_3h = septic_3h - control,
                          de_6h = septic_6h - control,
                          de_18h = septic_18h - control,
                          levels = designMatrix)

fit = glmQLFit(DGEobject, designMatrix, robust = T)
QFresults = list()
for (i in 1:dim(contrasts)[2]) {
  QFresults[[colnames(contrasts)[i]]] = glmQLFTest(fit, contrast = contrasts[,i])
}



# Volcanos ----------------------------------------------------------------

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

### FigInfection3hVolcano.pdf
qf.volcano(QFresults$de_3h, "3h septic")

### FigInfection6hVolcano.pdf
qf.volcano(QFresults$de_6h, "6h septic")

### FigInfection18hVolcano.pdf
qf.volcano(QFresults$de_18h, "18h septic")



# Heatmap -----------------------------------------------------------------

cond = c("control", "septic_3h", "septic_6h", "septic_18h")

hyper_gs = DGEobject$genes %>%
  dplyr::select(gene_id) %>%
  left_join(resting_H3K27me3 %>%
              dplyr::select(Geneid, class), 
            by = c(gene_id = "Geneid")) %>%
  group_by(class) %>%
  dplyr::summarize(m = n()) %>%
  dplyr::mutate(n = dim(DGEobject$genes)[1] - m)

rld_rows = DGEobject$genes %>%
  dplyr::select(gene_id) %>%
  left_join(resting_H3K27me3 %>%
              dplyr::select(Geneid, class), 
            by = c(gene_id = "Geneid")) 

rld_assay = SummarizedExperiment::assay(rld) %>%
  as.data.frame %>%
  dplyr::select(SampleDesctiptors$bam[SampleDesctiptors$condition %in% cond])


n_top = 300
nk = 5

metric = as.matrix(
  cbind(
    -log10(QFresults$de_3h$table$PValue),
    -log10(QFresults$de_6h$table$PValue)
  )
)
metric = rowmeans(metric)



row_select <- order(metric, decreasing = TRUE)[1:n_top]


rld_selected_counts = rld_assay[row_select,]
scaled_rld = rld_selected_counts %>% 
  t %>%
  scale %>%
  t

rld_cluster = DGEobject$genes %>%
  .[row_select,]

rld_cluster$cl = kmeans(scaled_rld, centers = nk, algorithm = 'MacQueen', iter.max = 50, nstart = 10)$cluster

heatmap_cols = SampleDesctiptors %>%
  filter(condition %in% cond) %>%
  dplyr::select(condition, bam) %>%
  mutate(condition = factor(condition, levels = cond)) %>%
  arrange(condition)

scaled_rld = scaled_rld %>%
  .[,as.character(heatmap_cols$bam)]


hm_rows = rld_rows[row_select,] %>%
  dplyr::mutate(Pc.M = case_when(class == "Pc-I" ~ "Pc-M",
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

### FigInfectionFullHeatmap.pdf

ht_list = Heatmap(matrix = scaled_rld,
                  column_order = 1 : dim(scaled_rld)[2],
                  column_split = heatmap_cols$condition,
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


t = draw(ht_list, row_split = rld_cluster$cl,
         show_row_dend = F)



### FigInfectionFullHeatmapTable.png

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

### FigInfectionFullHeatmapBoxplot.pdf

as.data.frame(scaled_rld) %>%
  bind_cols(rld_cluster %>% dplyr::select(gene_id, cl)) %>%
  pivot_longer(!gene_id:cl, names_to = "sample_id") %>%
  left_join(heatmap_cols, by = c("sample_id" = "bam")) %>%
  group_by(condition, gene_id, cl) %>%
  dplyr::summarise(mean_z = mean(value)) %>%
  dplyr::mutate(cl = factor(cl, levels = rev(as.numeric(t@ht_list[[1]]@row_title)))) %>%
  ggplot(aes(as.factor(cl), mean_z, fill = condition)) + 
  geom_boxplot(outlier.alpha = 0) +
  scale_fill_tableau(name = "Condition") +
  theme_bw() + ylab("mean z-score") + xlab("Cluster") +
  coord_flip()



# Enrichment --------------------------------------------------------------



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
    dplyr::mutate(p = phyper(q-1, m, n, k, lower.tail = F, log.p = F)) %>%
    dplyr::mutate(effect_size = (q/k)/(m/(m+n))) %>%
    filter(de != 0) %>%
    filter(class == "Pc-I") %>%
    arrange(p)
}


hyped = sapply(QFresults, simplify = F, hypMod, USE.NAMES = T)


temp.df = data.frame(-1, "Pc-I", 0, 4, 898, 6688, 1, 1)
colnames(temp.df) = colnames(hyped$de_18h)

hyped$de_18h = rbind(hyped$de_18h, 
                     temp.df)


y_coding = c(de_3h = "3h septic", de_6h = "6h septic", 
             de_18h = "18h septic")


# FigInfectionEnrichment.pdf

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



# PCA ---------------------------------------------------------------------


library(RColorBrewer)
library(plot3D)

pcaPlot$groups = factor(pcaPlot$groups)

t10 = ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value

pcaPlot = pcaPlot %>%
  filter(groups != "sterile") %>%
  mutate(groups = droplevels(groups))


## FigRNAiPCA.pdf

scatter3D(pcaPlot$PC1, pcaPlot$PC2, pcaPlot$PC3, 
          bgvar = as.numeric(as.factor(pcaPlot$groups)), 
          col = "black", bg = t10[as.numeric(pcaPlot$groups)],
          theta = 150, phi = 10, xlab = paste0("PC1: ", var.explained[1], "% var"),
          ylab = paste0("PC2: ", var.explained[2], "% var"), zlab = paste0("PC3: ", var.explained[3], "% var"), 
          bty ="b2", lwd = 1.5, cex = 2.5, colkey = F, alpha = .8, pch = 21)
legend("right",levels(pcaPlot$groups), 
       pt.bg = t10, pt.lwd = 1.2, cex = 1, pt.cex = 2.5, inset = c(0.07,0.07), pch = 21)





