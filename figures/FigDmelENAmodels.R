library(ggplot2)
library(plotly)
library(tidyverse)
library(openxlsx)
library(ComplexHeatmap)
library(ggthemes)
library(kableExtra)
library(circlize)
library(RColorBrewer)


# data import -------------------------------------------------------------



Paired = read.delim("data/DmelH3K27me3ENA/paired_end_samples.count", skip = 1)
Paired = Paired[,c(1,7:(dim(Paired)[2]-2))]
colnames(Paired)[2:dim(Paired)[2]] = sapply(strsplit(colnames(Paired)[2:dim(Paired)[2]], split = ".", fixed = T), function(x){x[length(x)-1]})

Single = read.delim("data/DmelH3K27me3ENA/single_end_samples.count", skip = 1)
Single = Single[,c(1,7:(dim(Single)[2]))]
colnames(Single)[2:dim(Single)[2]] = sapply(strsplit(colnames(Single)[2:dim(Single)[2]], split = ".", fixed = T), function(x){x[length(x)-1]})

AllCounts = cbind(Paired, Single[,2:dim(Single)[2]])

Sampels = read.xlsx("data/DmelH3K27me3ENA/MappedCountedSamples.xlsx")


Sample.selector = Sampels %>%
  filter(ip.type == "H3K27me3") %>%
  select(ID, Study.Accession, cell.type, Matched.input, File.accession) %>%
  left_join(Sampels %>%
              select(ID, File.accession),
            by = c("Matched.input" = "ID"),
            suffix = c(".chip", ".input")) %>%
  select(Study.Accession, cell.type, File.accession.chip, File.accession.input) %>%
  mutate(chip.read.count = colSums(AllCounts[File.accession.chip]),
         input.read.count = colSums(AllCounts[File.accession.input]))


Log2RatioTable = cbind(
  gene_id = Paired$Geneid,
  AllCounts[Sample.selector$File.accession.chip]) %>%
  pivot_longer(!gene_id, names_to = "File.accession.chip", values_to = "chip.counts") %>%
  left_join(Sample.selector,
            by = "File.accession.chip") %>%
  left_join(cbind(gene_id = Paired$Geneid,
                  AllCounts[Sample.selector$File.accession.input]) %>%
              pivot_longer(!gene_id, names_to = "File.accession.input", values_to = "input.counts"),
            by = c("File.accession.input", "gene_id")) %>%
  filter(chip.counts > 0 & input.counts > 0) %>%
  mutate(log2.ratio = log2(chip.counts) - log2(input.counts))

PlotLogRatios = Log2RatioTable %>%
  mutate(col = case_when(substr(cell.type, 1, 2) == "PS" ~ "parasegments",
                         cell.type %in% c("t1nb", "neuron", "OK107", "Tdc2", "R57C10") ~ "neuronal cells",
                         cell.type == "plasmat" ~ "plasmatocytes",
                         T ~ cell.type),
         study = case_when(Study.Accession == "Rob1me3" ~ "Our study",
                           T ~ paste0("Study #", as.numeric(factor(Study.Accession))))) %>%
  mutate(facet_labs = paste0(study, " - ", File.accession.chip))


ProjectDescriptors = read.xlsx("data/DmelH3K27me3ENA/StudyOverview.xlsx")
ProjectDescriptors = ProjectDescriptors[order(ProjectDescriptors$Study.Accession),]
ProjectDescriptors$Paper.Author = gsub("\r", "", ProjectDescriptors$Paper.Author, fixed = T)


# Cor heatmap -------------------------------------------------------------


matrix.for.cor = Log2RatioTable %>%
  pivot_wider(gene_id, names_from = File.accession.chip, 
              values_from = log2.ratio, values_fill = 0) %>%
  select(!gene_id) %>%
  as.matrix()

box_cols = ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value[1:6]
names(box_cols) = unique(ProjectDescriptors$Cell.types)


lgd_list = list(
  Legend(labels = names(box_cols), title = "cell\ntype", legend_gp = gpar(fill = box_cols)))

ha = rowAnnotation(foo = anno_empty(border = FALSE, 
                                    width = max_text_width(unlist(strsplit(ProjectDescriptors$Paper.Author, "\n", fixed = T))) + unit(4, "mm")))


### FigDmelENAmodelsCorHeatmap.pdf

t = Heatmap(cor(matrix.for.cor, method = "pearson"),
            col = colorRamp2(seq(-1,1, length.out = 7), rev(brewer.pal(7, "RdBu"))),
            row_split = Sample.selector$Study.Accession,
            column_split = Sample.selector$Study.Accession,
            border = T, show_row_names = F, show_column_names = T,
            row_title = NULL, column_title = NULL, show_row_dend = F,
            right_annotation = ha,
            name = "pearson\ncorrelation")
t2 = draw(t, annotation_legend_list = lgd_list)
for(i in 1:11) {
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, width = unit(2, "mm"), gp = gpar(
      fill = box_cols[ProjectDescriptors$Cell.types[unlist(t2@ht_list[[1]]@row_dend_slice)[i]]], col = NA), just = "left")
    grid.text(ProjectDescriptors$Paper.Author[unlist(t2@ht_list[[1]]@row_dend_slice)[i]], x = unit(4, "mm"), just = "left")
  })
}


# Fit models --------------------------------------------------------------

get_cor_spear = function(ids){
  if(length(ids) > 1){
    x = as.data.frame(matrix.for.cor)[ids]
    cor(x, method = "spearman")[1,2]
  }else{
    NA
  }
  
}

get_cor_pears = function(ids){
  if(length(ids) > 1){
    x = as.data.frame(matrix.for.cor)[ids]
    cor(x)[1,2]
  }else{
    NA
  }
  
}

Sample.selector.exclusion = Sample.selector %>%
  group_by(Study.Accession, cell.type) %>%
  mutate(replicates = n()) %>%
  mutate(rep_cor_spear = get_cor_spear(File.accession.chip),
         rep_cor_pears = get_cor_pears(File.accession.chip)) %>%
  ungroup() %>%
  mutate(study.celltype = paste0(Study.Accession, ".", cell.type))


Sample.selector.exclusion  = Sample.selector.exclusion %>%
  filter(replicates > 1) %>%
  filter(rep_cor_pears > 0.6)




source("EM/EMcalc.R")

em.model.list = list()
for (i in Sample.selector.exclusion$study.celltype) {
  S = as.matrix(as.data.frame(AllCounts)[
    Sample.selector.exclusion$File.accession.chip[
      Sample.selector.exclusion$study.celltype == i]])
  R = as.matrix(as.data.frame(AllCounts)[
    Sample.selector.exclusion$File.accession.input[
      Sample.selector.exclusion$study.celltype == i]])
  em.model.list[[i]] = BinomEMwrapperParallel(S+R, S, 3, ncores = 10)
}



EM.gene.df = data.table::rbindlist(
  lapply(1:9, 
         function(x){
           data.frame(gene_id = AllCounts$Geneid,
                      class = c("non-Pc", "Pc-M", "Pc-H")[as.numeric(em.model.list[[x]]$Group)],
                      study.celltype = names(em.model.list)[x],
                      total.chip.count = rowsums(as.matrix(as.data.frame(AllCounts)[
                        Sample.selector.exclusion$File.accession.chip[
                          Sample.selector.exclusion$study.celltype == names(em.model.list)[x]]])),
                      total.input.count = rowsums(as.matrix(as.data.frame(AllCounts)[
                        Sample.selector.exclusion$File.accession.input[
                          Sample.selector.exclusion$study.celltype == names(em.model.list)[x]]])))
         })
)


EM.gene.df = EM.gene.df %>%
  group_by(study.celltype) %>%
  filter(total.chip.count > 0 & total.input.count > 0) %>%
  mutate(total.log2.ratio = log2(total.chip.count/sum(total.chip.count)) - log2(total.input.count/sum(total.input.count)))



# Corr Plots --------------------------------------------------------------
library(gridExtra)
library(grid)


select.df = Sample.selector.exclusion %>%
  filter(replicates == 2 & rep_cor_spear > 0.6) %>%
  select(Study.Accession, cell.type, study.celltype) %>%
  distinct()

plt1 = list()

for (i in select.df$study.celltype) {
  temp.df = Log2RatioTable %>%
    mutate(study.celltype = paste0(Study.Accession, ".", cell.type)) %>%
    filter(study.celltype == i) %>%
    left_join(data.frame(
      gene_id = AllCounts$Geneid,
      gene_state = c("non-Pc", "Pc-M", "Pc-H")[as.numeric(em.model.list[[i]]$Group)]),
      by = "gene_id") %>%
    pivot_wider(id_cols = c(gene_id, gene_state, Study.Accession, cell.type),
                values_from = log2.ratio,
                names_from = File.accession.chip)
  
  colnames(temp.df)[5:6] = c("RepA", "RepB")
  
  plt1[[length(plt1) + 1]] = ggplot(temp.df, aes(RepA, fill = gene_state)) +
    geom_histogram(aes(y = 100*..count../sum(..count..)), bins = 100) + theme_bw() +
    scale_fill_tableau(palette = "Superfishel Stone") +
    ylab("% of Genes") + scale_y_continuous(breaks = c(0,2,4,6,8)) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(5,5,1,8.5), "points")) + guides(fill = "none") +
    ggtitle(paste0(temp.df$cell.type[1], ",  ", 
                   gsub("\n", " ", fixed = T,
                        ProjectDescriptors$Paper.Author[ProjectDescriptors$Study.Accession == temp.df$Study.Accession[1]])))
  
  plt1[[length(plt1) + 1]] = grid.rect(gp=gpar(col="white"))
  
  plt1[[length(plt1) + 1]] = ggplot(temp.df, aes(RepA, RepB, color = gene_state)) +
    geom_point(alpha = .4, shape = 16) + theme_bw() +
    scale_y_continuous(breaks = c(-8,-6,-4,-2,0,2,4)) +
    scale_color_tableau(palette = "Superfishel Stone") +
    theme(legend.position = c(.08, .75), 
          legend.background = element_rect(color = "black", size = .3),
          plot.margin = unit(c(1,5,5,5), "points")) +
    guides(color = "none")
  
  plt1[[length(plt1) + 1]] = ggplot(temp.df, aes(RepB, fill = gene_state)) +
    geom_histogram(aes(y = 100*..count../sum(..count..)), bins = 100) + theme_bw() +
    scale_fill_tableau(palette = "Superfishel Stone") +
    ylab("% of Genes") +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(1,5,5,1), "points")) + guides(fill = "none") +
    coord_flip()

}


plt1[[33]] = cowplot::get_legend(ggplot(temp.df, aes(RepA, RepB, color = gene_state)) +
  geom_point(shape = 16) + 
  scale_color_tableau(palette = "Superfishel Stone"))


lay = rbind(c(1,2,5,6,9,10),
            c(3,4,7,8,11,12),
            c(13,14,17,18,21,22),
            c(15,16,19,20,23,24),
            c(25,26,29,30,NA,NA),
            c(27,28,31,32,33,NA))

### FigDmelENAmodelsFullCorrPlot.pdf

grid.arrange(grobs = plt1, layout_matrix = lay, heights = rep(c(.3, .7), 3), widths = rep(c(.7, .2), 3))
