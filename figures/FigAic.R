library(ggplot2)
library(data.table)
library(tidyverse)
library(ggthemes)


chip_enrichment = function(chip, input){
  chip = chip/sum(chip)
  input = input/sum(input)
  (chip)/(input+chip)
}

log_ratio = function(chip, input){
  chip = (chip)/sum(chip)
  input = (input)/sum(input)
  return(log2(chip) - log2(input))
}


# GeneLevelH3K27me3 -------------------------------------------------------



resting_H3K27me3 = read.delim("data/fig1/H3K27me3ChIP3772.count", sep = "\t", skip = 1)
colnames(resting_H3K27me3)[7:10] = c("H3K27me3_A", "input_A", "H3K27me3_B", "input_B")

resting_H3K27me3 = resting_H3K27me3 %>% 
  dplyr::select(Geneid, H3K27me3_A:input_B) %>%
  dplyr::mutate(lr_A = log_ratio(H3K27me3_A, input_A),
         lr_B = log_ratio(H3K27me3_B, input_B)) %>%
  dplyr::mutate(h3k27me3_mean_lr = (lr_A + lr_B)/2)

plot_table = resting_H3K27me3 %>%
  dplyr::select(Geneid, h3k27me3_mean_lr)

lnl = c()

source("EM/EMcalc.R")
S = as.matrix(resting_H3K27me3[,c("H3K27me3_A","H3K27me3_B")])
N = as.matrix(resting_H3K27me3[,c("H3K27me3_A","H3K27me3_B")] + 
                resting_H3K27me3[,c("input_A","input_B")])

for (i in 1:10) {
  GeneLevelH3K27me3Cluster = BinomEMwrapperParallel(N, S, k=i, ncores = 10)
  lnl = c(lnl, GeneLevelH3K27me3Cluster$ModelLnLikelyhood)
  plot_table = cbind(plot_table, GeneLevelH3K27me3Cluster$Group)
  colnames(plot_table)[dim(plot_table)[2]] = paste0("p", i)
}

plot_table_melt = plot_table %>%
  pivot_longer(!Geneid:h3k27me3_mean_lr) %>%
  mutate(value = factor(value,
                           levels = as.character(1:10)),
         name = factor(paste0(substr(name, 2, 4), " classes"),
                        levels = paste0(1:10, " classes"))) %>%
  mutate(name = recode(name, "1 classes" = "1 class"))
  
lnl = data.frame(comp = 1:10, lnl = lnl, aic = -2*lnl + 2*1:10)

### FigAicGeneStateHist.pdf

ggplot(plot_table_melt, aes(h3k27me3_mean_lr, fill = value)) + geom_histogram(bins = 100) +
  facet_wrap(vars(name), ncol = 1, strip.position = "right") + theme_bw() + 
  scale_fill_tableau(palette = "Superfishel Stone", name = "component") + xlab("H3K27me3 log2 ratio") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())

### FigAicGeneStateAic.pdf

ggplot(lnl, aes(as.factor(comp), aic)) + geom_point(size = 2) + theme_bw() +
  geom_line(aes(group = 1 )) + xlab("Number of components") + ylab("Akaike Information Criterion") 



# ChromStateModel ---------------------------------------------------------

library(GenomicRanges)


## import 

chip_files = read.xlsx("data/fig1/SampleList.xlsx")
chip_file_selection = chip_files %>% 
  dplyr::select(Index:Input.ref, Replicate) %>%
  filter(Antibody != "Input") %>% 
  left_join(chip_files %>% 
              dplyr::select(Index, Bam.File), 
            by = c("Input.ref" = "Index"),
            suffix = c(".chip", ".input"))


genome = read.delim("/Users/streeck/Genomes/DmelBDGP6.91/chrNameLength.txt", header = F, stringsAsFactors = F)
genome = genome[1:7,]
gr = GRanges(genome[,1], IRanges(1, as.integer(genome[,2])))

gr_small = gr[1:2]

ncores = 10

path = "/Volumes/Promise_Leia/Thesis/ChIPReanalysis/ChIP_BAMs/"

ChIPBamCounts = parallel::mcmapply(bamsignals::bamProfile, bampath=paste(path, chip_file_selection$Bam.File.chip, sep = ""), MoreArgs =
                                     list(gr=gr_small, binsize = 200, paired.end="midpoint", tlenFilter = c(70, 500), 
                                          filteredFlag = 1024, verbose = F),
                                   mc.cores = ncores, SIMPLIFY = F)

InputBamCounts = parallel::mcmapply(bamsignals::bamProfile, bampath=paste(path, chip_file_selection$Bam.File.input, sep = ""), MoreArgs =
                                      list(gr=gr_small,binsize = 200, paired.end="midpoint", tlenFilter = c(70, 500), 
                                           filteredFlag = 1024, verbose = F),
                                    mc.cores = ncores, SIMPLIFY = F)

### number of ChIP sampels
ndim = length(chip_file_selection$Bam.File.chip)

### convert bamsignals output to a list
for(i in 1:ndim){ChIPBamCounts[[i]] = unlist(as.list(ChIPBamCounts[[i]]))}
for(i in 1:ndim){InputBamCounts[[i]] = unlist(as.list(InputBamCounts[[i]]))}

### convert bamsignal lists to matrices
ChIPBamMatrix = matrix(nrow = length(ChIPBamCounts[[1]]), ncol = ndim)
for(i in 1:ndim){ChIPBamMatrix[,i] = ChIPBamCounts[[i]]}
colnames(ChIPBamMatrix) = chip_file_selection$Bam.File.chip
InputBamMatrix = matrix(nrow = length(InputBamCounts[[1]]), ncol = ndim)
for(i in 1:ndim){InputBamMatrix[,i] = InputBamCounts[[i]]}
colnames(InputBamMatrix) = chip_file_selection$Bam.File.input

N = ChIPBamMatrix + InputBamMatrix

NanEliminator = rowsums(N == 0) == 0

## fit models


enrichments = log_ratio(ChIPBamMatrix, InputBamMatrix)
enrichments = enrichments[NanEliminator,]


mean_enrichments = matrix(nrow = sum(NanEliminator),
                          ncol = length(unique(chip_file_selection$Antibody)))
colnames(mean_enrichments) = unique(chip_file_selection$Antibody)

for(t in colnames(mean_enrichments)){
  bams = chip_file_selection$Bam.File.chip[chip_file_selection$Antibody == t]
  mean_enrichments[,t] = rowmeans(enrichments[,bams])
}

data_tables = list()
bin_lnl = c()
comp_size = list()
fits = list()

for (i in 4:14) {
  fit = BinomEMwrapperParallel(N[NanEliminator,], ChIPBamMatrix[NanEliminator,], k = i, ncores = ncores, maxiter = 200)
  fits[[i-3]] = fit
  bin_lnl = c(bin_lnl, fit$ModelLnLikelyhood)
  for (j in 1:i) {
    data_tables[[length(data_tables) + 1]] =
      data.frame(target = colnames(mean_enrichments),
                 mean_signal = apply(mean_enrichments[fit$Group == j,], 2, function(x){mean(x[!is.infinite(x) & !is.nan(x)])}),
                 mod_size = rep(paste0("model ", i), dim(mean_enrichments)[2]),
                 comp = rep(paste0(i, "_", j), dim(mean_enrichments)[2]))
  }
  comp_size[[i-3]] = data.frame(as.data.frame(table(fit$Group)), mod_size = paste0("model ", i))
}


df_bin_lnl = data_frame(comp = 4:14, lnl = bin_lnl, aic = -2*lnl + 2*4:14)

## heatmaps

hclust_levels = function(df, row_s, col_s, val_s){
  df_wide = df %>% dplyr::select(row_s, col_s, val_s) %>%
    pivot_wider(id_cols = row_s, names_from = col_s, values_from = val_s)
  df_wide = as.data.frame(df_wide)
  h_cluster = hclust(dist(df_wide[2:dim(df_wide)[2]]))
  return(as.character(df_wide[h_cluster$order,1]))
}

adjust_ggplot_facet_size = function(ggplot_object){
  # convert ggplot object to grob object
  gp <- ggplotGrob(p)
  
  # get gtable columns corresponding to the facets 
  facet.columns <- unique(gp$layout$l[grepl("panel", gp$layout$name)])
  
  # get gtable rows corresponding to the facets 
  facet.rows <- unique(gp$layout$b[grepl("panel", gp$layout$name)])
  
  # get the number of unique x-axis values per facet
  x.var <- sapply(ggplot_build(p)$layout$panel_scales_x,
                  function(l) length(l$range$range))
  
  # get the number of unique y-axis values per facet
  y.var <- sapply(ggplot_build(p)$layout$panel_scales_y,
                  function(l) length(l$range$range))
  
  # change the relative widths of the facet columns based on
  # how many unique x-axis values are in each facet
  gp$widths[facet.columns] <- gp$widths[facet.columns] * x.var
  
  # change the relative widths of the facet columns based on
  # how many unique x-axis values are in each facet
  gp$heights[facet.rows] <- gp$heights[facet.rows] * y.var
  
  return(gp)
}

states_table = rbindlist(data_tables)

states_table$target = factor(states_table$target, 
                             levels = hclust_levels(states_table, "target", "comp", "mean_signal"))

states_table$comp = factor(states_table$comp,
                           levels = hclust_levels(states_table, "comp", "target", "mean_signal"))

states_table$mod_size = factor(states_table$mod_size,
                               levels = paste0("model ", 4:14))

### FigAicChromStateHeatmap.pdf

p = ggplot(states_table, aes(target, comp, fill = mean_signal)) + geom_tile() + 
  theme_classic() + scale_fill_viridis_c(option = "magma", name = "mean\nlog2 ratio") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        axis.ticks = element_blank(), axis.line = element_blank(),
        axis.text = element_text(size = 10, face = "bold"),
        axis.text.y = element_blank()) +
  facet_wrap(vars(mod_size), ncol = 1, strip.position = "right", scales = "free_y") 

gp = adjust_ggplot_facet_size(p)
grid::grid.draw(gp)

### FigAicChromStateAIC.pdf

ggplot(df_bin_lnl, aes(as.factor(comp), aic)) + geom_point(size = 2) + theme_bw() +
  geom_line(aes(group = 1 )) + xlab("Number of components") + ylab("Akaike Information Criterion") 


comp_table = rbindlist(comp_size)

comp_table$Var1 = factor(comp_table$Var1, levels = 1:14)

comp_table$mod_size = factor(comp_table$mod_size,
                             levels = paste0("model ", 4:14))

### FigAicChromStateFractionOfBins.pdf

p = ggplot(comp_table, aes(Var1, Freq)) + geom_bar(stat="identity") + coord_flip() + 
  theme_bw() + ylab("Number of bins") + xlab("Component") +
  facet_wrap(vars(mod_size), ncol = 1, strip.position = "right", scales = "free_y") 

gp = adjust_ggplot_facet_size(p)
grid::grid.draw(gp)

