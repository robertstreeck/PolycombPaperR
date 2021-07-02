library(ggplot2)
library(plotly)
library(openxlsx)
library(Rfast)
library(ggplot2)
library(tidyverse)
library(rtracklayer)
library(ComplexHeatmap)
library(plyr)
library(circlize)
library(kableExtra)
library(ggthemes)
library(GGally)
library(edgeR)

source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

lower_ggpairs = function(data, mapping){
  ggplot(data, mapping) + 
    geom_point(alpha = .2, shape = 16) +
    geom_density_2d(bins = 8, binwidth = 1, color = "#90728f") +
    geom_smooth(method = "lm", color = "#9d983dbb", formula = y ~ x)
}


quantNorm = function(X, round = F){
  require(Rfast)
  X = as.matrix(X)
  ndim = dim(X)[2]
  OrderMat = matrix(nrow = dim(X)[1], ncol = ndim)
  for(i in 1:ndim){
    OrderMat[,i] = base::order(X[,i])
    X[,i] = X[OrderMat[,i],i]
  }
  if(round == F){
    RowMeans = Rfast::rowsums(X)/ndim
  }else{
    RowMeans = round(Rfast::rowsums(X)/ndim)
  }
  for(i in 1:ndim){
    X[OrderMat[,i],i] = RowMeans
  }
  X
}

GeneGFF = readGFF("/Users/streeck/Genomes/DmelBDGP6.91/Drosophila_melanogaster.BDGP6.91.gtf", filter = list(type = "gene"))


ChIPBamFiles = read.xlsx("data/fig2/ChIPSampleList.xlsx")
RestingK27me3 = read.delim("data/fig2/H3K27me3ChIP3772.count", skip = 1)[,c(1, 6:10)]
colnames(RestingK27me3)[3:6] = substr(colnames(RestingK27me3)[3:6], 14, 100)
colnames(RestingK27me3)[3:6] = paste(ChIPBamFiles$Antibody, ChIPBamFiles$Replicate, sep = "_")[match(colnames(RestingK27me3)[3:6], ChIPBamFiles$Bam.File)]
RestingK27me3$log.ratio.A = log2(RestingK27me3$H3K27me3_A) - log2(RestingK27me3$Input_A)
RestingK27me3$log.ratio.B = log2(RestingK27me3$H3K27me3_B) - log2(RestingK27me3$Input_B)

InfecedChIPBamFiles = read.xlsx("data/fig2/InfectionSampleList.xlsx")
InfectedK27me3 = read.delim("data/fig2/H3K27me3ChIP4045.count", skip = 1)[,c(1, 6:10)]
colnames(InfectedK27me3 )[3:6] = substr(colnames(InfectedK27me3 )[3:6], 14, 100)
colnames(InfectedK27me3 )[3:6] = paste(InfecedChIPBamFiles$Antibody, InfecedChIPBamFiles$Replicate, sep = "_")[match(colnames(InfectedK27me3 )[3:6], InfecedChIPBamFiles$Bam.File)]
InfectedK27me3$log.ratio.A = log2(InfectedK27me3$H3K27me3_A) - log2(InfectedK27me3$Input_A)
InfectedK27me3$log.ratio.B = log2(InfectedK27me3$H3K27me3_B) - log2(InfectedK27me3$Input_B)

TSSGFF = readGFF("/Users/streeck/Genomes/DmelBDGP6.91/Drosophila_melanogaster.BDGP6.91.gtf", filter = list(type = "transcript"))
TSSK27ac = read.delim("data/fig2/TSSK27acWithInfections.count", skip = 1)[,c(1,7:14)]
colnames(TSSK27ac) = c("TranscriptID", "RestingInputA", "RestingK27acA", "RestingInputB", "RestingK27acB", "InfInputA", "InfInputB", "InfK27acA", "InfK27acB")
TSSK27ac$GeneID = TSSGFF$gene_id[match(TSSK27ac$TranscriptID, TSSGFF$transcript_id)]
TSSK27acMaxTranscripts = ddply(TSSK27ac, .(GeneID), summarize, MaxTrasncriptID = TranscriptID[which.max((RestingK27acA + RestingK27acB + InfK27acA + InfK27acB)/(RestingInputA + RestingInputB + InfInputA + InfInputB))])
TSSK27ac = TSSK27ac[TSSK27ac$TranscriptID %in% TSSK27acMaxTranscripts$MaxTrasncriptID,]
TSSK27ac[,11:14] = log2(as.matrix(TSSK27ac[,c(3,5,8,9)])) - log2(as.matrix(TSSK27ac[,c(2,4,6,7)]))
colnames(TSSK27ac)[11:14] = paste0(rep(c("log.ratio.A", "log.ratio.B"), 2), 
                                   rep(c(".control", ".infection"), each = 2))

DiffChIP_data_prenorm = GeneGFF %>%
  filter(gene_biotype == "protein_coding") %>%
  filter(seqid %in% c("X", "2R", "2L", "3R", "3L")) %>%
  select(gene_id, gene_name) %>%
  left_join(RestingK27me3 %>%
              select(Geneid, log.ratio.A, log.ratio.B),
            by = c("gene_id" = "Geneid")) %>%
  left_join(InfectedK27me3 %>%
              select(Geneid, log.ratio.A, log.ratio.B),
            by = c("gene_id" = "Geneid"),
            suffix = c(".control", ".infection")) %>%
  left_join(TSSK27ac[,10:14],
            by = c("gene_id" = "GeneID"),
            suffix = c("H3K27me3", "H3K27ac")) %>%
  filter(is.finite(rowSums(.[3:10])))

### FigSubGeneLevelPreQuantileNormH3K27me3.pdf

DiffChIP_data_prenorm %>%
  .[,3:6] %>%
  ggpairs(lower = list(continuous = lower_ggpairs), title = "H3K27me3 before quantile normalization",
          xlab = "log2 ratio", ylab = "log2 ratio")

DiffChIP_data = DiffChIP_data_prenorm
DiffChIP_data[,3:6] = quantNorm(DiffChIP_data[,3:6])
DiffChIP_data[,7:10] = quantNorm(DiffChIP_data[,7:10])

### FigSubGeneLevelPostQuantileNormH3K27me3.pdf

DiffChIP_data %>%
  .[,3:6] %>%
  ggpairs(lower = list(continuous = lower_ggpairs), title = "H3K27me3 after quantile normalization",
          xlab = "log2 ratio", ylab = "log2 ratio")



load("data/fig1/resting_H3K27me3.Rdata")

PcM.Diff = DiffChIP_data %>%
  filter(gene_id %in% (resting_H3K27me3 %>%
                         filter(class == "Pc-I") %>%
                         filter(gene_biotype == "protein_coding") %>%
                         .$Geneid))



# Import Infection data ---------------------------------------------------




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


PcM.Diff = PcM.Diff %>%
  mutate(up3h = case_when((as.numeric(decideTests.DGELRT(QFresults$de_3h)) == 1)[match(PcM.Diff$gene_id, DGEobject$genes$gene_id)] ~ "up 3h",
                          T ~ "a"),
         up6h = case_when((as.numeric(decideTests.DGELRT(QFresults$de_6h)) == 1)[match(PcM.Diff$gene_id, DGEobject$genes$gene_id)] ~ "up 6h",
                          T ~ "a"),
         up18h = case_when((as.numeric(decideTests.DGELRT(QFresults$de_18h)) == 1)[match(PcM.Diff$gene_id, DGEobject$genes$gene_id)] ~ "up 18h",
                           T ~ "a"))




# Heatmaps ----------------------------------------------------------------


PcM.sum = PcM.Diff %>%
  select(up3h:up18h) %>%
  pivot_longer(up3h:up18h) %>%
  group_by(name, value) %>%
  dplyr::summarise(m = n()) %>%
  dplyr::mutate(n = dim(PcM.Diff)[1] - m) %>%
  filter(value != "a")

metric = rowVars(as.matrix(PcM.Diff[,3:6]))


### FigDiffChIPHeatmapH3K27me3Only.pdf


count.matrix = t(scale(t(PcM.Diff[,3:6])))[order(metric, decreasing = T),][1:200,]

lab.df = PcM.Diff[order(metric, decreasing = T),11:12][1:200,]

hm.clusters = kmeans(count.matrix, centers = 3, algorithm = 'MacQueen', iter.max = 50, nstart = 10)$cluster

hm_rows = data.frame(cl = hm.clusters) %>%
  left_join(cbind(lab.df, 
                  cl = hm.clusters) %>%
              pivot_longer(!cl, values_to = "reg") %>%
              group_by(cl, name, reg) %>%
              dplyr::summarise(q = n()) %>%
              group_by(cl, name) %>%
              dplyr::mutate(k = sum(q)) %>%
              filter(reg != "a")%>%
              ungroup() %>%
              left_join(PcM.sum, by = "name") %>%
              dplyr::mutate(p_val = phyper(q-1, m, n, k, lower.tail = F)) %>%
              mutate(p.adj = p.adjust(p_val)) %>%
              filter(p.adj <= 0.05) %>%
              mutate(p.adj = log10(p.adj)) %>%
              pivot_wider(cl, names_from = name, values_from = p.adj),
            by = "cl") %>%
  replace_na(list(up3h = 0, up6h = 0))

conditions = c("Resting H3K27me3", "Septic H3K27me3")

p_col_fun = colorRamp2(c(min(hm_rows[,2:3]), 0), c("#666666", "#FFFFFF"))


ht_list = Heatmap(matrix = count.matrix,
                  column_order = 1 : 4,
                  column_split = rep(conditions, each = 2),
                  show_column_names = FALSE,
                  col = colorRampPalette(rcartocolor::carto_pal(n = 7, name = "ArmyRose"))(10),
                  show_row_names = FALSE, name = "row\nz score")


ht_list = ht_list + Heatmap(lab.df[,1], 
                            name = "3h", 
                            col = c("white", "black"),
                            width = unit(.2, "cm"),
                            show_row_names = FALSE)

ht_list = ht_list + Heatmap(hm_rows[,2], 
                            name = "p3h",
                            col = p_col_fun,
                            width = unit(.2, "cm"),
                            show_column_names = FALSE)

ht_list = ht_list + Heatmap(lab.df[,2], 
                            name = "6h", 
                            col = c("white", "black"),
                            width = unit(.2, "cm"),
                            show_row_names = FALSE)


ht_list = ht_list + Heatmap(hm_rows[,3], 
                            name = "p6h",
                            col = p_col_fun,
                            width = unit(.2, "cm"),
                            show_column_names = FALSE)


t = draw(ht_list, row_split = hm.clusters, show_row_den = F)




cbind(PcM.Diff[order(metric, decreasing = T),11:12][1:200,], 
      cl = hm.clusters) %>%
  pivot_longer(!cl, values_to = "reg") %>%
  group_by(cl, name, reg) %>%
  dplyr::summarise(q = n()) %>%
  group_by(cl, name) %>%
  dplyr::mutate(k = sum(q)) %>%
  filter(reg != "a")%>%
  ungroup() %>%
  left_join(PcM.sum, by = "name") %>%
  dplyr::mutate(p_val = phyper(q-1, m, n, k, lower.tail = F, log.p = T)) %>%
  dplyr::mutate(effect_size = (q/k)/(m/(m+n))) %>%
  dplyr::mutate(p = exp(p_val)) %>%
  kable(caption = "p_vals", format.args = list(big.mark = ",")) %>% 
  kable_classic()



### FigDiffChIPHeatmapH3K27me3AndAc.pdf

count.matrix = cbind(t(scale(t(PcM.Diff[,3:6]))), t(scale(t(PcM.Diff[,7:10]))))[order(metric, decreasing = T),][1:200,]

lab.df = PcM.Diff[order(metric, decreasing = T),11:12][1:200,]

hm.clusters = kmeans(count.matrix, centers = 3, algorithm = 'MacQueen', iter.max = 50, nstart = 10)$cluster

hm_rows = data.frame(cl = hm.clusters) %>%
  left_join(cbind(lab.df, 
                  cl = hm.clusters) %>%
              pivot_longer(!cl, values_to = "reg") %>%
              group_by(cl, name, reg) %>%
              dplyr::summarise(q = n()) %>%
              group_by(cl, name) %>%
              dplyr::mutate(k = sum(q)) %>%
              filter(reg != "a")%>%
              ungroup() %>%
              left_join(PcM.sum, by = "name") %>%
              dplyr::mutate(p_val = phyper(q-1, m, n, k, lower.tail = F)) %>%
              mutate(p.adj = p.adjust(p_val)) %>%
              filter(p.adj <= 0.05) %>%
              mutate(p.adj = log10(p.adj)) %>%
              pivot_wider(cl, names_from = name, values_from = p.adj),
            by = "cl") %>%
  replace_na(list(up3h = 0, up6h = 0))



conditions = c("Resting H3K27me3", "Septic H3K27me3", "Resting H3K27ac", "Septic H3K27ac")

ht_list = Heatmap(matrix = count.matrix,
                  column_order = 1 : 8,
                  column_split = rep(conditions, each = 2),
                  show_column_names = FALSE,
                  col = colorRampPalette(rcartocolor::carto_pal(n = 7, name = "ArmyRose"))(10),
                  show_row_names = FALSE, name = "row\nz score")

p_col_fun = colorRamp2(c(min(hm_rows[,2:3]), 0), c("#672044", "#FFFFFF"))

ht_list = ht_list + Heatmap(lab.df[,1], 
                            name = "3h", 
                            col = c("white", "black"),
                            width = unit(.2, "cm"),
                            show_row_names = FALSE)

ht_list = ht_list + Heatmap(hm_rows[,2], 
                            name = "p3h",
                            col = p_col_fun,
                            width = unit(.2, "cm"),
                            show_column_names = FALSE)

ht_list = ht_list + Heatmap(lab.df[,2], 
                            name = "6h", 
                            col = c("white", "black"),
                            width = unit(.2, "cm"),
                            show_row_names = FALSE)


ht_list = ht_list + Heatmap(hm_rows[,3], 
                            name = "p6h",
                            col = p_col_fun,
                            width = unit(.2, "cm"),
                            show_column_names = FALSE)


t = draw(ht_list, row_split = hm.clusters, show_row_den = F)




cbind(PcM.Diff[order(metric, decreasing = T),11:12][1:200,], 
      cl = hm.clusters) %>%
  pivot_longer(!cl, values_to = "reg") %>%
  group_by(cl, name, reg) %>%
  dplyr::summarise(q = n()) %>%
  group_by(cl, name) %>%
  dplyr::mutate(k = sum(q)) %>%
  filter(reg != "a")%>%
  ungroup() %>%
  left_join(PcM.sum, by = "name") %>%
  dplyr::mutate(p_val = phyper(q-1, m, n, k, lower.tail = F, log.p = T)) %>%
  dplyr::mutate(effect_size = (q/k)/(m/(m+n))) %>%
  dplyr::mutate(p = exp(p_val)) %>%
  kable(caption = "p_vals", format.args = list(big.mark = ",")) %>% 
  kable_classic()




# boxplots ----------------------------------------------------------------



position_jitternudge <- function(jitter.width = NULL, jitter.height = 0,
                                 nudge.x = 0, nudge.y = 0, seed = NA) {
  if (!is.null(seed) && is.na(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
  }
  
  ggplot2::ggproto(NULL, PositionJitternudge,
                   jitter.width = jitter.width,
                   jitter.height = jitter.height,
                   nudge.x = nudge.x,
                   nudge.y = nudge.y,
                   seed = seed
  )
}

PositionJitternudge <- ggplot2::ggproto("PositionJitternudge", ggplot2::Position,
                                        jitter.width = NULL,
                                        jitter.height = NULL,
                                        nudge.x = NULL,
                                        nudge.y = NULL,
                                        
                                        required_aes = c("x", "y"),
                                        
                                        setup_params = function(self, data) {
                                          flipped_aes <- ggplot2::has_flipped_aes(data)
                                          data <- ggplot2::flip_data(data, flipped_aes)
                                          width <- self$jitter.width %||% (ggplot2::resolution(data$x, zero = FALSE) * 0.4)
                                          
                                          list(
                                            nudge.x = self$nudge.x,
                                            nudge.y = self$nudge.y,
                                            jitter.height = self$jitter.height,
                                            jitter.width = width / 2, #(ndodge + 2),
                                            seed = self$seed,
                                            flipped_aes = flipped_aes
                                          )
                                        },
                                        
                                        compute_panel = function(data, params, scales) {
                                          data <- ggplot2::flip_data(data, params$flipped_aes)
                                          
                                          trans_x <- if(params$jitter.width > 0) function(x) {jitter(x, amount = params$jitter.width) + params$nudge.x}
                                          trans_y <- if(params$jitter.height > 0) function(x) {jitter(x, amount = params$jitter.height)  + params$nudge.y}
                                          
                                          data <- ggplot2:::with_seed_null(params$seed, ggplot2::transform_position(data, trans_x, trans_y))
                                          ggplot2::flip_data(data, params$flipped_aes)
                                        }
)


raincloud_theme = theme(axis.title.x = element_blank(),
                        legend.position = "none", 
                        strip.background = element_rect(colour="black", fill="white"))

wilcox.p = function(x,y){wilcox.test(x,y, paired = T)$p.value}

boxplot.df = PcM.Diff %>%
  select(!up18h)

boxplot.df[,3:6] = t(scale(t(boxplot.df[,3:6]), scale = F))
boxplot.df[,7:10] = t(scale(t(boxplot.df[,7:10]), scale = F))

### H3K27me3

H3K27me3Box = boxplot.df %>%
  mutate(mean.z.H3K27me3.infected = (boxplot.df$log.ratio.A.infectionH3K27me3 + boxplot.df$log.ratio.B.infectionH3K27me3)/2,
         mean.z.H3K27me3.control = (boxplot.df$log.ratio.B.controlH3K27me3 + boxplot.df$log.ratio.A.controlH3K27me3)/2) %>%
  select(mean.z.H3K27me3.infected,  mean.z.H3K27me3.control, up3h, up6h) %>%
  pivot_longer(!up3h:up6h) %>%
  mutate(up3h = recode(up3h, a = "not de 3h"),
         up6h = recode(up6h, a = "not de 6h"),
         name = recode(name, mean.z.H3K27me3.infected = "infected", mean.z.H3K27me3.control = "control"))

### FigDiffChIPBoxplotH3K27me3_6hNoReplicate.pdf

ggplot(H3K27me3Box, aes(name, value, fill = name)) +
  geom_flat_violin(position = position_nudge(x = 0, y = 0), alpha = .8, scale = "width") +
  geom_point(aes(y = value, color = name), position = 
               position_jitternudge(jitter.width = .3, nudge.x = -.2),
             size = .5, alpha = 0.6) +
  geom_boxplot(width = .3, outlier.shape = NA, alpha = 0.5, position = 
                 position_nudge(x = -.2, y = 0), notch = T) +
  scale_y_continuous(name = "mean centered H3K27me3 signal") + 
  scale_color_tableau(palette = "Superfishel Stone") +
  scale_fill_tableau(palette = "Superfishel Stone") +
  theme_bw() + raincloud_theme +
  facet_wrap(vars(up6h), nrow = 1) +
  ggtitle("H3K27me3 change in Pc-M genes after septic injury")

### FigDiffChIPBoxplotH3K27me3_3hNoReplicate.pdf

ggplot(H3K27me3Box, aes(name, value, fill = name)) +
  geom_flat_violin(position = position_nudge(x = 0, y = 0), alpha = .8, scale = "width") +
  geom_point(aes(y = value, color = name), position = 
               position_jitternudge(jitter.width = .3, nudge.x = -.2),
             size = .5, alpha = 0.6) +
  geom_boxplot(width = .3, outlier.shape = NA, alpha = 0.5, position = 
                 position_nudge(x = -.2, y = 0), notch = T) +
  scale_y_continuous(name = "mean centered H3K27me3 signal") + 
  scale_color_tableau(palette = "Superfishel Stone") +
  scale_fill_tableau(palette = "Superfishel Stone") +
  theme_bw() + raincloud_theme +
  facet_wrap(vars(up3h), nrow = 1) +
  ggtitle("H3K27me3 change in Pc-M genes after septic injury")

boxplot.df %>%
  mutate(mean.z.H3K27me3.infected = (boxplot.df$log.ratio.A.infectionH3K27me3 + boxplot.df$log.ratio.B.infectionH3K27me3)/2,
         mean.z.H3K27me3.control = (boxplot.df$log.ratio.B.controlH3K27me3 + boxplot.df$log.ratio.A.controlH3K27me3)/2) %>%
  select(mean.z.H3K27me3.infected,  mean.z.H3K27me3.control, up3h, up6h) %>%
  pivot_longer(up3h:up6h) %>%
  group_by(name, value) %>%
  dplyr::summarise(wilcox.p = wilcox.p(mean.z.H3K27me3.infected, mean.z.H3K27me3.control),
                   t.test.p = t.test(mean.z.H3K27me3.infected, mean.z.H3K27me3.control, paired = T)$p.value)  %>%
  mutate(comparison = case_when(name == "up3h" & value == "a" ~ "not de 3h",
                                name == "up6h" & value == "a" ~ "not de 6h",
                                T ~ value)) %>%
  select(comparison, wilcox.p, t.test.p) %>%
  kable(caption = "p_vals", format.args = list(big.mark = ",")) %>% 
  kable_classic()


H3K27me3Boxv2 = boxplot.df %>%
  select(log.ratio.A.controlH3K27me3:log.ratio.B.infectionH3K27me3, up3h, up6h) %>%
  pivot_longer(!up3h:up6h)  %>%
  mutate(up3h = recode(up3h, a = "not de 3h"),
         up6h = recode(up6h, a = "not de 6h"),
         name = recode(name, 
                       log.ratio.B.infectionH3K27me3 = "infected.B", 
                       log.ratio.A.controlH3K27me3 = "control.A", 
                       log.ratio.A.infectionH3K27me3 = "infected.A", 
                       log.ratio.B.controlH3K27me3 = "control.B"))


### FigDiffChIPBoxplotH3K27me3_6hWithReplicate.pdf

ggplot(H3K27me3Boxv2, aes(name, value, fill = name)) +
  geom_flat_violin(position = position_nudge(x = 0, y = 0), alpha = .8, scale = "width") +
  geom_point(aes(y = value, color = name), position = 
               position_jitternudge(jitter.width = .3, nudge.x = -.2),
             size = .5, alpha = 0.6) +
  geom_boxplot(width = .3, outlier.shape = NA, alpha = 0.5, position = 
                 position_nudge(x = -.2, y = 0), notch = T) +
  scale_y_continuous(name = "mean centered H3K27me3 signal") + 
  scale_color_tableau(palette = "Superfishel Stone") +
  scale_fill_tableau(palette = "Superfishel Stone") +
  theme_bw() + raincloud_theme +
  facet_wrap(vars(up6h), nrow = 1) +
  ggtitle("H3K27me3 change in Pc-M genes after septic injury")

### FigDiffChIPBoxplotH3K27me3_3hWithReplicate.pdf

ggplot(H3K27me3Boxv2, aes(name, value, fill = name)) +
  geom_flat_violin(position = position_nudge(x = 0, y = 0), alpha = .8, scale = "width") +
  geom_point(aes(y = value, color = name), position = 
               position_jitternudge(jitter.width = .3, nudge.x = -.2),
             size = .5, alpha = 0.6) +
  geom_boxplot(width = .3, outlier.shape = NA, alpha = 0.5, position = 
                 position_nudge(x = -.2, y = 0), notch = T) +
  scale_y_continuous(name = "mean centered H3K27me3 signal") + 
  scale_color_tableau(palette = "Superfishel Stone") +
  scale_fill_tableau(palette = "Superfishel Stone") +
  theme_bw() + raincloud_theme +
  facet_wrap(vars(up3h), nrow = 1) +
  ggtitle("H3K27me3 change in Pc-M genes after septic injury")

boxplot.df %>%
  select(gene_id, log.ratio.A.controlH3K27me3:log.ratio.B.infectionH3K27me3, up3h, up6h) %>%
  pivot_longer(up3h:up6h) %>%
  pivot_longer(log.ratio.A.controlH3K27me3:log.ratio.B.infectionH3K27me3, names_to = "groups", values_to = "y") %>%
  group_by(name, value) %>%
  dplyr::summarise(friedman.p = friedman.test(y, groups, gene_id)$p.value) %>%
  mutate(comparison = case_when(name == "up3h" & value == "a" ~ "not de 3h",
                                name == "up6h" & value == "a" ~ "not de 6h",
                                T ~ value)) %>%
  select(comparison, friedman.p) %>%
  kable(caption = "p_vals", format.args = list(big.mark = ",")) %>% 
  kable_classic()


### H3K27ac


H3K27acBox = boxplot.df %>%
  mutate(mean.z.H3K27ac.infected = (boxplot.df$log.ratio.A.infectionH3K27ac + boxplot.df$log.ratio.B.infectionH3K27ac)/2,
         mean.z.H3K27ac.control = (boxplot.df$log.ratio.B.controlH3K27ac + boxplot.df$log.ratio.A.controlH3K27ac)/2) %>%
  select(mean.z.H3K27ac.infected,  mean.z.H3K27ac.control, up3h, up6h) %>%
  pivot_longer(!up3h:up6h) 

### FigDiffChIPBoxplotH3K27ac_6hNoReplicate.pdf

ggplot(H3K27acBox, aes(name, value, fill = name)) +
  geom_flat_violin(position = position_nudge(x = 0, y = 0), alpha = .8, scale = "width") +
  geom_point(aes(y = value, color = name), position = 
               position_jitternudge(jitter.width = .3, nudge.x = -.2),
             size = .5, alpha = 0.6) +
  geom_boxplot(width = .3, outlier.shape = NA, alpha = 0.5, position = 
                 position_nudge(x = -.2, y = 0), notch = T) +
  scale_y_continuous(name = "mean H3K27ac z-score") + 
  scale_color_tableau(palette = "Superfishel Stone") +
  scale_fill_tableau(palette = "Superfishel Stone") +
  theme_bw() + raincloud_theme +
  facet_wrap(vars(up6h), nrow = 1)

### FigDiffChIPBoxplotH3K27ac_3hNoReplicate.pdf

ggplot(H3K27acBox, aes(name, value, fill = name)) +
  geom_flat_violin(position = position_nudge(x = 0, y = 0), alpha = .8, scale = "width") +
  geom_point(aes(y = value, color = name), position = 
               position_jitternudge(jitter.width = .3, nudge.x = -.2),
             size = .5, alpha = 0.6) +
  geom_boxplot(width = .3, outlier.shape = NA, alpha = 0.5, position = 
                 position_nudge(x = -.2, y = 0), notch = T) +
  scale_y_continuous(name = "mean H3K27ac z-score") + 
  scale_color_tableau(palette = "Superfishel Stone") +
  scale_fill_tableau(palette = "Superfishel Stone") +
  theme_bw() + raincloud_theme +
  facet_wrap(vars(up3h), nrow = 1)


boxplot.df %>%
  mutate(mean.z.H3K27ac.infected = (boxplot.df$log.ratio.A.infectionH3K27ac + boxplot.df$log.ratio.B.infectionH3K27ac)/2,
         mean.z.H3K27ac.control = (boxplot.df$log.ratio.B.controlH3K27ac + boxplot.df$log.ratio.A.controlH3K27ac)/2) %>%
  select(mean.z.H3K27ac.infected,  mean.z.H3K27ac.control, up3h, up6h) %>%
  pivot_longer(up3h:up6h) %>%
  group_by(name, value) %>%
  dplyr::summarise(wilcox.p = wilcox.p(mean.z.H3K27ac.infected, mean.z.H3K27ac.control),
                   t.test.p = t.test(mean.z.H3K27ac.infected, mean.z.H3K27ac.control, paired = T)$p.value)  %>%
  mutate(comparison = case_when(name == "up3h" & value == "a" ~ "not de 3h",
                                name == "up6h" & value == "a" ~ "not de 6h",
                                T ~ value)) %>%
  select(comparison, wilcox.p, t.test.p) %>%
  kable(caption = "p_vals", format.args = list(big.mark = ",")) %>% 
  kable_classic()


H3K27acBoxv2 = boxplot.df %>%
  select(log.ratio.A.controlH3K27ac:log.ratio.B.infectionH3K27ac, up3h, up6h) %>%
  pivot_longer(!up3h:up6h)  %>%
  mutate(up3h = recode(up3h, a = "not de 3h"),
         up6h = recode(up6h, a = "not de 6h"),
         name = recode(name, 
                       log.ratio.B.infectionH3K27ac = "infected.B", 
                       log.ratio.A.controlH3K27ac = "control.A", 
                       log.ratio.A.infectionH3K27ac = "infected.A", 
                       log.ratio.B.controlH3K27ac = "control.B"))

### FigDiffChIPBoxplotH3K27ac_6hWithReplicate.pdf

ggplot(H3K27acBoxv2, aes(name, value, fill = name)) +
  geom_flat_violin(position = position_nudge(x = 0, y = 0), alpha = .8, scale = "width") +
  geom_point(aes(y = value, color = name), position = 
               position_jitternudge(jitter.width = .3, nudge.x = -.2),
             size = .5, alpha = 0.6) +
  geom_boxplot(width = .3, outlier.shape = NA, alpha = 0.5, position = 
                 position_nudge(x = -.2, y = 0), notch = T) +
  scale_y_continuous(name = "mean centered H3K27ac signal", limits = c(-2, 2)) + 
  scale_color_tableau(palette = "Superfishel Stone") +
  scale_fill_tableau(palette = "Superfishel Stone") +
  theme_bw() + raincloud_theme +
  facet_wrap(vars(up6h), nrow = 1) +
  ggtitle("H3K27ac change in Pc-M genes after septic injury")

### FigDiffChIPBoxplotH3K27ac_3hWithReplicate.pdf

ggplot(H3K27acBoxv2, aes(name, value, fill = name)) +
  geom_flat_violin(position = position_nudge(x = 0, y = 0), alpha = .8, scale = "width") +
  geom_point(aes(y = value, color = name), position = 
               position_jitternudge(jitter.width = .3, nudge.x = -.2),
             size = .5, alpha = 0.6) +
  geom_boxplot(width = .3, outlier.shape = NA, alpha = 0.5, position = 
                 position_nudge(x = -.2, y = 0), notch = T) +
  scale_y_continuous(name = "mean centered H3K27ac signal", limits = c(-2, 2)) + 
  scale_color_tableau(palette = "Superfishel Stone") +
  scale_fill_tableau(palette = "Superfishel Stone") +
  theme_bw() + raincloud_theme +
  facet_wrap(vars(up3h), nrow = 1) +
  ggtitle("H3K27me3 change in Pc-M genes after septic injury")

boxplot.df %>%
  select(gene_id, log.ratio.A.controlH3K27ac:log.ratio.B.infectionH3K27ac, up3h, up6h) %>%
  pivot_longer(up3h:up6h) %>%
  pivot_longer(log.ratio.A.controlH3K27ac:log.ratio.B.infectionH3K27ac, names_to = "groups", values_to = "y") %>%
  group_by(name, value) %>%
  dplyr::summarise(friedman.p = friedman.test(y, groups, gene_id)$p.value) %>%
  mutate(comparison = case_when(name == "up3h" & value == "a" ~ "not de 3h",
                                name == "up6h" & value == "a" ~ "not de 6h",
                                T ~ value)) %>%
  select(comparison, friedman.p) %>%
  kable(caption = "p_vals", format.args = list(big.mark = ",")) %>% 
  kable_classic()


# Correlation De and DiffChip ---------------------------------------------

CorrPlots = DiffChIP_data %>%
  filter(gene_id %in% (GeneGFF %>%
                         filter(gene_biotype == "protein_coding") %>%
                         filter(seqid %in% c("X", "2R", "2L", "3R", "3L")) %>%
                         .$gene_id)) %>%
  mutate(up3h = case_when((as.numeric(decideTests.DGELRT(QFresults$de_3h)) == 1)[match(gene_id, DGEobject$genes$gene_id)] ~ "up 3h",
                          T ~ ""),
         up6h = case_when((as.numeric(decideTests.DGELRT(QFresults$de_6h)) == 1)[match(gene_id, DGEobject$genes$gene_id)] ~ "up 6h",
                          T ~ "")) %>%
  filter(up3h == "up 3h" | up6h == "up 6h") %>%
  left_join(resting_H3K27me3 %>%
              select(class, Geneid),
            by = c("gene_id" = "Geneid")) %>%
  filter(class != "Pc-H") %>%
  mutate(H3K27me3.resting = (log.ratio.B.controlH3K27me3 + log.ratio.A.controlH3K27me3)/2,
         H3K27me3.change = ((log.ratio.A.infectionH3K27me3 + log.ratio.B.infectionH3K27me3) - 
                              (log.ratio.B.controlH3K27me3 + log.ratio.A.controlH3K27me3)),
         H3K27ac.change = ((log.ratio.A.infectionH3K27ac + log.ratio.B.infectionH3K27ac) - 
                             (log.ratio.B.controlH3K27ac + log.ratio.A.controlH3K27ac))) %>%
  select(gene_id, gene_name, up3h, up6h, class, H3K27me3.change, H3K27ac.change, H3K27me3.resting) %>%
  mutate(lfc.3h = QFresults$de_3h$table$logFC[match(gene_id, DGEobject$genes$gene_id)],
         lfc.6h = QFresults$de_6h$table$logFC[match(gene_id, DGEobject$genes$gene_id)])


### FigDiffChIPCorrToDe3h.pdf

CorrPlots %>%
  filter(up3h == "up 3h") %>%
  filter(lfc.3h > 1) %>%
  ggplot(aes(H3K27me3.change, lfc.3h,color = class)) +
  geom_point() + stat_smooth(method = "lm") +
  theme_bw() + ylab("RNA-seq LFC") + xlab("ChIP-seq LFC") + 
  ggtitle("Change Correlation at 3h Septic Injury") +
  theme(legend.position = c(.9, .88), 
        legend.background = element_rect(linetype = 2, size = 0.2, colour = 1))+
  scale_color_tableau(palette = "Superfishel Stone", name = "Gene state")


CorrPlots %>%
  filter(up3h == "up 3h") %>%
  select(class, H3K27me3.change, H3K27ac.change, lfc.3h) %>%
  pivot_longer(H3K27me3.change:H3K27ac.change, names_to = "modification", values_to = "h3k27.lfc") %>%
  dplyr::group_by(class, modification) %>%
  dplyr::summarise(broom::tidy(cor.test(h3k27.lfc, lfc.3h))) %>%
  kable(caption = "p_vals", format.args = list(big.mark = ",")) %>% 
  kable_classic()

### FigDiffChIPCorrToDe6h.pdf

CorrPlots %>%
  filter(up6h == "up 6h") %>%
  filter(lfc.6h > 0) %>%
  ggplot(aes(H3K27me3.change, lfc.6h, color = class)) +
  geom_point() + stat_smooth(method = "lm") +
  theme_bw() + ylab("RNA-seq LFC") + xlab("ChIP-seq LFC") + 
  ggtitle("Change Correlation at 6h Septic Injury") +
  theme(legend.position = c(.9, .88), 
        legend.background = element_rect(linetype = 2, size = 0.2, colour = 1))+
  scale_color_tableau(palette = "Superfishel Stone", name = "Gene state")


CorrPlots %>%
  filter(up6h == "up 6h") %>%
  select(class, H3K27me3.change, H3K27ac.change, lfc.6h) %>%
  pivot_longer(H3K27me3.change:H3K27ac.change, names_to = "modification", values_to = "h3k27.lfc") %>%
  dplyr::group_by(class, modification) %>%
  dplyr::summarise(broom::tidy(cor.test(h3k27.lfc, lfc.6h))) %>%
  kable(caption = "p_vals", format.args = list(big.mark = ",")) %>% 
  kable_classic()

