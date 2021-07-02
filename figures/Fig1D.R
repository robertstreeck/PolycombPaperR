library(tidyverse)
library(circlize)
library(ggthemes)
library(ComplexHeatmap)

load("data/fig1/SevenClassGenomeModel.Rdata")
load("data/fig1/resting_H3K27me3.Rdata")

class_dict = data.frame(group = 1:7,
                        group_name = c("EnhW", "Pc-I", "TEl", "EnhS", "Pc-H", "TSS", "Het"))

class_dict$group = as.character(class_dict$group)
GeneTable = GeneTable %>% left_join(class_dict, by = c("MajorityStateVote" = "group"))

GeneTableMerged = resting_H3K27me3 %>% select(Geneid, class) %>%
  right_join(GeneTable %>% select(gene_id, group_name), by = c("Geneid" = "gene_id")) %>%
  group_by(class, group_name) %>% summarise(count = n()) %>% filter(!is.na(group_name))

CircPlot = GeneTableMerged
colnames(CircPlot) = c("from", "to", "value")
CircPlot$from = factor(paste0("R_", CircPlot$from),
                       levels = c("R_non-Pc", "R_Pc-I", "R_Pc-H"))
CircPlot$to = factor(paste0("C_", CircPlot$to),
                     levels = c("C_EnhW", "C_EnhS", "C_TEl", "C_Pc-H", "C_Pc-I", "C_TSS", "C_Het"))

cols = c(ggthemes_data$tableau$`color-palettes`$regular$`Superfishel Stone`$value[1:3], 
         ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value[1:7])
names(cols) = c("R_non-Pc", "R_Pc-I", "R_Pc-H", 
                "C_TEl", "C_Pc-I", "C_Pc-H","C_TSS", "C_EnhS","C_EnhW", "C_Het")

plot.new()
circle_size = unit(1, "snpc")

pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))

circos.par(start.degree = 90)
cdm_res = chordDiagram(CircPlot, grid.col = cols,
                       order = c("R_Pc-H", "R_Pc-I", "R_non-Pc", "C_TSS", "C_EnhS", "C_TEl", 
                                 "C_EnhW", "C_Het", "C_Pc-I", "C_Pc-H"),
                       annotationTrack = c("grid", "axis"), annotationTrackHeight = mm_h(2), 
                       big.gap = 80, small.gap = 1,
                       preAllocateTracks = list(track.height = 0.1))

cex_adj = rep(.6, 11)
names(cex_adj) = names(cols) 

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim)+0.5, substr(si,3,10), sector.index = si, track.index = 1, 
              facing = "bending.inside", niceFacing = TRUE, col = "black", font = 2, cex = cex_adj[si])
}
