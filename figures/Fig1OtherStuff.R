### ChIP corr heatmap

###

library(ggplot2)
library(scales)
#### loading data set

# chromatin state assignment
load("data/fig1/resting_H3K27me3.Rdata")

# TF binding sites
load("/Users/streeck/Desktop/PolycombPaper/ExternalGeneSets/NarrowPeak_gene_list.Rdata")

# RNA-seq and other datasets
load("/Users/streeck/Desktop/PolycombPaper/ExternalGeneSets/Gene_sets_by_Alf.Rdata")



# functions ---------------------------------------------------------------


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


# Data analysis -----------------------------------------------------------

gene_set_df = resting_H3K27me3

for (i in names(Alfs_Gene_sets)) {
  gene_set_df[,i] = gene_set_df$Geneid %in% Alfs_Gene_sets[[i]]
}

for (i in names(Peak_gene_list)) {
  gene_set_df[,i] = gene_set_df$Geneid %in% Peak_gene_list[[i]]
}

bg_expectation = resting_H3K27me3 %>%
  dplyr::select(class) %>%
  group_by(class) %>%
  dplyr::summarise(count = n()) %>%
  ungroup() %>%
  dplyr::mutate(expected = count/sum(count))

gene_set_df2 = gene_set_df %>%
  dplyr::select(class:Stat92E_White) %>%
  pivot_longer(!class) %>%
  filter(value) %>%
  group_by(class, name) %>%
  dplyr::summarise(count = n()) %>%
  group_by(name) %>%
  dplyr::mutate(fraction = count/sum(count)) %>%
  left_join(bg_expectation %>%
              dplyr::select(class, expected),
            by = "class") %>%
  dplyr::mutate(odds = fraction/expected) %>%
  group_by(name) %>%
  filter(sum(count) >= 10)

y.groups = setNames(c(rep("EcR", 3), rep("dpp", 2), rep("Other", 2), "Immunity", rep("Development", 5), 
                      "Immunity", rep("Development", 2), "Jak-Stat",rep("Other", 5), "Immunity", "Other", 
                      "Jak-Stat", "Other", "Other",  "Immunity", "Other", "Jak-Stat", "Jak-Stat", rep("Immunity", 3),
                      "Other", "Immunity", "EcR", "Other", "Other", "Other", "Jak-Stat", "Jak-Stat"),
                    c(names(Alfs_Gene_sets), names(Peak_gene_list)))


gene_set_df2 = gene_set_df2 %>%
  dplyr::mutate(yfac = recode(name, !!!y.groups)) %>%
  filter(yfac != "dpp")

p = ggplot(gene_set_df2, aes(class, name, fill = odds)) + geom_tile() +
  scale_fill_distiller(palette = "Greens", direction = 1, limits = c(0, 2), oob = squish) +
  facet_wrap(vars(yfac), ncol = 1, scales = "free_y") +
  theme_bw()

gp = adjust_ggplot_facet_size(p)
grid::grid.draw(gp)


### heatmap
heatmap_table = as.matrix(gene_sets_wide_table[,2:11])
row.names(heatmap_table) = gene_sets_wide_table$sets
heatmap_table = t(scale(log(t(heatmap_table)+1)))
heatmap_table = reshape2::melt(data.frame(set = row.names(heatmap_table), heatmap_table, row.names = NULL))
heatmap_table$set = factor(as.character(heatmap_table$set), levels = names(c(Alfs_Gene_sets, Peak_gene_list)))
heatmap_table$variable = factor(as.character(heatmap_table$variable), levels = c("TEl", "TSS", "EnhS", "EnhW", "Het", "PcI", "PcH", "K27low", "K27med", "K27high"))
heatmap_table$x.facet = c("Chromatin State", "Gene State")[as.numeric(heatmap_table$variable %in% c("K27med", "K27low", "K27high"))+1]

p = ggplot(heatmap_table, aes(variable, set, fill = value)) + geom_tile() + scale_fill_distiller(palette = "RdBu", name = "z-score") +
  theme_classic() + theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = .8, hjust = .8)) +
  facet_grid(cols = vars(x.facet), scale = "free")
gp = adjust_ggplot_facet_size(p)
grid::grid.draw(gp)


heatmap_table2 = heatmap_table
heatmap.yfacets = data.frame(lable = unique(heatmap_table2$set))
heatmap.yfacets = data.frame(set = heatmap.yfacets[!grepl("Core Component", heatmap.yfacets$lable),])
heatmap.yfacets$y.facet = c("Jak-Stat", "Immunity", "Development", "EcR", "dpp", "Other")[c(4,4,4,5,5,6,6,2,3,3,3,3,3,2,3,3,1,1,1,2,2,2,6,2,4,6,6,6,1,1)]
heatmap_table2$y.facet = heatmap.yfacets$y.facet[match(heatmap_table2$set, heatmap.yfacets$set)]
heatmap_table2 = heatmap_table2[!grepl("Core Component", heatmap_table2$set),]

p = ggplot(heatmap_table2, aes(variable, set, fill = value)) + geom_tile() + scale_fill_distiller(palette = "RdBu", name = "z-score") +
  theme_classic() + theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = .8, hjust = .8)) +
  facet_grid(cols = vars(x.facet), rows = vars(y.facet), scale = "free")
gp = adjust_ggplot_facet_size(p)
grid::grid.draw(gp)