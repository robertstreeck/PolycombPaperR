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
  pivot_wider(name, names_from = class, values_from = count, values_fill = 0) %>%
  pivot_longer(!name, values_to = "count", names_to = "class") %>%
  group_by(name) %>%
  dplyr::mutate(fraction = count/sum(count)) %>%
  left_join(bg_expectation %>%
              dplyr::select(class, expected),
            by = "class") %>%
  dplyr::mutate(odds = fraction/expected) %>%
  group_by(name) %>%
  filter(sum(count) >= 10) %>%
  mutate(class = factor(class, 
                        levels = c("non-Pc", "Pc-I", "Pc-H")))

y.groups = setNames(c(rep("EcR", 3), rep("dpp", 2), rep("Other", 2), "Immunity", rep("Development", 5), 
                      "Immunity", rep("Development", 2), "Jak-Stat",rep("Other", 5), "Immunity", "Other", 
                      "Jak-Stat", "Other", "Other",  "Immunity", "Other", "Jak-Stat", "Jak-Stat", rep("Immunity", 3),
                      "Other", "Immunity", "EcR", "Other", "Other", "Other", "Jak-Stat", "Jak-Stat"),
                    c(names(Alfs_Gene_sets), names(Peak_gene_list)))


gene_set_df2 = gene_set_df2 %>%
  dplyr::mutate(yfac = recode(name, !!!y.groups)) %>%
  filter(yfac != "dpp")


### Fig1OtherStuffHeatmapOfGeneFunctions.pdf

p = ggplot(gene_set_df2, aes(class, name, fill = fraction)) + geom_tile() +
  scale_fill_gradient(low = "#FFFFFF", high = "#ef6f6a", 
                      limits = c(0, 1), oob = squish, name = "Fraction\nof genes") +
  facet_wrap(vars(yfac), ncol = 1, scales = "free_y", strip.position = "right") +
  theme_bw() + theme(axis.title = element_blank())

gp = adjust_ggplot_facet_size(p)
grid::grid.draw(gp)






