library(reshape2)
library(openxlsx)
library(plyr)
library(dplyr)
library(rcartocolor)


# Imports -----------------------------------------------------------------



read_single_cell_from_xlsx = function(con, col_selection = NULL){
  require(openxlsx)
  wb = loadWorkbook(con)
  sheet_names = sheets(wb)
  if(is.null(col_selection)){
    col_selection = 1:dim(readWorkbook(wb))[2]
  }
  long_table = readWorkbook(wb)
  long_table$cluster = sheet_names[1]
  for (i in 2:length(sheet_names)) {
    long_table = rbind(long_table, data.frame(readWorkbook(wb, sheet = i), cluster = sheet_names[i]))
  }
  return(long_table)
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


DRSC = read_single_cell_from_xlsx("Data/fig1/DRSC-allcombined.xlsx")
Sublemental = read_single_cell_from_xlsx("Data/fig1/elife-54818-supp3-v2.xlsx")

load("data/fig1/resting_H3K27me3.Rdata")
load("data/fig1/SevenClassGenomeModel.Rdata")

DRSC$gene_name = trimws(DRSC$gene_name)


log_ratio = function(chip, input){
  chip = (chip)/sum(chip)
  input = (input)/sum(input)
  return(log2(chip) - log2(input))
}



resting_H3K27me3 = resting_H3K27me3 %>% 
  select(Geneid, class, H3K27me3_A:input_B) %>%
  mutate(lr_A = log_ratio(H3K27me3_A, input_A),
         lr_B = log_ratio(H3K27me3_B, input_B)) %>%
  mutate(h3k27me3_mean_lr = (lr_A + lr_B)/2) %>%
  mutate(class = recode(class, "Pc-I" = "Pc-M")) %>%
  mutate(class = factor(class, levels = c("non-Pc", "Pc-M", "Pc-H")))


DRSC_ridge = DRSC %>% plyr::rename(c("sum" = "expression")) %>%
  left_join(resting_H3K27me3 %>% select(Geneid, class), by = c("gene_name" = "Geneid")) %>%
  filter(!is.na(class)) %>%
  group_by(cluster) %>%
  dplyr::mutate(cpm = 1e6*expression/sum(expression))

DRSC_expressed = DRSC_ridge %>% 
  filter(cpm >= 1) %>%
  select(gene_name) %>%
  group_by(gene_name) %>% 
  dplyr::summarise(n_expressed = n())

For_hyper = resting_H3K27me3 %>% select(Geneid, class) %>%
  left_join(DRSC_expressed, by = c("Geneid" = "gene_name")) %>%
  mutate_if(is.numeric, coalesce,0)

class_dict = data.frame(group = as.character(1:7),
                        group_name = c("EnhW", "Pc-I", "TEl", "EnhS", "Pc-H", "TSS", "Het"))

GeneTable = GeneTable %>% 
  left_join(resting_H3K27me3 %>% select(Geneid, h3k27me3_mean_lr, class), by = c("gene_id" = "Geneid")) %>%
  left_join(class_dict, by = c("MajorityStateVote" = "group"))


df = GeneTable %>% select(gene_id, MajorityStateVote) %>%
  left_join(class_dict, by = c("MajorityStateVote" = "group")) %>%
  select(gene_id, group_name) %>%
  left_join(DRSC_expressed, by = c("gene_id" = "gene_name")) %>%
  mutate_if(is.numeric, coalesce,0) %>%
  plyr::rename(c("gene_id" = "Geneid", "group_name" = "class")) %>%
  filter(!is.na(class))

For_hyper = rbind(data.frame(method = "H3K27me3\nEnrichment", For_hyper),
                  data.frame(method = "Modification\nstate", df))

### FigSingleCellHistograms.pdf

library(ggplot2)
library(ggridges)
library(ggthemes)
ggplot(DRSC_ridge, aes(cpm, fill=class, y = cluster, height = ..count..)) + 
  geom_density_ridges(stat = "binline", alpha = .5, bins = 100) + xlab("CPM") + ylab("scRNA-seq cell cluster") + 
  scale_fill_tableau(palette = "Superfishel Stone", name = "Gene State") + 
  theme_bw() + scale_x_log10() + coord_flip()


hyper = For_hyper %>% group_by(class, n_expressed, method) %>%
  dplyr::summarize(q = n()) %>%
  group_by(class, method) %>%
  dplyr::mutate(k = sum(q)) %>%
  group_by(n_expressed, method) %>%
  dplyr::mutate(m = sum(q)) %>%
  dplyr::mutate(n = dim(For_hyper)[1] - m) %>%
  dplyr::mutate(p_val = phyper(q-1, m, n, k, lower.tail = F, log.p = T)) %>%
  dplyr::mutate(effect_size = (q/k)/(m/(m+n)))

hyper$log10_p_val = -hyper$p_val/log(10)
hyper$y_fact = factor(as.character(hyper$n_expressed), levels = as.character(0:15))

hyper$class = factor(hyper$class, 
                     levels = c("non-Pc", "TSS", "TEl", "EnhS", "EnhW", "Het", "Pc-I", "Pc-M", "Pc-H"))

hyper$method = factor(hyper$method, 
                      levels = c("Modification\nstate", "H3K27me3\nEnrichment"))


### FigSingleCellEnrichment.pdf
require(scales)
pmax = 40
p = ggplot(hyper[hyper$log10_p_val > -log10(0.05),], aes(class, y_fact, fill = log10_p_val, size = effect_size)) +
  geom_point(shape = 21) + scale_y_discrete(name = "# of clusters expressed in") +
  scale_fill_viridis_c(option = "magma", limits = c(-log10(0.05), pmax), oob = squish,
                     breaks = c(10, 20, 30, 40), 
                     labels = c(expression(10^-10),expression(10^-20),expression(10^-30), expression("<"~10^-40))) + 
  facet_wrap(vars(method), nrow = 1, scales = "free_x") + theme_bw() + 
  guides(fill = guide_legend(order = 1, override.aes = list(size=6)),
         size = guide_legend(override.aes = list(fill="grey"))) +
  scale_size_continuous(name = "Effect\nsize", range = c(1.5,6)) +
  theme(axis.title.x = element_blank())

gp = adjust_ggplot_facet_size(p)
grid::grid.draw(gp)

