library(tidyverse)
library(ggplot2)
library(ggthemes)
library(openxlsx)
library(data.table)

lr_matrix = function(chip, input){
  chip = apply(chip + 1, 2, function(x){x/sum(x)})
  input = apply(input + 1, 2, function(x){x/sum(x)})
  return(log2(chip/input))
}

hclust_levels = function(df, row_s, col_s, val_s){
  df_wide = df %>% select(row_s, col_s, val_s) %>%
    pivot_wider(id_cols = row_s, names_from = col_s, values_from = val_s)
  df_wide = as.data.frame(df_wide)
  h_cluster = hclust(dist(df_wide[2:dim(df_wide)[2]]))
  return(as.character(df_wide[h_cluster$order,1]))
}

chip_files = read.xlsx("data/fig1/SampleList.xlsx")
chip_file_selection = chip_files %>% select(Index:Input.ref, Replicate) %>%
  filter(Antibody != "Input") %>% 
  left_join(chip_files %>% select(Index, Bam.File), by = c("Input.ref" = "Index"),
            suffix = c(".chip", ".input"))


load("data/fig1/SevenClassGenomeModel.Rdata")

heatmap_table = data.frame(
  group = c("EnhW", "Pc-I", "TEl", "EnhS", "Pc-H", "TSS", "Het")[as.numeric(multi_chip_fit$Group)],
  lr_matrix(multi_chip_fit$ChIP[multi_chip_fit$excluded,], multi_chip_fit$Input[multi_chip_fit$excluded,])) %>%
  pivot_longer(!group) %>%
  left_join(chip_files %>%
              mutate(name = gsub("-", ".", Bam.File))%>%
              select(name, Antibody),
            by = "name") %>%
  select(!name) %>%
  group_by(group, Antibody) %>%
  summarise(mean_lr = mean(value))


heatmap_table$Antibody = factor(heatmap_table$Antibody, 
                              levels = hclust_levels(heatmap_table, "Antibody", "group", "mean_lr"))

heatmap_table$group = factor(heatmap_table$group, 
                                levels = hclust_levels(heatmap_table, "group", "Antibody", "mean_lr"))

ggplot(heatmap_table, aes(Antibody, group, fill = mean_lr)) + geom_tile() +
  theme_classic() + viridis::scale_fill_viridis(name = "mean\nlog2 ratio", direction = 1, option = "viridis") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 30, hjust = .8, vjust = 1),
        axis.ticks = element_blank(), axis.line = element_blank())

  