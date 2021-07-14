library(tidyverse)
library(ggplot2)
library(ggthemes)

load(file = "2021-06-08-LargeOverwieTableForAlf.Rdata")

ViolinData = LargeOverviewTable %>%
  select(gene_id, H3K27me3.log_ratio, de_3h.de, de_6h.de, de_18h.de,
         eGFP_induction.de, hop_induction.de, PGRP_induction.de) %>%
  pivot_longer(!gene_id:H3K27me3.log_ratio, names_to = "contrast") %>%
  filter(!is.na(value)) %>%
  mutate(de.reduced = case_when(value == 1 ~ "Up",
                                T ~ "not Up")) %>%
  dplyr::mutate(xfact = case_when(grepl("induc", contrast) ~ "Switch Gal4",
                           T ~ "Septic injury")) %>%
  mutate(contrast = recode(contrast, "de_18h.de" = "18h", "de_6h.de" = "6h",
                       "de_3h.de" = "3h", "eGFP_induction.de" = "eGFP",
                       "hop_induction.de" = "hop", "PGRP_induction.de" = "PGRP-LC"))


### Fig2H3K27me3Violins.pdf

ggplot(ViolinData, aes(contrast, H3K27me3.log_ratio, fill = de.reduced)) + 
  geom_violin(width = 1) + 
  scale_fill_tableau(palette = "Tableau 10", name = "Genes") +
  theme_bw() + 
  geom_boxplot(aes(group = paste0(de.reduced, contrast)), 
               outlier.shape = NA, fill = "#FFFFFF00", width = .2,
               notch = T, position = position_dodge(width = 1)) +
  facet_wrap(vars(xfact), nrow = 1, scales = "free_x")
