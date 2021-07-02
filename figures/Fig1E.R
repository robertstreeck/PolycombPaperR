library(tidyverse)
library(ggplot2)
library(ggthemes)

load("data/fig1/resting_H3K27me3.Rdata")

tpm = function(x, len){
  x = (x+1)/len
  1e6*x/sum(x)
}

resting_RNAseq = read.delim("data/fig1/2426.count", sep = "\t", skip = 1)
colnames(resting_RNAseq)[7:9] = c("RepA", "RepB", "RepC")
resting_RNAseq = resting_RNAseq %>% select(Geneid, Length:RepC) %>%
  mutate(RepA_tpm = tpm(RepA, Length),
         RepB_tpm = tpm(RepB, Length),
         RepC_tpm = tpm(RepC, Length)) %>%
  mutate(mean_tpm = (RepA_tpm + RepB_tpm + RepC_tpm)/3) 

resting_RNAseq = resting_RNAseq %>%
  left_join(resting_H3K27me3 %>% 
              select(Geneid, gene_name, gene_biotype, class),
            by = "Geneid") %>% 
  filter(gene_biotype == "protein_coding")

resting_RNAseq$class = factor(resting_RNAseq$class, 
                              levels = c("non-Pc", "Pc-I", "Pc-H"))

source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")


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

ggplot(resting_RNAseq, aes(class, mean_tpm, fill = class)) +
  geom_flat_violin(position = position_nudge(x = 0, y = 0), alpha = .8, scale = "width") +
  geom_point(aes(y = mean_tpm, color = class), position = 
               position_jitternudge(jitter.width = .3, nudge.x = -.2),
             size = .5, alpha = 0.6) +
  geom_boxplot(width = .3, outlier.shape = NA, alpha = 0.5, position = 
                 position_nudge(x = -.2, y = 0)) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000), labels = c(0.1, 1, 10, 100, 1000), 
                name = "mean tags per million") + 
  scale_color_tableau(palette = "Superfishel Stone") +
  scale_fill_tableau(palette = "Superfishel Stone") +
  theme_bw() + raincloud_theme 
