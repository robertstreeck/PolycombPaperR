library(GenomicRanges)
library(rtracklayer)
library(fmsb)
library(tidyverse)
library(rcartocolor)

load("data/fig1/SevenClassGenomeModel.Rdata")

genome = read.delim("/Users/streeck/Genomes/DmelBDGP6.91/chrNameLength.txt", header = F, stringsAsFactors = F)
genome = genome[1:7,]
gr = GRanges(genome[,1], IRanges(1, as.integer(genome[,2])))

chr_lengths = gr@ranges@width
names(chr_lengths) = gr@seqnames

tiled_genome = tileGenome(chr_lengths, tilewidth = 200, cut.last.tile.in.chrom = T)
tiled_genome$cluster = NA
tiled_genome$cluster[multi_chip_fit$excluded] = multi_chip_fit$Group

map_to_clusters = function(DataSet){
  DataSet[,2] = as.numeric(DataSet[,2])
  DataSet[,3] = as.numeric(DataSet[,3])
  DataSet$Cluster = NA
  
  DataSet_GR = GRanges(DataSet$seqid, ranges = IRanges(DataSet$start, DataSet$end))
  
  for (i in 1:length(DataSet$seqid)) {
    DataSet$Cluster[i] = names(sort(-table(subsetByOverlaps(tiled_genome, DataSet_GR[i], minoverlap = 1)$cluster)))[1]
  }
  print(table(DataSet$Cluster))
  return(DataSet)
}


load("data/PREs/Ederle_PRE_HC_dm6.Rdata")
Ederle_PRE_HC_dm6 = map_to_clusters(Ederle_PRE_HC_dm6)
RadarPlot = data.frame(Set = "Ederle_PRE_HC", state = Ederle_PRE_HC_dm6$Cluster)

load("data/PREs/Kahn_PRE_dm6.Rdata")
Kahn_PRE_dm6 = map_to_clusters(Kahn_PRE_dm6)
RadarPlot = rbind(RadarPlot, data.frame(Set = "Kahn_PRE", state = Kahn_PRE_dm6$Cluster))

load("data/PREs/Schwartz_PRE_dm6.Rdata")
Schwartz_PRE_dm6 = map_to_clusters(Schwartz_PRE_dm6)
RadarPlot = rbind(RadarPlot, data.frame(Set = "Schwartz_PRE", state = Schwartz_PRE_dm6$Cluster))

RadarPlotSummary = RadarPlot %>%
  group_by(Set, state) %>%
  summarise(count = n()) %>%
  group_by(Set) %>%
  mutate(fraction = count/sum(count)) %>%
  mutate(state = c("EnhW", "Pc-I", "TEl", "EnhS", "Pc-H", "TSS", "Het")[as.numeric(state)]) %>%
  pivot_wider(state, values_from = fraction, names_from = Set, values_fill = 0)


missing.state = data.frame(state = c("TEl", "Het"),
                           Ederle_PRE_HC = c(0,0),
                           Kahn_PRE = c(0,0),
                           Schwartz_PRE = c(0,0))


RadarPlotSummary = rbind(RadarPlotSummary, missing.state)

RadarPlotSummaryTrans = t(RadarPlotSummary[,2:4])
colnames(RadarPlotSummaryTrans) = RadarPlotSummary$state

RadarPlotSummaryTrans = as.data.frame(rbind(matrix(1, nrow = 1, ncol = 7),
                              matrix(0, nrow = 1, ncol = 7),
                              RadarPlotSummaryTrans))


colors_border=rcartocolor::carto_pal(n = 4, "Bold")[1:3]
colors_in=rcartocolor::carto_pal(n = 4, "Bold")[1:3]


# Prepare title
mytitle <- c("Ederle et al.\n(high confidence)", "Kahn et al.", "Schwartz et al.")

plot.new()
# Split the screen in 3 parts
par(mar=rep(0.8,4))
par(mfrow=c(1,3))


## Fig1PREOverlabRadar.pdf

# Loop for each plot
for(i in 1:3){

  # Custom the radarChart !
  radarchart(RadarPlotSummaryTrans[c(1,2,i+2),], axistype=1, 
              
              #custom polygon
              pcol=paste0(colors_in, "D0")[i] , pfcol=paste0(colors_in, "80")[i] , plwd=4, plty=1 , 
              
              #custom the grid
              cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8,
              
              #custom labels
              vlcex=0.8
              
  )
  title(main = mytitle[i], adj=0.5, line = -1.5)
  text(x = 0, y = 0, c(157,200,170)[i])
}


