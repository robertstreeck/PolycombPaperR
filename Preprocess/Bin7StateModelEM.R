library(openxlsx)
library(GenomicRanges)
library(rtracklayer)

chip_files = read.xlsx("Data/fig1/SampleList.xlsx")
chip_file_selection = chip_files %>% select(Index:Input.ref, Replicate) %>%
  filter(Antibody != "Input") %>% 
  left_join(chip_files %>% select(Index, Bam.File), by = c("Input.ref" = "Index"),
            suffix = c(".chip", ".input"))


genome = read.delim("../../Genomes/DmelBDGP6.91/chrNameLength.txt", header = F, stringsAsFactors = F)
genome = genome[1:7,]
gr = GRanges(genome[,1], IRanges(1, as.integer(genome[,2])))


source("EM/EMdataIngestWrapper.R")

multi_chip_fit = BinEMClusterFromBAM(ChIPBamlist = chip_file_selection$Bam.File.chip, 
                                     InputBamList = chip_file_selection$Bam.File.input, 
                                     gr = gr, 
                                     path = "/Volumes/Promise_Leia/Thesis/ChIPReanalysis/ChIP_BAMs/", 
                                     target = "plasmatocytes", nmodel = 7, binsize = 200, ncores = 2)


chr_lengths = gr@ranges@width
names(chr_lengths) = gr@seqnames

tiled_genome = tileGenome(chr_lengths, tilewidth = 200, cut.last.tile.in.chrom = T)
tiled_genome$cluster = NA
tiled_genome$cluster[multi_chip_fit$excluded] = multi_chip_fit$Group

GenomeGFF = import.gff("/Users/streeck/Genomes/DmelBDGP6.91/Drosophila_melanogaster.BDGP6.91.gtf")
GenomeGFF = GenomeGFF[GenomeGFF$type == "exon"]
GeneTable = readGFF("/Users/streeck/Genomes/DmelBDGP6.91/Drosophila_melanogaster.BDGP6.91.gtf", filter = list(type = "gene"))
GeneTable$MajorityStateVote = NA
GeneTable = GeneTable[GeneTable$seqid %in% c("3R","3L","2R","X","2L","Y","4"),]

for(i in 1:length(GeneTable$gene_id)){
  j = names(sort(-table(subsetByOverlaps(tiled_genome, GenomeGFF[GenomeGFF$gene_id == GeneTable$gene_id[i]], minoverlap = 1)$cluster)))[1]
  if(length(j) == 0){
    j = NA
  }
  GeneTable$MajorityStateVote[i] = j
}




save(multi_chip_fit, GeneTable, file = "SevenClassGenomeModel.Rdata")

