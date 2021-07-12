library(rtracklayer)
library(edgeR)
library(knitr)
library(tidyverse)
library(kableExtra)


load("data/fig1/resting_H3K27me3.Rdata")
load("data/fig1/SevenClassGenomeModel.Rdata")

GFF = readGFF("/Users/streeck/Genomes/DmelBDGP6.91/Drosophila_melanogaster.BDGP6.91.gtf", filter = list(type = "gene"))


tpm = function(x, len){
  x = (x+1)/len
  1e6*x/sum(x)
}

resting_RNAseq = read.delim("data/fig1/2426.count", sep = "\t", skip = 1)
colnames(resting_RNAseq)[7:9] = c("RepA", "RepB", "RepC")
resting_RNAseq = resting_RNAseq %>%
  dplyr::select(Geneid, Length:RepC) %>%
  dplyr::mutate(RepA_tpm = tpm(RepA, Length),
         RepB_tpm = tpm(RepB, Length),
         RepC_tpm = tpm(RepC, Length)) %>%
  dplyr::mutate(mean_tpm = (RepA_tpm + RepB_tpm + RepC_tpm)/3) 


# load RNA-seq data -------------------------------------------------------


## switch

CountTable_switch = read.delim("data/fig2/unstranded.count", skip = 1)
CountTable_switch = CountTable_switch[,c(1,7:dim(CountTable_switch)[2])]

colnames(CountTable_switch)[2:dim(CountTable_switch)[2]] = sapply(
  strsplit(colnames(CountTable_switch)[2:dim(CountTable_switch)[2]], ".", fixed = T), 
  function(x) paste0(x[length(x)-2], "_", x[length(x)-1]))

CountTable_switch[,dim(CountTable_switch)[2]+1] = read.delim("data/fig2/unstrandedfeatureCounts.count", skip = 1)[7]
colnames(CountTable_switch)[dim(CountTable_switch)[2]] = "NG_26757_24"


SampleDesctiptors_switch = read.xlsx("data/fig2/2021-04-11MappedSamples.xlsx")
SampleDesctiptors_switch$sample_num = sapply(strsplit(SampleDesctiptors_switch$sample_id, "_", fixed = T), function(x) x[2])
SampleDesctiptors_switch$sample_id = gsub("-", "_", SampleDesctiptors_switch$sample_id, fixed = T)

## exclude samples 1 and 8 (look bad on PCA)
CountTable_switch = CountTable_switch[, !colnames(CountTable_switch) %in% c("NG_26757_1", "NG_26757_8")]
SampleDesctiptors_switch = SampleDesctiptors_switch %>%
  filter(!sample_num %in% c("1", "8"))

groups_switch = SampleDesctiptors_switch$condition
designMatrix_switch = model.matrix(~ 0 + groups_switch)
colnames(designMatrix_switch) = substr(colnames(designMatrix_switch), 14, 100)
colnames(designMatrix_switch) = gsub(" ", "_", colnames(designMatrix_switch), fixed = T)

DGEobject_switch = DGEList(CountTable_switch[,2:dim(CountTable_switch)[2]], group = groups_switch, genes = CountTable_switch[,1])
DGEobject_switch = DGEobject_switch[rowSums(cpm(DGEobject_switch)>1) >= max(table(groups_switch)), , keep.lib.sizes=FALSE]
DGEobject_switch$genes = GFF[match(DGEobject_switch$genes$genes, GFF$gene_id), c(9,10,12)]
DGEobject_switch = calcNormFactors(DGEobject_switch)
DGEobject_switch = estimateDisp(DGEobject_switch, design = designMatrix_switch, robust = T)

contrasts_switch = makeContrasts(eGFP_induction = EGFP_Ind - EGFP_control,
                                 hop_induction = hop_Ind - hop_control,
                                 PGRP_induction = PGRP_Ind - PGRP_control,
                                 levels = designMatrix_switch)

fit_switch = glmQLFit(DGEobject_switch, designMatrix_switch, robust = T)
QFresults_switch = list()
for (i in 1:dim(contrasts_switch)[2]) {
  QFresults_switch[[colnames(contrasts_switch)[i]]] = glmQLFTest(fit_switch, contrast = contrasts_switch[,i])
}

## infection


CountTableInfection = read.delim("data/fig2/OrRInfection.count", skip = 1)
CountTableInfection = CountTableInfection[,c(1,7:dim(CountTableInfection)[2])]

SampleDesctiptorsInfection = read.xlsx("data/fig2/OrRInfectionLibraries.xlsx")
colnames(CountTableInfection)[2:21] = sapply(strsplit(colnames(CountTableInfection)[2:21], "_", fixed = T), function(x){paste0("Lib_", x[2])})
SampleDesctiptorsInfection$bam = factor(SampleDesctiptorsInfection$bam, levels = colnames(CountTableInfection)[2:21])
SampleDesctiptorsInfection = SampleDesctiptorsInfection[order(SampleDesctiptorsInfection$bam),]

groupsInfection = SampleDesctiptorsInfection$condition
designMatrixInfection = model.matrix(~ 0 + groupsInfection)
colnames(designMatrixInfection) = substr(colnames(designMatrixInfection), 16, 100)
colnames(designMatrixInfection) = gsub(" ", "_", colnames(designMatrixInfection), fixed = T)

DGEobjectInfection = DGEList(CountTableInfection[,2:dim(CountTableInfection)[2]], group = groupsInfection, genes = CountTableInfection[,1])
DGEobjectInfection = DGEobjectInfection[rowSums(cpm(DGEobjectInfection)>1) >= max(table(groupsInfection)), , keep.lib.sizes=FALSE]
DGEobjectInfection$genes = GFF[match(DGEobjectInfection$genes$genes, GFF$gene_id), c(9,10,12)]
DGEobjectInfection = calcNormFactors(DGEobjectInfection)
DGEobjectInfection = estimateDisp(DGEobjectInfection, design = designMatrixInfection, robust = T)



## contrasts and fits
contrastsInfection = makeContrasts(de_3h = septic_3h - control,
                                   de_6h = septic_6h - control,
                                   de_18h = septic_18h - control,
                                   levels = designMatrixInfection)

fitInfection = glmQLFit(DGEobjectInfection, designMatrixInfection, robust = T)
QFresultsInfection = list()
for (i in 1:dim(contrastsInfection)[2]) {
  QFresultsInfection[[colnames(contrastsInfection)[i]]] = glmQLFTest(fitInfection, contrast = contrastsInfection[,i])
}


## RNAi import 

CountTableRNAi = read.delim("data/fig3/2021-03-20featureCounts.count", skip = 1)
CountTableRNAi = CountTableRNAi[,c(1,7:dim(CountTableRNAi)[2])]

colnames(CountTableRNAi)[2:dim(CountTableRNAi)[2]] = sapply(strsplit(colnames(CountTableRNAi)[2:dim(CountTableRNAi)[2]], ".", fixed = T), function(x) x[length(x)-1])


SampleDesctiptorsRNAi = read.xlsx("data/fig3/2021-03-21MappedSamples.xlsx")

groupsRNAi = SampleDesctiptorsRNAi$condition
designMatrixRNAi = model.matrix(~ 0 + groupsRNAi)
colnames(designMatrixRNAi) = substr(colnames(designMatrixRNAi), 11, 100)
colnames(designMatrixRNAi) = gsub(" ", "_", colnames(designMatrixRNAi), fixed = T)

DGEobjectRNAi = DGEList(CountTableRNAi[,2:dim(CountTableRNAi)[2]], group = groupsRNAi, genes = CountTableRNAi[,1])
DGEobjectRNAi = DGEobjectRNAi[rowSums(cpm(DGEobjectRNAi)>1) >= max(table(groupsRNAi)), , keep.lib.sizes=FALSE]
DGEobjectRNAi$genes = GFF[match(DGEobjectRNAi$genes$genes, GFF$gene_id), c(9,10,12)]
DGEobjectRNAi = calcNormFactors(DGEobjectRNAi)
DGEobjectRNAi = estimateDisp(DGEobjectRNAi, design = designMatrixRNAi, robust = T)

## contrasts and fits
contrastsRNAi = makeContrasts(luc = luc_septic - luc_untreated,
                              nej = nej_septic - nej_untreated,
                              Set1 = Set1_septic - Set1_untreated,
                              trr = trr_septic - trr_untreated,
                              trx = trx_septic - trx_untreated,
                              Utx = Utx_septic - Utx_untreated,
                              nej_baseline = nej_untreated - luc_untreated,
                              Set1_baseline =  Set1_untreated - luc_untreated,
                              trr_baseline = trr_untreated - luc_untreated,
                              trx_baseline = trx_untreated - luc_untreated,
                              Utx_baseline = Utx_untreated - luc_untreated,
                              nej_septic = nej_septic - luc_septic,
                              Set1_septic =  Set1_septic - luc_septic,
                              trr_septic = trr_septic - luc_septic,
                              trx_septic = trx_septic - luc_septic,
                              Utx_septic = Utx_septic - luc_septic,
                              nej_full = nej_septic - nej_untreated - (luc_septic - luc_untreated),
                              Set1_full = Set1_septic - Set1_untreated- (luc_septic - luc_untreated),
                              trr_full = trr_septic - trr_untreated- (luc_septic - luc_untreated),
                              trx_full = trx_septic - trx_untreated- (luc_septic - luc_untreated),
                              Utx_full = Utx_septic - Utx_untreated- (luc_septic - luc_untreated),
                              levels = designMatrixRNAi)

fitRNAi = glmQLFit(DGEobjectRNAi, designMatrixRNAi, robust = T)
QFresultsRNAi = list()
for (i in 1:dim(contrastsRNAi)[2]) {
  QFresultsRNAi[[colnames(contrastsRNAi)[i]]] = glmQLFTest(fitRNAi, contrast = contrastsRNAi[,i])
}

## Mosaics

CountTable2537 = read.delim("data/fig3/2537.count", skip = 1)
CountTable2537 = CountTable2537[,c(1,7:dim(CountTable2537)[2])]
CountTable2537 = CountTable2537[,c(T, grepl("NoDup", colnames(CountTable2537)[2:dim(CountTable2537)[2]]))]
colnames(CountTable2537)[2:dim(CountTable2537)[2]] = sapply(
  strsplit(colnames(CountTable2537)[2:dim(CountTable2537)[2]], "_", fixed = T), 
  function(x){paste0("Lib_2537_", x[2])})

CountTable3506 = read.delim("data/fig3/3506.count", skip = 1)
CountTable3506 = CountTable3506[,c(1,7:dim(CountTable3506)[2])]
colnames(CountTable3506)[2:dim(CountTable3506)[2]] = sapply(
  strsplit(colnames(CountTable3506)[2:dim(CountTable3506)[2]], "_", fixed = T), 
  function(x){paste0("Lib_3506_", x[2])})

CountTable = CountTable2537 %>%
  left_join(CountTable3506, by = "Geneid")

SampleDesctiptors = read.xlsx("data/fig3/MosaicSamples.xlsx") %>%
  mutate(bam = paste0("Lib_", Project, "_", Lib_id))
SampleDesctiptors$bam = factor(SampleDesctiptors$bam, levels = colnames(CountTable)[2:dim(CountTable)[2]])
SampleDesctiptors = SampleDesctiptors[order(SampleDesctiptors$bam),]

groups = SampleDesctiptors$genotype
designMatrix = model.matrix(~ 0 + groups)
colnames(designMatrix) = substr(colnames(designMatrix), 7, 100)
colnames(designMatrix) = gsub(" ", "_", colnames(designMatrix), fixed = T)

DGEobject = DGEList(CountTable[,2:dim(CountTable)[2]], group = groups, genes = CountTable[,1])
DGEobject = DGEobject[rowSums(cpm(DGEobject)>1) >= max(table(groups)), , keep.lib.sizes=FALSE]
DGEobject$genes = GFF[match(DGEobject$genes$genes, GFF$gene_id), c(9,10,12)]
DGEobject = calcNormFactors(DGEobject)
DGEobject = estimateDisp(DGEobject, design = designMatrix, robust = T)


## contrasts and fits
contrasts = makeContrasts(H3.K27R = H3.K27R - H3.wt,
                          Suz12 = Suz12.4 - Suz12.wt,
                          Ez = Ez.731 - Ez.wt,
                          levels = designMatrix)

fit = glmQLFit(DGEobject, designMatrix, robust = T)
QFresults = list()
for (i in 1:dim(contrasts)[2]) {
  QFresults[[colnames(contrasts)[i]]] = glmQLFTest(fit, contrast = contrasts[,i])
}



# combine data ------------------------------------------------------------

leftjoin.named.list = function(.data, qf.list){
  for (i in 1:length(qf.list)) {
    qf = qf.list[[i]]
    qf = data.frame(gene_id = qf$genes[,1], 
                    qf$table[,c(1,4)],
                    de = as.numeric(decideTests.DGELRT(qf)))
    colnames(qf)[2:4] = paste0(names(qf.list)[i], ".", colnames(qf)[2:4])
    .data = .data %>%
      left_join(qf, by = "gene_id")
  }
  return(.data)
}

leftjoin.cpm = function(.data, DGE){
  cpm.df = data.frame(gene_id = DGE$genes$gene_id,
                      cpmByGroup(DGE))
  colnames(cpm.df)[2:dim(cpm.df)[2]] = paste0("cpm.", colnames(cpm.df)[2:dim(cpm.df)[2]])
  return(.data %>%
    left_join(cpm.df, by = "gene_id"))
}


LargeOverviewTable = GFF %>%
  dplyr::select(seqid, strand, gene_id, gene_name, gene_biotype) %>%
  left_join(resting_RNAseq %>%
              dplyr::select(Geneid, mean_tpm),
            by = c("gene_id" = "Geneid")) %>%
  left_join(resting_H3K27me3 %>%
              mutate(H3K27me3.log_ratio = log2((H3K27me3_A + H3K27me3_B)/(input_A + input_B))) %>%
              dplyr::select(Geneid, class, H3K27me3.log_ratio),
            by = c("gene_id" = "Geneid")) %>%
  left_join(GeneTable %>%
              mutate(chrom_state = c("EnhW", "Pc-I", "TEl", "EnhS", "Pc-H", "TSS", "Het")[as.numeric(MajorityStateVote)]) %>%
              dplyr::select(gene_id, chrom_state),
            by = "gene_id") %>%
  leftjoin.cpm(DGEobjectInfection) %>%
  leftjoin.named.list(QFresultsInfection) %>%
  leftjoin.cpm(DGEobject_switch) %>%
  leftjoin.named.list(QFresults_switch) %>%
  leftjoin.cpm(DGEobject) %>%
  leftjoin.named.list(QFresults) %>%
  leftjoin.cpm(DGEobjectRNAi) %>%
  leftjoin.named.list(QFresultsRNAi)

save(LargeOverviewTable, file = "2021-06-08-LargeOverwieTableForAlf.Rdata")



