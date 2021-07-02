library(rtracklayer)
library(Rfast)

quantNorm = function(X, round = F){
  require(Rfast)
  X = as.matrix(X)
  ndim = dim(X)[2]
  OrderMat = matrix(nrow = dim(X)[1], ncol = ndim)
  for(i in 1:ndim){
    OrderMat[,i] = base::order(X[,i])
    X[,i] = X[OrderMat[,i],i]
  }
  if(round == F){
    RowMeans = Rfast::rowsums(X)/ndim
  }else{
    RowMeans = round(Rfast::rowsums(X)/ndim)
  }
  for(i in 1:ndim){
    X[OrderMat[,i],i] = RowMeans
  }
  X
}


genome = read.delim("/Users/streeck/Genomes/DmelBDGP6.91/chrNameLength.txt", header = F, stringsAsFactors = F)
genome = genome[1:7,]
gr = GRanges(genome[,1], IRanges(1, as.integer(genome[,2])))
TiledGenome = unlist(tile(gr, width = 50))

fileList = list(InfInputA = "data/fig2/CoverageBWs/4054_A_1_run573_ATCACGAT_S44_L008_NoDup.bw",
                InfInputB = "data/fig2/CoverageBWs/4054_A_2_run573_CGATGTAT_S45_L008_NoDup.bw",
                InfK27acA = "data/fig2/CoverageBWs/4054_A_3_run573_TTAGGCAT_S46_L008_NoDup.bw",
                InfK27acB = "data/fig2/CoverageBWs/4054_A_4_run573_TGACCAAT_S47_L008_NoDup.bw",
                InfK27me3A = "data/fig2/CoverageBWs/4054_A_5_run573_ACAGTGAT_S48_L008_NoDup.bw",
                InfK27me3B = "data/fig2/CoverageBWs/4054_A_6_run573_CAGATCAT_S49_L008_NoDup.bw",
                ContAcInputA = "data/fig2/CoverageBWs/mpimg_L10591_CH_098-802_S5_NoDup.bw",
                ContAcInputB = "data/fig2/CoverageBWs/mpimg_L10591_CH_098-804_S7_NoDup.bw",
                ContK27acA = "data/fig2/CoverageBWs/mpimg_L10591_CH_098-803_S6_NoDup.bw",
                ContK27acB = "data/fig2/CoverageBWs/mpimg_L10591_CH_098-805_S8_NoDup.bw",
                ContMeInputA = "data/fig2/CoverageBWs/K002000265_94327__NoDup.bw",
                ContMeInputB = "data/fig2/CoverageBWs/K002000265_94329__NoDup.bw",
                ContMeK27me3A = "data/fig2/CoverageBWs/K002000265_94326__NoDup.bw",
                ContMeK27me3B = "data/fig2/CoverageBWs/K002000265_94328__NoDup.bw")

BWlist = list()
for (i in 1:14) {
  BWlist[[i]] = import(fileList[[i]])
}


BWframe = matrix(ncol = 14, nrow = length(TiledGenome))
for (i in 1:14) {
  BWframe[,i] = BWlist[[i]]$score[findOverlaps(TiledGenome, BWlist[[i]], minoverlap = 25, select = "first")]
}


### log2Ratio Norm

S = BWframe[,c(13,14,9,10,5,6,3,4)] + 1
R = BWframe[,c(11,12,7,8,1,2,1,2)] + 1
S = S/matrix(colsums(S), nrow = dim(S)[1], ncol = dim(S)[2], byrow = T)
R = R/matrix(colsums(R), nrow = dim(S)[1], ncol = dim(S)[2], byrow = T)
E = S/R
E = log2(E)
E[is.nan(E)] = 0

library(GGally)

lower_ggpairs = function(data, mapping){
  ggplot(data, mapping) + 
    geom_point(alpha = .2, shape = 16) +
    geom_density_2d(bins = 8, binwidth = 1, color = "#90728f") +
    geom_smooth(method = "lm", color = "#9d983dbb", formula = y ~ x)
}


### FigSubBWPreQuantileNormH3K27me3.pdf

E %>%
  .[,c(1,2,5,6)] %>%
  .[sample(1:dim(E)[1], 5000),] %>%
  as.data.frame() %>%
  dplyr::rename(RestingRepA = "V1", RestingRepB = "V2", SepticgRepA = "V3", SepticRepB = "V4") %>%
  ggpairs(lower = list(continuous = lower_ggpairs), title = "H3K27me3 before quantile normalization",
          xlab = "log2 ratio", ylab = "log2 ratio")

### FigSubBWPreQuantileNormH3K27ac.pdf

E %>%
  .[,c(3,4,7,8)] %>%
  .[sample(1:dim(E)[1], 5000),] %>%
  as.data.frame() %>%
  dplyr::rename(RestingRepA = "V1", RestingRepB = "V2", SepticgRepA = "V3", SepticRepB = "V4") %>%
  ggpairs(lower = list(continuous = lower_ggpairs), title = "H3K27ac before quantile normalization",
          xlab = "log2 ratio", ylab = "log2 ratio")

E.norm = E
E.norm[,c(1,2,5,6)] = quantNorm(E[,c(1,2,5,6)])
E.norm[,c(3,4,7,8)] = quantNorm(E[,c(3,4,7,8)])
E.norm[is.nan(E.norm)] = 0
E.norm[is.na(E.norm)] = 0

### FigSubBWPostQuantileNormH3K27me3.pdf

E.norm %>%
  .[,c(1,2,5,6)] %>%
  .[sample(1:dim(E)[1], 5000),] %>%
  as.data.frame() %>%
  dplyr::rename(RestingRepA = "V1", RestingRepB = "V2", SepticgRepA = "V3", SepticRepB = "V4") %>%
  ggpairs(lower = list(continuous = lower_ggpairs), title = "H3K27me3 after quantile normalization",
          xlab = "log2 ratio", ylab = "log2 ratio")

### FigSubBWPostQuantileNormH3K27ac.pdf

E.norm %>%
  .[,c(3,4,7,8)] %>%
  .[sample(1:dim(E)[1], 5000),] %>%
  as.data.frame() %>%
  dplyr::rename(RestingRepA = "V1", RestingRepB = "V2", SepticgRepA = "V3", SepticRepB = "V4") %>%
  ggpairs(lower = list(continuous = lower_ggpairs), title = "H3K27ac after quantile normalization",
          xlab = "log2 ratio", ylab = "log2 ratio")


ex = TiledGenome
ex$score = rowsums(E[,1:2])/2
ex$score = ex$score - quantile(ex$score, probs = c(.1))[1]
ex@seqinfo@seqlengths = genome$V2
export.bw(ex, "data/fig2/nomrbw/RestingH3K27me3NormalizedLogRatio.bw")
ex$score = rowsums(E[,3:4])/2
ex$score = ex$score - quantile(ex$score, probs = c(.1))[1]
export.bw(ex, "data/fig2/nomrbw/RestingH3K27acNormalizedLogRatio.bw")
ex$score = rowsums(E[,5:6])/2
ex$score = ex$score - quantile(ex$score, probs = c(.1))[1]
export.bw(ex, "data/fig2/nomrbw/InfectedH3K27me3NormalizedLogRatio.bw")
ex$score = rowsums(E[,7:8])/2
ex$score = ex$score - quantile(ex$score, probs = c(.1))[1]
export.bw(ex, "data/fig2/nomrbw/InfectedH3K27acNormalizedLogRatio.bw")
ex$score = rowsums(E[,5:6]) - rowsums(E[,1:2])
export.bw(ex, "data/fig2/nomrbw/DiffH3K27me3NormalizedLogRatio.bw")
ex$score = rowsums(E[,7:8]) - rowsums(E[,3:4])
export.bw(ex, "data/fig2/nomrbw/DiffH3K27acNormalizedLogRatio.bw")

ex$score = E[,1]
ex$score = ex$score - quantile(ex$score, probs = c(.1))[1]
export.bw(ex, "data/fig2/nomrbw/RestingH3K27me3NormalizedLogRatioArep.bw")
ex$score = E[,3]
ex$score = ex$score - quantile(ex$score, probs = c(.1))[1]
export.bw(ex, "data/fig2/nomrbw/RestingH3K27acNormalizedLogRatioArep.bw")
ex$score = E[,5]
ex$score = ex$score - quantile(ex$score, probs = c(.1))[1]
export.bw(ex, "data/fig2/nomrbw/InfectedH3K27me3NormalizedLogRatioArep.bw")
ex$score = E[,7]
ex$score = ex$score - quantile(ex$score, probs = c(.1))[1]
export.bw(ex, "data/fig2/nomrbw/InfectedH3K27acNormalizedLogRatioArep.bw")
ex$score = E[,2]
ex$score = ex$score - quantile(ex$score, probs = c(.1))[1]
export.bw(ex, "data/fig2/nomrbw/RestingH3K27me3NormalizedLogRatioBrep.bw")
ex$score = E[,4]
ex$score = ex$score - quantile(ex$score, probs = c(.1))[1]
export.bw(ex, "data/fig2/nomrbw/RestingH3K27acNormalizedLogRatioBrep.bw")
ex$score = E[,6]
ex$score = ex$score - quantile(ex$score, probs = c(.1))[1]
export.bw(ex, "data/fig2/nomrbw/InfectedH3K27me3NormalizedLogRatioBrep.bw")
ex$score = E[,8]
ex$score = ex$score - quantile(ex$score, probs = c(.1))[1]
export.bw(ex, "data/fig2/nomrbw/InfectedH3K27acNormalizedLogRatioBrep.bw")





