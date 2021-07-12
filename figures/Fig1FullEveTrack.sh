python3.8 ../data/bws/SparK.py \
-pr chr2R:9958921-10033540 \
-dg Pal1 CG12133 eIF3j eve TER94 Pka-R2 \
-cf ../data/bws/RestingHemocyteK3772H3K27me3Mean.bdg ../data/bws/RestingHemocyteCH075H3K27me23Mean.bdg ../data/bws/RestingHemocyteCH105H3K9me3Mean.bdg \
../data/bws/RestingHemocyteCH075H3K4me1Mean.bdg ../data/bws/RestingHemocyteCH098H3K27acMean.bdg ../data/bws/RestingHemocyteCH105H3K9acMean.bdg \
../data/bws/RestingHemocyteCH075H3K4me3Mean.bdg ../data/bws/RestingHemocyteCH075H3K36me3Mean.bdg ../data/bws/RestingHemocyteCH105PolIIMean.bdg \
../data/bws/RestingHemocyteCH075H4K20me1Mean.bdg ../data/bws/RestingHemocyteCH075InputMean.bdg ../data/bws/HemocyteMerged2426Transcriptome.bdg \
-gtf /Users/streeck/Genomes/DmelBDGP6.91/Drosophila_melanogaster.BDGP6.91.gtf \
-gl H3K27me3 H3K27me23 H3K9me3 H3K4me1 H3K27ac H3K9ac H3K4me3 H3K36me3 RNA-PolII H4K20me1 Input polyA-RNA \
-bed ../data/fig1/MACS2Test/H3K27me3_infected_peaks_subset.broadPeak \
../data/fig1/MACS2Test/H3K27me3_infected_peaks.broadPeak \
-bedcol EF1414 FFBC00 \
-bedlab MACS2_high MACS2_all \
-o Fig1FullEveTrack
