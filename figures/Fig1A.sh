python3.8 ../data/bws/SparK.py \
-pr chr2R:9958921-10033540 \
-dg Pal1 CG12133 eIF3j eve TER94 Pka-R2 \
-cf ../data/bws/RestingHemocyteK3772H3K27me3Mean.bdg ../data/bws/RestingHemocyteCH105PolIIMean.bdg ../data/bws/HemocyteMerged2426Transcriptome.bdg \
-gtf /Users/streeck/Genomes/DmelBDGP6.91/Drosophila_melanogaster.BDGP6.91.gtf \
-gl H3K27me3 RNA-PolII polyA-RNA \
-bed ../data/fig1/MACS2Test/H3K27me3_infected_peaks_subset.broadPeak \
../data/fig1/MACS2Test/H3K27me3_infected_peaks.broadPeak \
-bedcol EF1414 FFBC00 \
-bedlab MACS2_high MACS2_all \
-o Fig1A
