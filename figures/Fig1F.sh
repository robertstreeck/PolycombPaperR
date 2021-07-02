#computeMatrix scale-regions -m 5000 -a 3000 -b 3000 -bs 50 --skipZeros --missingDataAsZero --metagene -p 18 -S RestingHemocyteCH105H3K9me3L2Ration.bw RestingHemocyteCH075H3K27me23L2Ration.bw RestingHemocyteK3772H3K27me3L2Ration.bw RestingHemocyteCH105PolIIL2Ration.bw RestingHemocyteCH075H3K36me3L2Ration.bw -R K27C1_Table.gtf K27C2_Table.gtf K27C3_Table.gtf -o GeneBody_K27C_ScaledRegions.gz;

Burg="#ffc6c4,#ee919b,#cc607d,#9e3963,#672044"
Peach="#fde0c5,#f9c098,#f59e72,#f17854,#eb4a40"
Teal="#d1eeea,#96d0d1,#68abb8,#45829b,#2a5674"



plotHeatmap -m GeneBody_K27C_ScaledRegions.gz -o FigS3_GeneBody_ScaledRegions_hm.pdf \
    --regionsLabel non-Pc Pc-I Pc-H --samplesLabel H3K9me3 H3K27me23 H3K27me3 PolII H3K36me3 \
    --colorList $Peach $Peach $Burg $Teal $Teal --zMin -1 -4 -4 0 -3 --zMax 4.5 0.5 3 4.5 2 \
    --whatToShow 'heatmap and colorbar' --heatmapWidth 5; 
plotProfile -m GeneBody_K27C_ScaledRegions.gz -o Fig1F_GeneBody_ScaledRegions_profile.pdf \
    --regionsLabel non-Pc Pc-I Pc-H --samplesLabel H3K9me3 H3K27me23 H3K27me3 PolII H3K36me3 \
    --colors "#88CCEE" "#CC6677" "#DDCC77";

#computeMatrix reference-point -a 3000 -b 3000 -bs 50 --skipZeros --missingDataAsZero --metagene -p 18 -S RestingHemocyteK3772H3K27me3L2Ration.bw RestingHemocyteCH075H3K4me1L2Ration.bw RestingHemocyteCH075H3K4me3L2Ration.bw RestingHemocyteCH105PolIIL2Ration.bw -R K27C1_Table.gtf K27C2_Table.gtf K27C3_Table.gtf -o critical_K27C_ReferencePoint.gz;


plotHeatmap -m critical_K27C_ReferencePoint.gz -o FigS3_ReferencePoint_hm.pdf \
    --regionsLabel non-Pc Pc-I Pc-H --samplesLabel H3K27me3 H3K4me1 H3K4me3 PolII \
    --colorList $Burg $Teal $Teal $Teal --zMin -4 -4 -4.2 0.5 --zMax 3.5 0 1 4.5 \
    --whatToShow 'heatmap and colorbar' --heatmapWidth 5;
plotProfile -m critical_K27C_ReferencePoint.gz -o Fig1F_ReferencePoint_profile.pdf \
    --regionsLabel non-Pc Pc-I Pc-H --samplesLabel H3K27me3 H3K4me1 H3K4me3 PolII \
    --colors "#88CCEE" "#CC6677" "#DDCC77";
