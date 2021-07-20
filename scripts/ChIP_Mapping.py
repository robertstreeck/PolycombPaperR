import os
############ Define a lot of paths
# directory in which the (STAR indexed) genome is
OSGenDir = "/Users/streeck/Genomes/DmelBDGP6.91/"
# Which genome to use in IGV foldes
IGVToolsGenome = "dm6"
# What directory the fastq.gz files are in
SourceDir = "/Users/streeck/Desktop/Seq4054/ChIP/ChIPfastq"
TargetDir = "/Users/streeck/Desktop/Seq4054/ChIP"
# Picard Command
Pic = "java -Xmx32G -jar /Users/streeck/picard-2.8.2/picard-2.8.2.jar"
#Fasta location
Fasta =  "/Users/streeck/Genomes/DmelBDGP6.91/Drosophila_melanogaster.BDGP6.dna.toplevel.fa"


######### Start scriptfile with general info
# opening target files for mapping and cleanup scripts
f = open("MappingScript.sh","w")
# modifying PATH to include necessary tools (might need adjustment)
f.write("export PATH=$PATH:/Users/streeck/IGVTools:/Users/streeck/UCSCTools:/Users/streeck/STAR-2.7.0e/bin/MacOSX_x86_64\n")

# empty Lists
SampleList = list()
# For files in SourceDir, take all files with .fastq.gz, remember files and sample type (e.g. library)
for file in os.listdir(SourceDir):
    if file.endswith(".fastq.gz"):
        SampleList.append(file[0:-16])
SampleSet = list(set(SampleList))

f.write("mkdir " + TargetDir + "/DupBAM\n")
f.write("mkdir " + TargetDir + "/BAM\n")
f.write("mkdir " + TargetDir + "/BW\n")
f.write("mkdir " + TargetDir + "/QCData\n")

# Map all files with Star and remove merged Fastqs
for sample in SampleSet:
    f.write("STAR --genomeDir " +OSGenDir + " --readFilesIn " + SourceDir + "/" + sample +  "_R1_001.fastq.gz " + SourceDir + "/" + sample +  "_R2_001.fastq.gz" + " --readFilesCommand gunzip -c --outSAMtype BAM Unsorted --clip3pNbases 30 --alignIntronMax 1 --outFileNamePrefix " + TargetDir + "/DupBAM/" + sample + " --runThreadN 6 \n")
    f.write("samtools sort -o" + TargetDir + "/DupBAM/" + sample + "sorted.bam " + TargetDir + "/DupBAM/" + sample + "Aligned.out.bam\n")
    f.write("samtools index " + TargetDir + "/DupBAM/" + sample +  "sorted.bam\n")
    f.write(Pic + " CollectAlignmentSummaryMetrics I=" + TargetDir + "/DupBAM/" + sample + "sorted.bam O=" + TargetDir + "/QCData/" + sample + "PreDubRemovalAlignmentSummaryMetrics.txt R=" + Fasta + " & ")
    f.write(Pic + " CollectInsertSizeMetrics I=" + TargetDir + "/DupBAM/" + sample + "sorted.bam O=" + TargetDir + "/QCData/" + sample + "PreDubRemovalInsertSizeMetrics.txt H=" + TargetDir + "/QCData/" + sample + "PreDubRemovalHist.pdf & ")
    f.write(Pic + " MarkDuplicates INPUT=" + TargetDir + "/DupBAM/" + sample + "sorted.bam OUTPUT=" + TargetDir + "/BAM/" + sample + "_NoDup.bam METRICS_FILE=" + TargetDir + "/QCData/" + sample + "MarkDuplicates_metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true TMP_DIR=/tmp REMOVE_DUPLICATES=true ASSUME_SORTED=true\n")
    f.write(Pic + " CollectAlignmentSummaryMetrics I=" + TargetDir + "/BAM/" + sample + "_NoDup.bam O=" + TargetDir + "/QCData/" + sample + "NoDupAlignmentSummaryMetrics.txt R=" + Fasta + " & ")
    f.write("igvtools count --minMapQuality 30 " + TargetDir + "/BAM/" + sample + "_NoDup.bam stdout " + IGVToolsGenome + " | wigToBigWig -clip stdin /Users/streeck/IGVTools/genomes/" + IGVToolsGenome + ".chrom.sizes " + TargetDir + "/BW/" + sample + "_NoDup.bw \n")

f.close()

