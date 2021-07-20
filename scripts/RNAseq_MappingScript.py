#!/usr/bin/env python

#%%
import os
import sys
import pandas as pd
import argparse
from multiprocessing import Pool
from datetime import datetime


#%%
if __name__ == '__main__':
    print("started mapping script")
    # setting up some variables
    ## paths to some tools that might not be in the env
    fastqc = "/Users/streeck/FastQC.app/Contents/MacOS/fastqc"
    STAR = "/Users/streeck/STAR-2.7.6a/bin/MacOSX_x86_64/STAR"
    Pic = "java -Xmx32G -jar /Users/streeck/picard-2.8.2/picard-2.8.2.jar"
    IGV = "/Users/streeck/IGVTools"
    UCSC = "/Users/streeck/UCSCTools"
    featureCounts = "/Users/streeck/subread-1.6.0-MacOSX-x86_64/bin/featureCounts"

    ## parse arguments for analysis
    parser = argparse.ArgumentParser(description = 'Script for mapping of RNA-seq data')
    parser.add_argument('--fastq', help = 'Source folder for fastq files', dest = 'fastq_directory', required = True)
    parser.add_argument('--xlsx', help = 'xlsx table of files and sample correspondance', dest = "xlsx", required = True)
    parser.add_argument('-o', '-out', help = 'directory to write the results to', dest = 'out', required = True)
    parser.add_argument('--not_paired', help = 'is the data paired end?', dest = 'paired', action = 'store_false')
    parser.add_argument('-n, --nCores', help = 'maximum number of cores to use', dest = 'nCores', default = 6)

    ## Add defaults for these
    parser.add_argument('--ref_path', help = 'path to the reference genome', dest = 'ref_genome', default = '/Users/streeck/Genomes/DmelBDGP6.91/')
    parser.add_argument('--genome', help = 'IGV tools genome name', dest = 'IGVToolsGenome', default = 'dm6')
    parser.add_argument('--gtf', help = 'GTF file path', dest = 'GTF', default = "/Users/streeck/Genomes/DmelBDGP6.91/Drosophila_melanogaster.BDGP6.91.gtf")
    parser.add_argument('--RSeQC_bed', help = 'Bed file to be used with RSeQC', dest = 'RefBed', default = "/Users/streeck/Genomes/DmelBDGP6.91/UCSC2R.bed")

    args = parser.parse_args(sys.argv[1:])
    out = args.out.rstrip("/")
    fastq_directory = args.fastq_directory.rstrip("/") + "/"
    paired = args.paired
    xlsx = args.xlsx
    nCores = args.nCores
    ref_genome = args.ref_genome
    IGVToolsGenome = args.IGVToolsGenome
    GTF = args.GTF
    RefBed = args.RefBed
    ChromSizes = "{}/genomes/{}.chrom.sizes".format(IGV, IGVToolsGenome)

    print("Starting mapping run on samples from {}".format(fastq_directory))

    ## generate some directories
    os.system("mkdir -p " + out + "/fastqc")
    os.system("mkdir -p " + out + "/BAM")
    os.system("mkdir -p " + out + "/temp")
    os.system("mkdir -p " + out + "/BW")
    os.system("mkdir -p " + out + "/QC")


    # generate a good list of samples
    ## Read xlsx to pd and make another frame with columns: sample_id, fwd_fastq, [rv_fastq], condition
    source_table = pd.read_excel(xlsx)

    ## check columns
    assert "fwd_fastq" in source_table, "missing column called fwd_fastq containing fastq files from fwd reads"
    if paired:
        assert "rv_fastq" in source_table, "missing column called rv_fastq containing fastq files from rv reads"
    assert "sample_id" in source_table, "missing column called sample_id containing ids for every library/RNA-sample sequenced"
    assert "condition" in source_table, "missing column called condition containing biological condition for that sample"
    print("found {} fastq (pairs) in {}".format(len(source_table.index), xlsx))

    ## generate a data frame grouped by samples to use for mapping
    select_columns = ["fwd_fastq", "rv_fastq", "sample_id", "condition"] if paired else ["fwd_fastq", "sample_id", "condition"]
    mapping_df = source_table[select_columns]
    mapping_df = mapping_df.groupby(by = ["sample_id", "condition"]).agg(lambda x: " ".join(x))
    mapping_df = mapping_df.reset_index()

    print("found {} samples in {}".format(len(mapping_df.index), xlsx))

    # Process to final BAMs and BWs using STAR, samtools, picardTools, USCS and IGV
    print("starting mapping")
    bamlist = []
    for row in mapping_df.iterrows():

        ## if paired cat fw and rv files and map using STAR
        if paired:
            os.system("cat {} > {}/temp/fwd_cat.fastq.gz" \
            .format(" ".join(fastq_directory + i for i in row[1]["fwd_fastq"].split()), out))
            os.system("cat {} > {}/temp/rv_cat.fastq.gz" \
            .format(" ".join(fastq_directory + i for i in row[1]["rv_fastq"].split()), out))
            os.system("{} --genomeDir {} --readFilesIn {}/temp/fwd_cat.fastq.gz {}/temp/rv_cat.fastq.gz" \
            " --readFilesCommand gunzip -c --outSAMtype BAM Unsorted --alignIntronMin 12"
            " --outFileNamePrefix {}/temp/{} --runThreadN {}"
            .format(STAR, ref_genome, out, out, out, row[1]["sample_id"].replace(" ", "_"), nCores))

        ## if not paired cat fw and map using STAR
        else:
            os.system("cat {} > {}/temp/fwd_cat.fastq.gz" \
            .format(" ".join(fastq_directory + i for i in row[1]["fwd_fastq"].split()), out))
            os.system("{} --genomeDir {} --readFilesIn {}/temp/fwd_cat.fastq.gz " \
            " --readFilesCommand gunzip -c --outSAMtype BAM Unsorted --alignIntronMin 12"
            " --outFileNamePrefix {}/temp/{} --runThreadN {}"
            .format(STAR, ref_genome, out, out, row[1]["sample_id"].replace(" ", "_"), nCores))

        ## sort and index bam using samtools
        os.system("samtools sort -o {}/temp/{}_Aligned.sortedByCoord.out.bam" \
        " -@ {} {}/temp/{}Aligned.out.bam"
        .format(out, row[1]["sample_id"].replace(" ", "_"), nCores, \
        out, row[1]["sample_id"].replace(" ", "_")))
        os.system("samtools index {}/temp/{}_Aligned.sortedByCoord.out.bam" \
        .format(out, row[1]["sample_id"].replace(" ", "_")))

        ## if paired remove duplicate reads using picardTools
        if paired:
            os.system("{} MarkDuplicates INPUT={}/temp/{}_Aligned.sortedByCoord.out.bam " \
            "OUTPUT={}/BAM/{}.bam METRICS_FILE={}/QC/{}_MarkDuplicates_metrics.txt" \
            " OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true TMP_DIR=/tmp REMOVE_DUPLICATES=true ASSUME_SORTED=true\n"
            .format(Pic, out, row[1]["sample_id"].replace(" ", "_"), out, row[1]["sample_id"].replace(" ", "_"), \
            out, row[1]["sample_id"].replace(" ", "_")))
            os.system("mv {}/BAM/{}.bai {}/BAM/{}.bam.bai" \
            .format(out, row[1]["sample_id"].replace(" ", "_"), out, row[1]["sample_id"].replace(" ", "_")))            

        ## is not paired move and rename bam
        else:
            os.system("mv {}/temp/{}_Aligned.sortedByCoord.out.bam {}/BAM/{}.bam" \
            .format(out, row[1]["sample_id"].replace(" ", "_"), out, row[1]["sample_id"].replace(" ", "_")))
            os.system("mv {}/temp/{}_Aligned.sortedByCoord.out.bam.bai {}/BAM/{}.bam.bai" \
            .format(out, row[1]["sample_id"].replace(" ", "_"), out, row[1]["sample_id"].replace(" ", "_")))

        ## make bigWigs from final BAMs
        os.system("{}/igvtools count --minMapQuality 30 {}/BAM/{}.bam stdout {} " \
        "| {}/wigToBigWig -clip stdin {} {}/BW/{}.bw" \
        .format(IGV, out, row[1]["sample_id"].replace(" ", "_"), IGVToolsGenome, \
        UCSC, ChromSizes, out, row[1]["sample_id"].replace(" ", "_")))

        ## append the name for the final BAM
        bamlist.append(row[1]["sample_id"].replace(" ", "_"))

    # Add columns to table of samples
    mapping_df["BAMs"] = [i + ".bam" for i in bamlist]
    mapping_df["BWs"] = [i + ".bw" for i in bamlist]

    # Count reads/fragments
    print("Counting reads using featureCounts")
    if paired:
        os.system("{} -p -s 2 -T {} -t exon -g gene_id -a {} -o {}/{}featureCounts.count {}" \
        .format(featureCounts, nCores, GTF, out, datetime.today().strftime('%Y-%m-%d'), " ".join("{}/BAM/{}.bam".format(out, i) for i in bamlist)))
    else:
        os.system("{} -s 2 -T {} -t exon -g gene_id -a {} -o {}/{}featureCounts.count {}" \
        .format(featureCounts, nCores, GTF, out, datetime.today().strftime('%Y-%m-%d'), " ".join("{}/BAM/{}.bam".format(out, i) for i in bamlist)))

    # Run RSeQC and fastqc
    print("generateing RSeQC reports")
    RSeQC = ["geneBody_coverage.py -r {} -i {} -o {}/QC/Combined"\
    .format(RefBed, ",".join("{}/BAM/{}.bam".format(out, i) for i in bamlist), out)]
    for bam in bamlist:
        RSeQC.append("read_distribution.py -i {}/BAM/{}.bam -r {} > {}/QC/{}_bam_stats.txt".format(out, bam, RefBed, out, bam))
        RSeQC.append("read_duplication.py -i {}/BAM/{}.bam -o {}/QC/{}".format(out, bam, out, bam, RefBed))
        RSeQC.append("junction_saturation.py -i {}/BAM/{}.bam -o {}/QC/{} -m 20 -r {}".format(out, bam, out, bam, RefBed))
        RSeQC.append("infer_experiment.py -i {}/BAM/{}.bam -r {} > {}/QC/{}_infer_experiment.txt ".format(out, bam, RefBed, out, bam))
        RSeQC.append("bam_stat.py -i {}/BAM/{}.bam > {}/QC/{}_bam_stats.txt".format(out, bam, out, bam))
        if paired:
            RSeQC.append("inner_distance.py -i {}/BAM/{}.bam -o {}/QC/{} -r {}".format(out, bam, out, bam, RefBed))
    for file in source_table["fwd_fastq"]:
        RSeQC.append("{} {}/{} --out {}".format(fastqc, fastq_directory, file, out + "/fastqc/"))
    if paired:
        for file in source_table["rv_fastq"]:
            RSeQC.append("{} {}/{} --out {}".format(fastqc, fastq_directory, file, out + "/fastqc/"))
    pool = Pool(int(nCores))                         # Create a multiprocessing Pool
    pool.map_async(os.system, RSeQC)
    pool.close()
    pool.join()

    #finish up with multiqc
    os.system("mv {}/{}featureCounts.count.summary {}/QC/{}featureCounts.count.summary".format(out, datetime.today().strftime('%Y-%m-%d'), out, datetime.today().strftime('%Y-%m-%d')))
    os.system("mv {}/temp/*.out {}/QC/.".format(out, out))
    os.system("rm {}/temp".format(out))
    os.system("multiqc {}/QC/ {}/fastqc/ -o {}".format(out, out, out))
    mapping_df.to_excel("{}/{}MappedSamples.xlsx".format(out, datetime.today().strftime('%Y-%m-%d')))
    print("all done")

