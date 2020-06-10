#########################################
# Snakemake pipeline for RNA-Seq analysis
#########################################


###########
# Libraries
###########
import pandas as pd
import glob

###############
# Configuration
###############

configfile: "config.yaml" # where to find parameters
WORKING_DIR = config["working_dir"]
RESULT_DIR = config["result_dir"]

# fetch URL to transcriptome multi fasta from configfile
genome_url = config["refs"]["genome"]
transcriptome_gtf_url= config["refs"]["transcriptome_gtf"]

functional_annotation = config["refs"]["annotation"]


########################
# Samples and conditions
########################

# read the tabulated separated table containing the sample, condition and fastq file information∂DE
units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)

# create lists containing the sample names and conditions
SAMPLES = units.index.get_level_values('sample').unique().tolist()
#SAMPLES = ["SRR212121","SRR212122"]
samples = pd.read_csv(config["units"], dtype=str,index_col=0,sep="\t")
#CONDITIONS = list(pd.read_table(config["units"])["condition"])
samplefile = config["units"]

#paired = dict(zip(list(units["sample"]),list(units["paired"])))
#print(paired)

###########################
# Input functions for rules
###########################

def sample_is_single_end(sample):
    """This function checks if raeds are single or paired end"""
    if os.stat(sample+"_2.fastq.gz").st_size < 100:
        return True
    else:
        return False

def get_fastq(wildcards):
    """ This function checks if the sample has paired end or single end reads
    and returns 1 or 2 names of the fastq files """
    if sample_is_single_end(wildcards.sample):
        return WORKING_DIR + wildcards.sample + "_1.fastq"
    else:
        return [WORKING_DIR + wildcards.sample + "_1.fastq", WORKING_DIR + wildcards.sample + "_2.fastq"]

def get_trimmed(wildcards):
    """ This function checks if sample is paired end or single end
    and returns 1 or 2 names of the trimmed fastq files """
    if sample_is_single_end(wildcards.sample):
        return WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz"
    else:
        return [WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz", WORKING_DIR + "trimmed/" + wildcards.sample + "_R2_trimmed.fq.gz"]


#################
# Desired outputs
#################

rule all:
    input:
        fw  = expand(WORKING_DIR + "fastq/{sample}_1.fastq",sample = SAMPLES),
        
        #bam = expand(WORKING_DIR + "mapped/{sample}.bam", sample = SAMPLES),
        #RESULT_DIR + "counts.txt",
        fq1 = expand(WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",sample = SAMPLES),
        fq2 = expand(WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz",sample = SAMPLES),

    message:
        "Job done! Removing temporary directory"

################
##### Rules ####
################

###########################################
# Download genome, annotation and SRA-files
###########################################

rule get_genome_fasta:
    output:
        WORKING_DIR + "genome/genome.fasta"
    message:
        "downloading the required genomic fasta file"
    #conda:
    #    "envs/wget.yaml"
    shell:
        "wget -O {output} {genome_url}"

rule get_SRR_files:
    output:
        fw = WORKING_DIR + "fastq/{sample}_1.fastq",
        rev= WORKING_DIR + "fastq/{sample}_2.fastq"
    params:
       SRA = "{sample}",
       DIR = "fastq"
    message:
        "using fastq-dump to download SRA data files."
    #conda:
    #    "envs/wget.yaml"
    shell:
        "touch {output.rev}; fastq-dump --split-files {params.SRA} -O {params.DIR}"       

#rule get_transcriptome_gtf:
#    output:
#        WORKING_DIR + "genome/ref_transcriptome.gff"
#    message:
#        "downloading required transcriptome gtf file"
#    #conda:
#    #    "envs/wget.yaml"
#    shell:
#        "wget -O {output} {transcriptome_gtf_url}"


##################################
# Fastp
##################################

rule fastp:
    input:
        fw = WORKING_DIR + "fastq/{sample}_1.fastq",
        rev= WORKING_DIR + "fastq/{sample}_2.fastq"
    output:
        fq1  = WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",
        fq2  = WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz",
        html = RESULT_DIR + "fastp/{sample}.html"
    message:"trimming {wildcards.sample} reads"
    threads: 10
    log:
        RESULT_DIR + "fastp/{sample}.log.txt"
    params:
        paired = "{lambda wildcards:paired[wildcards.SAMPLES]}",
        sampleName = "{sample}",
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"]
    run:
        if sample_is_single_end(params.sampleName):
            shell("fastp --thread {threads}  --html {output.html} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --in1 {input.fw} --out1 {output} \
            2> {log}; \
            touch {output.fq2}")
        else:
            shell("fastp --thread {threads}  --html {output.html} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --detect_adapter_for_pe \
            --in1 {input.fw} --in2 {input.rev} --out1 {output.fq1} --out2 {output.fq2}; \
            2> {log}")


#########################
# RNA-Seq read alignement
#########################

rule index:
    input:
        WORKING_DIR + "genome/genome.fasta"
    output:
        [WORKING_DIR + "genome/genome." + str(i) + ".ht2" for i in range(1,9)]
    message:
        "indexing genome"
    params:
        WORKING_DIR + "genome/genome"
    threads: 10
    shell:
        "hisat2-build -p {threads} {input} {params} --quiet"

rule hisat_mapping:
    input:
        get_trimmed,
        indexFiles = [WORKING_DIR + "genome/genome." + str(i) + ".ht2" for i in range(1,9)]
    output:
        bams  = WORKING_DIR + "mapped/{sample}.bam",
        sum   = RESULT_DIR + "logs/{sample}_sum.txt",
        met   = RESULT_DIR + "logs/{sample}_met.txt"
    params:
        indexName = str(WORKING_DIR + "genome/genome"),
        sampleName = "{sample}"
    # conda:
    #     "envs/hisat_mapping.yaml"
    message:
        "mapping reads to genome to bam files."
    threads: 10
    run:
        if sample_is_single_end(params.sampleName):
            shell("hisat2 -p {threads} --summary-file {output.sum} --met-file {output.met} -x {params.indexName} \
            -U {input} | samtools view -Sb -F 4 -o {output.bams}")
        else:
            shell("hisat2 -p {threads} --summary-file {output.sum} --met-file {output.met} -x {params.indexName} \
            -1 {input[0]} -2 {input[1]} | samtools view -Sb -F 4 -o {output.bams}")


######################
# Get raw counts table
######################

#rule create_counts_table:
#    input:
#        bams = expand(WORKING_DIR + "mapped/{sample}.bam", sample = SAMPLES),
#        gff  = WORKING_DIR + "genome/ref_transcriptome.gff"
#    output:
#        RESULT_DIR + "counts.txt"
    # conda:
    #     "envs/subread.yaml"
#    shell:
#        "featureCounts -O -t mRNA -g ID -F 'gtf' -a {input.gff} -o {output} {input.bams}"
