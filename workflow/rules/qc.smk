#!/bin/bash snakemake

## This group of rules aim to perform Quality Control.
## -----------
## 1. FastQC
## 2. Cutadapt
## 3. multiqc

rule fastqc:
    input:
        os.path.join(rawpath, "{sample}_{read}.fastq")
    output:
        html = "{outpath}/{project}/qc/fastqc/{sample}_{read}_fastqc.html",
        zip = "{outpath}/{project}/qc/fastqc/{sample}_{read}_fastqc.zip"
    params:
        extra = "--quiet"
    threads: 1
    resources: 
        mem_mb = 1024
    wrapper:
        "v2.6.0/bio/fastqc"

rule fastp_pe:
    input:
        sample = [os.path.join(rawpath, "{sample}_1.fastq"),
                    os.path.join(rawpath, "{sample}_2.fastq")]
    output:
        trimmed = ["{outpath}/{project}/qc/fastp/{sample}_1.trimmed.fastq",
            "{outpath}/{project}/qc/fastp/{sample}_2.trimmed.fastq"],
        
        unpaired1 = "{outpath}/{project}/qc/fastp/{sample}.u1.fastq",
        unpaired2 = "{outpath}/{project}/qc/fastp/{sample}.u2.fastq",
        
        merged = "{outpath}/{project}/qc/fastp/{sample}.merged.fastq",
        failed = "{outpath}/{project}/qc/fastp/{sample}.failed.fastq",
        html = "{outpath}/{project}/qc/fastp/{sample}.html",
        json = "{outpath}/{project}/qc/fastp/{sample}.json"
    params:
        extra="--merge"
    threads: 4  # set desired number of threads here
    wrapper:
        "v2.6.0/bio/fastp"


rule fastqc_trimmed:
    input:
        os.path.join(outpath, "{project}/qc/fastp/{sample}_{read}.trimmed.fastq")
    output:
        html = "{outpath}/{project}/qc/fastp/{sample}_{read}.trimmed_fastqc.html",
        zip = "{outpath}/{project}/qc/fastp/{sample}_{read}.trimmed_fastqc.zip"
    params:
        extra = "--quiet"
    threads: 1
    resources: 
        mem_mb = 1024
    wrapper:
        "v2.6.0/bio/fastqc"

rule multiqc:
    # input:
    #     expand("{outpath}/{project}/qc/fastqc/{sample}.html",
    #         outpath=outpath, project=project, sample=basenames)
    output:
        "{outpath}/{project}/qc/trimmed/multiqc_report.html",
        "{outpath}/{project}/qc/untrimmed/multiqc_report.html"
    params:
        in_raw = "{outpath}/{project}/qc/fastqc",
        out_raw = "{outpath}/{project}/qc/untrimmed",
        in_trimmed = "{outpath}/{project}/qc/fastp",
        out_trimmed = "{outpath}/{project}/qc/trimmed",
    conda: "envs/multiqc.yaml"
    shell:
        '''
            #!/bin/bash multiqc

            multiqc {params.in_raw} --outdir {params.out_raw}
            multiqc {params.in_trimmed} --outdir {params.out_trimmed}
        '''