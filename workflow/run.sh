#!/bin/bash snakemake

conda activate methylseq

rules=""

analysis=$1
tool=$2

if [ "$analysis" = "qc" ]
then
    rules="fastqc fastp_pe fastqc_trimmed multiqc"
elif [ "$analysis" = "alignment" ]
then
    if [ "$tool" = "bwa" ]
    then
        rules="align_bwa samtools_sort_bwa samtools_index_bwa markduplicates_bam_bwa verify_bam_id_bwa"
    else
        rules="align_abismal abismal_sort abismal_index abismal_filter verify_bam_id_abismal"
    fi
else
    if [ "$tool" = "dnmtools" ]
    then
        rules="bsrate base_methylation collpase_meth_counts global_meth_stats hypomethylated_regions hypermethylated_regions allele_methylation"
    else
        # rules=""
        echo "Unfortunately MethylDackel is supported as of now :("
        exit 1
    fi
fi

snakemake --cores $cores \
    --use-conda --conda-frontend conda \
    --allowed-rules $rules \
    --latency-wait 2800