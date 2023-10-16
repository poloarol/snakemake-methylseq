#!/bin/bash snakemake

# conda activate methylseq

rules=""

analysis=$1
tool=$2

if [ "$analysis" = "all" ]
then
    rules="all"
else
    if [ "$analysis" = "qc" ]
    then
        echo "Preparing to launch snakemake-methylseq pipeline to run Quality Control using fastqc, fastp and multqc"
        rules="fastqc fastp_pe fastqc_trimmed multiqc"
    elif [ "$analysis" = "alignment" ]
    then
        if [ "$tool" = "bwa" ]
        then
            echo "Preparing to launch snakemake-methylseq pipeline to run alignment using bwa-meth, samtools and picard"
            rules="align_bwa samtools_sort_bwa samtools_index_bwa markduplicates_bam_bwa verify_bam_id_bwa"
        else
            echo "Preparing to launch snakemake-methylseq pipeline to run alignment using abismal and samtools"
            rules="align_abismal abismal_sort abismal_index abismal_filter verify_bam_id_abismal"
        fi
    else
        if [ "$tool" = "dnmtools" ]
        then
            echo "Preparing to launch snakemake-methylseq pipeline to run methylation analysis using DNMTools"
            rules="bsrate base_methylation collpase_meth_counts global_meth_stats hypomethylated_regions hypermethylated_regions allele_methylation"
        else
            # rules=""
            echo "Unfortunately MethylDackel is not fully supported :("
            exit 1
        fi
    fi
fi

snakemake -n \
    --use-conda --conda-frontend conda \
    --allowed-rules $rules