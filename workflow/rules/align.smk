#!/bin/bash snakemake

## This group of rules consist of programs necessary for alignment
# 1. BWA-METH
# 2. samtools
# 3. Picard
# 4. VerifyBamID2

threads: int = config['threads']

rule align_bwa:
    input:
        read1 = os.path.join(outpath, "{project}/qc/fastp/{sample}_1.trimmed.fastq"),
        read2 = os.path.join(outpath, "{project}/qc/fastp/{sample}_2.trimmed.fastq")
    output:
        "{outpath}/{project}/alignment/bwa/unsorted/{sample}.sam"
    params:
        reference = config['ref']['genome']
    threads: 16
    conda: "envs/bwameth.yaml"
    shell:
        '''
        #!/bin/bash bwa-meth

        bwameth.py --threads {threads} \
            --reference {params.reference} \
            {input.read1} {input.read2} > {output}

        '''


rule align_abismal:
    input:
        read1 = os.path.join(outpath, "{project}/qc/fastp/{sample}_1.trimmed.fastq"),
        read2 = os.path.join(outpath, "{project}/qc/fastp/{sample}_2.trimmed.fastq")
    output:
        bam = "{outpath}/{project}/alignment/abismal/{sample}.bam",
        # filtered_bam = "{outpath}/{project}/alignment/abismal/{sample}.filtered.bam"
        stats = "{outpath}/{project}/alignment/abismal/stats/{sample}.metrics.yaml"
    params:
        reference = config['ref']['genome']
    threads: 16
    conda: "envs/dnmtools.yaml"
    shell:
        '''
        #!/bin/bash abismal

        dnmtools abismal -threads {threads} \
            -genome {params.reference} \
            -stats {output.stats} \
            -outfile {output.bam} \
            {input.read1} {input.read2}
        '''

rule abismal_sort:
    input:
        "{outpath}/{project}/alignment/abismal/{sample}.bam",
    output:
        "{outpath}/{project}/alignment/abismal/{sample}.sorted.bam"
    params:
        extra="-m 4G",
    threads: 8
    wrapper:
        "v2.6.0/bio/samtools/sort"


rule abismal_index:
    input:
        "{outpath}/{project}/alignment/abismal/{sample}.sorted.bam"
    output:
        "{outpath}/{project}/alignment/abismal/{sample}.sorted.bai"
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v2.6.0/bio/samtools/index"


rule abismal_filter:
    input:
        bam = "{outpath}/{project}/alignment/abismal/{sample}.sorted.bam"
    output:
        bam = "{outpath}/{project}/alignment/abismal/{sample}.sorted.filtered.bam",
        stats = "{outpath}/{project}/alignment/abismal/stats/{sample}.filtered.metrics.yaml",
    params:
        reference = config['ref']['genome']
    threads: 16
    conda: "envs/dnmtools.yaml"
    shell:
        '''
        #!/bin/bash abismal

        dnmtools uniq -S {output.stats} {input.bam} {output.bam}
        '''

rule samtools_sort_bwa:
    input:
        "{outpath}/{project}/alignment/bwa/unsorted/{sample}.sam",
    output:
        "{outpath}/{project}/alignment/bwa/sorted/{sample}.bam"
    params:
        extra="-m 4G",
    threads: 8
    wrapper:
        "v2.6.0/bio/samtools/sort"


rule samtools_index_bwa:
    input:
        "{outpath}/{project}/alignment/bwa/sorted/{sample}.bam"
    output:
        "{outpath}/{project}/alignment/bwa/sorted/{sample}.bai"
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v2.6.0/bio/samtools/index"


rule markduplicates_bam_bwa:
    input:
        bams = "{outpath}/{project}/alignment/bwa/sorted/{sample}.bam",
    # optional to specify a list of BAMs; this has the same effect
    # of marking duplicates on separate read groups for a sample
    # and then merging
    output:
        bam = "{outpath}/{project}/alignment/bwa/sorted/picard/{sample}.bam",
        metrics = "{outpath}/{project}/alignment/bwa/sorted/picard/{sample}.metrics.txt",
    params:
        extra="--REMOVE_DUPLICATES true --CREATE_INDEX true",
    # optional specification of memory usage of the JVM that snakemake 
    # will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission 
    # as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024,
    wrapper:
        "v2.6.0/bio/picard/markduplicates"


rule verify_bam_id_abismal:
    input:
        bam = "{outpath}/{project}/alignment/abismal/{sample}.sorted.filtered.bam",
        ref = config['ref']['genome'],
        # optional - this can be used to specify custom resource files if
        # necessary (if using GRCh37 or GRCh38 instead simply specify
        # params.genome_build="38", for example)
        # N.B. if svd_mu={prefix}.mu, then {prefix}.bed, {prefix}.UD, and
        # {prefix}.V must also exist
        svd_mu=config['ref']['svd_mu'],
    output:
        selfsm = "{outpath}/{project}/alignment/abismal/verify_bam_id/{sample}.selfsm",
        ancestry = "{outpath}/{project}/alignment/abismal/verify_bam_id/{sample}.ancestry",
    params:
        # optional - see note for input.svd_mu
        # current choices are {37,38}
        genome_build="38",
    wrapper:
        "v2.6.0/bio/verifybamid/verifybamid2"


rule verify_bam_id_bwa:
    input:
        bam = "{outpath}/{project}/alignment/bwa/sorted/picard/{sample}.bam",
        ref = config['ref']['genome'],
        # optional - this can be used to specify custom resource files if
        # necessary (if using GRCh37 or GRCh38 instead simply specify
        # params.genome_build="38", for example)
        # N.B. if svd_mu={prefix}.mu, then {prefix}.bed, {prefix}.UD, and
        # {prefix}.V must also exist
        svd_mu=config['ref']['svd_mu'],
    output:
        selfsm = "{outpath}/{project}/alignment/bwa/sorted/verify_bam_id/{sample}.selfsm",
        ancestry = "{outpath}/{project}/alignment/bwa/sorted/verify_bam_id/{sample}.ancestry",
    params:
        # optional - see note for input.svd_mu
        # current choices are {37,38}
        genome_build="38",
    wrapper:
        "v2.6.0/bio/verifybamid/verifybamid2"