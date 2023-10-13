#!/bin/bash snakemake

## This group rules consist of programs necessary for alignment
# 1. BWA-METH
# 2. samtools
# 3. Picard
# 4. VerifyBamID2


rule alignment:
    input:
        read1 = ...,
        read2 = ...
    output:
        ...
    params:
        reference = ...
        prefix = ...
    threads: 16
    shell:
        '''
        #!/bin/bash bwa-meth

        bwameth.py --threads {threads} \
            --reference {params.reference} \
            {input.read1} {input.read2} > {params.prefix}

        '''


rule samtools_view:
    input:
        ...,
    output:
        bam = ...,
        idx = ...,
    log:
        "{sample}.log",
    params:
        extra="",  # optional params string
        region="",  # optional region string
    threads: 2
    wrapper:
        "v2.6.0/bio/samtools/view"


rule markduplicates_bam:
    input:
        bams=...,
    # optional to specify a list of BAMs; this has the same effect
    # of marking duplicates on separate read groups for a sample
    # and then merging
    output:
        bam=...,
        metrics=...,
    params:
        extra="--REMOVE_DUPLICATES true",
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


rule verify_bam_id:
    input:
        bam=...,
        # optional - this can be used to specify custom resource files if
        # necessary (if using GRCh37 or GRCh38 instead simply specify
        # params.genome_build="38", for example)
        # N.B. if svd_mu={prefix}.mu, then {prefix}.bed, {prefix}.UD, and
        # {prefix}.V must also exist
        svd_mu="ref.vcf.mu",
    output:
        selfsm=...,
        ancestry=...,
    params:
        # optional - see note for input.svd_mu
        # current choices are {37,38}
        genome_build="38",
    wrapper:
        "v2.6.0/bio/verifybamid/verifybamid2"