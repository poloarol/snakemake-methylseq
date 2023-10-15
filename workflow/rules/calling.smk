#!/bin/bash snakemake

## This group of rules consist of programs necessary for methylation calling
# 1. DNMTools


rule bsrate:
    input:
        "{outpath}/{project}/alignment/abismal/{sample}.sorted.filtered.bam"
    output:
        "{outpath}/{project}/methyl/stats/conversion_rate/{sample}.bsrate"
    params:
        # chromosomes = ..., # Do -ve control, +ve control, chrM and 3, 20, 15
        reference = config['ref']['genome']
    conda: "envs/dnmtools.yaml"
    threads: config['threads']
    shell:
        '''
        #!/bin/bash
        # Create a temp file for each chromosone and merge using bs-merge
        dnmtools bsrate -t {threads} -c {params.reference} -o {output} {input}
        '''

rule base_methylation:
    input:
        "{outpath}/{project}/alignment/abismal/{sample}.sorted.filtered.bam"
    output:
        "{outpath}/{project}/methyl/stats/counts/{sample}_single_base.meth"
    params:
        reference = config['ref']['genome']
    conda: "envs/dnmtools.yaml"
    threads: config['threads']
    shell:
        '''
        #!/bin/bash
        dnmtools counts -t {threads} -c {params.reference} -o {output} {input}
        '''

rule collapse_meth_counts:
    input:
        "{outpath}/{project}/methyl/stats/counts/{sample}_single_base.meth"
    output:
        "{outpath}/{project}/methyl/stats/counts/{sample}_symmetric.meth"
    params:
        reference = config['ref']['genome']
    conda: "envs/dnmtools.yaml"
    threads: config['threads']
    shell:
        '''
        #!/bin/bash
        dnmtools sym -o {output} {input}
        '''

rule global_meth_stats:
    input:
        "{outpath}/{project}/methyl/stats/counts/{sample}_single_base.meth"
    output:
        "{outpath}/{project}/methyl/stats/counts/{sample}_global.meth"
    params:
        reference = config['ref']['genome']
    conda: "envs/dnmtools.yaml"
    threads: config['threads']
    shell:
        '''
        #!/bin/bash
        dnmtools levels -o {output} {input}
        '''

rule hypomethylated_regions:
    input:
        "{outpath}/{project}/methyl/stats/counts/{sample}_symmetric.meth"
    output:
        "{outpath}/{project}/methyl/methylome/{sample}.hmr"
    conda: "envs/dnmtools.yaml"
    threads: config['threads']
    shell:
        '''
        #!/bin/bash
        dnmtools hmr -o {output} {input}
        '''


rule hypermethylated_regions:
    input:
        "{outpath}/{project}/methyl/stats/counts/{sample}_symmetric.meth"
    output:
        "{outpath}/{project}/methyl/methylome/{sample}.hypermr"
    conda: "envs/dnmtools.yaml"
    threads: config['threads']
    shell:
        '''
        #!/bin/bash
        dnmtools hypermr -o {output} {input}
        '''


rule allele_methylation:
    input:
        "{outpath}/{project}/alignment/abismal/{sample}.sorted.filtered.bam"
    output:
        epiread = "{outpath}/{project}/methyl/methylome/{sample}.epiread",
        sorted_epiread = "{outpath}/{project}/methyl/methylome/{sample}.sorted.epiread",
    params:
        reference = config['ref']['genome']
    conda: "envs/dnmtools.yaml"
    threads: config['threads']
    shell:
        '''
        #!/bin/bash

        dnmtools states -c {params.reference} -o {output.epiread} {input}
        sort -k1,1 -k2,2g {output.epiread} -o {output.sorted_epiread}
        '''


rule entropy:
    input:
        epiread = "{outpath}/{project}/methyl/methylome/{sample}.sorted.epiread"
    output:
        meth = "{outpath}/{project}/methyl/methylome/{sample}.entropy.meth"
    params:
        reference = config['ref']['genome']
    conda: "envs/dnmtools.yaml"
    threads: config['threads']
    shell:
        '''
        #!/bin/bash

        dnmtools entropy -w 5 -v -o {output} {params.reference} {input}
        '''


rule avg_meth_level_region:
    input:
        "{outpath}/{project}/methyl/stats/counts/{sample}_symmetric.meth"
    output:
        "{outpath}/{project}/methyl/methylome/{sample}.avg.bed"
    params:
        bed = config['ref']['region']
    conda: "envs/dnmtools.yaml"
    threads: config['threads']
    shell:
        '''
        #!/bin/bash
        dnmtools roi -o {output} {params.bed} {input}
        '''


# rule hydroxymethylation_estimation:
#     input:
#         bs_counts = ...,
#         ox_counts = ...,
#         tab_counts = ...
#     output:
#         ...
#     conda: "envs/dnmtools.yaml"
#     threads: config['threads']
#     shell:
#         '''
#         #!/bin/bash

#         dnmtools mlml -u {input.bs_counts} \
#             -m {input.ox_counts} -o {output}
#         '''