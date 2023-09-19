# =================================================================================================
#     Dependencies
# =================================================================================================


configfile: "config.yaml"

replicates = list(range(0, config["replicates"]))
optima_index = list(range(0, config["optima_qty"]))

## this rule runs a python script that will generate the bed file containing the contributing loci and their effect sizes based on values of beta dn alpha
## the bed file will be then used to annotate a vcf file that will be used by SliM to run the simulations


rule all:
    input:
        expand(
            'results/arq_{allele_freq}_{pi}_{beta}/{selection}/optima{optima_index}/subp{replicates}_tree_output.trees',
            allele_freq=config['allele_freq'],
            pi=config["pi"],
            beta=config["beta"],
            selection=config["selection"],
            optima_index=optima_index,
            replicates=replicates
        ),

rule build_population_for_sim:
    input:
        og_tree_offset=config["og_tree_offset"],
        og_vcf_offset=config["og_vcf_offset"],
    output:
        tree_seq_causalloci="results/arq_{allele_freq}_{pi}_{beta}/tree_seq_causalloci.trees",
        loci_effectsize="results/arq_{allele_freq}_{pi}_{beta}/loci_effectsize.csv",
        phenotypes="results/arq_{allele_freq}_{pi}_{beta}/phenotypes.csv",
        optima_values="results/arq_{allele_freq}_{pi}_{beta}/optima_values.txt",
        variance_values="results/arq_{allele_freq}_{pi}_{beta}/variance_values.txt",
    params:
        optima_qty=config["optima_qty"],
        pi=lambda wildcards: str(wildcards.pi),
        beta=lambda wildcards: str(wildcards.beta),
        allele_freq=lambda wildcards: str(wildcards.allele_freq),
        lowfreq=config["lowfreq"],
        mediumfreq=config["mediumfreq"],
        highfreq=config["highfreq"],
        monogen=config["monogen"],
        fivepoly=config["fivepoly"],
        twentypoly=config["twentypoly"],
        onehpoly=config["onehpoly"],
        lowbeta=config["lowbeta"],
        highbeta=config["highbeta"],
    resources:
        mem_mb=30720,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/build_population_for_sim.py"

rule run_slim_simulation:
    input:
        tree_seq_causalloci="results/arq_{allele_freq}_{pi}_{beta}/tree_seq_causalloci.trees",
        optima_values="results/arq_{allele_freq}_{pi}_{beta}/optima_values.txt",
        variance_values="results/arq_{allele_freq}_{pi}_{beta}/variance_values.txt",
    output: 
        output_tree="results/arq_{allele_freq}_{pi}_{beta}/{selection}/optima{optima_index}/subp{replicates}_tree_output.trees",
        output_pop_size="results/arq_{allele_freq}_{pi}_{beta}/{selection}/optima{optima_index}/subp{replicates}_pop_size.txt",
    params:
        optima_index=lambda wildcards: str(wildcards.optima_index),        
        selection=lambda wildcards: str(wildcards.selection),
    resources:
        mem_mb=20480,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/slim.sh"