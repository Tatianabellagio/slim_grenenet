# =================================================================================================
#     Dependencies
# =================================================================================================


configfile: "config.yaml"

replicates = list(range(0, config["replicates"]))
optima = list(range(0, config["optima"]))

## this rule runs a python script that will generate the bed file containing the contributing loci and their effect sizes based on values of beta dn alpha
## the bed file will be then used to annotate a vcf file that will be used by SliM to run the simulations


rule all:
    input:
        expand(
            "results/arq_{allele_freq}_{pi}_{beta}_{selection}/optima{optima}/subp{replicates}_vcf_output_vcf",
            allele_freq=config['allele_freq'],
            pi=config["pi"],
            beta=config["beta"],
            selection=config["selection"],
            optima=optima,
            replicates=replicates,    
        ),

rule build_population_for_sim:
    input:
        og_tree_offset=config["og_tree_offset"],
    output:
        tree_seq_causalloci="results/arq_{allele_freq}_{pi}_{beta}_{selection}/tree_seq_causalloci.trees",
        bed_sc="results/arq_{allele_freq}_{pi}_{beta}_{selection}/loci_effectsize.csv",
        phenotypes="results/arq_{allele_freq}_{pi}_{beta}_{selection}/phenotypes.csv",
        optima_values="results/arq_{allele_freq}_{pi}_{beta}_{selection}/optima_values.txt",
        variance_values="results/arq_{allele_freq}_{pi}_{beta}_{selection}/variance_values.txt",
    params:
        n_optima=config["optima"],
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
        strongsel=config["strongsel"],
        moderatesel=config["moderatesel"],
        lowsel=config["lowsel"],

    resources:
        mem_mb=30720,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/sc_optima_calc.py"

rule run_slim_simulation:
    input:
        tree_seq_causalloci="results/arq_{allele_freq}_{pi}_{beta}_{selection}/tree_seq_causalloci.trees",
        optima_values="results/arq_{allele_freq}_{pi}_{beta}_{selection}/optima_values.txt",
    output: 
        "results/arq_{allele_freq}_{pi}_{beta}_{selection}/optima{optima}/subp{replicates}_tree_output.trees",
    params:
        optima=lambda wildcards: str(wildcards.optima),
    resources:
        mem_mb=51200,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/slim.sh"


rule tree_postprocessing:
    input:
        expand(
            "results/arq_{{allele_freq}}_{{pi}}_{{beta}}_{{selection}}/optima{optima}/subp{replicates}_tree_output.trees",
            optima=optima,    
            replicates=replicates,    
        ),

    output:
        "results/arq_{allele_freq}_{pi}_{beta}_{selection}/optima{optima}/subp{replicates}_vcf_output_vcf",
    params:
        pi=lambda wildcards: str(wildcards.pi),
        beta=lambda wildcards: str(wildcards.beta),
        allele_freq=lambda wildcards: str(wildcards.allele_freq),
    resources:
        mem_mb=30720,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/allele_freq_pool_sizes.py"