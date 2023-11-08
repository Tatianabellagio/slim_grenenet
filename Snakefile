# =================================================================================================
#     Dependencies
# =================================================================================================


configfile: "config.yaml"


## this rule runs a python script that will generate the bed file containing the contributing loci and their effect sizes based on values of dn alpha
## the bed file will be then used to annotate a vcf file that will be used by SliM to run the simulations


rule all:
    input:
        expand(
            'results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/lmm/lmm_nopc_results10env.csv',
            allele_freq=config['allele_freq'],
            pi=config["pi"],
            selection=config["selection"],
            heritability=config["heritability"],
            replicates_arq=config["replicates_arq"],
        ),
        expand(
            'results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/lmm/lmm_pc_results10env.csv',
            allele_freq=config['allele_freq'],
            pi=config["pi"],
            selection=config["selection"],
            heritability=config["heritability"],
            replicates_arq=config["replicates_arq"],
        ),

rule run_lmm_wpc:
    input:
        env_sites="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/lmm/env_variable10env.csv",
        p_norm="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/allele_freq_norm10env.csv",
        pop_structure ="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/lmm/pop_structure10env.csv",
    output:
        lmm_results ="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/lmm/lmm_pc_results10env.csv",
    benchmark:
        "benchmarks/lmm10env/arq_{allele_freq}_{pi}_{replicates_arq}_{heritability}_{selection}.txt"
    resources:
        mem_mb=61440,
    threads: 20,
    conda:
        "envs/r.yaml"
    script:
        "scripts/run_lmm_lmer_wpc.R"

rule run_lmm_nopc:
    input:
        env_sites="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/lmm/env_variable10env.csv",
        p_norm="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/allele_freq_norm10env.csv",
        pop_structure ="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/lmm/pop_structure10env.csv",
    output:
        lmm_results ="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/lmm/lmm_nopc_results10env.csv",
    benchmark:
        "benchmarks/lmm_nopc10env/arq_{allele_freq}_{pi}_{replicates_arq}_{heritability}_{selection}.txt"
    resources:
        mem_mb=61440 ,
    threads: 20,
    conda:
        "envs/r.yaml"
    script:
        "scripts/run_lmm_lmer_nopc.R"

rule create_famfile_gwa:
    input:
        ecotype_counts_file="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/ecotype_counts.csv",
        fam_file_input=config["fam_file"]
    output:
        fam_file_ouput = "results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/gwa/geno.fam",
    params:
        optima=lambda wildcards: str(wildcards.optima),
    resources:
        mem_mb=15350,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/generate_fam_files.py"

rule run_gwa:
    input:
        fam_file="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/gwa/geno.fam",
        bed_file=config['bed_file'],
        bim_file=config['bim_file'],
        kinship=config['kinship']
    output:
        output_gwas="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/gwa/output/results_nmaf.assoc.txt",
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/gwa_nmaf/arq_{allele_freq}_{pi}_{replicates_arq}_{heritability}_{selection}_optima{optima}.txt"
    conda:
        "envs/gwas.yaml"
    script:
        "scripts/gwas_gemma.sh"



