# =================================================================================================
#     Dependencies
# =================================================================================================


configfile: "config.yaml"


## this rule runs a python script that will generate the bed file containing the contributing loci and their effect sizes based on values of dn alpha
## the bed file will be then used to annotate a vcf file that will be used by SliM to run the simulations


rule all:
    input:
        expand(
            "results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_pop_size_early.txt",
            pi=config["pi"],
            selection=config["selection"],
            heritability=config["heritability"],
            replicates_arq=config["replicates_arq"],
            optima=config["optima"],
            replicates_sim=config["replicates_sim"],    
        ),


rule build_population_for_sim:
    input:
        og_tree_offset=config["og_tree_offset"],
        og_vcf_offset=config["og_vcf_offset"],
        polygenicity_params=config["polygenicity_params"],
    output:
        tree_seq_causalloci="results/arq_pi{pi}_{replicates_arq}/tree_seq_causalloci.trees",
        loci_effectsize="results/arq_pi{pi}_{replicates_arq}/loci_effectsize.csv",
        phenotypes="results/arq_pi{pi}_{replicates_arq}/phenotypes.csv",
    params:
        pi=lambda wildcards: str(wildcards.pi),
        replicates_arq=lambda wildcards: str(wildcards.replicates_arq),
        beta=config["beta"],
    resources:
        mem_mb=30720,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/build_population_for_sim.py"

rule run_slim_simulation:
    input:
        tree_seq_causalloci="results/arq_pi{pi}_{replicates_arq}/tree_seq_causalloci.trees",
    output: 
        #output_tree_gen4=temp("results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_tree_output_gen4.trees"),
        #output_tree_gen10="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_tree_output_gen10.trees",
        output_pop_size_early="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_pop_size_early.txt",
        output_va="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_va.txt",
        output_mfitness="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_mfitness.txt",
        output_vfitness="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_vfitness.txt",
        output_mpheno="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_mpheno.txt",
        output_vpheno="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_vpheno.txt",
        output_new_optimum="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_new_optimum.txt",
        output_adj_variance="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_adj_variance.txt",
        output_maxphenotype="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_maxphenotype.txt",
        output_minphenotype="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_minphenotype.txt",

    resources:
        mem_mb=40960,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/slim.sh"

rule tree_postprocessing:
    input:
        og_tree_offset=config["og_tree_offset"],
        mapper_ids=config['mapper_realid_metadataid'],
        output_sim_tree="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_tree_output_gen4.trees",
    output:
        output_vcf=temp("results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_vcfgen4_output.vcf"),
    resources:
        mem_mb=30720,
        limit_space=1,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/tree_postprocessing.py"

rule calc_ecotype_counts:
    input:
        nonhet_pos=config['nonhet_pos'],
        og_vcf_offset=config["og_vcf_offset"],
        ecotypes_grenenet=config['ecotypes_grenenet'],
        output_vcf_offset="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_vcfgen4_output.vcf",
    output:
        ecotype_counts="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_ecotype_counts10env.csv",
    resources:
        mem_mb=30720,
    threads: 20,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/calc_ecotype_counts.py"

rule gen_allele_freq:
    input:
        pos_vcf_og_offset=config['pos_vcf_og_offset'],
        output_vcf = expand(
            "results/arq_pi{{pi}}_{{replicates_arq}}/{{heritability}}/{{selection}}/optima{optima}/subp{replicates_sim}_vcfgen4_output.vcf",
            optima=config["optima"],
            replicates_sim=config["replicates_sim"],    
        ),
    output:
        allele_counts ="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/allele_counts10env.csv",
        allele_freq ="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/allele_freq10env.csv",
    resources:
        mem_mb=30720,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/allele_freq_calc.py"

rule prep_lmm:
    input:
        allele_freq_founder_offset=config['allele_freq_founder_offset'],
        pc_founders=config['pc_founders'],
        ecotype_counts ="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/ecotype_counts10env.csv",
        allele_freq ="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/allele_freq10env.csv",
    output:
        env_variable ="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/lmm/env_variable10env.csv",
        allele_freq_norm="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/allele_freq_norm10env.csv",
        pop_structure ="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/lmm/pop_structure10env.csv",
    resources:
        mem_mb=30720,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/prep_lmm.py"

rule run_lmm_wpc:
    input:
        env_sites="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/lmm/env_variable10env.csv",
        p_norm="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/allele_freq_norm10env.csv",
        pop_structure ="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/lmm/pop_structure10env.csv",
    output:
        lmm_results ="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/lmm/lmm_pc_results10env.csv",
    benchmark:
        "benchmarks/lmm10env/arq_pi{pi}_{replicates_arq}_{heritability}_{selection}.txt"
    resources:
        mem_mb=61440,
    threads: 20,
    conda:
        "envs/r.yaml"
    script:
        "scripts/run_lmm_lmer_wpc.R"

rule run_lmm_nopc:
    input:
        env_sites="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/lmm/env_variable10env.csv",
        p_norm="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/allele_freq_norm10env.csv",
        pop_structure ="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/lmm/pop_structure10env.csv",
    output:
        lmm_results ="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/lmm/lmm_nopc_results10env.csv",
    benchmark:
        "benchmarks/lmm_nopc10env/arq_pi{pi}_{replicates_arq}_{heritability}_{selection}.txt"
    resources:
        mem_mb=61440 ,
    threads: 20,
    conda:
        "envs/r.yaml"
    script:
        "scripts/run_lmm_lmer_nopc.R"

rule create_famfile_gwa:
    input:
        ecotype_counts_file="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/ecotype_counts.csv",
        fam_file_input=config["fam_file"]
    output:
        fam_file_ouput = "results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/gwa/geno.fam",
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
        fam_file="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/gwa/geno.fam",
        bed_file=config['bed_file'],
        bim_file=config['bim_file'],
        kinship=config['kinship']
    output:
        output_gwas="results/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/gwa/output/results_nmaf.assoc.txt",
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/gwa_nmaf/arq_pi{pi}_{replicates_arq}_{heritability}_{selection}_optima{optima}.txt"
    conda:
        "envs/gwas.yaml"
    script:
        "scripts/gwas_gemma.sh"



