# =================================================================================================
#     Dependencies
# =================================================================================================


configfile: "config.yaml"


replicates = list(range(0, config["replicates"]))
optima = list(range(0, config["optima"]))
generation = list(range(2,7))

## this rule runs a python script that will generate the bed file containing the contributing loci and their effect sizes based on values of beta dn alpha
## the bed file will be then used to annotate a vcf file that will be used by SliM to run the simulations


rule all:
    input:
       expand(
           "results/arq_{allele_freq}_{pi}_{beta}/lfmm_generation{generation}/manhattan_mean.png",
           allele_freq=config['allele_freq'],
           pi=config["pi"],
           beta=config["beta"],
           generation=generation,
       ),
        expand(
            "results/arq_{allele_freq}_{pi}_{beta}/gwas_generation{generation}/optima{optima}/output/results.assoc.txt",
            allele_freq=config['allele_freq'],
            pi=config["pi"],
            beta=config["beta"],
            generation=generation,
            optima=optima,
        ),



rule run_python_script:
    input:
        base_vcf=config["base_vcf"],
        optimum_ecotypes=config["optimum_ecotypes"],
        base_fasta=config["base_fasta"],
    output:
        bed_sc="results/arq_{allele_freq}_{pi}_{beta}/selection_coef.bed",
        slim="results/arq_{allele_freq}_{pi}_{beta}/optima_files/optima_slim.txt",
    params:
        n_optima=config["optima"],
        pi=lambda wildcards: str(wildcards.pi),
        beta=lambda wildcards: str(wildcards.beta),
        allele_freq=lambda wildcards: str(wildcards.allele_freq),
        chr_number=config["chr_number"],
        n_ecotypes=config["n_ecotypes"],
        lowfreq=config["lowfreq"],
        mediumfreq=config["mediumfreq"],
        highfreq=config["highfreq"],
	    monogen=config["monogen"],
        fivepoly=config["fivepoly"],
        tenpoly=config["tenpoly"],
        fifthpoly=config["fifthpoly"],
        twentypoly=config["twentypoly"]
    resources:
        mem_mb=30720,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/sc_optima_calc.py"


rule run_bash_script:
    input:
        vcf_file=config["base_vcf"],
        bed_sc="results/arq_{allele_freq}_{pi}_{beta}/selection_coef.bed",
    output:
        "results/arq_{allele_freq}_{pi}_{beta}/annotated.vcf",
    threads: 10
    log:
        out="log/arq_{allele_freq}_{pi}_{beta}_run_bash_script_stdout.log",
        err="log/arq_{allele_freq}_{pi}_{beta}_run_bash_script_stderr.err",
    resources:
        mem_mb=30720,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/annotate.sh"


rule run_slim_script:
    input:
        fasta=config["base_fasta"],
        vcf="results/arq_{allele_freq}_{pi}_{beta}/annotated.vcf",
        optima_slim="results/arq_{allele_freq}_{pi}_{beta}/optima_files/optima_slim.txt",
    output: 
        "results/arq_{allele_freq}_{pi}_{beta}/optima{optima}/subp{replicates}_slim_output.txt",
    params:
        optima=lambda wildcards: str(wildcards.optima),
        initial_pop=config["initial_pop"],
    resources:
        mem_mb=51200,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/slim.sh"


rule create_allele_freq_pool_size_tables:
    input:
        expand(
            "results/arq_{{allele_freq}}_{{pi}}_{{beta}}/optima{optima}/subp{replicates}_slim_output.txt",
            optima=optima,    
            replicates=replicates,    
        ),
        og_allele_freq=config["og_allele_freq"],

    output:
        expand(
            ["results/arq_{{allele_freq}}_{{pi}}_{{beta}}/allele_freq_generation{generation}.csv",
            "results/arq_{{allele_freq}}_{{pi}}_{{beta}}/deltap_norm_generation{generation}.csv",
            "results/arq_{{allele_freq}}_{{pi}}_{{beta}}/pop_sizes_generation{generation}.csv"],
            generation=generation,  
        ),
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


rule gen_lfmm_files:
    input:
        allele_freq_table = "results/arq_{allele_freq}_{pi}_{beta}/allele_freq_generation{generation}.csv",
        pool_sizes = "results/arq_{allele_freq}_{pi}_{beta}/pop_sizes_generation{generation}.csv",
        optimas = "results/arq_{allele_freq}_{pi}_{beta}/optima_files/optima_slim.txt"
    output:
        gen = "results/arq_{allele_freq}_{pi}_{beta}/lfmm_generation{generation}/geno_mean.csv",
        env = "results/arq_{allele_freq}_{pi}_{beta}/lfmm_generation{generation}/env_mean.csv",
        num_components = "results/arq_{allele_freq}_{pi}_{beta}/lfmm_generation{generation}/num_components_mean.txt",
    params:
        n_replicates=config["replicates"],
    resources:
        mem_mb=15350,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/gen_lfmm_files.py"


rule run_lfmm:
    input:
        gen = "results/arq_{allele_freq}_{pi}_{beta}/lfmm_generation{generation}/geno_mean.csv",
        env = "results/arq_{allele_freq}_{pi}_{beta}/lfmm_generation{generation}/env_mean.csv",
        num_components = "results/arq_{allele_freq}_{pi}_{beta}/lfmm_generation{generation}/num_components_mean.txt"
    output:
        p_values_lfmm = "results/arq_{allele_freq}_{pi}_{beta}/lfmm_generation{generation}/p_values_mean.csv",
        qq_plot = "results/arq_{allele_freq}_{pi}_{beta}/lfmm_generation{generation}/qq_plot_mean.png"
    resources:
        mem_mb=30720,
    conda:
        "envs/r.yaml"
    script:
        "scripts/run_lfmm.R"

rule results_lfmm:
    input:
        allele_freq_table = "results/arq_{allele_freq}_{pi}_{beta}/allele_freq_generation{generation}.csv",
        effect_sizes = "results/arq_{allele_freq}_{pi}_{beta}/selection_coef.bed",
        p_values_lfmm = "results/arq_{allele_freq}_{pi}_{beta}/lfmm_generation{generation}/p_values_mean.csv"
    output:
        plot_lfmm = "results/arq_{allele_freq}_{pi}_{beta}/lfmm_generation{generation}/manhattan_mean.png",
        p_values_vs_causal = "results/arq_{allele_freq}_{pi}_{beta}/lfmm_generation{generation}/causal_loci_vs_flmm_mean.csv"
    resources:
        mem_mb=15350,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/plotting_lfmm.py"

rule create_pheno_og:
    input:
        vcf_file=config["base_vcf"],
        bed_sc="results/arq_{allele_freq}_{pi}_{beta}/selection_coef.bed",
    output:
        output_pheno="results/arq_{allele_freq}_{pi}_{beta}/pheno_og.csv",
    resources:
        mem_mb=15350,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/pheno_og.py"  

rule create_pheno_gwas:
    input:
        allele_counts = "results/arq_{allele_freq}_{pi}_{beta}/allele_freq_generation{generation}.csv",
        pheno_og="results/arq_{allele_freq}_{pi}_{beta}/pheno_og.csv",
        og_famfile=config["base_fam"]
    output:
        "results/arq_{allele_freq}_{pi}_{beta}/gwas_generation{generation}/optima{optima}/geno.fam",
    params:
        pi=lambda wildcards: str(wildcards.pi),
        beta=lambda wildcards: str(wildcards.beta),
        allele_freq=lambda wildcards: str(wildcards.allele_freq),
        optima=lambda wildcards: str(wildcards.optima),
        generation=lambda wildcards: str(wildcards.generation),
    resources:
        mem_mb=15350,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/generate_fam_files.py"

rule run_gwa:
    input:
        fam_file="results/arq_{allele_freq}_{pi}_{beta}/gwas_generation{generation}/optima{optima}/geno.fam",
        
        bed_file=config['bed_file'],
        bim_file=config['bim_file'],
        kinship=config['kinship']
    output:
        output_gwas="results/arq_{allele_freq}_{pi}_{beta}/gwas_generation{generation}/optima{optima}/output/results.assoc.txt",
    resources:
        mem_mb=30720,
    conda:
        "envs/gwas.yaml"
    script:
        "scripts/gwas_gemma.sh"


#        "results/arq_pi{{pi}}_beta{{beta}}/optima{{optima}}/pop_sizes.csv",
#    log:
#        out="log/arq_pi{pi}_beta{beta}_optima{optima}_stdout.log",
#        err="log/arq_pi{pi}_beta{beta}_optima{optima}_stderr.err",

# rule ld_pruning:
#     input:
#         vcf_file = config['base_vcf'],
#         maf= config[maf]
#         corr_ld= config[corr_ld]
#     output:
#         "data/chr5_grenenet_filteredmaf.prune.in"
# resources:
#     mem_mb=61440
#     log:
#         out = "log/pi{pi}_beta{beta}_run_bash_script_stdout.log",
#         err = "log/pi{pi}_beta{beta}_run_bash_script_stderr.err"
#     conda:
#         'envs/base_env.yaml'
#     script:
#         "ld_pruning.sh"

# rule apply_ld_pruning:
#     input:
#         allele_counts_table = 'results/arq_pi{pi}_beta{beta}/simulation_results/allele_counts.csv',
#         snps_to_keep = 'data/chr5_grenenet_filteredmaf.prune.in'
#     output:
#         allele_counts_table_filtered = 'results/arq_pi{pi}_beta{beta}/simulation_results/allele_counts_filtered.csv',
#     resources:
#         mem_mb=61440
#     log:
#         out = "log/pi{pi}_beta{beta}_run_bash_script_stdout.log",
#         err = "log/pi{pi}_beta{beta}_run_bash_script_stderr.err"
#     conda:
#         'envs/base_env.yaml'
#     script:
#         "apply_ld_pruning.sh"


## for output vcf
# expand(
#    "results/arq_pi{{pi}}_beta{{beta}}/vcf_files/optima{{optima}}/subp{replicates}.vcf",
#    replicates = range(0,12)
# )


## backlog
## this was the begining thing


# TAIR10_chr_all.fas | awk '/^>Chr5 /{p=1}p' > chr5.fasta
# ## then in the python script this file will be modified so it is int he slim format

# ##### phase 1 #######################################
# ## filter ch5 and grenenet ecotypes from the big file

# bcftools view 1001_genomes_snps_missing0.8_merged_imputed_biallelic_named.vcf.gz --regions 5 -o chr5.vcf.gz

# bcftools view -S ecotypes_grenenet.txt chr5.vcf.gz -o chr5_grenenet_nofilt.vcf

## AND NOW BECASUE IF HAVE FILTERED THE SAMPLES I ENED TO ALSO FILTER FOR ONLY VARAIBLE SITES!

# bcftools view --min-ac=1 --no-update chr5_grenenet_nofilt.vcf > chr5_grenenet.vcf


#bcftools view --min-ac=1 --max-ac=449 --no-update chr5_grenenet_nofilt.vcf > chr5_grenenet.vcf
#bcftools view --min-ac=1 --max-ac=449 --no-update chr5_grenenet.vcf > chr5_grenenet_new.vcf

## original vcf
# vcf_og = allel.read_vcf('/home/tbellagio/snakemake_test/data/chr5_grenenet.vcf')
# genotypes = allel.GenotypeArray(vcf_og['calldata/GT'])

# alt_allele_count = genotypes.sum(axis=2).sum(axis=1)

# alt_allele_freq = alt_allele_count / 450

# pos_og = vcf_og['variants/POS']
# og_allele_freq = pd.DataFrame({'pos':pos_og, 'a_freq':alt_allele_freq})
# og_allele_freq['pos'] = og_allele_freq['pos'].astype(int)
# og_allele_freq.to_csv('og_allele_freq.csv')


## maybe add the conversion fo the files from vcf to bim gam with plink adn the calculation fo the kinship with gemma 