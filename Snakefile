# =================================================================================================
#     Dependencies
# =================================================================================================


configfile: "config.yaml"


## this rule runs a python script that will generate the bed file containing the contributing loci and their effect sizes based on values of dn alpha
## the bed file will be then used to annotate a vcf file that will be used by SliM to run the simulations


rule all:
    input:
        expand(
            'results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/allele_freq.csv',
            allele_freq=config['allele_freq'],
            pi=config["pi"],
            selection=config["selection"],
            heritability=config["heritability"],
            replicates_arq=config["replicates_arq"],
        ),

rule build_population_for_sim:
    input:
        og_tree_offset=config["og_tree_offset"],
        og_vcf_offset=config["og_vcf_offset"],
    output:
        tree_seq_causalloci="results/arq_{allele_freq}_{pi}_{replicates_arq}/tree_seq_causalloci.trees",
        loci_effectsize="results/arq_{allele_freq}_{pi}_{replicates_arq}/loci_effectsize.csv",
        phenotypes="results/arq_{allele_freq}_{pi}_{replicates_arq}/phenotypes.csv",
    params:
        pi=lambda wildcards: str(wildcards.pi),
        replicates_arq=lambda wildcards: str(wildcards.replicates_arq),
        allele_freq=lambda wildcards: str(wildcards.allele_freq),
        lowfreq=config["lowfreq"],
        mediumfreq=config["mediumfreq"],
        highfreq=config["highfreq"],
        monogen=config["monogen"],
        fivepoly=config["fivepoly"],
        twentypoly=config["twentypoly"],
        onehpoly=config["onehpoly"],
        beta=config["beta"],
    resources:
        mem_mb=30720,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/build_population_for_sim.py"

rule run_slim_simulation:
    input:
        tree_seq_causalloci="results/arq_{allele_freq}_{pi}_{replicates_arq}/tree_seq_causalloci.trees",
    output: 
        output_tree_gen4="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_tree_output_gen4.trees",
        output_tree_gen10="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_tree_output_gen10.trees",
        output_pop_size="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_pop_size.txt",
        output_va="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_va.txt",
        output_vpheno="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_vpheno.txt",
        output_mfitness="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_mfitness.txt",
        output_vfitness="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_vfitness.txt",

    params:
        optima=lambda wildcards: str(wildcards.optima),        
        selection=lambda wildcards: str(wildcards.selection),
        heritability=lambda wildcards: str(wildcards.heritability),

    resources:
        mem_mb=40960,
    conda:
        "envs/base_env.yaml"
    shell:
        "scripts/slim.sh {input} {params} {output}"

rule tree_postprocessing:
    input:
        og_tree_offset=config["og_tree_offset"],
        mapper_ids=config['mapper_realid_metadataid'],
        output_sim_tree="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_tree_output_gen4.trees",
    output:
        output_sim_tree_wm ="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_tree_output_gen4_wm.trees",
        output_vcf ="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_vcfgen4_output.vcf",
    resources:
        mem_mb=30720,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/tree_postprocessing.py"

rule fix_positions_vcf:
    input:
        output_vcf ="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_vcfgen4_output.vcf",
    output: 
        output_vcf_fixpos ="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_vcfgen4_output_rp.vcf",
    resources:
        mem_mb=10240,
    conda:
        "envs/base_env.yaml"
    shell:
        "scripts/fix_positions.sh {input} {output}"

rule gen_allele_freq:
    input:
        pos_vcf_og=config['pos_vcf_og'],
        output_vcf_fixpos = expand(
            "results/arq_{{allele_freq}}_{{pi}}_{{replicates_arq}}/{{heritability}}/{{selection}}/optima{optima}/subp{replicates_sim}_vcfgen4_output_rp.vcf",
            optima=config["optima"],
            replicates_sim=config["replicates_sim"],    
        ),
    output:
        allele_freq ="results/arq_{allele_freq}_{pi}_{replicates_arq}/{heritability}/{selection}/allele_freq.csv",
    resources:
        mem_mb=30720,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/allele_freq_calc.py"
