import pandas as pd
import tskit
import allel
import random
import numpy as np
import tsinfer
import pyslim

og_vcf_offset = snakemake.input['og_vcf_offset'] 
og_tree_offset = snakemake.input['og_tree_offset'] 
beta = int(snakemake.params['beta'])

polygenicity_params_file = snakemake.input['polygenicty_params'] 

pi_option =  snakemake.params['pi']
pi = pd.read_csv(polygenicity_params_file,header=None, usecols=[int(pi_option)]).values[0][0]

#get the actual values
output_tree_seq_causalloci = snakemake.output["tree_seq_causalloci"]
output_loci_effectsize = snakemake.output["loci_effectsize"]
output_phenotypes = snakemake.output["phenotypes"]

def calc_pos_sc(alt_al_per_pos, pos, n_ecotypes, pi, beta):
    # calculate the total number of alternative alleles at each position
    alt_al_count = alt_al_per_pos.sum(axis=1)
    # create a dataframe with the positions and the frequency of the alt allele at each posiiton 
    alelle_dist = pd.DataFrame({'alt_al_count':alt_al_count, 'pos':pos})
    alelle_dist['alt_al_freq'] = alelle_dist['alt_al_count'] / (n_ecotypes*2)
    # randomly sample the causal snps based on pi
    alelle_dist = alelle_dist.sample(pi)
    # randomly draw the effect sizes 
    sc = np.random.normal(0, beta, pi)
    alelle_dist['sc'] = sc
    return alelle_dist

def calc_phenotypes(pos,pos_sc, alt_al_per_pos):
    # create a mask including only the contributing loci 
    mask_positions = pd.Series(pos).isin(pos_sc['pos'])
    # use this mask to mask the genotype matrix
    alt_al_per_pos_selected_sites = alt_al_per_pos[mask_positions]
    # calculate the phenotypes by multiplying the genotype matrix by the effect sizes at each position 
    phenotypes = []
    for i in range(alt_al_per_pos_selected_sites.shape[1]):
        gen_effectsize = np.multiply(alt_al_per_pos_selected_sites[:, i] , pos_sc['sc'])
        phenotypes.append(gen_effectsize.sum())
    return phenotypes

def keep_only_causal_sites_and_mutations(og_tree_offset, pos_sc):
    ts = tskit.load(og_tree_offset)

    # dumpt the tables from the tree
    tables = ts.dump_tables()

    ## extract all teh sites from the og tree
    complete_sites = pd.Series(tables.sites.position)

    # create a mask to filter only the ones present in the selected sites (causal loci)
    mask_delete_sites = complete_sites.isin(pos_sc['pos'])

    ## replace the table only with the causal sites, and same for mutation tables
    tables.sites.replace_with(tables.sites[mask_delete_sites])
    tables.mutations.replace_with(tables.mutations[mask_delete_sites])
    ## extract the new site index
    tables.mutations.site = np.array(range(0, len(tables.mutations))).astype('int32')

    ## ge tthe positions and sc in teh right order 
    pos_table = pd.Series(tables.sites.position).reset_index()
    right_order_pos = pos_sc.merge(pos_table, left_on='pos',right_on =0).sort_values('index')

    ## create the tree to then modify it 
    new_ts = tables.tree_sequence()
    tables = new_ts.dump_tables()

    ## chance the ancestral state to empty or slim will complain
    tables.sites.clear()
    for s in new_ts.sites():
        tables.sites.append(s.replace(ancestral_state=""))

    ## add the selection coefficient and the rigth emtadata fro slim 
    tables.mutations.clear()
    for k, (m, sc) in enumerate(zip(new_ts.mutations(), right_order_pos['sc'])):
        mm = pyslim.default_slim_metadata('mutation_list_entry')
        mm['selection_coeff'] = sc
        tables.mutations.append(
            m.replace(derived_state=str(k), metadata={'mutation_list': [mm]}))
        
    return tables.tree_sequence()

# read the original vcf with 'offset' positions so that all the chromosomes are 'like one long chromosome'
vcf_og = allel.read_vcf(og_vcf_offset, fields=["calldata/GT", 'variants/POS' , 'samples'])
# extract the genotypes matrix from the vcf
geno_og = vcf_og["calldata/GT"]
# extract the sample names to calculate the number of different ectoypes present in the vcf file
samples = vcf_og['samples']
n_ecotypes = len(vcf_og['samples'])
# extract the position of all the snps
pos = vcf_og['variants/POS']
# convert the genotype matrix of 0 and 1 into an alternative alleles matrix with 0,1 and 2 based on the number of copies 
# in each position 
alt_al_per_pos = geno_og.sum(axis=2) 

# get the causal loci
pos_sc = calc_pos_sc(alt_al_per_pos, pos, n_ecotypes , pi, beta)
# calculate the phenotypes for each individual
phenotypes = calc_phenotypes(pos,pos_sc, alt_al_per_pos)

## save phenotypes and details about causal loci 
pd.Series(phenotypes).to_csv(output_phenotypes)
pos_sc.to_csv(output_loci_effectsize)

### filter tree to only keep causal mutations 
pre_slim_tree = keep_only_causal_sites_and_mutations(og_tree_offset, pos_sc)

## save tree
pre_slim_tree.dump(output_tree_seq_causalloci)