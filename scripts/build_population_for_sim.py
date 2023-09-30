import pandas as pd
import tskit
import allel
import random
import numpy as np
import tsinfer
import pyslim

og_vcf_offset = snakemake.input['og_vcf_offset'] 
og_tree_offset = snakemake.input['og_tree_offset'] 

pi_option =  snakemake.params['pi']
pi = int(snakemake.params[pi_option])

beta = int(snakemake.params['beta'])

alelle_freq_option = snakemake.params['allele_freq']
allele_freq = snakemake.params[alelle_freq_option]
lower_bound = float(allele_freq[0])
upper_bound = float(allele_freq[1])


#get the actual values
output_tree_seq_causalloci = snakemake.output["tree_seq_causalloci"]
output_loci_effectsize = snakemake.output["loci_effectsize"]
output_phenotypes = snakemake.output["phenotypes"]

def calc_pos_sc(alt_al_per_pos, pos, n_ecotypes, allele_freq, pi, beta):
    alt_al_count = alt_al_per_pos.sum(axis=1)
    alelle_dist = pd.DataFrame({'alt_al_count':alt_al_count, 'pos':pos})
    alelle_dist['alt_al_freq'] = alelle_dist['alt_al_count'] / (n_ecotypes*2)
    sim_freq_pos = alelle_dist[(alelle_dist['alt_al_freq'] < upper_bound) & (alelle_dist['alt_al_freq'] >= lower_bound)]['pos']
    selected_sites = sim_freq_pos.sample(pi).values
    print(selected_sites)
    sc = np.random.normal(0, beta, pi)
    pos_sc = pd.DataFrame({'pos': selected_sites, 'sc': sc})
    return pos_sc

def calc_phenotypes_st(pos,pos_sc, alt_al_per_pos):
    mask_positions = pd.Series(pos).isin(pos_sc['pos'])
    alt_al_per_pos_selected_sites = alt_al_per_pos[mask_positions]
    phenotypes = []
    for i in range(alt_al_per_pos_selected_sites.shape[1]):
        gen_effectsize = np.multiply(alt_al_per_pos_selected_sites[:, i] , pos_sc['sc'])
        phenotypes.append(gen_effectsize.sum())
    ## calcualte mean and std of phenotypes to standarize them 
    pheno_mean = np.array(phenotypes).mean()
    pheno_std = np.array(phenotypes).std()

    phenotypes_st = (phenotypes - pheno_mean)/ pheno_std
    return phenotypes_st

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


## for this im gonna use the og vcf file wth the offset to be able to map the positions correctly 
vcf_og = allel.read_vcf(og_vcf_offset, fields=["calldata/GT", 'variants/POS' , 'samples'])
geno_og = vcf_og["calldata/GT"]
samples = vcf_og['samples']
pos = vcf_og['variants/POS']

n_ecotypes = len(vcf_og['samples'])
alt_al_per_pos = geno_og.sum(axis=2) 

pos_sc = calc_pos_sc(alt_al_per_pos, pos, n_ecotypes, allele_freq, pi, beta)

phenotypes = calc_phenotypes_st(pos,pos_sc, alt_al_per_pos)

## save 

pd.Series(phenotypes).to_csv(output_phenotypes)
pos_sc.to_csv(output_loci_effectsize)

### filter tree

pre_slim_tree = keep_only_causal_sites_and_mutations(og_tree_offset, pos_sc)

## save tree

pre_slim_tree.dump(output_tree_seq_causalloci)