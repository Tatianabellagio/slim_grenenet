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
print(pi)

beta_option = snakemake.params['beta']
beta = int(snakemake.params[beta_option])
print(beta)

alelle_freq_option = snakemake.params['allele_freq']
allele_freq = snakemake.params[alelle_freq_option]
lower_bound = float(allele_freq[0])
upper_bound = float(allele_freq[1])


#get the actual values
optima_qty = str(snakemake.params['optima_qty']) 

output_tree_seq_causalloci = snakemake.output["tree_seq_causalloci"]
output_loci_effectsize = snakemake.output["loci_effectsize"]
output_phenotypes = snakemake.output["phenotypes"]
output_optima_values = snakemake.output["optima_values"]
output_variance_values = snakemake.output["variance_values"]

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

def calc_phenotypes(pos,pos_sc, alt_al_per_pos):
    mask_positions = pd.Series(pos).isin(pos_sc['pos'])
    alt_al_per_pos_selected_sites = alt_al_per_pos[mask_positions]
    phenotypes = []
    for i in range(alt_al_per_pos_selected_sites.shape[1]):
        gen_effectsize = np.multiply(alt_al_per_pos_selected_sites[:, i] , pos_sc['sc'])
        phenotypes.append(gen_effectsize.sum())
    return phenotypes

def calc_optima(phenotypes):
    max_pheno = max(phenotypes)
    min_pheno = min(phenotypes)

    length = max_pheno - min_pheno
    step = length/(int(optima_qty) - 1)
    optima = [round(min_pheno + i * step, 4) for i in range(0, int(optima_qty))]
    return optima

def calc_variances(phenotypes, optima):
    range_pheno =  max(phenotypes) - min(phenotypes)
    dist_between_env = range_pheno / len(optima)
    ## Strong selection, fitness 0 in the adyacent environemnt
    sd1 = dist_between_env / 3
    variance1 = sd1**2
    ## moderate selection, fitness 0 in the other extremee environemnt 
    sd2 = (dist_between_env * 4) / 3 ## 3 sd will be in between 4  environments 
    variance2 = sd2**2
    ## weak selection, half fitness in the other extreme environment 
    sd3 = (dist_between_env * 8) / 3
    variance3 = sd3**2
    variances = [sd1, sd2, sd3]
    return variances

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

phenotypes = calc_phenotypes(pos,pos_sc, alt_al_per_pos)

optima = calc_optima(phenotypes)

variances = calc_variances(phenotypes, optima)

variances

pd.DataFrame(index =  ['strongsel','moderatesel','lowsel'],data = { 'var': variances}).to_csv('sel_var.csv')

pd.DataFrame(data = {'selection': ['strongsel','moderatesel','lowsel'], 'var': variances}).to_csv('sel_var.csv')

## save 

pd.Series(phenotypes).to_csv(output_phenotypes)
pos_sc.to_csv(output_loci_effectsize)

with open(output_optima_values, 'w') as file:
    for element in optima:
        file.write(str(element) + '\n')  # Write element followed by a newline

with open(output_variance_values, 'w') as file:
    for element in variances:
        file.write(str(round(element,4)) + '\n')  # Write element followed by a newline

### filter tree

pre_slim_tree = keep_only_causal_sites_and_mutations(og_tree_offset, pos_sc)

## save tree

pre_slim_tree.dump(output_tree_seq_causalloci)