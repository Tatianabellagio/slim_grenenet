print('reading')
import pandas as pd
import tskit
import allel
import random
import numpy as np
import tsinfer
import pyslim
import os

og_tree_offset = snakemake.input['og_tree_offset'] 
mapper_realid_metadataid = snakemake.input['mapper_ids'] 
output_sim_tree = snakemake.input['output_sim_tree'] 
output_sim_tree_wm = snakemake.output['output_sim_tree_wm'] 
output_vcf = snakemake.output['output_vcf'] 

def overlap_neutral_mut (ts_new, ts, mapper_realid_metadataid):
    ## extract surviving ndoes and comapre them to our old ndoes to place mtuations in the right place
    surviving_nodes = []
    for i in ts_new.tables.nodes:
        surviving_nodes.append(i.metadata['slim_id'])
    ## new nodes id and the ids i gave them in the past
    new_mapper = pd.DataFrame({'new_ids': range(0, len(ts_new.tables.nodes)), 'my_ids_metadata':surviving_nodes})
    ## map old nodes with new nodes
    mapper_lost_nodes = new_mapper.merge(mapper_realid_metadataid, left_on = 'my_ids_metadata', right_on = 'my_ids_metadata', how= 'right')

    ## create a mask to only keep from the old nodes the ones that survived the simulation
    mask = mapper_lost_nodes['new_ids'].notna()

    tables_og = ts.dump_tables()

    ## now filter old tables only based on surviving nodes 
    tables_og.nodes.replace_with(tables_og.nodes[mask])

    ## now filter mutation table based on the surviving nodes, for that, extract the nodes 
    old_nodes = tables_og.mutations.node

    old_nodes = pd.Series(old_nodes)

    old_nodes.name = 'old_nodes'

    ## create a dataframe relating the new and old nodes
    replace_oldbynew_nodes = pd.merge(old_nodes, mapper_lost_nodes, left_on ='old_nodes', right_on = 'real_id', how= 'left')

    ## create a mask to filter out all the mutations than has been lost 
    mask_mutations_lost = replace_oldbynew_nodes['new_ids'].notna()

    ## filter out mutations that has been lost 
    table_mutations = tables_og.mutations[mask_mutations_lost]

    ## now replace the old nodes ids by the new nodes ids with the mapper
    ids_to_replace = replace_oldbynew_nodes.dropna()['new_ids']
    table_mutations.node = np.array(ids_to_replace.astype('int32'))

    ## and jsut set the sites from 0 to the length of mutation table 
    table_mutations.site = np.array(range(0, len(table_mutations))).astype('int32')

    ## apply the same filter from the mutations table to the sites table 
    table_sites = tables_og.sites[mask_mutations_lost]  

    ## now replace all this filter old tables in the new tree seq! 
    new_tables = ts_new.dump_tables()

    new_tables.mutations.replace_with(table_mutations)

    new_tables.sites.replace_with(table_sites)

    ## make sure to compure mutations parents
    new_tables.compute_mutation_parents()

    ## create tree seq based on tables
    tree_nm = new_tables.tree_sequence()

    return tree_nm.simplify()

def convert_tree_to_vcf (tree,name_vcf):
    # create a vcf file from the treeseq 
    with open(name_vcf, 'w') as file:
        # Pass the file object as the output parameter
        tree.write_vcf(output=file)

#import the old tree
ts_old = tskit.load(og_tree_offset)
#import mapper old nodes to new nodes
mapper_realid_metadataid = pd.read_csv(mapper_realid_metadataid)

if os.path.exists(output_sim_tree) and os.path.getsize(output_sim_tree) <= 1:
    print('empty_tree')
    with open(output_vcf, "w"):
        pass  # Create an empty vcf file 
    with open(output_sim_tree_wm, "w"):
        pass  # Create an empty tree file 
elif os.path.exists(output_sim_tree) and os.path.getsize(output_sim_tree) > 1:
    ts_new = tskit.load(output_sim_tree)
    ts_nm = overlap_neutral_mut(ts_new, ts_old, mapper_realid_metadataid)
    ts_nm.dump(output_sim_tree_wm)
    convert_tree_to_vcf(ts_nm, output_vcf)