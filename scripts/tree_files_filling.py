import os

filename = snakemake.input['pheno_file'] 


output_gen3_trees = filename.replace('_st_phenov.txt', '_tree_output_gen3.trees')
output_gen10_trees = filename.replace('_st_phenov.txt', '_tree_output_gen10.trees')

def create_empty_file_if_not_exists(filename):
    if not os.path.exists(filename):
        # Create an empty file
        with open(filename, 'w') as file:
            pass

create_empty_file_if_not_exists(output_gen3_trees)
create_empty_file_if_not_exists(output_gen10_trees)