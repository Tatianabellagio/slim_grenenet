import glob
import pandas as pd
import os
import tskit

# Use the glob module to search for files in subfolders
pattern = os.path.join(path, '**', 'allele_freq.csv')
trees = glob.glob(pattern, recursive=True)

path = '/home/tbellagio/scratch/slim_grenenet/'


# Use the glob module to search for files in subfolders
pattern = os.path.join(path, '**', '*_tree_output.trees')
trees = glob.glob(pattern, recursive=True)



pop_size = dict()
for i in trees:
    print(i)
    if os.path.getsize(i) > 1:
        tree = tskit.load(i)
        pop_s = len(tree.tables.individuals)
        pop_size[i] = pop_s
    else:
        pop_size[i] = 0

pd.DataFrame([pop_size]).T.to_csv('pop_sizes.csv')