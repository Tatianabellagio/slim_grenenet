import pandas as pd
import glob
import os

pattern = os.path.join('results', '**', 'ecotype_counts.csv')
ecotype_counts_files = glob.glob(pattern, recursive=True)
print(ecotype_counts_files)
for i in ecotype_counts_files:
    ec = pd.read_csv(i,nrows=1)
    if ec.shape[1] < 100:
        print(ec.shape)
