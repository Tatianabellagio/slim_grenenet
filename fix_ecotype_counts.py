import os
import pandas as pd
import glob
pattern = os.path.join('results', '**', 'ecotype_counts.csv')
ecotype_counts_files = glob.glob(pattern, recursive=True)

to_fix = []
for i in ecotype_counts_files:
    df = pd.read_csv(i,  nrows =1).drop(['Unnamed: 0', 'ecotype'],axis=1)
    length_col_names = [len(i) for i in df.columns]
    if 7 in length_col_names:
        to_fix.append(i)
        print(length_col_names)


for df_path in to_fix:
    df = pd.read_csv(df_path).drop('Unnamed: 0',axis=1)
    right_col_name = []
    for i in df.columns:
        if i == 'ecotype':
            right_col_name.append(i)
        elif len(i)< 10:
            i = 'optima' + i 
            right_col_name.append(i)
        else:
            i = 'optima' + i.split('optima')[1].split('_vcf')[0].replace('/', '_')
            right_col_name.append(i)
    df.columns = right_col_name
    print(df.head(5))
    df.to_csv(df_path)
