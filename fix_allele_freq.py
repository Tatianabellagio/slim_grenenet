import numpy as np
import pandas as pd
df_wmc = ['results/arq_lowfreq_fivepoly_5/lowh/estrongsel/allele_freq.csv',
'results/arq_highfreq_fivepoly_5/highh/exstrongsel/allele_freq.csv',
'results/arq_highfreq_fivepoly_5/mediumh/estrongsel/allele_freq.csv',
'results/arq_mediumfreq_onehpoly_1/highh/exstrongsel/allele_freq.csv',
'results/arq_highfreq_monogen_5/mediumh/estrongsel/allele_freq.csv',
'results/arq_highfreq_fivepoly_3/highh/exstrongsel/allele_freq.csv',
'results/arq_mediumfreq_onehpoly_4/highh/exstrongsel/allele_freq.csv',
'results/arq_lowfreq_twentypoly_2/mediumh/estrongsel/allele_freq.csv',
'results/arq_highfreq_fivepoly_1/mediumh/estrongsel/allele_freq.csv',
'results/arq_mediumfreq_fivepoly_5/mediumh/estrongsel/allele_freq.csv',
'results/arq_lowfreq_twentypoly_4/lowh/estrongsel/allele_freq.csv',
'results/arq_lowfreq_fivepoly_1/mediumh/estrongsel/allele_freq.csv',
'results/arq_highfreq_monogen_1/highh/strongsel/allele_freq.csv',
'results/arq_lowfreq_onehpoly_2/lowh/estrongsel/allele_freq.csv',
'results/arq_highfreq_onehpoly_3/lowh/estrongsel/allele_freq.csv',
'results/arq_mediumfreq_onehpoly_2/highh/exstrongsel/allele_freq.csv']

for df_af in df_wmc:
    df_ac = df_af.replace('allele_freq', 'allele_counts')
    
    allele_counts = pd.read_csv(df_ac,nrows = 5)
    allele_freq = pd.read_csv(df_af)
    missing_col = set(allele_counts.columns) - set(allele_freq.columns) 
    print(missing_col)
    for i in list(missing_col):
        allele_freq[i] = np.nan
        
    allele_freq.to_csv(df_af,index=None)
