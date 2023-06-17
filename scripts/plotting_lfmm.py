import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

## snakemake 

input_freq_counts = snakemake.input['allele_freq_table'] 
input_effect_sizes = snakemake.input['effect_sizes'] 
input_p_values = snakemake.input['p_values_lfmm'] 

output_plot = snakemake.output['plot_lfmm'] 
output_table = snakemake.output['p_values_vs_causal'] 

## snakemake 

print('got files')

allele_freq = pd.read_csv(input_freq_counts, usecols=['pos'])

effect_sizes = pd.read_csv(input_effect_sizes, sep = '\t', header=None,usecols=[2,3])



### check if pvalues fiel is too small is because most populations deid so there is not eoguh data to run the model 
if os.path.getsize(input_p_values) < 10:
    ## create fake files for snakemake 
    plt.figure()
    plt.savefig(output_plot)
    empty_df = pd.DataFrame()
    empty_df.to_csv(output_table, index=False)
else:
    pvalues = pd.read_csv(input_p_values)
    effect_sizes.columns=['pos', 'effect_size']

    pvalues.columns = ['pvalues']

    pvalues_pos = pd.concat([pvalues, allele_freq],axis=1)

    pvalues_pos = pd.merge(pvalues_pos, effect_sizes,how='left')

    pvalues_pos.to_csv(output_table, index=False)

    pvalues_pos["causal"] = pvalues_pos['effect_size']!=0

    pvalues_pos['pvalues'] = -np.log10(pvalues_pos['pvalues'])

    causal_pos = pvalues_pos[pvalues_pos['effect_size']!=0]
    
    plt.figure(figsize=(12, 4))
    plt.scatter(pvalues_pos['pos'], pvalues_pos['pvalues'], color='grey', s=2)
    plt.xlabel("SNP")
    plt.ylabel("-Log P")
    plt.plot(causal_pos['pos'], causal_pos['pvalues'], ".", markersize=7, color='red')
    plt.savefig(output_plot)
    plt.close()  # Close the figure without displaying it