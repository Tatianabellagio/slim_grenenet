import pandas as pd
##params 
number_of_env = 32
grenenet_1001g_ecotypes = 'ecotypes_grenenet_1001g.txt'
path_worldclim_ecotypesdata = 'worldclim_ecotypesdata.csv'
## import all the ecotypes that are in the 1001g and grenenet 
ecotypes_grenenet_1001g = pd.read_csv(grenenet_1001g_ecotypes).columns.astype(int)
## get the mean temp of their place of origin 
ecotypes_temp = pd.read_csv(path_worldclim_ecotypesdata, usecols=['ecotypeid', 'bio1'])
#BIO1 = Annual Mean Temperature
## filter by the ones coming from the 1001 genomes
ecotypes_temp = ecotypes_temp[ecotypes_temp['ecotypeid'].isin(ecotypes_grenenet_1001g)]
### get the min and max values just for completeness of the range and then sample the rest 
# Find the minimum and maximum elements
ecot_max = ecotypes_temp[ecotypes_temp['bio1'] == ecotypes_temp['bio1'].max()]
ecot_min = ecotypes_temp[ecotypes_temp['bio1'] == ecotypes_temp['bio1'].min()]
# Remove the minimum and maximum elements from the dataframe 
ecotypes_temp = ecotypes_temp[~ecotypes_temp['bio1'].isin([ecotypes_temp['bio1'].max(),ecotypes_temp['bio1'].min()])].copy()

## and then subsample the rest 
samples = ecotypes_temp.sample(number_of_env-2)
## concat all 
samples = pd.concat([samples,ecot_max,ecot_min],axis=0).reset_index(drop=True)

samples.to_csv('optimum_ecotypes.csv')