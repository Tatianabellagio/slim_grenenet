import pandas as pd
import numpy as np
import os
import glob

def arrange_dataset(files):
    joint = pd.DataFrame()
    for i in files:
        one_file =pd.read_csv(i,header=None)
        one_file.columns = [i]
        joint = pd.concat([joint, one_file],axis=1)
    joint = joint.T
    #joint.index = joint.index.str.split('/subp').str[0]
    #joint_meansubp = joint.reset_index().groupby('index')[joint.columns].mean()
    return joint

path2 = '/home/tbellagi/bigscratch/slim_grenenet/'

# Use the glob module to search for files in subfolders
pattern = os.path.join(path2, '**', '*_pop_size.txt')
pop_size_files = glob.glob(pattern, recursive=True)
pop_size = arrange_dataset(pop_size_files)
pop_size.to_csv('to_transfer/pop_size_oct10.csv')

pattern = os.path.join(path2, '**', '*_va.txt')
va_files = glob.glob(pattern, recursive=True)
va = arrange_dataset(va_files)
va.to_csv('to_transfer/va_oct10.csv')

# Use the glob module to search for files in subfolders
pattern = os.path.join(path2, '**', '*_mfitness.txt')
mfitness_files = glob.glob(pattern, recursive=True)
mfitness = arrange_dataset(mfitness_files)
mfitness.to_csv('to_transfer/mfitnes_oct10.csv')

# Use the glob module to search for files in subfolders
pattern = os.path.join(path2, '**', '*_vfitness.txt')
vfitness_files = glob.glob(pattern, recursive=True)
vfitness = arrange_dataset(vfitness_files)
vfitness.to_csv('to_transfer/vfitnes_oct10.csv')

pattern = os.path.join(path2, '**', '*_vpheno.txt')
vpheno_files = glob.glob(pattern, recursive=True)
vpheno = arrange_dataset(vpheno_files)
vpheno.to_csv('to_transfer/vpheno_oct10.csv')
