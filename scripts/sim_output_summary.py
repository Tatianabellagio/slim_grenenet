import pandas as pd
import numpy as np
import os
import glob

path2 = '/home/tbellagi/bigscratch/slim_grenenet/results/'

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

def generatefile(file_type, pat, final_name):
    pattern = os.path.join(path2, pat, file_type)
    files = glob.glob(pattern, recursive=True)
    dataset = arrange_dataset(files)
    dataset.to_csv(final_name)
    return None

generatefile('*_pop_size.txt', '**/','to_transfer/pop_size_oct13.csv')
generatefile('*_va.txt', '**/','to_transfer/va_oct13.csv')
generatefile('*_mfitness.txt','**/','to_transfer/mfitnes_oct13.csv')
generatefile('*_vfitness.txt', '**/', 'to_transfer/vfitnes_oct13.csv')
generatefile('*_vpheno.txt', '**/', 'to_transfer/vpheno_oct13.csv')
