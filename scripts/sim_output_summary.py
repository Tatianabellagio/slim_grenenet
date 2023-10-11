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

def generatefile(file_type, part1_pat, part2_pat, final_name):
    pattern = os.path.join(path2, part1_pat, file_type)
    files1 = glob.glob(pattern, recursive=True)
    pattern = os.path.join(path2, part2_pat, file_type)
    files2 = glob.glob(pattern, recursive=True)
    files = files1 + files2
    dataset = arrange_dataset(files)
    dataset.to_csv(final_name)  
    return None


generatefile('*_pop_size.txt', '**/**/estrongsel/**/','**/**/vstrongsel/**/', 'to_transfer/pop_size_oct11.csv')
generatefile('*_va.txt', '**/**/estrongsel/**/','**/**/vstrongsel/**/', 'to_transfer/va_oct11.csv')
generatefile('*_mfitness.txt', '**/**/estrongsel/**/','**/**/vstrongsel/**/', 'to_transfer/mfitnes_oct11.csv')
generatefile('*_vfitness.txt', '**/**/estrongsel/**/','**/**/vstrongsel/**/', 'to_transfer/vfitnes_oct11.csv')
generatefile('*_vpheno.txt', '**/**/estrongsel/**/','**/**/vstrongsel/**/', 'to_transfer/vpheno_oct11.csv')
