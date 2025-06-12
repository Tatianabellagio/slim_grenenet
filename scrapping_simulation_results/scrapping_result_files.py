import pandas as pd
import os
import glob
import multiprocessing

def generate_dataset(args):
    file_type, pat, final_name, path2 = args
    pattern = os.path.join(path2, pat, file_type)
    print(pattern)


    files = glob.glob(pattern, recursive=True)

    joint = dict()
    for i in files:
        one_file = pd.read_csv(i, header=None)
        joint[i] = one_file
    joint = pd.concat(joint, axis=1).T

    joint.to_csv(final_name)
    return None

if __name__ == "__main__":
    # Define the base path for your datasets 
    path = '/global/scratch/users/tbellg/slim_grenenet/results/'
    # Prepare arguments for each dataset generation
    tasks = [
       ('*_pop_size_early.txt', '**/', 'es_dep_af_drift/pop_size_early.csv', path),
       ('*_va.txt', '**/', 'es_dep_af_drift/va_ve.csv', path),
       ('*_mfitness.txt', '**/', 'es_dep_af_drift/mfitnes.csv', path),
       ('*_vfitness.txt', '**/', 'es_dep_af_drift/vfitnes.csv', path),
       ('*_mpheno.txt', '**/', 'es_dep_af_drift/mpheno.csv', path),
       ('*_vpheno.txt', '**/', 'es_dep_af_drift/vpheno.csv', path),
       #('*_adj_variance.txt', '**/', 'adj_variance.csv', path),
       #('*_new_optimum.txt', '**/', '_new_optimum.csv', path),
       ('*_minphenotype.txt', '**/', 'es_dep_af_drift/minphenotype.csv', path),
       ('*_maxphenotype.txt', '**/', 'es_dep_af_drift/maxphenotype.csv', path),
    ]

    # Number of CPU cores to use (adjust as needed)
    num_cores = 8

    # Create a pool of worker processes
    with multiprocessing.Pool(processes=num_cores) as pool:
        # Use the pool to parallelize the processing of tasks
        pool.map(generate_dataset, tasks)