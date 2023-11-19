import subprocess

# List of batch sizes to run
batch_sizes = [1,2,3,4]

# Loop through batch sizes
for batch_size in batch_sizes:
    command = f"snakemake --rerun-triggers mtime --use-conda --conda-frontend mamba --profile profiles/slurm/ --scheduler greedy --batch all={batch_size}/4"
    subprocess.run(command, shell=True)
