import subprocess

# List of batch sizes to run
batch_sizes = [1,2,3,4,5,6, 7,8,9,10]

# Loop through batch sizes
for batch_size in batch_sizes:
    command = f"snakemake --rerun-triggers mtime --profile profiles/slurm/ --latency-wait 20 --batch all={batch_size}/10"
    subprocess.run(command, shell=True)
