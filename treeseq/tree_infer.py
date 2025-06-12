import sys
import tsinfer
import time
import os

input_filename = sys.argv[1]
output_filename = sys.argv[2]

print(input_filename)

start_time = time.time()

samplefile = tsinfer.load(input_filename) ## 
ts = tsinfer.infer(samplefile,num_threads=20)
ts.dump(output_filename)

end_time = time.time()
elapsed_time = end_time - start_time

# Extract the base name of the input file without the extension
input_basename = os.path.splitext(os.path.basename(input_filename))[0]
time_filename = f"{input_basename}_runtime.txt"

with open(time_filename, "w") as time_file:
    time_file.write(f"Elapsed time: {elapsed_time:.2f} seconds\n")
