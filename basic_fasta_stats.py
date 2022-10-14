import os
import pandas as pd
import argparse
import re
import subprocess
parser = argparse.ArgumentParser()
parser.add_argument("--input_fasta_dir", "-i", default = 'fdggfgg', help = "directory containing fasta files")
parser.add_argument('--output_dir', '-o', default = 'quast_outputs', help = 'directory to put the quast qc outputs')
args = parser.parse_args()

if not os.path.isdir(args.output_dir):
    os.makedirs(args.output_dir)
fastas = [os.path.join(args.input_fasta_dir, i) for i in os.listdir(args.input_fasta_dir) if re.match('.*.fasta', i)]

step_size = 1000
# version of quast used: QUAST v5.0.0

for i in range(0, len(fastas), step_size):
    idx_start = i
    idx_end = min(i+step_size, len(fastas))
    fasta_subset = fastas[idx_start:idx_end]
    fasta_string = ' '.join(fasta_subset)
    quast_cmd = f'quast.py -t 10 -o quast_{idx_start}_{idx_end} -m 0 --no-plots --no-html {fasta_string}'
    sbatch_cmd = f"sbatch -D {args.output_dir} -c 10 --mem=8G -J quast_on_ecoli_{idx_start}_{idx_end} -p covid --wrap=\'{quast_cmd}\'"
    process  = subprocess.Popen(sbatch_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, env=os.environ, encoding='utf-8')
    process_out, process_err = process.communicate()
    print(process_out, process_err)
