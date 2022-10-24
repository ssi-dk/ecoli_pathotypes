import os
import pandas as pd
import argparse
import re
import subprocess
parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", "-i", default = 'fdggfgg', help = "directory containing fastas")
parser.add_argument('--kma_executable', '-k', default = 'kma', help = 'location of the kma executable')
parser.add_argument("--kma_db", '-db', default='kma_test_db/genes_as_db', help = 'which indexed kma db to use')
parser.add_argument('--output_dir', '-o', default = 'quast_outputs', help = 'output dir for putting kma output')
parser.add_argument('--chunk_size', '-chunk', default = 1000, help = 'chunk_size', type=int)
args = parser.parse_args()

if not os.path.isdir(args.output_dir):
    os.makedirs(args.output_dir)

fastas = [os.path.join(args.input_dir, i) for i in os.listdir(args.input_dir) if re.match('.*.fasta', i)]
chunk_size = args.chunk_size
chunks = list(range(0, len(fastas)+chunk_size, chunk_size))

#python run_kma_with_genes_as_db.py -i fasta_files -k $repos/kma/kma -db kma_test_db/genes_as_db -o test_1000_samples -j 0,10

for i in range(len(chunks)-1):
    idxs = [chunks[i], min(chunks[i+1], len(fastas))]
    idxs_str = ','.join([str(i) for i in idxs])
    cmd = f'python run_kma_with_genes_as_db.py -i {args.input_dir} -k {args.kma_executable} -db {args.kma_db} -o {args.output_dir} -j {idxs_str}'
    slurm_cmd = f"sbatch -D . -c 2 --mem=8G -J kma_on_fastas_{idxs[0]}_{idxs[1]} -p covid --wrap=\'{cmd}\'"
    print(slurm_cmd)
    process = subprocess.Popen(slurm_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, env=os.environ, encoding='utf-8')
    process_out, process_err = process.communicate()
    print(process_out, process_err)

