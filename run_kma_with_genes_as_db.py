import os
import pandas as pd
import argparse
import re
import subprocess
parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", "-i", default = 'fdggfgg', help = "directory containing fastas")
parser.add_argument("--fasta_idxs", "-j", default = None, help = 'indexes to slice fasta files, comma separated with two values')
parser.add_argument('--kma_executable', '-k', default = 'kma', help = 'location of the kma executable')
parser.add_argument("--kma_db", '-db', default='kma_test_db/genes_as_db', help = 'which indexed kma db to use')
parser.add_argument('--output_dir', '-o', default = 'quast_outputs', help = 'output dir for putting kma output')
args = parser.parse_args()

if not os.path.isdir(args.output_dir):
    os.makedirs(args.output_dir)

fastas = [os.path.join(args.input_dir, i) for i in os.listdir(args.input_dir) if re.match('.*.fasta', i)]

if args.fasta_idxs != None:
    idxs = [int(i) for i in args.fasta_idxs.split(',')]
    fastas = fastas[idxs[0]:idxs[1]]
print(len(fastas))
for i in fastas:
    file_name = os.path.basename(i)
    cmd = f'{args.kma_executable} -i {i} -o {args.output_dir}/{file_name} -t_db {args.kma_db} -ca  -nf -nc'
    #cmd = f'{args.kma_executable} -i {i} -o {args.output_dir}/{file_name} -t_db {args.kma_db} -nf -nc'
    #print(cmd)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, env=os.environ, encoding='utf-8')
    process_out, process_err = process.communicate()
    #print(process_out, process_err)
