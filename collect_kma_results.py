import os
import pandas as pd
import argparse
import re
import subprocess
parser = argparse.ArgumentParser()
parser.add_argument("--kma_outputs", "-i", default = 'fdggfgg', help = "directory containing kma_results sub directories")
parser.add_argument('--output_tsv', '-o', default = 'quast_outputs', help = 'file containing concatenated kma tsv')
args = parser.parse_args()

files = [os.path.join(args.kma_outputs, i) for i in os.listdir(args.kma_outputs) if re.match('.*res', i)]

def read_kma_results(file_):
    tsv = pd.read_csv(file_, sep='\t')
    tsv['fasta_name'] = os.path.basename(file_)
    return tsv
tsvs = [read_kma_results(i) for i in files]
concat_tsv = pd.concat(tsvs, axis=0)
concat_tsv.to_csv(args.output_tsv, sep='\t', index=False)
