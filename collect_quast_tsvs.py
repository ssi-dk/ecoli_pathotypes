import os
import pandas as pd
import argparse
import re
import subprocess
parser = argparse.ArgumentParser()
parser.add_argument("--input_quast_dir", "-i", default = 'fdggfgg', help = "directory containing quast sub directories")
parser.add_argument('--output_tsv', '-o', default = 'quast_outputs', help = 'file containing concatenated quast tsvs')
args = parser.parse_args()

dirs = [os.path.join(args.input_quast_dir, i) for i in os.listdir(args.input_quast_dir) if re.match('quast.*', i)]

tsvs = [pd.read_csv(f'{i}/transposed_report.tsv', sep='\t') for i in dirs]
concat_tsv = pd.concat(tsvs, axis=0)
concat_tsv.to_csv(args.output_tsv, sep='\t', index=False)
