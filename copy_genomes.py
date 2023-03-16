import pandas as pd
import os
import re
import functools
import yaml
import traceback
import shutil
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", "-i", default = 'fdggfgg', help = "directory containing fastas")
parser.add_argument('--meta_tsv', '-tsv', default = 'dsfdf', help = 'emtadata file containing an Assembly barcode field')
parser.add_argument('--output_dir', '-o', default = 'some_output', help = 'output dir containing the copied files')
args = parser.parse_args()

meta = pd.read_csv(args.meta_tsv, sep='\t')
outdir = args.output_dir
if not os.path.isdir(outdir):
    os.makedirs(outdir)

assemblies = meta['Assembly barcode'].tolist()
#assemblies = assemblies[0:10]
for i in assemblies:
    fasta_path = os.path.join(args.input_dir, i + '.fasta')
    out_path = os.path.join(outdir, i + '.fasta')
    #print(fasta_path, out_path)
    shutil.copyfile(fasta_path, out_path)
