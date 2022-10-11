from urllib.request import urlopen
from urllib.error import HTTPError
import urllib
import base64
import json
import os
import pandas as pd
import argparse
import logging

parser = argparse.ArgumentParser()
parser.add_argument("--input_fasta_links", "-i", default = 'fdggfgg', help = "tsv file that has a field called download_fasta_link")
parser.add_argument('--output_dir', '-o', default = 'fasta_files', help = 'output directory to put fasta files into')
args = parser.parse_args()
fasta_table = pd.read_csv(args.input_fasta_links, sep='\t')
# You must have a valid API Token
API_TOKEN = os.getenv('ENTEROBASE_API_TOKEN', None)

SERVER_ADDRESS = 'https://enterobase.warwick.ac.uk'
DATABASE = 'ecoli'

def __create_request(request_str):
    base64string = base64.b64encode('{0}: '.format(API_TOKEN).encode('utf-8'))
    headers = {"Authorization": "Basic {0}".format(base64string.decode())}
    request = urllib.request.Request(request_str, None, headers)
    return request

OUTDIR = args.output_dir
if not os.path.exists(OUTDIR):
    os.mkdir(OUTDIR)

fasta_links = fasta_table['download_fasta_link'].tolist() # turn the columns to lists because paranoid of any type of indexing issues with pandas
assembly_barcodes = fasta_table['assembly_barcode'].tolist()

for i,j in zip(fasta_links, assembly_barcodes):
    try:
        response = urlopen(__create_request(i))
        with open(os.path.join(f'{OUTDIR}/{j}.fasta'), 'w') as out_ass:
            out_ass.write(response.read().decode())
    except HTTPError as Response_error:
        #print(Response_error)
        logging.error('%d %s. <%s>\n Reason: %s' %(Response_error.code,
                                              Response_error.reason,
                                              Response_error.geturl(),
                                              Response_error.read()))
