from urllib.request import urlopen
from urllib.error import HTTPError
import urllib
import base64
import json
import os
import functools
# You must have a valid API Token
API_TOKEN = os.getenv('ENTEROBASE_API_TOKEN', None)
SERVER_ADDRESS = 'https://enterobase.warwick.ac.uk'
DATABASE = 'ecoli'

def __create_request(request_str):
    base64string = base64.b64encode('{0}: '.format(API_TOKEN).encode('utf-8'))
    headers = {"Authorization": "Basic {0}".format(base64string.decode())}
    request = urllib.request.Request(request_str, None, headers)
    return request

outdir = 'ppp'
if not os.path.exists(outdir):
    os.mkdir(outdir)
address = f'{SERVER_ADDRESS}/api/v2.0/{DATABASE}/straindata?assembly_status=Assembled&limit=1000&only_fields=strain_barcode,download_fasta_link&offset=0'
curr_address = address
max_iter_test = 1000
strain_dicts = []

while curr_address != None and max_iter_test > 0: 
    try:
        response = urlopen(__create_request(curr_address))
        data = json.load(response)
        strain_dicts.append(data['straindata'])
        curr_address = data['links']['paging'].get('next')
        max_iter_test -= 1
        print(max_iter_test)

    except HTTPError as Response_error:
        logging.error('%d %s. <%s>\n Reason: %s' %(Response_error.code,
                                                   Response_error.reason,
                                                   Response_error.geturl(),
                                                   Response_error.read()))
combined_dict = functools.reduce(lambda x,y: {**x, **y}, strain_dicts)
combined_json = json.dumps(combined_dict, indent=4)
with open("enterobase_meta_with_links.json", 'w') as output:
            output.write(combined_json)
