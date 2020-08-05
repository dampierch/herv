'''
this script queries the gdc legacy archive via the search and retrieve api and
returns msi_status object (from files endpoint on legacy)
-- get uuids of xml files with the msi annotations from legacy server
-- download each xml file
-- parse xml files to extract msi annotations for each subject

script should be called from within gdc_ann_make, which itself should be called
as part of snakemake pipeline
-- usage: snakemake setup_tcga
'''


import io
import json
import os
import pandas as pd
import requests
import re
import subprocess
import glob
import xml.etree.ElementTree as ET


modname = 'gdc_req_legacy'


def set_filters():
    '''
    set filters for gdc legacy files endpoint search
    -- json format
    -- for files.data_type, values for MSI status are 'Auxiliary test' and
    'Microsatellite instability'
    -- here use 'Auxiliary test' per TCGAbiolinks examples
    '''
    filters = {
        'op':'and',
        'content':[
            {'op':'or',
            'content':[
                {'op':'in',
                'content':{
                    'field':'cases.project.project_id',
                    'value':'TCGA-COAD'
                }
                },
                {'op':'in',
                'content':{
                    'field':'cases.project.project_id',
                    'value':'TCGA-READ'
                }
                }
            ]
            },
            {'op':'and',
            'content':[
                {'op':'in',
                'content':{
                    'field':'files.data_category',
                    'value':'Other'
                }
                },
                {'op':'in',
                'content':{
                    'field':'files.data_type',
                    'value':'Auxiliary test'
                }
                },
                {'op':'in',
                'content':{
                    'field':'files.access',
                    'value':'open'
                }
                }
            ]
            }
        ]
    }
    filters = json.dumps(filters)
    return filters


def set_fields():
    '''
    set fields for extraction from endpoint
    '''
    fields = [
        'file_name',
        'file_id',
        'md5sum',
        'file_size',
        'state'
    ]
    fields = ','.join(fields)
    return fields


def set_params(filters,fields):
    '''
    set parameters for https get request to endpoint
    -- set size parameter empirically to a level greater than number of target
    cases to get all records at once
    '''
    params = {
        'filters': filters,
        'fields': fields,
        'format': 'TSV',
        'size': '1500'
    }
    return params


def get_results(endpoint,params):
    '''
    given an endpoint and parameters, execute https GET request for xml file_id
    entities and build results dataframe with msi results
    '''
    response = requests.get(endpoint, params=params)
    object = io.StringIO(response.content.decode('utf-8'))
    results = pd.read_table(object)
    return results


def download_xml_uuid(files_res,dest):
    '''
    download xml files one at a time by uuid
    '''
    file_count = 0
    for uuid in files_res.id:
        cmd = ' '.join(['gdc-client download',uuid,'-d',dest])
        subprocess.call(cmd, shell=True)
        print(' '.join([uuid,'downloaded']))
        file_count = file_count + 1
    print(' '.join([str(file_count),'files downloaded']))


def download_xml_manifest(files_res,dest):
    '''
    -- create manifest object
    -- write manifest to file
    -- use manifest for bulk download
    '''
    select = ['file_id', 'file_name', 'md5sum', 'file_size', 'state']
    manifest = files_res[select]
    manifest.columns = ['id', 'filename', 'md5', 'size', 'state']
    manifest = manifest.sort_values(by=['id'])
    out_file = dest + 'manifest.tsv'
    manifest.to_csv(out_file, sep='\t', index=False)
    cmd = ' '.join(['gdc-client download','-m',out_file,'-d',dest])
    subprocess.call(cmd, shell=True)
    print('manifest downloaded')


def parse_xml(files_res,dest):
    '''
    parse xml files to extract msi status
    '''
    msi_dict = {}
    msi_dict['subject_id'] = []
    msi_dict['msi_status'] = []
    tag1 = 'mononucleotide_and_dinucleotide_marker_panel_analysis_status'
    tag2 = 'mononucleotide_marker_panel_analysis_status'
    file_count = 0
    for uuid in files_res.id:
        pattern = dest + uuid + '/*.xml'
        fn = glob.glob(pattern)[0]
        tree = ET.parse(fn)
        for elem in tree.getiterator():
            if 'bcr_patient_barcode' in elem.tag:
                subject_id = elem.text
            if tag1 in elem.tag and elem.text != None:
                msi_status = elem.text
            elif tag2 in elem.tag and elem.text != None:
                msi_status = elem.text
        msi_dict['subject_id'].append(subject_id)
        msi_dict['msi_status'].append(msi_status)
        file_count = file_count + 1
    print(' '.join([str(file_count),'files parsed']))
    msi_res = pd.DataFrame.from_dict(msi_dict)
    return msi_res


def check_outpath(out_path):
    '''
    check for presence of absence of out_path and make directory if absent
    '''
    l = out_path.strip('/').split('/')
    d = ''
    for e in l:
        d = d + '/' + e
        if os.path.exists(d):
            print(d,'present')
        else:
            print(d,'absent')
            print('making',d,'now')
            os.mkdir(d)


def main():
    endpoint = 'https://api.gdc.cancer.gov/legacy/files/'
    filters = set_filters()
    fields = set_fields()
    params = set_params(filters, fields)
    files_res = get_results(endpoint, params)
    dest = os.environ['ann_dir'] + 'tcga/msi/'
    check_outpath(dest)
    download_xml_manifest(files_res, dest)
    msi_res = parse_xml(files_res, dest)
    return msi_res


if __name__ == '__main__':
    print('This script is not meant to be run as main. See usage statment:')
    print('usage: snakemake setup_tcga')
else:
    msi_res = main()
