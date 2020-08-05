'''
-- this script creates manifest and pheno objects from gdc query results
and writes them to file for future use in download or parsing
-- the key field offered by cases endpoint is the subject_id (submitter_id)
-- script should be called as part of snakemake pipeline
-- usage: snakemake setup_tcga
'''


import io
import json
import os
import pandas as pd


from gdc_req_files import files_res
from gdc_req_cases import cases_res
from gdc_req_legacy import msi_res


modname = 'gdc_ann_make'


def merge_results(files_res,cases_res,msi_res):
    '''
    -- add cases_res to files_res matching on case_id, which also adds
    subject_id
    -- add msi_res to merged results matching on subject_id
    -- merge function performs inner join by default
    -- merge will perform left join with how='left' option
    '''
    results = files_res.merge(cases_res, on='case_id')
    results = results.merge(msi_res, how='left', on='subject_id')
    return results


def make_manifest(gdc_ann,ann_dir):
    '''
    -- create manifest object from gdc_ann and write to file
    -- manifest fields: id, filename, md5, size, state
    '''
    select = ['file_id', 'file_name', 'md5sum', 'file_size', 'state']
    manifest = gdc_ann[select]
    manifest.columns = ['id', 'filename', 'md5', 'size', 'state']
    manifest = manifest.sort_values(by=['id'])
    manifest.to_csv(ann_dir + 'manifest.tsv', sep='\t', index=False)


def make_pheno(gdc_ann,ann_dir):
    '''
    create pheno annotations object from gdc_ann and write to file
    '''
    select = [
        'file_id', 'file_name', 'sample_id', 'case_id', 'submitter_id',
        'subject_id', 'project_id', 'disease_type', 'primary_site',
        'tissue_type', 'stage', 'msi_status', 'birth_year', 'vital_status',
        'age_at_index', 'days_to_death', 'height', 'weight', 'bmi', 'race',
        'sex', 'hospital', 'is_ffpe', 'read_format', 'total_seq'
        ]
    pheno = gdc_ann[select]
    pheno = pheno.sort_values(by=['file_id'])
    pheno.to_csv(ann_dir + 'ann.tsv', sep='\t', index=False)


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
    '''
    make a manifest and phenotype annotation for dataset of interest
    '''
    ann_dir = os.environ['ann_dir'] + 'tcga/'
    gdc_ann = merge_results(files_res,cases_res,msi_res)
    check_outpath(ann_dir)
    make_manifest(gdc_ann,ann_dir)
    make_pheno(gdc_ann,ann_dir)


if __name__ == '__main__':
    main()
else:
    print('functions loaded for', modname)
