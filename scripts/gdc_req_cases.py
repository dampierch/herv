'''
this script queries the gdc archive via the search and retrieve api and returns
cases_res object (results from cases endpoint query)

script should be called from within gdc_ann_make, which itself should be called
as part of snakemake pipeline
-- usage: snakemake setup_tcga
'''


import io
import json
import pandas as pd
import requests


modname = 'gdc_req_cases'


def set_filters(exp_strat):
    '''
    set filters to target relevant units on endpoint
    -- json format
    '''
    filters = {
        'op':'and',
        'content':[
            {'op':'or',
            'content':[
                {'op':'in',
                'content':{
                    'field':'project.project_id',
                    'value':'TCGA-COAD'
                }
                },
                {'op':'in',
                'content':{
                    'field':'project.project_id',
                    'value':'TCGA-READ'
                }
                }
            ]
            },
            {'op':'in',
            'content':{
                'field':'primary_site',
                'value':['Colon', 'Rectum', 'Rectosigmoid junction']
            }
            },
            {'op':'in',
            'content':{
                'field':'files.experimental_strategy',
                'value':exp_strat
            }
            },
            {'op':'in',
            'content':{
                'field':'files.data_category',
                'value':'Sequencing Reads'
            }
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
        'submitter_id',  ## good
        'case_id',  ## good
        'disease_type',  ## good
        'diagnoses.tumor_stage'  ## good
        # 'sample_ids' ## ok

        # "submitter_analyte_ids", #?
        # "analyte_ids", #?
        # "portion_ids", #?
        # "submitter_portion_ids", #?

        # 'diagnoses.tumor_grade',  ## 'not reported'
        # 'exposures.alcohol_history',  ## 'Not Reported'
        # 'exposures.alcohol_intensity',  ## empty
        # 'exposures.cigarettes_per_day',  ## empty
        # 'exposures.years_smoked'  ## empty
        # "days_to_index", # empty
        # 'samples.biospecimen_anatomic_site', # empty
        # 'samples.distance_normal_to_tumor', # empty
        # 'samples.annotations.case_id', # empty
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
        'size': '1000'
    }
    return params


def get_results(endpoint,params):
    '''
    given an endpoint and parameters, execute https GET request and build
    cases_res dataframe with subject info
    '''
    response = requests.get(endpoint, params=params)
    object = io.StringIO(response.content.decode('utf-8'))
    results = pd.read_table(object)
    return results


def fix_format(cases_res):
    '''
    format cases_res
    -- submitter_id should be subject_id for clarity
    '''
    column_dict = {
        'submitter_id': 'subject_id',
        'diagnoses.0.tumor_stage': 'stage',
    }
    column_select = ['subject_id', 'case_id', 'disease_type', 'stage']
    cases_res = cases_res.rename(column_dict, axis='columns')
    cases_res = cases_res[column_select]
    return cases_res


def main():
    '''
    get all case_ids from TCGA-COAD and READ
    -- 630 unique case_id entries w/o filters for experimental_strategy and
    data_category
    -- 620 unique case_id entries (and 695 files) w/ filters
    '''
    endpoint = 'https://api.gdc.cancer.gov/cases'
    exp_strat = 'RNA-Seq'
    filters = set_filters(exp_strat)
    fields = set_fields()
    params = set_params(filters, fields)
    cases_res = get_results(endpoint, params)
    cases_res = fix_format(cases_res)
    return cases_res


if __name__ == '__main__':
    print('This script is not meant to be run as main. See usage statment:')
    print('usage: snakemake setup_tcga')
else:
    cases_res = main()
