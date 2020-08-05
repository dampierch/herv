'''
this script queries the gdc archive via the search and retrieve api and returns
files_res object (results from files endpoint query)

script should be called from within gdc_ann_make, which itself should be called
as part of snakemake pipeline
-- usage: snakemake setup_tcga
'''


import io
import json
import pandas as pd
import requests
import re


modname = 'gdc_req_files'


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
            {'op':'in',
            'content':{
                'field':'cases.primary_site',
                'value':['Colon', 'Rectum', 'Rectosigmoid junction']
            }
            },
            {'op':'in',
            'content':{
                'field':'experimental_strategy',
                'value':exp_strat
            }
            },
            {'op':'in',
            'content':{
                'field':'data_category',
                'value':'Sequencing Reads'
            }
            }
        ]
    }
    filters = json.dumps(filters)
    return filters


def set_fields(field_type):
    '''
    set fields for extraction from endpoint
    '''
    if field_type == 'default':
        fields = [
            'file_name',
            'file_id',
            'file_size',
            'md5sum',
            'state',
            'submitter_id',
            'experimental_strategy',  ## 'WXS', 'RNA-Seq'
            'data_category',  ## 'Sequencing Reads'
            'data_format',  ## 'BAM'
            'cases.case_id',
            'cases.project.project_id',
            'cases.primary_site',
            'cases.demographic.gender',
            'cases.demographic.race',
            'cases.demographic.year_of_birth',
            'cases.demographic.vital_status',
            'cases.demographic.days_to_death',
            'cases.demographic.age_at_index',
            'cases.samples.sample_type',  ## tumor vs normal
            'cases.samples.is_ffpe',
            'cases.samples.sample_id',
            'cases.tissue_source_site.name',
            'cases.exposures.weight',
            'cases.exposures.bmi',
            'cases.exposures.height'

            ## nested nature makes more suitable for special lookup
            #'analysis.metadata.read_groups.read_group_qcs.total_sequences',

            ## nested nature makes more suitable for invididual lookup
            # 'analysis.metadata.read_groups.target_capture_kit_name',
            # 'analysis.metadata.read_groups.target_capture_kit',
            # 'analysis.metadata.read_groups.target_capture_kit_target_region',
            # 'analysis.metadata.read_groups.target_capture_kit_catalog_number',
            # 'analysis.metadata.read_groups.target_capture_kit_vendor'

            ## useless fields
            # 'cases.demographic.ethnicity', # less useful than race; 'not hispanic or latino'
            # 'platform', # generic Illumina
            # 'data_type', # 'Aligned Reads'
            # 'annotations.entity_id', # same as cases.case_id
            # 'annotations.case_id', # same as cases.case_id
            # 'msi_score', # empty
            # 'msi_status', # empty
            # 'proportion_reads_mapped', # empty
            # 'proportion_targets_no_coverage' # empty
            # 'cases.lost_to_followup', # empty
            # 'cases.samples.annotations.entity_id', # empty
            # 'cases.follow_ups.molecular_tests.mismatch_repair_mutation', # empty
            # 'cases.follow_ups.molecular_tests.ploidy', # empty
            # 'analysis.metadata.read_groups.target_capture_kit_version', # empty
            # 'analysis.metadata.read_groups.RIN', # empty
            # 'analysis.metadata.read_groups.instrument_model', # empty
            # 'analysis.metadata.read_groups.read_group_qcs.per_tile_sequence_quality', # PASS/FAIL/WARN
            # 'analysis.metadata.read_groups.read_group_qcs.per_base_sequence_content', # PASS/FAIL/WARN
            # 'analysis.metadata.read_groups.read_group_qcs.per_base_sequence_quality', # PASS/FAIL/WARN
            # 'analysis.metadata.read_groups.read_group_qcs.per_sequence_gc_content', # PASS/FAIL/WARN
            # 'analysis.metadata.read_groups.read_group_qcs.per_sequence_quality_score', # PASS/FAIL/WARN
            # 'analysis.metadata.read_groups.read_group_qcs.workflow_type', # empty
            # 'analysis.metadata.read_groups.read_group_qcs.workflow_version', # empty
            # 'analysis.metadata.read_groups.read_group_qcs.fastq_name', # empty
            # 'analysis.metadata.read_groups.read_group_qcs.basic_statistics', # empty
            # 'downstream_analyses.output_files.proportion_reads_mapped', # empty
            # 'downstream_analyses.output_files.proportion_targets_no_coverage', # empty
            # 'archive.archive_id', # empty
            # 'proportion_coverage_30X', # empty
            # 'average_base_quality', # empty
            # 'mean_coverage', # empty
            # 'cases.sample_ids', # emtpy
            # 'cases.samples.tissue_type', # typically 'not reported'
            # 'cases.samples.distance_normal_to_tumor', # empty
            # 'cases.samples.biospecimen_anatomic_site', # empty
            # 'cases.exposures.pack_years_smoked', # empty
            # 'cases.exposures.alcohol_history', # typically 'not reported'
        ]
    elif field_type == 'depth':
        fields = [
            'file_id',
            'analysis.metadata.read_groups.read_group_qcs.total_sequences'
        ]
    else:
        sys.exit('must enter valid field_type (i.e. default, depth)')
    fields = ','.join(fields)
    return fields


def set_params(filters,fields):
    '''
    set parameters for https get request to endpoint
    -- set size parameter empirically to a level greater than number of target
    files to get all records at once
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
    given an endpoint and parameters, execute https GET request and build
    files_res dataframe with file info
    '''
    response = requests.get(endpoint, params=params)
    object = io.StringIO(response.content.decode('utf-8'))
    results = pd.read_table(object)
    return results


def parse_depth(read_depths):
    '''
    given read_depths dataframe, return total_seq dataframe with aggregate read
    depth for each file_id as well as read format (i.e. paired vs single)
    -- regex identifies all possible seq groups, takes first of paired-end
    seq counts
    '''
    read_depths['read_format'] = pd.Series()
    se_idx = read_depths[
        'analysis.metadata.read_groups.0.read_group_qcs.1.total_sequences'].isna()
    read_depths.loc[se_idx,'read_format'] = 'single'
    read_depths.loc[-se_idx,'read_format'] = 'paired'
    cols = list(read_depths.columns)
    pattern = re.compile(
        r'analysis.metadata.read_groups.\d*.read_group_qcs.0.total_sequences'
        )
    read_depth_sum = read_depths[
        list(filter(pattern.match, cols))].sum(axis=1).astype(int)
    frame = {
        'file_id': read_depths['file_id'],
        'read_format': read_depths['read_format'],
        'total_seq': read_depth_sum
        }
    total_seq = pd.DataFrame(frame)
    return total_seq


def fix_format(files_res):
    '''
    format files_res
    '''
    column_dict = {
        'cases.0.samples.0.is_ffpe': 'is_ffpe',
        'cases.0.exposures.0.height': 'height',
        'cases.0.tissue_source_site.name': 'hospital',
        'cases.0.project.project_id': 'project_id',
        'cases.0.demographic.year_of_birth': 'birth_year',
        'id': 'uuid',
        'cases.0.samples.0.sample_type': 'tissue_type',
        'cases.0.demographic.days_to_death': 'days_to_death',
        'cases.0.demographic.vital_status': 'vital_status',
        'cases.0.samples.0.sample_id': 'sample_id',
        'cases.0.exposures.0.weight': 'weight',
        'cases.0.demographic.race': 'race',
        'cases.0.demographic.age_at_index': 'age_at_index',
        'cases.0.demographic.gender': 'sex',
        'cases.0.case_id': 'case_id',
        'cases.0.primary_site': 'primary_site',
        'cases.0.exposures.0.bmi': 'bmi'
    }
    column_order = [
        'file_id', 'uuid', 'file_name', 'sample_id', 'case_id', 'submitter_id',
        'project_id', 'primary_site', 'tissue_type', 'birth_year',
        'vital_status', 'age_at_index', 'days_to_death', 'height', 'weight',
        'bmi', 'race', 'sex', 'hospital', 'is_ffpe', 'experimental_strategy',
        'data_format', 'data_category', 'read_format', 'total_seq', 'file_size',
        'md5sum', 'state'
        ]
    files_res = files_res.rename(column_dict, axis='columns')
    files_res = files_res[column_order]
    return files_res


def main():
    ## set endpoint
    endpoint = 'https://api.gdc.cancer.gov/files/'
    exp_strat = 'RNA-Seq'
    ## get baseline files uuids and other info
    filters = set_filters(exp_strat)
    fields = set_fields('default')
    params = set_params(filters, fields)
    files_res = get_results(endpoint, params)
    ## add sequencing depths
    fields = set_fields('depth')
    params = set_params(filters, fields)
    read_depths = get_results(endpoint, params)
    total_sequences = parse_depth(read_depths)
    files_res = files_res.merge(total_sequences, on='file_id')
    ## fix formatting
    files_res = fix_format(files_res)
    files_res = files_res.drop_duplicates('file_id')
    ## return results object
    return files_res


if __name__ == '__main__':
    print('This script is not meant to be run as main. See usage statment:')
    print('usage: snakemake setup_tcga')
else:
    files_res = main()
