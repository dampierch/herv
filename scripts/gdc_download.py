'''
this script downloads sequencing files from gdc with gdc-client using file uuid
instead of manifest due to time contraints on front end
-- executes download on front end
-- script should be called as part of snakemake pipeline
-- used from within snakemake, so print statements not part of interactive
testing or module loading should include flush=True
-- usage: snakemake download_tcga
'''


import subprocess
import pandas as pd
import os
import argparse
import glob
import logging


modname = 'gdc_download'


def read_inputs(in_man, in_ann):
    '''
    reads targets
    '''
    man = pd.read_csv(in_man, sep='\t')
    ann = pd.read_csv(in_ann, sep='\t')
    return man, ann


def unique_subjects(df):
    '''
    set unique subject tag for given dataframe
    '''
    unq = pd.DataFrame(df['subject_id'].value_counts() == 1)
    unq['unq'] = unq['subject_id']
    unq['subject_id'] = unq.index
    df = df.merge(unq, on='subject_id')
    return df


def sample_preference(df):
    '''
    set download tag with the following priority list:
    -- unique subject
    -- paired read format for replicates
    -- max read count for replicates
    '''
    df['download'] = pd.Series()
    fids = []
    for i in df.index:
        if df.loc[i,'unq']:
            df.loc[i,'download'] = True
        else:
            df.loc[i,'download'] = False
            id = df.loc[i,'subject_id']
            idx = df['subject_id'] == id
            subdf = df.loc[idx,['file_id','subject_id','read_format','total_seq']]
            if any(subdf['read_format'] == 'paired'):
                idx = subdf['read_format'] == 'paired'
                pedf = subdf.loc[idx,:]
                idx = pedf['total_seq'] == pedf['total_seq'].max()
                fids.append(list(pedf.loc[idx,'file_id'])[0])
            else:
                idx = subdf['total_seq'] == subdf['total_seq'].max()
                fids.append(list(subdf.loc[idx,'file_id'])[0])
    for i in df.index:
        if not df.loc[i,'download']:
            if df.loc[i,'file_id'] in set(fids):
                df.loc[i,'download'] = True
    return df


def set_selection(man, ann):
    '''
    identifies subset of files to download based on following criteria:
    -- not from ffpe specimen
    -- from Adenomas and Adenocarcinomas
    -- from Primary Tumor OR Solid Tissue Normal
    -- from unique subject within tissue_type OR
        -- from replicated subject BUT
            -- from paired-end sample with highest number of reads within
            tissue_type when paired-end sample present OR
            -- from single-end sample with highest number of reads within
            tissue_type when paired-end sample absent
    '''
    ## ffpe filter
    idx = ann['is_ffpe'] == False
    ann = ann.loc[idx,:]
    ## adenocarcinoma filter
    idx = ann['disease_type'] == 'Adenomas and Adenocarcinomas'
    ann = ann.loc[idx,:]
    ## tissue type filter
    idx = ann['tissue_type'] == 'Primary Tumor'
    pri = ann.loc[idx,:]
    idx = ann['tissue_type'] == 'Solid Tissue Normal'
    nrm = ann.loc[idx,:]
    ## unique subject tag
    pri = unique_subjects(pri)
    nrm = unique_subjects(nrm)
    ## read format and max depth tag
    pri = sample_preference(pri)
    nrm = sample_preference(nrm)
    ann = pd.concat([pri, nrm])
    ## unique subject, read format, and max depth filter
    idx = ann['download'] == True
    ann = ann.loc[idx,:]
    ann = ann.drop(columns=['unq','download'])
    idx = man['id'].isin(ann['file_id'])
    man = man.loc[idx,:]
    return man, ann


def write_outputs(df1, df2, fn1, fn2):
    df1.to_csv(fn1, sep='\t')
    df2.to_csv(fn2, sep='\t')


def get_uuids(ann, read_format):
    '''
    assemble uuids (ann file_ids, man ids) for download into list by read_format
    '''
    idx = ann['read_format'] == read_format
    uuids = list(ann.loc[idx,'file_id'])
    return uuids


def parse_files(in_man, in_ann, out_man, out_ann, read_format):
    '''
    -- parse full gdc annotation and manifest
    -- write select gdc annotation and manifest for download
    -- return uuids for download
    '''
    man, ann = read_inputs(in_man, in_ann)
    man, ann = set_selection(man, ann)
    write_outputs(man, ann, out_man, out_ann)
    uuids = get_uuids(ann, read_format)
    return uuids


def check_bams(uuids, dest):
    '''
    -- compares list of files_ids for download to file_ids downloaded already
    -- returns list of file_ids to be downloaded
    '''
    bams = glob.glob(dest + '*/*.bam')
    uuids_done = [i.split('/')[-2] for i in bams if i.split('/')[-2] in uuids]
    uuids_todo = [i for i in uuids if i not in uuids_done]
    return uuids_todo


def download_files(dest, token, uuids_todo, live):
    '''
    download files specified in uuids_todo on frontend
    '''
    logging.info('download attempt start')
    file_count = 0
    targets = uuids_todo
    for uuid in targets:
        cmd = ' '.join(['gdc-client download', uuid, '-d', dest, '-t', token])
        if live:
            subprocess.call(cmd, shell=True)
        else:
            print(cmd, flush=True)
        file_count = file_count + 1
        logging.info('%s attempted', uuid)
    logging.info('%s files attempted', file_count)
    logging.info('download attempt end')


def main(args):
    '''
    parse uuids and download selected set
    '''
    logging.basicConfig(
        filename='.'.join([modname,'log']),
        format='%(levelname)s:%(asctime)s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        level=logging.INFO
    )
    in_man = os.environ['ann_dir'] + 'tcga/' + 'manifest.tsv'
    in_ann = os.environ['ann_dir'] + 'tcga/' + 'ann.tsv'
    out_man = os.environ['ann_dir'] + 'tcga/' + 'manifest_select.tsv'
    out_ann = os.environ['ann_dir'] + 'tcga/' + 'ann_select.tsv'
    read_format = args.read_format
    uuids = parse_files(in_man, in_ann, out_man, out_ann, read_format)
    dest = os.environ['download_dir']
    if args.token == 'gdc-user-token':
        token = glob.glob(dest + args.token + '*')
    else:
        token = args.token
    uuids_todo = check_bams(uuids, dest)
    live = args.run_type
    download_files(dest, token, uuids_todo, live)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='gdc_download')
    parser.add_argument(
        '--read_format', help='paired or single', action='store', dest='read_format'
    )
    parser.add_argument(
        '--token', help='gdc access token', action='store', dest='token'
    )
    parser.add_argument(
        '--live', help='execute a live run', action='store_true',
        dest='run_type'
    )
    parser.add_argument(
        '--test', help='execute a test run', action='store_false',
        dest='run_type'
    )
    parser.set_defaults(
        read_format='paired', token='gdc-user-token', run_type=False
    )
    args = parser.parse_args()
    main(args)
else:
    print('functions loaded for', modname)
