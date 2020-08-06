'''
this script downloads sequencing files from gdc with gdc-client using file uuid
instead of manifest due to time contraints on front end
-- executes download on front end
-- script should be called as part of snakemake pipeline
-- used from within snakemake, so print statements not part of interactive
testing or module loading should include flush=True
-- usage: snakemake download_tcga

-- 695 uniq file_id
-- 686 uniq sample_id
-- 620 uniq case_id, subject_id
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


def set_download(man, ann):
    '''
    identifies subset of files to download based on following criteria:
    -- not from ffpe specimen
    -- from unique subject OR
    -- from replicated subject but is
        -- paired-end sample with highest number of reads OR if no paired-end
        -- single-end sample with highest number of reads
    '''
    ## ffpe filter
    idx = ann['is_ffpe'] == False
    ann = ann.loc[idx,:]
    ## unique tag, includes unique subject filter
    idx = pd.DataFrame(ann['subject_id'].value_counts() == 1)
    idx['unq'] = idx['subject_id']
    idx['subject_id'] = idx.index
    ann = ann.merge(idx, on='subject_id')
    ## download tag, includes max read count filter for replicates
    ann['download'] = pd.Series()
    fids = []
    for i in ann.index:
        if ann.loc[i,'unq']:
            ann.loc[i,'download'] = True
        else:
            ann.loc[i,'download'] = False
            id = ann.loc[i,'subject_id']
            idx = ann['subject_id'] == id
            df = ann.loc[idx,['file_id','subject_id','read_format','total_seq']]
            if any(df['read_format'] == 'paired'):
                idx = df['read_format'] == 'paired'
                df_pe = df.loc[idx,:]
                idx = df_pe['total_seq'] == df_pe['total_seq'].max()
                fids.append(list(df_pe.loc[idx,'file_id'])[0])
            else:
                idx = df['total_seq'] == df['total_seq'].max()
                fids.append(list(df.loc[idx,'file_id'])[0])
    for i in ann.index:
        if not ann.loc[i,'download']:
            if ann.loc[i,'file_id'] in set(fids):
                ann.loc[i,'download'] = True
    idx = ann['download'] == True
    ann = ann.loc[idx,:]
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
    man, ann = set_download(man, ann)
    write_outputs(man, ann, out_man, out_ann)
    uuids = get_uuids(ann, read_format)
    return uuids


def check_bams(uuids, dest):
    '''
    -- compares list of files_ids for download to file_ids downloaded already
    -- returns list of file_ids to be downloaded
    '''
    bams = glob.glob(dest + '*/*.bam')
    uuids_done = [i.split('/')[-2] for i in bams if i in uuids]
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
