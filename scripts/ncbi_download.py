'''
this script downloads sequencing files from ncbi with sratoolkit using file
accession numbers
-- executes download on front end
-- script should be called as part of snakemake pipeline
-- used from within snakemake, so print statements not part of interactive
testing or module loading should include flush=True
-- usage: snakemake download_ncbi
'''


import subprocess
import pandas as pd
import os
import argparse
import glob
import logging


modname = 'ncbi_download'


def read_inputs(a, b, c):
    '''
    read target tables and merge them into a single dataframe
    '''
    dfa = pd.read_csv(a, sep='\t')
    dfb = pd.read_csv(b, sep='\t')
    dfc = pd.read_csv(c, sep='\t')
    df = pd.concat([dfa, dfb, dfc], ignore_index=True)
    return df


def get_parameters(df):
    '''
    assemble accession number, source database, and seq format for each sample
    '''
    pars = {}
    for i in df.index:
        if df.loc[i,'data'] == 'gdc':
            continue
        acc = df.loc[i,'dirName']
        pars[acc] = {}
        pars[acc]['accession'] = df.loc[i,'dirName']
        pars[acc]['database'] = df.loc[i,'data']
        pars[acc]['format'] = df.loc[i,'format']
    return pars


def check_fastq(accs, dest):
    '''
    -- compare list of accession numbers for download to file names already
    downloaded
    -- return list of accession numbers to be downloaded
    '''
    fq = glob.glob(dest + '*.fastq.gz')
    l = [i.split('/')[-1][:-len('.fastq.gz')].split('_')[0] for i in fq]
    acc_done = [i for i in set(l) if i in accs]
    acc_todo = [i for i in accs if i not in acc_done]
    return acc_todo


def download_cmd(d, threads, key):
    '''
    construct fasterq-dump download command using dictionary of parameters
    '''
    cmd = 'fasterq-dump'
    if d['database'] == 'sradbg':
        cmd = ' '.join([cmd, '--ngc', key])
    thr =  ' '.join(['--threads', str(threads)])
    opx = '--skip-technical --details'
    if d['format'] == 'paired':
        opx = ' '.join([opx, '--split-files'])
    tmp = ' '.join(['--temp', os.environ['download_dir']])
    odr = ' '.join(['--outdir', os.environ['download_dir']])
    ofi = ' '.join(['--outfile', '.'.join([d['accession'], 'fastq'])])
    cmd = ' '.join([cmd, thr, opx, tmp, odr, ofi, d['accession']])
    return cmd


def get_files(a, b, c, dest, threads, key, live):
    '''
    download target fastq files by accession number
    '''
    df = read_inputs(a, b, c)
    pars = get_parameters(df)
    acc_todo = check_fastq(pars.keys(), dest)
    for e in acc_todo:
        logging.info('Working on accession %s', e)
        cmd = download_cmd(pars[e], threads, key)
        if pars[e]['format'] == 'paired':
            logging.info('Downloading paired-end files. This may take a while ...')
        else:
            logging.info('Downloading single-end file. This may take a while ...')
        if live:
            subprocess.call(cmd, shell=True)
            logging.info('Download command sent!')
        else:
            print(cmd, flush=True)
        logging.info('Done with accession %s', e)


def main(args):
    '''
    get fastq files from NCBI using accession numbers from input tables
    '''
    logging.basicConfig(
        filename='.'.join([modname,'log']),
        format='%(levelname)s:%(asctime)s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        level=logging.INFO
    )
    a = args.tsv_a
    b = args.tsv_b
    c = args.tsv_c
    dest = os.environ['download_dir']
    threads = 8
    key = args.key
    live = args.run_type
    get_files(a, b, c, dest, threads, key, live)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='ncbi_download')
    parser.add_argument(
        '--tsv_a', help='path to cohort A TSV', action='store', dest='tsv_a'
    )
    parser.add_argument(
        '--tsv_b', help='path to cohort B TSV', action='store', dest='tsv_b'
    )
    parser.add_argument(
        '--tsv_c', help='path to cohort C TSV', action='store', dest='tsv_c'
    )
    parser.add_argument(
        '--key', help='path to dbGaP access key', action='store', dest='key'
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
        tsv_a='/scratch/chd5n/herv/todo_A.tsv',
        tsv_b='/scratch/chd5n/herv/todo_B.tsv',
        tsv_c='/scratch/chd5n/herv/todo_C.tsv',
        key='',
        run_type=False
    )
    args = parser.parse_args()
    main(args)
else:
    print('functions loaded for', modname)
