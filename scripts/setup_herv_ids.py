'''
this script assembles a table mapping transcript id to HERV name for select
HERV transcripts from Attig et al.
-- executes on front end
-- script should be called as part of snakemake pipeline
-- used from within snakemake, so print statements not part of interactive
testing or module loading should include flush=True
-- usage: snakemake setup_dge
'''


import os
import gzip
import pandas as pd
import argparse

from setup_herv_bed import make_att_dct, gtf_to_dct


modname = 'setup_herv_ids'


def load_fasta_ids(fasta):
    '''
    read herv transcript ids quantified with salmon from fasta into dict
    '''
    with open(fasta, 'r') as infile:
        d = {}
        for inline in infile:
            if inline[0] == '>':
                k = inline.strip('>\n').split('::')[0]
                v = inline.strip('>\n')
                d[k] = {}
                d[k]['ID'] = k
                d[k]['tx_id'] = v
    return d


def map_ids(d):
    '''
    make dataframe with mapping of herv transcript ids to herv names
    '''
    fp = os.environ['attig_dir']
    fn = 'attig_flat_transcriptome.gtf.gz'
    target = ''.join([fp, fn])
    txdct, exdct = gtf_to_dct(target)
    dnm = {k: txdct[k][0]['repeatIDs'].split('%3b') for k in txdct if k in d}
    tx_id = []
    herv_id = []
    for k in dnm:
        tx_id.append(d[k]['tx_id'])
        l = [e for e in dnm[k] if 'HERV' in e]
        if len(l) == 1:
            herv_id.append(l[0])
        else:
            herv_id.append(';'.join(l))
    dmp = {'tx_id': tx_id, 'herv_id': herv_id}
    df = pd.DataFrame(dmp)
    return df


def write_map(df):
    fp = os.environ['genmod_dir']
    fn = 'herv_ids.tsv'
    target = ''.join([fp, fn])
    df.to_csv(target, sep='\t', index=False)


def main(args):
    fasta = args.fasta
    d = load_fasta_ids(fasta)
    df = map_ids(d)
    write_map(df)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='setup_herv_ids')
    parser.add_argument(
        '--fasta', help='path to herv fasta', action='store', dest='fasta'
    )
    parser.set_defaults(
        fasta='/scratch/chd5n/herv/genmod/herv_transcripts.fa'
    )
    args = parser.parse_args()
    main(args)
else:
    print('functions loaded for', modname)
