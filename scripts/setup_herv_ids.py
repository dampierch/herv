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


def load_ids(bed):
    '''
    read herv transcript ids quantified with salmon into list
    '''
    with open(bed, 'r') as infile:
        l = []
        for inline in infile:
            l.append(inline.strip('\n').split('\t')[3])
    return l


def make_id_df(l):
    '''
    make dataframe with mapping of herv transcript ids to herv names
    '''
    fp = os.environ['attig_dir']
    fn = 'attig_flat_transcriptome.gtf.gz'
    target = ''.join([fp, fn])
    txdct, exdct = gtf_to_dct(target)
    d = {k: txdct[k] for k in txdct if k in l}
    d = {k: d[k][0]['repeatIDs'].split('%3b') for k in d}
    tx_id = []
    herv_id = []
    for k in d:
        tx_id.append(k)
        l = [e for e in d[k] if 'HERV' in e]
        if len(l) == 1:
            herv_id.append(l[0])
        else:
            herv_id.append(';'.join(l))
    d = {'tx_id': tx_id, 'herv_id': herv_id}
    df = pd.DataFrame(d)
    return df


def main(args):
    bed = args.bed
    l = load_ids(bed)
    df = make_id_df(l)
    fp = os.environ['genmod_dir']
    fn = 'herv_ids.tsv'
    target = ''.join([fp, fn])
    df.to_csv(target, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='setup_herv_ids')
    parser.add_argument(
        '--bed', help='path to herv bed', action='store', dest='bed'
    )
    parser.set_defaults(
        bed='',
    )
    args = parser.parse_args()
    main(args)
else:
    print('functions loaded for', modname)
