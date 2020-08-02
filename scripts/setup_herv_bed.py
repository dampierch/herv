'''
-- module built for python3

-- parse attig ERE transcript annotation from GTF/GFF format into BED12 format,
output should have one unannotated HERV transcript per line and exclude introns
-- this requires not selecting for HERV transripts directly, as they include
introns, but selecting for exons associated with a transcript and then grouping
by parent transcript into BED blocks

-- select LTR-overlapping transcripts from Attig, where LTRderived=TRUE
    -- evidently exon features are considered LTRderived independent of their
    parent and sibling features, but transcript (i.e. parent) features are
    considered LTRderived if any part (i.e. any of their child exons) overlap
    with LTR
-- select transcripts identified as HERVs
-- select transcripts that do not overlap known genes
-- extract exons for selected LTR-overlapping transcripts
-- ouput coordinates for aggregated exons (i.e. intron-less transcripts) in
BED12 format

-- requires config.sh to set environmental variables
-- for use in snakemake workflow
-- usage: snakemake setup_herv_bed.py

-- input specifications:
    -- GTF: \t delimited, start index is 1-based
    -- columns: <seqname>\t<source>\t<feature>\t
                <start>\t<end>\t<score>\t
                <strand>\t<frame>\t[attributes;more]\t[comments]
    -- BED12: \t delimited, start index is 0-based
    -- columns: <chrom>\t<chromStart>\t<chromEnd>\t
                <name>\t<score>\t<strand>\t
                <thickStart>\t<thickEnd>\t<itemRgb>\t
                <blockCount>\t<blockSizes>\t<blockStarts>
'''


import os
import gzip
import re


modname = 'setup_herv_bed'


def make_att_dct(s):
    '''
    -- split attributes string into elements and then split elements into
    key-value pairs
    -- make dict from key-value pairs
    '''
    l = [e.split('=') for e in s.strip(';').split(';')]
    k = [e[0] for e in l]
    v = [e[1] for e in l]
    d = dict(list(zip(k,v)))
    return d


def gtf_to_dct(target):
    '''
    -- parse attig herv annotation GTF into dicts of transcripts and exons
    -- calculate feature length based on start and end positions; GTF/GFF
    coordinates are 1-based, so add 1 position to start,end interval
    '''
    infields = ['seqname', 'source', 'feature', 'start', 'end', 'score',
        'strand', 'frame', 'attributes']
    txdct = {}
    exdct = {}
    with gzip.open(target, 'rt') as infile:
        for inline in infile:
            if inline[0] == '#':
                continue
            else:
                indata = inline.strip('\n').split('\t')
                indct = dict(list(zip(infields, indata)))
                indct['length'] = int(indct['end']) - int(indct['start']) + 1
                d = make_att_dct(indct['attributes'])
                for k in d:
                    if k in indct:
                        print('key already in dictionary\n')
                    else:
                        indct[k] = d[k]
                if indct['feature'] == 'transcript':
                    if indct['ID'] in txdct:
                        txdct[indct['ID']].append(indct)
                    else:
                        txdct[indct['ID']] = []
                        txdct[indct['ID']].append(indct)
                elif indct['feature'] == 'exon':
                    if indct['Parent'] in exdct:
                        exdct[indct['Parent']].append(indct)
                    else:
                        exdct[indct['Parent']] = []
                        exdct[indct['Parent']].append(indct)
    return txdct,exdct


def select_tx(d):
    '''
    -- takes dict of transcripts as input, selects subset of transcripts that
    match criteria, outputs list of transcript IDs (i.e. keys)
    -- criteria are i. overlap with LTR, ii. identified as HERV, iii. no overlap
    with any known gene
    '''
    l = []
    for tx in d:
        if d[tx][0]['LTRderived'] == 'TRUE':
            if re.search('HERV', d[tx][0]['repeatIDs']):
                if d[tx][0]['geneID'] == '':
                    l.append(tx)
    return l


def dct_to_bed(l,d):
    '''
    -- takes list of selected transcripts and dict of all transcripts along with
    child exons; for each selected transcript, selects child exons, assembles
    BED12 record for selected transcript from child exons
    -- output is in unpredictable order
    '''
    bed = ''
    for tx in l:
        ex_sorted = sorted(d[tx], key=lambda ex: int(ex['start']))
        chrom = ex_sorted[0]['seqname']
        chromStart = str(int(ex_sorted[0]['start']) - 1)
        chromEnd = ex_sorted[-1]['end']
        name = tx
        score = '.'
        strand = ex_sorted[0]['strand']
        thickStart = '.'
        thickEnd = '.'
        itemRgb = '.'
        blockCount = str(len(ex_sorted))
        blockSizes = ','.join([str(ex['length']) for ex in ex_sorted])
        blockStarts = ','.join(
            [str(int(ex['start']) - int(ex_sorted[0]['start']))
                for ex in ex_sorted]
            )
        txout = '\t'.join(
            [chrom, chromStart, chromEnd, name, score, strand, thickStart,
                thickEnd, itemRgb, blockCount, blockSizes, blockStarts]
            )
        bed = bed + txout + '\n'
    return bed


def write_bed(bed,target):
    '''
    takes bed string and target file name, writes string to target file
    '''
    with open(target,'w') as outfile:
        outfile.write(bed)


def dct_to_bed_write(l,d,target):
    '''
    -- takes list of selected transcripts and dict of all transcripts along with
    child exons; for each selected transcript, selects child exons, assembles
    BED12 record for selected transcript from child exons, writes ouput in BED12
    '''
    bed = ''
    with open(target,'w') as outfile:
        for tx in l:
            ex_sorted = sorted(d[tx], key=lambda ex: int(ex['start']))
            chrom = ex_sorted[0]['seqname']
            chromStart = str(int(ex_sorted[0]['start']) - 1)
            chromEnd = ex_sorted[-1]['end']
            name = tx
            score = '.'
            strand = ex_sorted[0]['strand']
            thickStart = '.'
            thickEnd = '.'
            itemRgb = '.'
            blockCount = str(len(ex_sorted))
            blockSizes = ','.join([str(ex['length']) for ex in ex_sorted])
            blockStarts = ','.join(
                [str(int(ex['start']) - int(ex_sorted[0]['start']))
                    for ex in ex_sorted]
                )
            txout = '\t'.join(
                [chrom, chromStart, chromEnd, name, score, strand, thickStart,
                    thickEnd, itemRgb, blockCount, blockSizes, blockStarts]
                )
            outfile.write(txout + '\n')


def messages(k,target=''):
    '''
    messages to keep user informed
    '''
    d = {
        'startup': ' '.join(['Launching',modname,'module']),
        'in': ' '.join(['Parsing input from',target]),
        'out': ' '.join(['Writing output to',target]),
        'end': ' '.join(['Module',modname,'done!']),
        'exit': 'Cleaning up and exiting ...'
        }
    print(d[k])


def main():
    '''
    -- converts Attig ERE GTF into BED for unannotated HERVs
    -- dct_to_bed(txl,exdct) and write_bed(bed,target) create output in
    unpredictable order, dct_to_bed_write(txl,exdct,target) creates output
    line by line in order of GTF lines
    '''
    messages('startup')
    fp = os.environ['attig_dir']
    fn = 'attig_flat_transcriptome.gtf.gz'
    target = fp + fn
    messages('in',target)
    txdct,exdct = gtf_to_dct(target)
    txl = select_tx(txdct)
    fp = os.environ['genmod_dir']
    fn = 'herv_transcriptome.bed'
    target = fp + fn
    messages('out',target)
    dct_to_bed_write(txl,exdct,target)
    messages('end')
    messages('exit')


if __name__ == '__main__':
    main()
else:
    print('functions loaded from', modname)
