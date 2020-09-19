'''
this script assembles a pheno file for a given cohort based on annotations from
dampier et al. field effect study and updated tcga annotations
-- executes on front end
-- script should be called as part of snakemake pipeline
-- used from within snakemake, so print statements not part of interactive
testing or module loading should include flush=True
-- usage: snakemake setup_dge
'''


import os
import pandas as pd
import argparse


modname = 'pheno'


def fix_columns(df, dftcga):
    '''
    1. drop following superfluous columns from A, B, and C:
    -- sampId, dist_cm, daysToDeath, wt_kg, ht_cm, rin, sequencer, platform
    -- percDup, percGc, seqLen, rdProc, rdMap, percMap
    2. add msi_status to A, B, C
    '''
    l1 = ['dirName', 'projId', 'subId', 'sampType', 'sex', 'race', 'tStage',
            'ageAtDiagDays', 'sampSite', 'rnaMethod', 'format', 'study', 'data']
    l2 = ['uuid', 'proj_id', 'subject_id', 'phenotype', 'sex', 'race', 'stage',
            'age_days', 'subsite', 'rna_method', 'format', 'study', 'source']
    l3 = ['uuid', 'proj_id', 'subject_id', 'phenotype', 'sex', 'race', 'stage',
            'age_days', 'subsite', 'msi_status', 'rna_method', 'format',
            'study', 'source']
    df = df[l1]
    df.columns = l2
    dftcga = dftcga.query('tissue_type == "Primary Tumor"')
    l = ['subject_id', 'msi_status']
    dftcga = dftcga[l]
    df = df.merge(dftcga, how='left', on='subject_id')
    df = df[l3]
    return df


def add_tcga(dftarget, dfref, dftcga, id_in, id_ref):
    '''
    -- for given subject, takes id of file to include in target cohort and id of
    file to exclude from reference cohort, copies subject information from
    reference cohort, updates file information for file to include based on tcga
    annotation, then adds information for file to target cohort and returns
    updated target cohort
    '''
    q = '=='.join(['uuid', '"%s"' % id_ref])
    df = dfref.query(q)
    targetinfo = df.to_dict('records')[0]
    q = '=='.join(['file_id', '"%s"' % id_in])
    df = dftcga.query(q)
    d = {'tissue_type': {'Primary Tumor': 'CRC', 'Solid Tissue Normal': 'NAT'}}
    df = df.replace(d)
    tcgainfo = df.to_dict('records')[0]
    targetinfo['uuid'] = tcgainfo['file_id']
    targetinfo['phenotype'] = tcgainfo['tissue_type']
    targetinfo['format'] = tcgainfo['read_format']
    df = pd.DataFrame.from_dict(targetinfo, orient='index')
    df = df.transpose()
    dftarget = pd.concat([dftarget, df], ignore_index=True)
    return dftarget


def fix_rows_A(dfa, dfc, dftcga):
    '''
    add rows:
    -- 18625fe4-3c19-45d9-9d7c-a295fbf83f2e (TCGA-A6-2684)
    drop rows:
    -- b7e78979-c177-427d-a627-b5f241f36232 (TCGA-A6-2677)
    -- ff6d6688-c19c-4a7f-8058-4d1bc0249d83 (TCGA-A6-2672)
    '''
    id_in = '18625fe4-3c19-45d9-9d7c-a295fbf83f2e'
    id_ref = 'ad1cc824-7bbc-4883-96da-55e51d1b6f0d'
    dfa = add_tcga(dfa, dfc, dftcga, id_in, id_ref)
    id1 = 'b7e78979-c177-427d-a627-b5f241f36232'
    id2 = 'ff6d6688-c19c-4a7f-8058-4d1bc0249d83'
    cond1 = '!='.join(['uuid', '"%s"' % id1])
    cond2 = '!='.join(['uuid', '"%s"' % id2])
    q = '&'.join([cond1, cond2])
    dfa = dfa.query(q)
    return dfa


def fix_rows_C(dfc, dfa, dftcga):
    '''
    add rows:
    -- 90832632-cf57-463b-9d08-c76975066f56 (TCGA-A6-2677)
    -- f08dc7f4-3cc3-4743-a84e-d586d74af8d1 (TCGA-A6-2672)
    drop rows:
    -- ad1cc824-7bbc-4883-96da-55e51d1b6f0d (TCGA-A6-2684)
    '''
    id_in = '90832632-cf57-463b-9d08-c76975066f56'
    id_ref = 'b7e78979-c177-427d-a627-b5f241f36232'
    dfc = add_tcga(dfc, dfa, dftcga, id_in, id_ref)
    id_in = 'f08dc7f4-3cc3-4743-a84e-d586d74af8d1'
    id_ref = 'ff6d6688-c19c-4a7f-8058-4d1bc0249d83'
    dfc = add_tcga(dfc, dfa, dftcga, id_in, id_ref)
    id = 'ad1cc824-7bbc-4883-96da-55e51d1b6f0d'
    q = '!='.join(['uuid', '"%s"' % id])
    dfc = dfc.query(q)
    return dfc


def read_inputs(cohort, tcga, a, b, c):
    '''
    start with field effect annotations, modify columns for herv study, revise
    cohorts A and C based on updated tcga annotations
    '''
    dftcga = pd.read_csv(tcga, sep='\t')
    dfa = pd.read_csv(a, sep='\t')
    dfb = pd.read_csv(b, sep='\t')
    dfc = pd.read_csv(c, sep='\t')
    dfa = fix_columns(dfa, dftcga)
    dfb = fix_columns(dfb, dftcga)
    dfc = fix_columns(dfc, dftcga)
    if cohort == 'A':
        dfa = fix_rows_A(dfa, dfc, dftcga)
        dfmain = dfa
    elif cohort == 'B':
        dfmain = dfb
    elif cohort == 'C':
        dfc = fix_rows_C(dfc, dfa, dftcga)
        dfmain = dfc
    else:
        raise Exception('Unrecognized cohort:', cohort)
    return dftcga, dfa, dfb, dfc, dfmain


def set_paths(df):
    '''
    set paths for tximport
    '''
    qp = os.environ['salqnt_dir']
    df['path'] = pd.Series()
    for i in df.index:
        if df.loc[i,'source'] != 'gdc':
            sd = df.loc[i,'uuid']
        else:
            if df.loc[i,'format'] == 'paired':
                if df.loc[i,'phenotype'] == 'NAT':
                    sd = '-'.join([df.loc[i,'subject_id'], 'Normal'])
                else:
                    sd = '-'.join([df.loc[i,'subject_id'], 'Tumor'])
            else:
                sd = df.loc[i,'subject_id']
        path = ''.join([qp, sd, '/quant.sf'])
        if os.path.exists(path):
            df.loc[i,'path'] = path
        else:
            raise Exception(path, 'does not exist')
    return df


def write_pheno(df, fn):
    fp = os.environ['ann_dir']
    target = ''.join([fp, fn])
    df.to_csv(target, sep='\t', index=False)


def main(args):
    '''
    load annotations and revise cohorts, set paths for each cohort, write pheno
    file for each cohort with field for tximport path
    '''
    cohort = args.cohort
    tcga = args.tcga_ann
    a = args.tsv_a
    b = args.tsv_b
    c = args.tsv_c
    dftcga, dfa, dfb, dfc, dfmain = read_inputs(cohort, tcga, a, b, c)
    df = set_paths(dfmain)
    fn = 'pheno_%s.tsv' % cohort
    write_pheno(df, fn)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='pheno')
    parser.add_argument(
        '--cohort', help='A, B, or C', action='store', dest='cohort'
    )
    parser.add_argument(
        '--tcga_ann', help='path to TCGA annotations', action='store', dest='tcga_ann'
    )
    parser.add_argument(
        '--tsv_a', help='path to cohort A TSV', action='store', dest='tsv_a'
    )
    parser.add_argument(
        '--tsv_b', help='path to cohort B TSV', action='store', dest='tsv_b'
    )
    parser.add_argument(
        '--tsv_c', help='path to cohort C TSV', action='store', dest='tsv_c'
    )
    parser.set_defaults(
        cohort='A',
        tcga_ann='/scratch/chd5n/herv/annotations/tcga/ann_select.tsv',
        tsv_a='/scratch/chd5n/field-effect/annotations/cohort_A.tsv',
        tsv_b='/scratch/chd5n/field-effect/annotations/cohort_B.tsv',
        tsv_c='/scratch/chd5n/field-effect/annotations/cohort_C.tsv'
    )
    args = parser.parse_args()
    main(args)
else:
    print('functions loaded for', modname)
