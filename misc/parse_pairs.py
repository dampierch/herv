'''
-- draft script to find all tumor-normal sample pairs in field effect meta
   analysis
-- only samples included in analysis cohorts, which were purposefully pruned
   to exclude subject replicates, were downloaded for herv project
-- need to re-download all excluded tumor and normal samples for a proper
   paired analysis
'''

import pandas as pd
import os

def load_data(target):
    df = pd.read_csv(target, sep='\t')
    return df

def filter_data(df, col, val=0):
    if col == 'subId':
        if val == 0:
            dff = df.loc[df.duplicated(subset=[col], keep=False), :]
        else:
            print('no filter\n')
            dff = df
    elif col in ['format']:
        if val != 0:
            dff = df.loc[df.loc[:, col] == val, :]
        else:
            print('no filter\n')
            dff = df
    elif col == 'sampType':
        if val == 'HLT':
            print('filtering out', val, '\n')
            dff = df.loc[df.loc[:, col] != val, :]
        else:
            print('no filter\n')
            dff = df
    else:
        print('no filter\n')
        dff = df
    return dff

def parse_pairs(df):
    l = list()
    for i in df.loc[:, 'subId'].unique():
        idx = df.loc[:, 'subId'] == i
        df1 = df.loc[idx, :]
        if df1.shape[0] == 2:
            if 'NAT' in list(df1.sampType) and 'CRC' in list(df1.sampType):
                l.append(list(df1.sampId))
            else:
                print('no tumor normal pair')
        else:
            if 'NAT' in list(df1.sampType) and 'CRC' in list(df1.sampType):
                l2 = list()
                df2 = df1.loc[df1.loc[:, 'sampType'] == 'NAT', :]
                l2.append(list(df2.loc[df2.rdMap == max(df2.rdMap), 'sampId'])[0])
                df2 = df1.loc[df1.loc[:, 'sampType'] == 'CRC', :]
                l2.append(list(df2.loc[df2.rdMap == max(df2.rdMap), 'sampId'])[0])
                l.append(l2)
            else:
                print('no tumor normal pair')
    return l 

def select_pairs(df, l):
    i6i = pd.Int64Index([], dtype='int64')
    for e in l:
        i6i = i6i.append(df.index[df.sampId == e[0]])
        i6i = i6i.append(df.index[df.sampId == e[1]])
    df1 = df.loc[i6i, :]
    return df1

def set_dir(df):
    p = '/scratch/chd5n/herv/salmon/quant/'
    for i in df.index:
        if df.loc[i, 'data'] == 'gdc':
            if df.loc[i, 'sampType'] == 'NAT':
                df.loc[i, 'd'] = '-'.join([df.loc[i, 'subId'], 'Normal'])
            elif df.loc[i, 'sampType'] == 'CRC':
                df.loc[i, 'd'] = '-'.join([df.loc[i, 'subId'], 'Tumor'])
        else:
            df.loc[i, 'd'] = df.loc[i, 'dirName']
    df.loc[:, 'dir'] = df.d.apply(lambda x: ''.join([p, x, '/']))
    df.loc[:, 'exist'] = df.dir.apply(lambda x: os.path.exists(x))
    return df

def main():
    target = '/'.join([os.environ["HOME"], 'metaPheno.tsv'])
    df = load_data(target)
    df = filter_data(df, 'sampType', 'HLT')
    dfm = filter_data(df, 'subId')
    dfp = filter_data(dfm, 'format', val='paired')
    dfs = filter_data(dfm, 'format', val='single')
    dfpm = filter_data(dfp, 'subId')
    dfsm = filter_data(dfs, 'subId')
    lp = parse_pairs(dfpm)
    ls = parse_pairs(dfsm)
    dfp = select_pairs(df, lp)
    dfs = select_pairs(df, ls)
    df = pd.concat([dfp, dfs])
    df = set_dir(df)
    return df

if __name__ == '__main__':
    df = main()
else:
    print('functions loaded')
