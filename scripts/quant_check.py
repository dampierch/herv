'''
-- module built for python3
-- identify samples from field-effect study missing from herv quants
-- for use in snakemake workflow
-- usage: snakemake quant_check
'''


import os
import argparse


modname = 'quant_check'


def get_info(target):
    '''
    read sample information from field-effect study
    '''
    infodict = {}
    keys = ['dirName', 'projId', 'subId', 'sampId', 'format', 'study', 'data']
    with open(target,'r') as infile:
        for inline in infile:
            if inline.strip('\n').split('\t')[0] == 'dirName':
                columns = inline.strip('\n').split('\t')
            else:
                data_list = inline.strip('\n').split('\t')
                data_dict = dict(list(zip(columns,data_list)))
                if data_dict['format'] == 'paired':
                    if data_dict['data'] == 'gdc':
                        if data_dict['sampType'] == 'NAT':
                            infodict['-'.join([data_dict['subId'], 'Normal'])] =  {key: data_dict[key] for key in keys}
                        else:
                            infodict['-'.join([data_dict['subId'], 'Tumor'])] = {key: data_dict[key] for key in keys}
                    else:
                        infodict[data_dict['dirName']] = {key: data_dict[key] for key in keys}
                else:
                    if data_dict['data'] == 'gdc':
                        infodict[data_dict['subId']] = {key: data_dict[key] for key in keys}
                    else:
                        infodict[data_dict['dirName']] = {key: data_dict[key] for key in keys}
    return infodict


def write_missing(target,info,done):
    '''
    write sample info for missing samples
    '''
    l = [info[k] for k in info if k not in done]
    if len(l) > 0:
        with open(target,'w') as outfile:
            outfile.write(''.join(['\t'.join(l[0].keys()),'\n']))
            for d in l:
                outfile.write(''.join(['\t'.join(d.values()),'\n']))
    else:
        print('Nothing to write to', target)
        with open(target,'w') as outfile:
            outfile.write('')


def main(args):
    a = args.tsv_a
    b = args.tsv_b
    c = args.tsv_c
    d = ''.join([os.environ['work_dir'], 'done.txt'])
    ainfo = get_info(a)
    binfo = get_info(b)
    cinfo = get_info(c)
    dlist = []
    with open(d,'r') as infile:
        for inline in infile:
            dlist.append(inline.strip('\n'))
    target = ''.join([os.environ['work_dir'], 'todo_A.tsv'])
    write_missing(target,ainfo,dlist)
    target = ''.join([os.environ['work_dir'], 'todo_B.tsv'])
    write_missing(target,binfo,dlist)
    target = ''.join([os.environ['work_dir'], 'todo_C.tsv'])
    write_missing(target,cinfo,dlist)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='quant_check')
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
        tsv_a='/scratch/chd5n/field-effect/annotations/cohort_A.tsv',
        tsv_b='/scratch/chd5n/field-effect/annotations/cohort_B.tsv',
        tsv_c='/scratch/chd5n/field-effect/annotations/cohort_C.tsv'
    )
    args = parser.parse_args()
    main(args)
else:
    print('functions loaded for', modname)
