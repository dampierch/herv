'''
-- module built for python3
-- identify bam files for conversion to fastq
-- submit jobs for biobambam conversion on rivanna
-- requires config.sh to set environmental variables
-- for use in snakemake workflow
-- usage: snakemake convert_tcga
'''


import os
import glob
import subprocess
import argparse
import math


modname = 'gdc_biobambam'


def run_mod(**par):
    '''
    function to execute submission of sbatch job associated with this module
    '''
    os.chdir(par['command_dir'])
    fp = par['out_path']
    fn = '.'.join([modname, 'idx'])
    with open(fp + fn, 'w') as outfile:
        for k in par['file_dict']:
            outfile.write('%s::%s\n' % (k, par['file_dict'][k]))
    output = '_'.join([modname, '%a.out'])
    shell_args = ' '.join(
            [
                par['in_path'], par['ext'], par['read_format'],
                str(par['batch_size']), par['out_path']
            ]
        )
    cmd = ' '.join(
            [
                'sbatch', '--output=' + output,
                '--array=1-' + str(par['array_max']),
                '.'.join([modname,'sh']), shell_args
            ]
        )
    if par['live'] == True:
        subprocess.call(cmd, shell=True)
    else:
        print(cmd)


def get_ids(ann, read_format):
    '''
    load selected file_ids and subject_ids into dict for processing
    '''
    ids = {}
    line_num = 0
    with open(ann, 'r') as in_file:
        for in_line in in_file:
            if line_num == 0:
                fields = in_line.strip('\n').split('\t')
            else:
                d = dict(list(zip(fields, in_line.strip('\n').split('\t'))))
                if d['read_format'] == read_format:
                    k = '-'.join([d['subject_id'],d['tissue_type'].split()[-1]])
                    ids[k] = d['file_id']
            line_num = line_num + 1
    return ids


def parse_filesets(
        ids, in_path, infile_ext, read_format, out_path, outfile_ext,
        batch_size
    ):
    '''
    determine which files to process and how many batches are required
    -- return files to process in dict with subject_id as key
    '''
    fs_in = glob.glob(in_path + '*/*' + infile_ext)
    d_in = {}
    for f in fs_in:
        if f[len(in_path):].split('/')[0] in ids.values():
            d_in[f[len(in_path):].split('/')[0]] = f[len(in_path):]
    fs_out = glob.glob(out_path + '*' + outfile_ext)
    l_out = []
    for f in fs_out:
        if f[len(out_path):-len(outfile_ext)].split('_')[0] in ids:
            l_out.append(f[len(out_path):-len(outfile_ext)].split('_')[0])
    filenum = len(d_in) - len(set(l_out))  ## unprocessed minus processed
    array_max = math.ceil(filenum / batch_size)
    l_out = [d_in[ids[id]] for id in set(l_out)]
    fs_done = [f for f in l_out]
    fs_todo = [d_in[f] for f in d_in if d_in[f] not in fs_done]
    ids = {v: k for k, v in ids.items()}  ## invert ids dict to write todo list
    d_todo = {}
    for f in fs_todo:
        d_todo[ids[f.split('/')[0]]] = ''.join([in_path, f])
    return array_max, d_todo


def parse_params(
        command_dir, in_path, infile_ext, read_format, batch_size,
        out_path, array_max, d_todo, run_type
    ):
    '''
    set parameters
    '''
    params = {
        'command_dir': command_dir,
        'in_path': in_path,
        'ext': infile_ext,
        'read_format': read_format,
        'batch_size': batch_size,
        'out_path': out_path,
        'array_max': array_max,
        'file_dict': d_todo,
        'live': run_type
    }
    return params


def set_source(modname, source):
    '''
    set input path and file extension depending on module and source
    '''
    input_path_dict = {
        'gdc_biobambam': {
            'tcga': os.environ['download_dir']
        }
    }
    input_ext_dict = {
        'gdc_biobambam': {
            'tcga': '.bam'
        }
    }
    return input_path_dict[modname][source], input_ext_dict[modname][source]


def set_dest(modname):
    '''
    set output path and file extension depending on module
    '''
    out_path_dict = {
        'gdc_biobambam': os.environ['fastq_dir']
    }
    out_ext_dict = {
        'gdc_biobambam': '.fastq.gz'
    }
    return out_path_dict[modname], out_ext_dict[modname]


def check_outpath(out_path, live):
    '''
    check for presence of absence of out_path and make directory if absent
    '''
    l = out_path.strip('/').split('/')
    d = ''
    for e in l:
        d = d + '/' + e
        if os.path.exists(d):
            print(d,'present')
        else:
            print(d,'absent')
            print('making',d,'now')
            if live:
                os.mkdir(d)


def main(args):
    '''
    function to organize execution of tasks across batches
    '''
    command_dir = os.environ['command_dir']
    work_dir = os.environ['work_dir']
    in_path, infile_ext = set_source(modname, 'tcga')
    out_path, outfile_ext = set_dest(modname)
    read_format = args.read_format
    batch_size = int(args.batch_size)
    run_type = args.run_type
    ann = os.environ['ann_dir'] + 'tcga/' + 'ann_select.tsv'
    ids = get_ids(ann, read_format)
    check_outpath(out_path, run_type)
    array_max, d_todo = parse_filesets(
            ids, in_path, infile_ext, read_format, out_path, outfile_ext,
            batch_size
        )
    params = parse_params(
            command_dir, in_path, infile_ext, read_format, batch_size,
            out_path, array_max, d_todo, run_type
        )
    run_mod(**params)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='biobambam')
    parser.add_argument(
        '--read_format', help='paired or single reads', action='store',
        dest='read_format'
    )
    parser.add_argument(
        '--batch_size', help='number of samples to process in single batch',
        action='store', dest='batch_size'
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
        read_format='paired', batch_size=10, run_type=False
    )
    args = parser.parse_args()
    main(args)
else:
    print('functions loaded for', modname)
