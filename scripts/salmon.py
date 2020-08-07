'''
-- module built for python3
-- identify fastq files for input to salmon
-- submit jobs for salmon quantification on rivanna
-- requires config.sh to set environmental variables
-- for use in snakemake workflow
-- usage: snakemake salmon
'''


import os
import glob
import subprocess
import argparse
import math


modname = 'salmon'


def run_mod(**par):
    '''
    function to execute submission of sbatch job associated with this module
    '''
    os.chdir(par['command_dir'])
    fp = par['out_path']
    fn = '.'.join([par['source'], modname, 'idx'])
    with open(fp + fn, 'w') as outfile:
        for f in par['file_list']:
            outfile.write('%s\n' % f)
    output = '_'.join([modname, par['source'], '%a.out'])
    shell_args = ' '.join(
            [
                par['source'], par['in_path'], par['ext'], par['read_format'],
                str(par['batch_size']), par['idx'], par['out_path']
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


def parse_filesets(
        source, in_path, infile_ext, read_format, out_path, outfile_ext,
        batch_size
    ):
    '''
    determine which files to process and how many batches are required
    '''
    fs_in = glob.glob(in_path + '*' + infile_ext)
    if read_format == 'paired':
        l_in = [f[len(in_path):-len(infile_ext)].split('_')[0] for f in fs_in]
    else:
        l_in = [f[len(in_path):-len(infile_ext)] for f in fs_in]
    fs_out = glob.glob(out_path + '*/*' + outfile_ext)
    l_out = [f.split('/')[-2] for f in fs_out if f.split('/')[-2] in l_in]
    if read_format == 'paired':
        if source == 'barcuva':
            r1 = [id + '_R1_val_1' for id in l_out]
            r2 = [id + '_R2_val_2' for id in l_out]
        else:
            r1 = [id + '_1' for id in l_out]
            r2 = [id + '_2' for id in l_out]
        l_out = r1 + r2  ## expand single-file output for pairs
    filenum = len(fs_in) - len(l_out)  ## unprocessed minus processed
    if read_format == 'paired':
        array_max = math.ceil((filenum/2) / batch_size)
    else:
        array_max = math.ceil(filenum / batch_size)
    fs_done = [in_path + f + infile_ext for f in l_out]
    fs_todo = [f for f in fs_in if f not in fs_done]
    return array_max, fs_todo


def parse_params(
        command_dir, source, in_path, infile_ext, read_format, batch_size, idx,
        out_path, array_max, fs_todo, run_type
    ):
    '''
    set parameters
    '''
    params = {
        'command_dir': command_dir,
        'source': source,
        'in_path': in_path,
        'ext': infile_ext,
        'read_format': read_format,
        'batch_size': batch_size,
        'idx': idx,
        'out_path': out_path,
        'array_max': array_max,
        'file_list': fs_todo,
        'live': run_type
    }
    return params


def set_source(modname, source):
    '''
    set input path and file extension depending on module and source
    '''
    input_path_dict = {
        'salmon': {
            'barcuva': os.environ['barcuva_fqdir'],
            'sra': os.environ['fastq_dir'],
            'gtex': os.environ['fastq_dir'],
            'tcga': os.environ['fastq_dir']
        }
    }
    input_ext_dict = {
        'salmon': {
            'barcuva': '.fq.gz',
            'sra': '.fastq.gz',
            'gtex': '.fastq.gz',
            'tcga': '.fastq.gz'
        }
    }
    return input_path_dict[modname][source], input_ext_dict[modname][source]


def set_dest(modname):
    '''
    set output path and file extension depending on module
    '''
    out_path_dict = {
        'salmon': os.environ['salqnt_dir']
    }
    out_ext_dict = {
        'salmon': 'quant.sf'
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
    source = args.source
    in_path, infile_ext = set_source(modname, source)
    out_path, outfile_ext = set_dest(modname)
    batch_size = int(args.batch_size)
    run_type = args.run_type
    check_outpath(out_path, run_type)
    read_format = args.read_format
    array_max, fs_todo = parse_filesets(
            source, in_path, infile_ext, read_format, out_path, outfile_ext,
            batch_size
        )
    params = parse_params(
            command_dir, source, in_path, infile_ext, read_format, batch_size,
            args.idx, out_path, array_max, fs_todo, run_type
        )
    run_mod(**params)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='salmon')
    parser.add_argument(
        '--source', help='source of fastq files', action='store',
        dest='source'
    )
    parser.add_argument(
        '--read_format', help='paired or single reads', action='store',
        dest='read_format'
    )
    parser.add_argument(
        '--batch_size', help='number of samples to process in single batch',
        action='store', dest='batch_size'
    )
    parser.add_argument(
        '--idx', help='path to salmon index files', action='store',
        dest='idx'
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
        source='sra', read_format='paired', batch_size=10, run_type=False
    )
    args = parser.parse_args()
    main(args)
else:
    print('functions loaded for', modname)
