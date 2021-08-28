import subprocess as sp

from helpers import timestamp
from constants import seq2vec_path

def count_kmers(k, input, input_path, write_path, out_path, err_file):
    try:
        n, frag = input
        input_file = f'{input_path}/{n}'
        write_file = f"{write_path}/{n}"

        args = [
            seq2vec_path,
            '-f', str(input_file),
            '-o', str(write_file),
            '-t', '8',
            '-k', str(k),
        ]    
        proc = sp.run(args,  stdout=sp.PIPE, stderr=sp.PIPE)
    except Exception as err:
        with open(err_file, 'a') as fout:
            fout.write(f"{timestamp()} Error using seq2vec for file {n} {err}\n")
