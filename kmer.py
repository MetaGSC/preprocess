import subprocess as sp

from progress_bar import create_progress_bars, update_progress_bar
from helpers import timestamp, create_kmer_files
from constants import *
import os

def count_kmers(k, frag_path, write_path, progress_bar):
    try:
        for filename in os.listdir(frag_path):
            name = filename.split(".")[0]
            frag_file = f"{frag_path}/{filename}"
            write_file = f"{write_path}/{name}"
            args = [
                seq2vec_path,
                '-f', str(frag_file),
                '-o', str(write_file),
                '-t', '8',
                '-k', str(k),
            ]    
            proc = sp.run(args, stdout=sp.PIPE, stderr=sp.PIPE)

            update_progress_bar(progress_bar, batch_size)

    except Exception as err:
        with open(err_file, 'a') as fout:
            fout.write(f"{timestamp()} Error using seq2vec for file {n} {err}\n")

if __name__ == "__main__":
    k = 7
    plas_bar, chrom_bar, ex_plas_bar = create_progress_bars()
    create_kmer_files()
    count_kmers(k, plas_write_path, plas_7mer_write_path, plas_bar)
    count_kmers(k, chrom_write_path, chrom_7mer_write_path, chrom_bar)
    count_kmers(k, extra_plasmid_write_path, ex_plas_7mer_write_path, ex_plas_bar)
