import subprocess as sp

from progress_bar import *
from helpers import timestamp, create_kmer_files
from constants import *
import os

def count_kmers(k, frag_path, write_path, pb_desc):
    try:
        progress_bar = create_progress_bar(pb_desc)
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

        close_progress_bar(progress_bar)
    except Exception as err:
        with open(err_file, 'a') as fout:
            fout.write(f"{timestamp()} Error using seq2vec for file {name} {err}\n")

if __name__ == "__main__":
    k = 7
    create_kmer_files()
    count_kmers(k, plas_write_path, plas_7mer_write_path, plas_bar_desc)
    count_kmers(k, chrom_write_path, chrom_7mer_write_path, chrom_bar_desc)
    count_kmers(k, extra_plasmid_write_path, ex_plas_7mer_write_path, ex_plas_bar_desc)
