## Single threaded circular.py is faster than the multiprocessing version here

import subprocess as sp
import concurrent.futures
from itertools import repeat
from Bio import SeqIO
import os
import math

from progress_bar import *
from helpers import timestamp
from constants import *
from helpers import create_circ_dirs, find_mean, delete_dir_if_exist

def circ_worker(filename, frag_path, split_path, out_path, write_path, pb_i, file_count):
    name = filename.split(".")[0]
    try:
        batch = []
        for record in SeqIO.parse(f"{frag_path}/{filename}", 'fasta'):
            n = int(record.id)
            seq_mid = int(len(record.seq)/2)
            seq_a = str(record.seq[:seq_mid])
            seq_b = str(record.seq[seq_mid:])

            seq_a_file = f'{split_path}/{n}_a.fasta'
            seq_b_file = f'{split_path}/{n}_b.fasta'
            with open(seq_a_file, mode='w+') as fout:
                fout.write(f'>{n}_a\n')
                fout.write(seq_a + '\n ')
            
            with open(seq_b_file, mode='w+') as fout:
                fout.write(f'>{n}_b\n')
                fout.write(seq_b + '\n ')

            out_file = f"{out_path}/{n}_temp"
            
            cmd = [
            nucmer_path,
            '-f',  # only forward strand
            '-l', '40',  # increase min match length to 40 bp
            '-p', out_file,
            seq_a_file,
            seq_b_file
            ]
            proc = sp.run(
                cmd,
                stdout=sp.PIPE,
                stderr=sp.PIPE,
                universal_newlines=True,
            )
            alignment_a_mean, alignment_b_mean, mismatches_mean, count = 0, 0, 0, 0
            if(proc.returncode == 0):
                has_match = False
                with open(f"{out_file}.delta") as fout:
                    for line in fout:
                        line = line.rstrip()
                        if(line[0] == '>'):
                            has_match = True
                        elif(has_match):
                            cols = line.split(' ')
                            if(len(cols) == 7):
                                start_b = int(cols[0])
                                end_b = int(cols[1])
                                start_a = int(cols[2])
                                end_a = int(cols[3])
                                mismatches = int(cols[4])
                                alignment_a = end_a - start_a + 1
                                alignment_b = end_b - start_b + 1
                                mismatches_mean += mismatches
                                alignment_a_mean += alignment_a
                                alignment_b_mean += alignment_b
                                count+=1
            
            find_mean(alignment_a_mean, count)
            find_mean(alignment_b_mean, count)
            find_mean(mismatches_mean, count)

            data = f"{n}\t{alignment_a_mean}\t{alignment_b_mean}\t{mismatches_mean}\t{count}\n"
            batch.append(data)

            if(pb_i == 1): progress_bar = plas_bar
            elif(pb_i == 2): progress_bar = chrom_bar
            else: progress_bar = ex_plas_bar
            if(n%10 == 0):
                update_progress_bar(progress_bar, 10*min(file_count, thread_count))

            delete_dir_if_exist(seq_a_file)
            delete_dir_if_exist(seq_b_file)
            delete_dir_if_exist(f"{out_file}.delta")
            delete_dir_if_exist(f"{out_file}.ntref")
            delete_dir_if_exist(f"{out_file}.mgaps")

        with open(f"{write_path}/{name}", "w") as f:
                f.writelines(batch)

    except Exception as err:
        with open(err_file, 'a') as fout:
            fout.write(f"{timestamp()} Error reading file {filename}: {err}\n")

def circularity(frag_path, split_path, out_path, write_path, progress_bar):

    file_list = os.listdir(frag_path)
    file_count = len(file_list)
    executor = concurrent.futures.ProcessPoolExecutor(thread_count)
    futures = [
        # executor.submit(circ_worker, filename, frag_path, split_path, out_path, write_path, progress_bar) for filename in file_list
        executor.submit(circ_worker, filename, frag_path, split_path, out_path, write_path, progress_bar, file_count) for filename in file_list
        ]
    concurrent.futures.wait(futures)
    executor.shutdown(wait=True)
    


if __name__ == "__main__":
    create_circ_dirs()
    plas_bar = create_progress_bar(plas_bar_desc)
    circularity(plas_write_path, plas_frag_split_path, plas_circ_out_path, plas_circ_write_path, 1)
    close_progress_bar(plas_bar)
    chrom_bar = create_progress_bar(chrom_bar_desc)
    circularity(chrom_write_path, chrom_frag_split_path, chrom_circ_out_path, chrom_circ_write_path, 2)
    close_progress_bar(chrom_bar)
    ex_plas_bar = create_progress_bar(ex_plas_bar_desc)
    circularity(extra_plasmid_write_path, ex_plas_frag_split_path, ex_plas_circ_out_path, ex_plas_circ_write_path, 3)
    close_progress_bar(ex_plas_bar)
