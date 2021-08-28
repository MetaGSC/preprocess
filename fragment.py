import concurrent.futures
from random import randint
from Bio import SeqIO
import gzip
import os
import os.path
from math import log

from progress_bar import *
from constants import *
from helpers import create_fragment_dirs, timestamp
from write_file_id import write_file_id
    
def plas_frag_generator(dir_path, frag_len, coverage, err_file):
    for filename in os.listdir(dir_path):
        try:
            for record in SeqIO.parse(f"{dir_path}/{filename}", 'fasta'):
                length = len(record.seq)
                for _ in range(int(length*coverage/frag_len)):
                    if(length<frag_len):
                        rand_i = 0
                    else:
                        rand_i = randint(0, length-frag_len)
                    yield {"id":record.description, "seq":str(record.seq)[rand_i:rand_i+min(frag_len, length)]}
        except Exception as err:
            with open(err_file, 'a') as fout:
                fout.write(f"{timestamp()} Error reading plasmid file {filename}: {err}\n")

def chrom_frag_generator(dir_path, frag_len, coverage, err_file):
    for filename in os.listdir(dir_path):
        try:
            file_comp = filename.split("_")
            name = file_comp[0]+"_"+file_comp[1]
            with gzip.open(f"{dir_path}/{filename}", "rt") as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    length = len(record.seq)
                    if("plasmid" in record.description): 
                        type = "plasmid"
                    else:
                        type = "chromosome"
                    for _ in range(int(len(record.seq)*coverage/frag_len)):
                        if(length<frag_len):
                            rand_i = 0
                        else:
                            rand_i = randint(0, length-frag_len)
                        yield {
                            "id":record.description, "name":name, 
                            "seq":str(record.seq)[rand_i:rand_i+min(frag_len,length)], "type": type,
                            }
        except Exception as err:
            with open(err_file, 'a') as fout:
                fout.write(f"{timestamp()} Error reading chromosome file {filename}: {err}, {len(record)}\n")

def write_plas_frags(batch_i, input_arr, frag_path, err_file):
    try:
        with open(f'{frag_path}/{batch_i}.fasta', 'w+') as fout:
            for n, frag in input_arr:
                fout.write(f'>{n} {frag["id"]}\n{frag["seq"]}\n')
        write_file_id(batch_i, input_arr, plas_target_path, err_file)
    except Exception as err:
        with open(err_file, 'a') as fout:
            fout.write(f"{timestamp()} Error writing plasmid file for batch {batch_i}: {err}\n")

def write_chrom_frags(batch_i, input_arr, frag_path, err_file):
    try:
        with open(f'{frag_path}/{batch_i}.fasta', 'w+') as fout:
            for n, frag in input_arr:
                fout.write(f'>{n} {frag["name"]} {frag["id"]}\n{frag["seq"]}\n')

    except Exception as err:
        with open(err_file, 'a') as fout:
            fout.write(f"{timestamp()} Error writing plasmid file for batch {batch_i}: {err}\n")

def process():
    
    ## create file structure
    create_fragment_dirs()

    ## read and write plasmid/ chromosome files
    plas_frag_gen = plas_frag_generator(plas_db_path, frag_len, coverage, err_file)
    chrom_frag_gen = chrom_frag_generator(chrom_db_path, frag_len, coverage, err_file)
    chrom_frag_gen_c = filter(lambda x: x["type"] == "chromosome", chrom_frag_gen)
    chrom_frag_gen_p = filter(lambda x: x["type"] == "plasmid", chrom_frag_gen)

    plas_bar = create_progress_bar(plas_bar_desc)
    batch_i = 0
    for input in enumerate(plas_frag_gen):
        n, _ = input
        if(n % batch_size == 0):
            input_arr = []
        input_arr.append(input)
        if ((n+1)%batch_size == 0):
            write_plas_frags(batch_i, input_arr, plas_write_path, err_file)
            write_file_id(batch_i, input_arr, plas_target_path, err_file)
            batch_i+=1
        if(n%10 == 0):
            update_progress_bar(plas_bar, 10)

    if not ((n+1)%batch_size == 0):
        write_plas_frags(batch_i, input_arr, plas_write_path, err_file)
        write_file_id(batch_i, input_arr, plas_target_path, err_file)
    close_progress_bar(plas_bar)

    chrom_bar = create_progress_bar(chrom_bar_desc)
    batch_i = 0
    for input in enumerate(chrom_frag_gen_c):
        n, _ = input
        if(n % batch_size == 0):
            input_arr = []
        input_arr.append(input)
        if ((n+1)%batch_size == 0):
            write_chrom_frags(batch_i, input_arr, chrom_write_path, err_file)
            write_file_id(batch_i, input_arr, chrom_target_path, err_file)
            batch_i+=1
        if(n%10 == 0):
            update_progress_bar(chrom_bar, 10)
        
    if not ((n+1)%batch_size == 0):
        write_chrom_frags(batch_i, input_arr, chrom_write_path, err_file)
        write_file_id(batch_i, input_arr, chrom_target_path, err_file)
    close_progress_bar(chrom_bar)

    ex_plas_bar = create_progress_bar(ex_plas_bar_desc)
    batch_i = 0
    for input in enumerate(chrom_frag_gen_p):
        n, _ = input
        if(n % batch_size == 0):
            input_arr = []
        input_arr.append(input)
        if ((n+1)%batch_size == 0):
            write_chrom_frags(batch_i, input_arr, extra_plasmid_write_path, err_file)
            write_file_id(batch_i, input_arr, ex_plas_target_path, err_file)
            batch_i+=1
        if(n%10 == 0):
            update_progress_bar(ex_plas_bar, 10)

    if not ((n+1)%batch_size == 0):
        write_chrom_frags(batch_i, input_arr, extra_plasmid_write_path, err_file)
        write_file_id(batch_i, input_arr, ex_plas_target_path, err_file)
    close_progress_bar(ex_plas_bar)

if __name__ == "__main__":
    process()
