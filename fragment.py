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
    
def plas_frag_generator(dir_path, frag_len, max_frag):
    frag_count = 0
    try:
        while(True):
            for filename in os.listdir(dir_path):
                for record in SeqIO.parse(f"{dir_path}/{filename}", 'fasta'):
                    length = len(record.seq)
                    if(length<frag_len):
                        rand_i = 0
                    else:
                        rand_i = randint(0, length-frag_len)
                    yield {"id":record.description, "seq":str(record.seq)[rand_i:rand_i+min(frag_len, length)]}
                    frag_count+=1
                    if(frag_count>=max_frag):
                        return
    except Exception as err:
        with open(err_file, 'a') as fout:
            fout.write(f"{timestamp()} Error reading plasmid file {filename}: {err}\n")

def chrom_frag_generator(dir_path, frag_len, max_frag, type):
    frag_count = 0
    try:
        while(True):
            for filename in os.listdir(dir_path):
                file_comp = filename.split("_")
                name = file_comp[0]+"_"+file_comp[1]
                with gzip.open(f"{dir_path}/{filename}", "rt") as handle:
                    for record in SeqIO.parse(handle, 'fasta'):
                        length = len(record.seq)
                        if((type== "chromosome" and "plasmid" in record.description)
                            or (type== "plasmid" and "chromosome" in record.description)):
                            continue
                        if(length<frag_len):
                            rand_i = 0
                        else:
                            rand_i = randint(0, length-frag_len)
                        yield {
                            "id":record.description, "name":name, 
                            "seq":str(record.seq)[rand_i:rand_i+min(frag_len,length)],
                            }
                        frag_count+=1
                        if(frag_count>=max_frag):
                            return
    except Exception as err:
        with open(err_file, 'a') as fout:
            fout.write(f"{timestamp()} Error reading chromosome file {filename}: {err}, {len(record)}\n")

def write_plas_frags(batch_i, input_arr, frag_path, err_file):
    try:
        with open(f'{frag_path}/{batch_i}.fasta', 'w+') as fout:
            for n, frag in input_arr:
                fout.write(f'>{n} {frag["id"]}\n{frag["seq"]}\n')
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
    plas_frag_gen = plas_frag_generator(plas_db_path, frag_len, plas_max)
    chrom_frag_gen_c = chrom_frag_generator(chrom_db_path, frag_len, chrom_max, "chromosome")
    chrom_frag_gen_p = chrom_frag_generator(chrom_db_path, frag_len, ex_plas_max, "plasmid")

    # chrom_frag_gen_c = filter(lambda x: x["type"] == "chromosome", chrom_frag_gen)
    # chrom_frag_gen_p = filter(lambda x: x["type"] == "plasmid", chrom_frag_gen)

    plas_bar = create_progress_bar(plas_bar_desc)
    batch_i = 0
    n, input_arr = 0, []
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
    n, input_arr = 0, []
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
    n, input_arr = 0, []
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
