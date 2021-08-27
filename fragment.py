import concurrent.futures
from random import randint
from Bio import SeqIO
import gzip
import os
import os.path
from math import log

from progress_bar import create_progress_bars, update_progress_bar
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
        with open(f'{frag_path}/{batch_i}', 'w+') as fout:
            for n, frag in input_arr:
                fout.write(f'>{n} {frag["id"]}\n{frag["seq"]}\n')
        write_file_id(batch_i, input_arr, plas_target_path, err_file)
    except Exception as err:
        with open(err_file, 'a') as fout:
            fout.write(f"{timestamp()} Error writing plasmid file for batch {batch_i}: {err}\n")

def write_chrom_frags(batch_i, input_arr, frag_path, err_file):
    try:
        with open(f'{frag_path}/{batch_i}', 'w+') as fout:
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

    plas_bar, chrom_bar, ex_plas_bar = create_progress_bars()

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

if __name__ == "__main__":
    process()


# def plasmid_worker(
#     input, k, frag_path, 
#     target_path,
#     kmer_write_path, kmer_out_path,
#     circ_write_path, 
#     db_path, inc_out_path, inc_write_path,
#     rrna_out_path, rrna_write_path,
#     orit_out_path, orit_write_path,
#     progress_bar, err_file, log_file):

#     n, frag = input

#     write_plas_frags(input, frag_path, err_file)
#     write_file_id(input, target_path, err_file)
#     count_kmers(k, input, frag_path, kmer_write_path, kmer_out_path, err_file)
#     circularity(input, circ_write_path, err_file)
#     inc_factor(input, db_path, frag_path, inc_out_path, inc_write_path, err_file)
#     rrna_search(input, frag_path, rrna_out_path, rrna_write_path, db_path, err_file)
#     orit_search(input, frag_path, db_path, orit_out_path, orit_write_path, err_file)

#     if(n%10 == 0):
#         update_progress_bar(progress_bar, 10)

#     if(int(log(n, 2)) == log(n, 2)):
#         with open(log_file, 'a') as fout:
#             fout.write(f"{timestamp()} Completed writing plasmid file {n}\n")

# def chrom_worker(
#     input, k, frag_path,
#     target_path,
#     kmer_write_path, kmer_out_path,
#     circ_write_path,
#     db_path, inc_out_path, inc_write_path,
#     rrna_out_path, rrna_write_path,
#     orit_out_path, orit_write_path,
#     progress_bar, err_file, log_file):

#     n, frag = input

#     write_chrom_frags(input, frag_path, err_file)
#     write_file_id(input, target_path, err_file)
#     count_kmers(k, input, frag_path, kmer_write_path, kmer_out_path, err_file)
#     circularity(input, circ_write_path, err_file)
#     inc_factor(input, db_path, frag_path, inc_out_path, inc_write_path, err_file)
#     rrna_search(input, frag_path, rrna_out_path, rrna_write_path, db_path, err_file)
#     orit_search(input, frag_path, db_path, orit_out_path, orit_write_path, err_file)

#     if(n%10 == 0):
#         update_progress_bar(progress_bar, 10)

#     if(int(log(n, 2)) == log(n, 2)):
#         with open(log_file, 'a') as fout:
#             fout.write(f"{timestamp()} Completed writing chromosome file {n}\n")

# def process():
    
#     ## create file structure
#     create_fragment_dirs()

#     ## read and write plasmid/ chromosome files
#     plas_frag_gen = plas_frag_generator(plas_db_path, frag_len, plas_frag_split_path, plas_circ_out_path, coverage, err_file)
#     chrom_frag_gen = chrom_frag_generator(chrom_db_path, frag_len, chrom_frag_split_path, chrom_circ_out_path, coverage, err_file)
#     chrom_frag_gen_c = filter(lambda x: x["type"] == "chromosome", chrom_frag_gen)
#     chrom_frag_gen_p = filter(lambda x: x["type"] == "plasmid", chrom_frag_gen)

#     plas_bar, chrom_bar, ex_plas_bar = create_progress_bars()

#     batch_i = 0
#     for input in enumerate(plas_frag_gen):
#         n, frag = input
#         if(n % batch_size == 0):
#             input_arr = []
#         input_arr.append(input)
#         if ((n+1)%batch_size == 0):
#             write_plas_frags(batch_i, input_arr, plas_write_path, err_file)
#             write_file_id(batch_i, input_arr, plas_target_path, err_file)
#             batch_i+=1
#     if not ((n+1)%batch_size == 0):
#         write_plas_frags(batch_i, input_arr, plas_write_path, err_file)
#         write_file_id(batch_i, input_arr, plas_target_path, err_file)

#     batch_i = 0
#     for input in enumerate(chrom_frag_gen_c):
#         n, frag = input
#         if(n % batch_size == 0):
#             input_arr = []
#         input_arr.append(input)
#         if ((n+1)%batch_size == 0):
#             write_chrom_frags(batch_i, input_arr, chrom_write_path, err_file)
#             write_file_id(batch_i, input_arr, chrom_target_path, err_file)
#             batch_i+=1
#     if not ((n+1)%batch_size == 0):
#         write_chrom_frags(batch_i, input_arr, chrom_write_path, err_file)
#         write_file_id(batch_i, input_arr, chrom_target_path, err_file)
    
#     batch_i = 0
#     for input in enumerate(chrom_frag_gen_p):
#         n, frag = input
#         if(n % batch_size == 0):
#             input_arr = []
#         input_arr.append(input)
#         if ((n+1)%batch_size == 0):
#             write_chrom_frags(batch_i, input_arr, extra_plasmid_write_path, err_file)
#             write_file_id(batch_i, input_arr, ex_plas_target_path, err_file)
#             batch_i+=1
#     if not ((n+1)%batch_size == 0):
#         write_chrom_frags(batch_i, input_arr, extra_plasmid_write_path, err_file)
#         write_file_id(batch_i, input_arr, ex_plas_target_path, err_file)


            
    # with concurrent.futures.ThreadPoolExecutor(max_workers=thread_count) as executor:
    #     try:
    #         executor.map(
    #             plasmid_worker, enumerate(plas_frag_gen), repeat(k), repeat(plas_write_path),
    #             repeat(plas_target_path),
    #             repeat(plas_7mer_write_path), repeat(plas_7mer_out_path),
    #             repeat(plas_circ_write_path),
    #             repeat(db_path), repeat(plas_inc_out_path), repeat(plas_inc_write_path),
    #             repeat(plas_rrna_out_path), repeat(plas_rrna_write_path),
    #             repeat(plas_orit_out_path), repeat(plas_orit_write_path),
    #             repeat(plas_bar),
    #             repeat(err_file), repeat(log_file))
    #     except Exception as err:
    #         with open(err_file, 'a') as fout:
    #             fout.write(f"{timestamp()} Concurrency Error: {err}\n")

    #     try:
    #         executor.map(
    #             chrom_worker, enumerate(chrom_frag_gen_c), repeat(k), repeat(chrom_write_path),
    #             repeat(chrom_target_path),
    #             repeat(chrom_7mer_write_path), repeat(chrom_7mer_out_path),
    #             repeat(chrom_circ_write_path),
    #             repeat(db_path), repeat(chrom_inc_out_path), repeat(chrom_inc_write_path),
    #             repeat(chrom_rrna_out_path), repeat(chrom_rrna_write_path),
    #             repeat(chrom_orit_out_path), repeat(chrom_orit_write_path),
    #             repeat(chrom_bar),
    #             repeat(err_file), repeat(log_file))
    #     except Exception as err:
    #         with open(err_file, 'a') as fout:
    #             fout.write(f"{timestamp()} Concurrency Error: {err}\n")

    #     try:
    #         executor.map(chrom_worker, enumerate(chrom_frag_gen_p), repeat(k), repeat(extra_plasmid_write_path), 
    #         repeat(ex_plas_target_path),
    #         repeat(ex_plas_7mer_write_path), repeat(ex_plas_7mer_out_path),
    #         repeat(ex_plas_circ_write_path),
    #         repeat(db_path), repeat(chrom_inc_out_path), repeat(chrom_inc_write_path),
    #         repeat(ex_plas_rrna_out_path), repeat(ex_plas_rrna_write_path),
    #         repeat(ex_plas_orit_out_path), repeat(ex_plas_orit_write_path),
    #         repeat(ex_plas_bar),
    #         repeat(err_file), repeat(log_file))
    #     except Exception as err:
    #         with open(err_file, 'a') as fout:
    #             fout.write(f"{timestamp()} Concurrency Error: {err}\n")

# if __name__ == "__main__":
#     process()
