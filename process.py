import concurrent.futures
from random import randint
from Bio import SeqIO
import gzip
import os
import os.path
from itertools import repeat
import subprocess as sp
from pathlib import Path
from math import log

from progress_bar import create_progress_bars, update_progress_bar
from constants import *
from helpers import create_file_structure, timestamp
from circular import check_circularity, circularity
from kmer import count_kmers
from inc_fac import inc_factor
from rrna import rrna_search
from orit import orit_search
from specialgene import special_gene_search
from write_file_id import write_file_id
    
def plas_frag_generator(dir_path, frag_len, split_path, out_path, coverage, err_file):
    for filename in os.listdir(dir_path):
        try:
            for record in SeqIO.parse(f"{dir_path}/{filename}", 'fasta'):
                circular = check_circularity(record, split_path, out_path)
                length = len(record.seq)
                for _ in range(int(length*coverage/frag_len)):
                    if(length<frag_len):
                        rand_i = 0
                    else:
                        rand_i = randint(0, length-frag_len)
                    yield {"id":record.description, "seq":str(record.seq)[rand_i:rand_i+min(frag_len, length)], "circular": circular}
        except Exception as err:
            with open(err_file, 'a') as fout:
                fout.write(f"{timestamp()} Error reading plasmid file {filename}: {err}\n")

def chrom_frag_generator(dir_path, frag_len, split_path, out_path, coverage, err_file):
    for filename in os.listdir(dir_path):
        try:
            file_comp = filename.split("_")
            name = file_comp[0]+"_"+file_comp[1]
            with gzip.open(f"{dir_path}/{filename}", "rt") as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    circular = check_circularity(record, split_path, out_path)
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
                            "id":record.description, "name":name, "seq":str(record.seq)[rand_i:rand_i+min(frag_len,length)], 
                            "type": type, "circular": circular
                            }
        except Exception as err:
            with open(err_file, 'a') as fout:
                fout.write(f"{timestamp()} Error reading chromosome file {filename}: {err}, {len(record)}\n")

def write_plas_frags(input, path, txt_path, err_file, log_file):
    n, frag = input
    try:
        if (not os.path.isfile(f'{path}/{n}.gz')):
            with open(f'{path}/{n}', 'w+') as fout:
                fout.write(f'>{n} {frag["id"]}\n{frag["seq"]}\n')
            with open(f'{txt_path}/{n}', 'w') as fout:
                fout.write(f'>{n} {frag["id"]}\n{frag["seq"]}\n')
            
    except Exception as err:
        with open(err_file, 'a') as fout:
            fout.write(f"{timestamp()} Error writing plasmid file {n}: {err}\n")

def write_chrom_frags(input, path, txt_path, err_file, log_file):
    n, frag = input
    try:
        if (not os.path.isfile(f'{path}/{n}')):
            with open(f'{path}/{n}', 'w+') as fout:
                fout.write(f'>{n} {frag["name"]} {frag["id"]}\n{frag["seq"]}\n')
            with open(f'{txt_path}/{n}', 'w') as fout:
                fout.write(f'>{n} {frag["name"]} {frag["id"]}\n{frag["seq"]}\n')
    except Exception as err:
        with open(err_file, 'a') as fout:
            fout.write(f"{timestamp()} Error writing chromosome file {n}: {err}\n")

def plasmid_worker(
    input, k, frag_path, frag_txt_path, 
    target_path,
    kmer_write_path, kmer_out_path,
    circ_write_path, 
    db_path, inc_out_path, inc_write_path,
    rrna_out_path, rrna_write_path,
    orit_out_path, orit_write_path,
    mob_out_path, rep_out_path, con_out_path, mob_write_path, rep_write_path, con_write_path,
    progress_bar, err_file, log_file):

    n, frag = input

    write_plas_frags(input, frag_path, frag_txt_path, err_file, log_file)
    write_file_id(input, target_path, err_file)
    count_kmers(k, input, frag_txt_path, kmer_write_path, kmer_out_path, err_file)
    circularity(input, circ_write_path, err_file)
    inc_factor(input, db_path, frag_txt_path, inc_out_path, inc_write_path, err_file)
    rrna_search(input, frag_txt_path, rrna_out_path, rrna_write_path, db_path, err_file)
    orit_search(input, frag_txt_path, db_path, orit_out_path, orit_write_path, err_file)
    # special_gene_search(input,frag_txt_path, mob_out_path, rep_out_path, con_out_path, mob_write_path, rep_write_path, con_write_path, db_path, err_file)

    if(n%10 == 0):
        update_progress_bar(progress_bar, 10)

    if(int(log(n, 2)) == log(n, 2)):
        with open(log_file, 'a') as fout:
            fout.write(f"{timestamp()} Completed writing plasmid file {n}\n")

def chrom_worker(
    input, k, frag_path, frag_txt_path,
    target_path,
    kmer_write_path, kmer_out_path,
    circ_write_path,
    db_path, inc_out_path, inc_write_path,
    rrna_out_path, rrna_write_path,
    orit_out_path, orit_write_path,
    mob_out_path, rep_out_path, con_out_path, mob_wchrom_mob_out_pathrite_path, rep_write_path, con_write_path,
    progress_bar, err_file, log_file):

    n, frag = input

    write_chrom_frags(input, frag_path, frag_txt_path, err_file, log_file)
    write_file_id(input, target_path, err_file)
    count_kmers(k, input, frag_txt_path, kmer_write_path, kmer_out_path, err_file)
    circularity(input, circ_write_path, err_file)
    inc_factor(input, db_path, frag_txt_path, inc_out_path, inc_write_path, err_file)
    rrna_search(input, frag_txt_path, rrna_out_path, rrna_write_path, db_path, err_file)
    orit_search(input, frag_txt_path, db_path, orit_out_path, orit_write_path, err_file)
    # special_gene_search(input,frag_txt_path, mob_out_path, rep_out_path, con_out_path, mob_write_path, rep_write_path, con_write_path, db_path, err_file)

    if(n%10 == 0):
        update_progress_bar(progress_bar, 10)

    if(int(log(n, 2)) == log(n, 2)):
        with open(log_file, 'a') as fout:
            fout.write(f"{timestamp()} Completed writing chromosome file {n}\n")

def process():
    
    ## create file structure
    create_file_structure()

    ## read and write plasmid/ chromosome files
    plas_frag_gen = plas_frag_generator(plas_db_path, frag_len, plas_frag_split_path, plas_circ_out_path, coverage, err_file)
    chrom_frag_gen = chrom_frag_generator(chrom_db_path, frag_len, chrom_frag_split_path, chrom_circ_out_path, coverage, err_file)

    chrom_frag_gen_c = filter(lambda x: x["type"] == "chromosome", chrom_frag_gen)
    chrom_frag_gen_p = filter(lambda x: x["type"] == "plasmid", chrom_frag_gen)

    plas_bar, chrom_bar, ex_plas_bar = create_progress_bars()

    with concurrent.futures.ThreadPoolExecutor(max_workers=thread_count) as executor:
        try:
            executor.map(
                plasmid_worker, enumerate(plas_frag_gen), repeat(k), repeat(plas_write_path), repeat(plas_txt_write_path),
                repeat(plas_target_path),
                repeat(plas_7mer_write_path), repeat(plas_7mer_out_path),
                repeat(plas_circ_write_path),
                repeat(db_path), repeat(plas_inc_out_path), repeat(plas_inc_write_path),
                repeat(plas_rrna_out_path), repeat(plas_rrna_write_path),
                repeat(plas_orit_out_path), repeat(plas_orit_write_path),
                repeat(plas_mob_out_path), repeat(plas_rep_out_path), repeat(plas_con_out_path), 
                repeat(plas_mob_write_path), repeat(plas_rep_write_path), repeat(plas_con_write_path),
                repeat(plas_bar),
                repeat(err_file), repeat(log_file))
        except Exception as err:
            with open(err_file, 'a') as fout:
                fout.write(f"{timestamp()} Concurrency Error: {err}\n")

        try:
            executor.map(
                chrom_worker, enumerate(chrom_frag_gen_c), repeat(k), repeat(chrom_write_path), repeat(chrom_txt_write_path), 
                repeat(chrom_target_path),
                repeat(chrom_7mer_write_path), repeat(chrom_7mer_out_path),
                repeat(chrom_circ_write_path),
                repeat(db_path), repeat(chrom_inc_out_path), repeat(chrom_inc_write_path),
                repeat(chrom_rrna_out_path), repeat(chrom_rrna_write_path),
                repeat(chrom_orit_out_path), repeat(chrom_orit_write_path),
                repeat(chrom_mob_out_path), repeat(chrom_rep_out_path), repeat(chrom_con_out_path), 
                repeat(chrom_mob_write_path), repeat(chrom_rep_write_path), repeat(chrom_con_write_path),
                repeat(chrom_bar),
                repeat(err_file), repeat(log_file))
        except Exception as err:
            with open(err_file, 'a') as fout:
                fout.write(f"{timestamp()} Concurrency Error: {err}\n")

        try:
            executor.map(chrom_worker, enumerate(chrom_frag_gen_p), repeat(k), repeat(extra_plasmid_write_path), repeat(extra_plasmid_txt_write_path), 
            repeat(ex_plas_target_path),
            repeat(ex_plas_7mer_write_path), repeat(ex_plas_7mer_out_path),
            repeat(ex_plas_circ_write_path),
            repeat(db_path), repeat(chrom_inc_out_path), repeat(chrom_inc_write_path),
            repeat(ex_plas_rrna_out_path), repeat(ex_plas_rrna_write_path),
            repeat(ex_plas_orit_out_path), repeat(ex_plas_orit_write_path),
            repeat(ex_plas_mob_out_path), repeat(ex_plas_rep_out_path), repeat(ex_plas_con_out_path), 
            repeat(ex_plas_mob_write_path), repeat(ex_plas_rep_write_path), repeat(ex_plas_con_write_path),
            repeat(ex_plas_bar),
            repeat(err_file), repeat(log_file))
        except Exception as err:
            with open(err_file, 'a') as fout:
                fout.write(f"{timestamp()} Concurrency Error: {err}\n")

if __name__ == "__main__":
    process()
