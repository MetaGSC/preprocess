import subprocess as sp
import os
from Bio import SeqIO

from helpers import timestamp, create_gccontent_dirs
from progress_bar import *
from constants import *

def count_bases(seq):
    sq_ln = len(seq)

    dir = {}
    dir['a'] = (seq.count("a") + seq.count("A"))/sq_ln
    dir['c'] = (seq.count("c") + seq.count("C"))/sq_ln
    dir['g'] = (seq.count("g") + seq.count("G"))/sq_ln
    dir['t'] = (seq.count("t") + seq.count("T"))/sq_ln
    dir['gc'] = dir["g"] + dir["c"]

    return dir


def gc_content(frag_path, write_path, pb_desc):
    try:
        progress_bar = create_progress_bar(pb_desc)
        for filename in os.listdir(frag_path):
            name = filename.split(".")[0]
            write_file = f'{write_path}/{name}'
            match_array = []
            for record in SeqIO.parse(f"{frag_path}/{filename}", 'fasta'):
                n = int(record.id)
                percentages = count_bases(record.seq)
                match_array.append(f"{n}\t{percentages['a']}\t{percentages['c']}\t{percentages['g']}\t{percentages['t']}\t{percentages['gc']}\n")
                

            with open(f"{write_file}", 'w+') as fout:
                fout.writelines(match_array)
                
            update_progress_bar(progress_bar, batch_size)
    
        close_progress_bar(progress_bar)
       
    except Exception as err:
        with open(err_file, 'a') as fout:
            fout.write(f"{timestamp()} Error calculating gc content {err}\n")

if __name__ == "__main__":
    create_gccontent_dirs()
    gc_content(plas_write_path, plas_gccontent_write_path, plas_bar_desc)
    gc_content(chrom_write_path, chrom_gccontent_write_path, chrom_bar_desc)
    gc_content(extra_plasmid_write_path, ex_plas_gccontent_write_path, ex_plas_bar_desc)
