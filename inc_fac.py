import subprocess as sp
import os
from Bio import SeqIO

from helpers import timestamp, create_inc_factor_dirs
from progress_bar import *
from constants import *

def inc_factor(frag_path, db_path, out_path, write_path, pb_desc):
    try:
        progress_bar = create_progress_bar(pb_desc)
        for filename in os.listdir(frag_path):
            name = filename.split(".")[0]
            input_file = f'{frag_path}/{name}.fasta'
            out_file = f"{out_path}/_temp.inc.blast.out"
            write_file = f'{write_path}/{name}'
            db_file = f"{db_path}/inc-types.fasta"

            cmd = [
                blastn_path,
                '-query', db_file,
                '-subject', input_file,
                '-num_threads', str(thread_count),    
                '-perc_identity', '90',
                '-culling_limit', '1',
                '-outfmt', '6 qseqid sseqid pident length bitscore',
                '-out', out_file,
            ]
            proc = sp.run(
                cmd,
                stdout=sp.PIPE,
                stderr=sp.PIPE,
                universal_newlines=True
            )

            matches = {}
            if proc.returncode == 0:
                with open(out_file, 'r') as fh:
                    for line in fh:
                        line = line.rstrip()
                        cols = line.split('\t')
                        n = int(cols[1])
                        identity = float(cols[2])
                        length = int(cols[3])
                        bitscore = float(cols[4])

                        if n in matches.keys():
                            matches[n]['identity'] = (matches[n]['identity']*matches[n]['count'] + identity) / (matches[n]['count'] + 1)
                            matches[n]['length'] = (matches[n]['length']*matches[n]['count'] + length) / (matches[n]['count'] + 1)
                            matches[n]['bitscore'] = (matches[n]['bitscore']*matches[n]['count'] + bitscore) / (matches[n]['count'] + 1)
                            matches[n]['count'] += 1
                        else:
                            matches[n] = {
                                'identity': identity,
                                'length': length,
                                'bitscore': bitscore,
                                'count': 1
                            }

            match_array = []
            for record in SeqIO.parse(f"{frag_path}/{filename}", 'fasta'):
                n = int(record.id)
                if n in matches.keys():
                    match = matches[n]
                    match_array.append(f"{n}\t{match['identity']}\t{match['length']}\t{match['bitscore']}\t{match['count']}\n")
                else:
                    match_array.append(f"{n}\t0\t0\t0\t0\n")
                

            with open(f"{write_file}", 'w+') as fout:
                fout.writelines(match_array)
                
            update_progress_bar(progress_bar, batch_size)
    
        close_progress_bar(progress_bar)
       
    except Exception as err:
        with open(err_file, 'a') as fout:
            fout.write(f"{timestamp()} Error calculating inc. factor for file:{name} {err}\n")

if __name__ == "__main__":
    create_inc_factor_dirs()
    inc_factor(plas_write_path, db_path, plas_inc_out_path, plas_inc_write_path, plas_bar_desc)
    inc_factor(chrom_write_path, db_path, chrom_inc_out_path, chrom_inc_write_path, chrom_bar_desc)
    inc_factor(extra_plasmid_write_path, db_path, ex_plas_inc_out_path, ex_plas_inc_write_path, ex_plas_bar_desc)
