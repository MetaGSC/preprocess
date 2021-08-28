import subprocess as sp
import os
import csv
from Bio import SeqIO

from helpers import timestamp, create_rrna_dirs
from constants import *
from progress_bar import create_progress_bars, update_progress_bar

def rrna_search(frag_path, db_path, out_path, write_path, progress_bar):
    try:
        for filename in os.listdir(frag_path):
            name = filename.split(".")[0]
            input_file = f'{frag_path}/{filename}'
            out_file = f"{out_path}/_temp.rrna.cmscan.tsv"
            write_file = f'{write_path}/{name}'
            db_file = f'{db_path}/rRNA'

            cmd = [
                cmscan_path,
                '--noali',
                '--cut_tc',
                '--cpu', str(thread_count),
                '--tblout', str(out_file),
                db_file,
                input_file,
            ]

            proc = sp.run(
                cmd,
                stdout=sp.PIPE,
                stderr=sp.PIPE,
                universal_newlines=True
            )

            matches = {}
            if proc.returncode == 0:
                tsvfile = open(out_file)
                fh = csv.reader(tsvfile, delimiter="\t")
                for line in fh:
                    if(line[0][0] != '#'):
                        cols = line[0].strip().split()
                        n = int(cols[2])
                        length = abs(int(cols[8])- int(cols[7]))
                        bitscore = float(cols[14])

                        if n in matches.keys():
                            matches[n]['length'] = (matches[n]['length']*matches[n]['count'] + length) / (matches[n]['count'] + 1)
                            matches[n]['bitscore'] = (matches[n]['bitscore']*matches[n]['count'] + bitscore) / (matches[n]['count'] + 1)
                            matches[n]['count'] += 1
                        else:
                            matches[n] = {
                                'length': length,
                                'bitscore': bitscore,
                                'count': 1
                            }

            match_array = []
            for record in SeqIO.parse(f"{frag_path}/{filename}", 'fasta'):
                n = int(record.id)
                if n in matches.keys():
                    match = matches[n]
                    match_array.append(f"{n}\t{match['length']}\t{match['bitscore']}\t{match['count']}\n")
                else:
                    match_array.append(f"{n}\t0\t0\t0\n")
                
            with open(f"{write_file}", 'w+') as fout:
                fout.writelines(match_array)
            update_progress_bar(progress_bar, batch_size)

    except Exception as err:
        with open(err_file, 'a') as fout:
            fout.write(f"{timestamp()} Error calculating rrna factor for file:{n} {err}\n")

if __name__ == "__main__":
    k = 7
    plas_bar, chrom_bar, ex_plas_bar = create_progress_bars()
    create_rrna_dirs()
    rrna_search(plas_write_path, db_path, plas_rrna_out_path, plas_rrna_write_path, plas_bar)
    rrna_search(chrom_write_path, db_path, chrom_rrna_out_path, chrom_rrna_write_path, chrom_bar)
    rrna_search(extra_plasmid_write_path, db_path, ex_plas_rrna_out_path, ex_plas_rrna_write_path, ex_plas_bar)
