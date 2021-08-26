import gzip
import subprocess as sp

from helpers import timestamp
from constants import nucmer_path

def check_circularity(record, split_path, out_path):
    id = record.id
    seq_mid = int(len(record.seq)/2)
    seq_a = str(record.seq[:seq_mid])
    seq_b = str(record.seq[seq_mid:])

    seq_a_path = f'{split_path}/{id}_a'
    seq_b_path = f'{split_path}/{id}_b'
    with open(seq_a_path, mode='w+') as fout:
        fout.write(f'>{id}_a\n')
        fout.write(seq_a + '\n ')
    
    with open(seq_b_path, mode='w+') as fout:
        fout.write(f'>{id}_b\n')
        fout.write(seq_b + '\n ')

    out_file = f"{out_path}/{id}"

    cmd = [
    nucmer_path,
    '-f',  # only forward strand
    '-l', '40',  # increase min match length to 40 bp
    '-p', out_file,
    seq_a_path,
    seq_b_path
    ]
    proc = sp.run(
        cmd,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True,
    )
    circular = 0
    has_match = False
    if(proc.returncode == 0):
        with open(f"{out_file}.delta", ) as fout:
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
                        if(alignment_a == alignment_b
                                and alignment_a > 100
                                and (mismatches / alignment_a) < 0.05
                                and end_b == len(seq_b)
                                and start_a == 1):
                            circular = 1
    return circular

def circularity(input, write_path, err_file):
    try:
        n, frag = input
        with gzip.open(f"{write_path}/{n}.gz", mode='wt') as fout:
            fout.write(str(frag["circular"]))
    except Exception as err:
        with open(err_file, 'a') as fout:
            fout.write(f"{timestamp()} Error using nucmer for file {n} {err}\n")
