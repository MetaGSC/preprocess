import gzip
import subprocess as sp
import os

from helpers import find_mean, timestamp
from constants import blastn_path

def inc_factor(input, db_path, input_path, out_path, result_path, err_file):
    try:
        n, frag = input
        input_file = f'{input_path}/{n}'
        out_file = f"{out_path}/{n}.inc.blast.out"

        db_file = f"{db_path}/inc-types.fasta"

        cmd = [
            blastn_path,
            '-query', db_file,
            '-subject', input_file,
            '-num_threads', '1',    
            '-perc_identity', '90',
            '-culling_limit', '1',
            '-outfmt', '6 qseqid sstart send sstrand pident qcovs bitscore',
            '-out', out_file,
        ]
        proc = sp.run(
            cmd,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            universal_newlines=True
        )
        bitscore_mean, coverage_mean, length_mean, count = 0, 0, 0, 0
        if(proc.returncode == 0 and not os.stat(out_file).st_size == 0):
            bitscores = []
            coverages = []
            lengths = []
            hit_poses = []
            count = 0

            with out_file.open() as fh:
                for line in fh:
                    cols = line.rstrip().split('\t')

                    start = int(cols[1])
                    end = int(cols[2])
                    length = abs(end - start) + 1

                    strand = '+' if cols[3] == 'plus' else '-',
                    coverage = float(cols[5]) / 100
                    bitscore = float(cols[6])

                    hit_pos = end if strand == '+' else start
                    
                    if (hit_pos in hit_poses):
                        former_hit = hit_poses.index(hit_pos)
                        if (bitscores[former_hit] < bitscore):
                            bitscores[former_hit] = bitscore
                            coverages[former_hit] = coverage
                            lengths[former_hit] = length
                            count+=1
                        else:
                            bitscores.append(bitscore)
                            coverages.append(coverage)
                            lengths.append(length)
                            count+=1
                    
            bitscore_mean = find_mean(bitscores)
            coverage_mean = find_mean(coverages)
            length_mean = find_mean(lengths)
        with gzip.open(f"{result_path}/{n}.gz", mode='wt') as fout:
            fout.write(f"{bitscore_mean}\t{coverage_mean}\t{length_mean}\t{count}\n")
    except Exception as err:
        with open(err_file, 'a') as fout:
            fout.write(f"{timestamp()} Error calculating inc. factor for file:{n} {err}\n")
