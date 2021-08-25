import subprocess as sp
import gzip
import csv

from helpers import find_mean, timestamp

def orit_search(input, input_path, db_path, out_path, write_path, err_file):
    try:
        n, frag = input
        input_file = f'{input_path}/{n}'
        db_file = f'{db_path}/orit'
        out_file = f'{out_path}/{n}.orit.blast.out'
        write_file = f'{write_path}/{n}'

        cmd = [
            'blastn',
            '-query', str(input_file),
            '-db', str(db_file),
            '-num_threads', '1',
            '-culling_limit', '1',
            '-perc_identity', '90',
            '-evalue', '1E-5',
            '-outfmt', '6 sseqid qstart qend sstart send slen length nident',
            '-out', str(out_file)
        ]

        proc = sp.run(
            cmd,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            universal_newlines=True
        )
        count = 0
        lengths = []
        coverages = []
        identities = []
        
        if proc.returncode == 0:
            with open(out_file, 'r') as fh:
                for line in fh:
                    line = line.rstrip()
                    cols = line.split('\t')
                    
                    lengths.append(int(cols[5]))
                    coverages.append(float(cols[6]) / int(cols[5]))
                    identities.append(float(cols[7]) / float(cols[6]))
                    count+=1
                    print(n, count) # remove later
                    
        length_mean = find_mean(lengths)
        coverage_mean = find_mean(coverages)
        identity_mean = find_mean(identities)

        with gzip.open(f"{write_file}.gz", mode='wt') as fout:
            fout.write(f"{identity_mean}\t{coverage_mean}\t{length_mean}\t{count}\n")

    except Exception as err:
        with open(err_file, 'a') as fout:
            fout.write(f"{timestamp()} Error calculating oriT factor for file:{n} {err}\n")
