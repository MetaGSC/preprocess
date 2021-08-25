import subprocess as sp
import gzip
import csv

from helpers import find_mean

def rRNA_search(input, input_path, out_path, write_path, db_path, err_file):
    try:
        n, frag = input
        input_file = f'{input_path}/{n}'
        out_file = f"{out_path}/{n}.rrna.cmscan.tsv"
        db_file = f"{db_path}/rRNA"

        cmd = [
            'cmscan',
            '--noali',
            '--cut_tc',
            '--cpu', '1',
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
        tsvfile = open(out_file)
        fh = csv.reader(tsvfile, delimiter="\t")
        lengths = []
        bitscores = []
        count = 0

        for line in fh:
            if(line[0][0] != '#'):
                cols = line[0].strip().split()
                lengths.append(abs(int(cols[8])- int(cols[7])))
                bitscores.append(float(cols[14]))
                count+=1
        tsvfile.close()
        
        bitscore_mean = find_mean(bitscores)
        length_mean = find_mean(lengths)
        with gzip.open(f"{write_path}/{n}.gz", mode='wt') as fout:
            fout.write(f"{bitscore_mean}\t{length_mean}\t{count}\n")
    except Exception as err:
        print("Error calculating rrna factor for file:"+n+" "+str(err))
        with open(err_file, 'a') as fout:
            fout.write("Error calculating rrna factor for file:"+n+" "+str(err)+"\n")