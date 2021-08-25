import gzip
import subprocess as sp

def count_kmers(k, input, input_path, write_path, out_path, err_file):
    try:
        n, frag = input
        input_file = f'{input_path}/{n}'
        kmer_out_file = f'{out_path}/{n}'
        write_file = f"{write_path}/{n}"

        # with gzip.open(f'{input_path}.gz', 'rb') as fin:
        #     with open(temp_in, 'wb') as fout:
        #         fout.write(fin.read())
        args = [
            'seq2vec',
            '-f', str(input_file),
            '-o', str(kmer_out_file),
            '-t', '8',
            '-k', str(k),
        ]    
        proc = sp.run(args,  stdout=sp.PIPE, stderr=sp.PIPE)
        with gzip.open(f'{write_file}.gz', 'wb') as fout:
            with open(kmer_out_file, 'rb') as fin:
                fout.write(fin.read())
    except Exception as err:
        print("Error using seq2vec for file:"+input_file+" "+str(err))
        with open(err_file, 'a') as fout:
            fout.write("Error using seq2vec for file "+input_file+" "+str(err)+"\n")