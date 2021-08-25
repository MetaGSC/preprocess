import subprocess as sp
from helpers import timestamp

def special_gene_search(
    input, input_path, mob_out_path, rep_out_path, con_out_path, mob_write_path, 
    rep_write_path, con_write_path, db_path, err_file
    ):
    try:
        n, frag = input

        input_file = f'{input_path}/{n}'

        mob_db_file = f"{db_path}/mobilization"
        rep_db_file = f"{db_path}/replication"
        con_db_file = f"{db_path}/conjugation"

        mob_out_file = f"{mob_out_path}/{n}.mob.hmm.out"
        rep_out_file = f"{rep_out_path}/{n}.rep.hmm.out"
        con_out_file = f"{con_out_path}/{n}.con.hmm.out"

        mob_write_file = f"{mob_write_path}/{n}"
        rep_write_file = f"{rep_write_path}/{n}"
        con_write_file = f"{con_write_path}/{n}"

        mob_cmd = [
            'hmmsearch',
            '--noali',
            '--cpu', '1',
            '-E', '1E-10',
            '--tblout', str(mob_out_file),
            str(mob_db_file),
            str(input_file)
        ] 
        proc = sp.run(
            mob_cmd,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            universal_newlines=True
        )
    
    except Exception as err:
        with open(err_file, 'a') as fout:
            fout.write(f"{timestamp()} Error calculating rrna factor for file:{n} {err}\n")
