import os

from constants import *

def find_mean(data):
  if len(data) == 0:
    return 0
  return float(sum(data))/len(data)

def timestamp():
    import datetime
    return datetime.datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")

def create_dir_if_needed(path):
    if not os.path.exists(path):
        os.makedirs(path)

def create_file_structure():

    ## create dirs if they don't exist
    # fragments
    create_dir_if_needed(plas_write_path)
    create_dir_if_needed(chrom_write_path)
    create_dir_if_needed(extra_plasmid_write_path)

    # frag txt
    create_dir_if_needed(plas_txt_write_path)
    create_dir_if_needed(chrom_txt_write_path)
    create_dir_if_needed(extra_plasmid_txt_write_path)

    # kmer files
    create_dir_if_needed(plas_7mer_write_path)
    create_dir_if_needed(chrom_7mer_write_path)
    create_dir_if_needed(ex_plas_7mer_write_path)

    create_dir_if_needed(plas_7mer_out_path)
    create_dir_if_needed(chrom_7mer_out_path)
    create_dir_if_needed(ex_plas_7mer_out_path)

    # circularity files
    create_dir_if_needed(plas_circ_write_path)
    create_dir_if_needed(chrom_circ_write_path)
    create_dir_if_needed(ex_plas_circ_write_path)

    create_dir_if_needed(plas_circ_out_path)
    create_dir_if_needed(chrom_circ_out_path)
    create_dir_if_needed(ex_plas_circ_out_path)

    create_dir_if_needed(plas_frag_split_path)
    create_dir_if_needed(chrom_frag_split_path)
    create_dir_if_needed(ex_plas_frag_split_path)

    # inc factor files
    create_dir_if_needed(plas_inc_out_path)
    create_dir_if_needed(chrom_inc_out_path)
    create_dir_if_needed(ex_plas_inc_out_path)

    create_dir_if_needed(plas_inc_write_path)
    create_dir_if_needed(chrom_inc_write_path)
    create_dir_if_needed(ex_plas_inc_write_path)

    # rrna files
    create_dir_if_needed(plas_rrna_out_path)
    create_dir_if_needed(chrom_rrna_out_path)
    create_dir_if_needed(ex_plas_rrna_out_path)

    create_dir_if_needed(plas_rrna_write_path)
    create_dir_if_needed(chrom_rrna_write_path)
    create_dir_if_needed(ex_plas_rrna_write_path)

    # orit files
    create_dir_if_needed(plas_orit_out_path)
    create_dir_if_needed(chrom_orit_out_path)
    create_dir_if_needed(ex_plas_orit_out_path)
    
    create_dir_if_needed(plas_orit_write_path)
    create_dir_if_needed(chrom_orit_write_path)
    create_dir_if_needed(ex_plas_orit_write_path)

    # mobilization files
    create_dir_if_needed(plas_mob_out_path)
    create_dir_if_needed(chrom_mob_out_path)
    create_dir_if_needed(ex_plas_mob_out_path)

    create_dir_if_needed(plas_mob_write_path)
    create_dir_if_needed(chrom_mob_write_path)
    create_dir_if_needed(ex_plas_mob_write_path)

    # replication files
    create_dir_if_needed(plas_rep_out_path)
    create_dir_if_needed(chrom_rep_out_path)
    create_dir_if_needed(ex_plas_rep_out_path)

    create_dir_if_needed(plas_rep_write_path)
    create_dir_if_needed(chrom_rep_write_path)
    create_dir_if_needed(ex_plas_rep_write_path)

    # conjugation files
    create_dir_if_needed(plas_con_out_path)
    create_dir_if_needed(chrom_con_out_path)
    create_dir_if_needed(ex_plas_con_out_path)

    create_dir_if_needed(plas_con_write_path)
    create_dir_if_needed(chrom_con_write_path)
    create_dir_if_needed(ex_plas_con_write_path)
