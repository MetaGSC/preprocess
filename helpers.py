import os
from constants import *

def find_mean(data):
  if len(data) == 0:
    return 0
  return float(sum(data))/len(data)

def timestamp():
    import datetime
    return datetime.datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")

def createDirIfNeeded(path):
    if not os.path.exists(path):
        os.makedirs(path)

def createFileStructure():
        ## create dirs if they don't exist
    # fragments
    createDirIfNeeded(plas_write_path)
    createDirIfNeeded(chrom_write_path)
    createDirIfNeeded(extra_plasmid_write_path)

    # frag txt
    createDirIfNeeded(plas_txt_write_path)
    createDirIfNeeded(chrom_txt_write_path)
    createDirIfNeeded(extra_plasmid_txt_write_path)

    # kmer files
    createDirIfNeeded(plas_7mer_write_path)
    createDirIfNeeded(chrom_7mer_write_path)
    createDirIfNeeded(ex_plas_7mer_write_path)

    createDirIfNeeded(plas_7mer_out_path)
    createDirIfNeeded(chrom_7mer_out_path)
    createDirIfNeeded(ex_plas_7mer_out_path)

    # circularity files
    createDirIfNeeded(plas_circ_write_path)
    createDirIfNeeded(chrom_circ_write_path)
    createDirIfNeeded(ex_plas_circ_write_path)

    createDirIfNeeded(plas_circ_out_path)
    createDirIfNeeded(chrom_circ_out_path)
    createDirIfNeeded(ex_plas_circ_out_path)

    createDirIfNeeded(plas_frag_split_path)
    createDirIfNeeded(chrom_frag_split_path)
    createDirIfNeeded(ex_plas_frag_split_path)

    # inc factor files
    createDirIfNeeded(plas_inc_out_path)
    createDirIfNeeded(chrom_inc_out_path)
    createDirIfNeeded(ex_plas_inc_out_path)

    createDirIfNeeded(plas_inc_write_path)
    createDirIfNeeded(chrom_inc_write_path)
    createDirIfNeeded(ex_plas_inc_write_path)

    # rrna files
    createDirIfNeeded(plas_rrna_out_path)
    createDirIfNeeded(chrom_rrna_out_path)
    createDirIfNeeded(ex_plas_rrna_out_path)

    createDirIfNeeded(plas_rrna_write_path)
    createDirIfNeeded(chrom_rrna_write_path)
    createDirIfNeeded(ex_plas_rrna_write_path)

    # orit files
    createDirIfNeeded(plas_orit_out_path)
    createDirIfNeeded(chrom_orit_out_path)
    createDirIfNeeded(ex_plas_orit_out_path)
    
    createDirIfNeeded(plas_orit_write_path)
    createDirIfNeeded(chrom_orit_write_path)
    createDirIfNeeded(ex_plas_orit_write_path)

    # mobilization files
    createDirIfNeeded(plas_mob_out_path)
    createDirIfNeeded(chrom_mob_out_path)
    createDirIfNeeded(ex_plas_mob_out_path)

    createDirIfNeeded(plas_mob_write_path)
    createDirIfNeeded(chrom_mob_write_path)
    createDirIfNeeded(ex_plas_mob_write_path)

    # replication files
    createDirIfNeeded(plas_rep_out_path)
    createDirIfNeeded(chrom_rep_out_path)
    createDirIfNeeded(ex_plas_rep_out_path)

    createDirIfNeeded(plas_rep_write_path)
    createDirIfNeeded(chrom_rep_write_path)
    createDirIfNeeded(ex_plas_rep_write_path)

    # conjugation files
    createDirIfNeeded(plas_con_out_path)
    createDirIfNeeded(chrom_con_out_path)
    createDirIfNeeded(ex_plas_con_out_path)

    createDirIfNeeded(plas_con_write_path)
    createDirIfNeeded(chrom_con_write_path)
    createDirIfNeeded(ex_plas_con_write_path)
