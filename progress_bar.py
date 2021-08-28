from time import sleep
from tqdm import tqdm

def create_progress_bars():
    plas_bar = tqdm(desc = "plasmid      ")
    chrom_bar = tqdm(desc = "chromosome   ")
    ex_plas_bar = tqdm(desc = "extra plasmid")

    return plas_bar, chrom_bar, ex_plas_bar

def update_progress_bar(pbar,i):
    pbar.update(i)
