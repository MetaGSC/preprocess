from time import sleep
from tqdm import tqdm

def create_progress_bar(desc):
    return tqdm(desc = desc)

def update_progress_bar(pbar,i):
    pbar.update(i)

def close_progress_bar(pbar):
    pbar.close()