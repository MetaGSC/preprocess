# README

## Prerequisites

Following tools should be installed.

-   **Seq2Vec** - used for kmer count
    
-   **Infernal** - used for calculate the rRNA availability
    
-   **Mummer** - used to calculate circularity
    
-   **Blast+** - used for OriT sequence availability and Incompatibility sequence availability

-   **Biopython** - used to read fasta files

-   **tqdm** - used to display runtime progress
 
### Infernal (rRNA) - [github link](https://github.com/EddyRivasLab/infernal)
 
```
wget eddylab.org/infernal/infernal-1.1.4.tar.gz
```
```
tar xf infernal-1.1.4.tar.gz
```
```    
cd infernal-1.1.4
```
```
./configure
```
```
make
```
```
make install
```
```
cd easel
```
```
make install
```
```
cmscan
```
```
chmod 777  '/infernal-1.1.4/src/cmscan.c'
```

### MUMMER (Circularity) - [documentation link](http://mummer.sourceforge.net/manual/#installation)

```
wget https://downloads.sourceforge.net/project/mummer/mummer/3.23/MUMmer3.23.tar.gz
```
```
tar -xvf /MUMmer3.23.tar.gz
```
```   
cd /MUMmer3.23/
```
```
sudo apt-get install csh
```
```
make check
```
``` 
make install
```
```
./nucmer -h
```
 
### BLAST+ (OriT and Incomp) - [documentation link](https://www.ncbi.nlm.nih.gov/books/NBK569861/)
 
```
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz
```
```
tar zxvpf ncbi-blast-2.7.1+-x64-linux.tar.gz
```

### BIOPYTHON
```
conda install -c conda-forge biopython
```

### TQDM
```
conda install -c conda-forge tqdm
```

> [OPTIONAL] Add the paths to the executables to your path for tools `seq2vec`, `blastn` in BLAST+ `cmscan` in Infernal and `nucmer` in MUMMER.

## Usage

1. Provide the correct paths to the chromosome and plasmid databases in the constants file `plas_db_path` `chrom_db_path`.
```
plas_db_path = "../databases/DNA-ML_FYP_2021/plasmid_refs"
chrom_db_path = '../databases/DNA-ML_FYP_2021/bacterial_references'
```

2. Provide the `db_path` - the biomarker db path in the constants file.

```
db_path = "../../references/biomarkerdbs"
```

3. If you haven't added the paths to the executables to your path, add them to the constants file.
``` python
cmscan_path = "/path/to/executable"
blastn_path = "/path/to/executable"
seq2vec_path = "/path/to/executable"
nucmer_path = "/path/to/executable"
```

3. Run the file `process.py`
```
python process.py
```