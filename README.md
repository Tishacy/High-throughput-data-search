# High-throughput-data-search

## Introduction

This repo is doing data preparation for SUNXING's graduate design.  All source code are contained in `src` folder. 

- `dataset`: a directory to store the original dataset as well as the code generated files.
- `img`: a directory to store images used in `Searching Data.ipynb`, in order to describe the workflow.
- `Searching Data.ipynb`: a jupyter notebook to describe the workflow.  
- `processData.py`: process data in dataset. This would generate two `txt` files in each own dataset: 
  - `OTUs_need_search.txt`: store the strain names and their OTUs that need search.
  - `fasta_need_search.txt`: store the strain names and their gene sequences that need search.
- `searchSeqID.py`: search the gene sequences stored in `fasta_need_search.txt` on [Ribosomal Database Project (RDP)](http://rdp.cme.msu.edu/seqmatch/seqmatch_intro.jsp), get sequence IDs that have match scores greater than 0.95, and store the sequence IDs in `seq_id_searched.txt` in dataset.
- `extractData.py`: extract sequence IDs from `seq_id_searched.txt`, search the genus of each `sequence ID`, and store results in `best_match_strain.xlsx`.

## Usage

1. Download the repository and unzip it.

2. Get into `src` by typing `cd src` in terminal or just open the `src` folder.

3. Preprocess data

   ```python
   python processData.py
   ```

   Type the dataset name (names of folders in `dataset`), and return it.

4. Search the sequence IDs

   ```python
   python searchSeqID.py
   ```

   Type the dataset name and return it.

5. Extract data to get genus names of sequence IDs

   ```python
   python extractData.py
   ```

   Type the dataset name and return it.

## Requirements

Using `pip` to install all the dependencies required.

```bash
pip install -r requirements.txt
```
