# High-throughput-data-search

This repo is doing data preparation for SUNXING's graduate design.  

- `dataset`: a directory to store the original dataset as well as the code generated files.
- `images`: a directory to store images used in `Searching Data.ipynb`, in order to describe the workflow.
- `Searching Data.ipynb`: a jupyter notebook to describe the whole workflow.  
- `processData.py`: processment of data in dataset. This would generate two json files in each own dataset:
  - `OTUs_need_search.json`: used to store the strain names and their OTUs that need search.
  - `fasta_need_search.json`: used to store the strain names and their gene sequences that need search.
- `searchSeqID.py`: searching the gene sequences stored in `fasta_need_search.json` on [rdp](http://rdp.cme.msu.edu/seqmatch/seqmatch_intro.jsp), getting sequence IDs that have match scores greater than 0.95, and store the sequence IDs in `seq_id_searched.txt` in dataset.

