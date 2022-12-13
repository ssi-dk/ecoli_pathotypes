# A bunch of stuff for downloading and checking stuff on a bunch of ecoli genomes from enterobase
The scripts work together based on the assumption that stuff is being done in a particular project directory. The R scripts take this directory as a command line argument. It's also assumed at the [slurm](https://slurm.schedmd.com) scheduler is used for sending jobs on the cluster.
# Dependencies
We used [quast](http://bioinf.spbau.ru/quast) for checking basic statistics on the assemblies, [kma](https://bitbucket.org/genomicepidemiology/kma/src/master/) for checking matches to specific genes.
Python dependencies: pandas  
R dependencies: tidyverse, Biostrings, jsonlite  
Additionally, these scripts assume we're using the slurm job manager.  
# Downloads
First, metadata containing assemblystats, cgmlst, wgmlst phylotypes and achtman mlst was downloaded off enterobase manually, in our case such that the the release date was at most 22-09-27, and named as following:  
'meta_with_assemblystats.tsv'  
'meta_with_cgmlst.tsv'  
meta_with_wgmlst.tsv  
meta_with_phylotypes.tsv  
wrangle_ecoli_clean.R was used for combining these metadata files and filtering them to get a list of desired samples (e.g. must have collection year, source).  
download_links.py was used for downloading assembly links.  
generate_download.R is a script that calls download_fastas.py on a subset of assembly links and does this using slurm.  
The provided genes used for classifying pathotypes were indexed with kma and then send_kma_jobs.py is used which calls run_kma_with_genes_as_db.py and uses kma to compute basic qc statistics on the assemblies:  
python send_kma_jobs.py -i fasta_files -k /path/to/kma_binary -db /path/to/indexed_genes -o concatenated_genes_kma_results -chunk 10000  
The resulting tsvs produced by kma for each assembly are then concatenated using collect_kma_results.py and the QC results are appended to the metadata using get_final_meta.R  
Then finally generate_gene_matrix.R is used to both generate a binary gene matrix for the samples but also call pathotypes depending on the presence of certain genes in the assemblies