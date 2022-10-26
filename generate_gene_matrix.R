library(tidyverse)
library(magrittr)
library(Biostrings)

args = commandArgs(trailingOnly=T)
PROJECT_DIR = args[1] # directory containing files and subfolders, must end with a /

seqs <- readDNAStringSet(paste0(PROJECT_DIR, "Ecoli_VF_v2.txt"))

# some of the sequences in the file are duplicated, we will check if these are the same sequences
duplicated_names <- names(seqs)[names(seqs) %>% duplicated]
check_if_dups_equal <- duplicated_names %>% lapply(function(x){
  inds <- which(names(seqs)==x)
  seq1 <- seqs[inds[1]]
  seq2 <- seqs[inds[2]]
  return(as.character(seq1) == as.character(seq2))
}) %>% unlist()

print(paste0('are all duplicated sequence ids also duplicate sequence strings? ', all(check_if_dups_equal)))

unique_gene_names <- unique(names(seqs))

meta_final <- read_tsv(paste0(PROJECT_DIR, 'joined_meta_filtered_no_missing_with_quast_results.tsv'), guess_max = 100000)

kma_res <- read_tsv('/Volumes/data/MPV/projects/ecoli_pathotypes/kma_results_concatenated.tsv')
kma_res %<>% mutate('assembly_barcode'=gsub('.fasta.res', '', fasta_name))


# split the kma results by assembly so we can create gene vectors etc.
kma_res_by_assembly <- kma_res %$% split(., assembly_barcode)

# we will get a binary gene vector for the genes for each assembly based on matches in the kma result submatrices
get_gene_vector <- function(kma.result, ecoli.seq.names){
  binary.vec <- as.integer(ecoli.seq.names %in% kma.result$'#Template')
  names(binary.vec) <- ecoli.seq.names
  return(binary.vec)
}

gene_vectors <- kma_res_by_assembly %>% lapply(get_gene_vector, unique_gene_names) %>% bind_rows(.id='assembly_barcode')

# some samples are missing from the kma results because no matches, we add zero vectors for these
missing_assemblies <- setdiff(meta_final$`Assembly barcode`, kma_res$assembly_barcode)
missing_mat <- matrix(0, length(missing_assemblies), length(unique_gene_names))
colnames(missing_mat) <- unique_gene_names
missing_mat <- bind_cols(assembly_barcode=missing_assemblies, missing_mat)

final_df <- bind_rows(gene_vectors, missing_mat)

write_tsv(final_df, paste0(PROJECT_DIR, 'assembly_gene_matrix.tsv'))
