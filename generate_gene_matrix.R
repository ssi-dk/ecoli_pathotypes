library(tidyverse)
library(magrittr)
library(Biostrings)
library(data.table)
library(cowplot)

args = commandArgs(trailingOnly=T)
PROJECT_DIR = args[1] # directory containing files and subfolders, must end with a /

seqs <- readDNAStringSet(paste0(PROJECT_DIR, "concatenated_genes.fasta"))

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

kma_res <- read_tsv(paste0(PROJECT_DIR, 'concatenated_genes_kma_results.tsv'))
kma_res %<>% mutate('assembly_barcode'=gsub('.fasta.res', '', fasta_name))
# also filter based on template identity and template coverage
kma_res %<>% filter(Template_Identity>=90, Template_Coverage>=90)

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

assembly_gene_matrix <- bind_rows(gene_vectors, missing_mat)

# save the gene matrix with all genes from the fasta file
write_tsv(assembly_gene_matrix, paste0(PROJECT_DIR, 'assembly_gene_matrix.tsv'))

# generate classifications and append them to the metadata file
# generate a vector of genes supersets (i.e. not variants) that should be checked for matches in the kma results
# note that some capitalizations have been changed/added according to what exists in the concatenated_genes.fasta file
gene_super_sets <- c('chuA', 'fyuA', 'vat', 'yfcV', # UPEC
                     'pap', 'sfa', 'foc', 'afa', 'Afa', 'Dra', 'iutA', 'kpsMII', # EXPEC
                     'eltIAB', 'eltIIAB', 'estap', 'estah', 'estb', # ETEC
                     'aggR', 'aggA', 'aafA', 'agg3A', 'agg4A', 'agg5A', 'C719-09_CS22', # EAEC
                     'ipaH', # EIEC
                     'eae', # AEEC
                     'stx1', 'stx2', # STEC
                     'bfpA', 'eae', # EPEC
                     'iutA','hlyF', 'iss', 'iroN', 'ompT', # APEC
                     'afa','Dra','bma', 'Daa' # DAEC
                     ) %>% unique()
gene_bases <- kma_res$`#Template` %>% gsub('kpsMIII.*', 'not_interested_in_these', .) %>% str_extract(paste0(gene_super_sets, collapse='|')) # will grab the first match if several
kma_res %<>% mutate(gene_bases=gene_bases)
# check the ones where the gene base is NA
na_kma_res <- kma_res[is.na(gene_bases), ]
na_kma_res$`#Template` %>% unique


start_time <- Sys.time()
super_sets_per_sample <- assembly_gene_matrix %>% select(-assembly_barcode) %>% as.matrix %>% apply(1, function(gene.vector){
  gene.vector.hits <- gene.vector[gene.vector > 0] %>% names()
  gene.bases <- gene.vector.hits %>% gsub('kpsMIII.*', 'not_interested_in_these', .) %>% str_extract(paste0(gene_super_sets, collapse='|')) %>% unique()
  return(gene.bases)
})
end_time <- Sys.time()
end_time - start_time

check_kma_results_by_some_index_in_superset_result <- function(x){
  return(kma_res %>% filter(assembly_barcode==assembly_gene_matrix$assembly_barcode[x]) %>% select(-c(Depth, q_value, p_value, fasta_name)))
}


# which samples have no matches
rows_in_assembly_matrix_no_match <- super_sets_per_sample %>% lapply(function(x){length(x)<1}) %>% unlist() %>% which()
assemblies_no_matches <- assembly_gene_matrix$assembly_barcode[rows_in_assembly_matrix_no_match]

# check concordance of NA in the matches fetched with both the whole template vector and applying it over the gene matrix
super_sets_per_sample %>% unlist %>% is.na() %>% sum()
is.na(gene_bases) %>% sum() # difference between these two is caused by the uniqueification during the processing of gene vectors

# looks like it's right, start making the function for classifying the types
assign_pathotype <- function(vector.of.present.genes){ # vector of present genes (not variants) for a given assembly
  # we just enumerate the genes again because some of the conditions are not super easy to 'logic' out without checking for presence/absence
  init.match <- c()
  sfa.foc.matches <- any(c('sfa', 'foc') %in% vector.of.present.genes) %>% sum() # same gene so much be condensed
  afa.dra.matches <- any(c('afa', 'Afa', 'Dra') %in% vector.of.present.genes) %>% sum() # same gene, also one is named Afa instead of afa
  ETEC <- any(c('eltIAB', 'eltIIAB', 'estap', 'estah', 'estb') %chin% vector.of.present.genes)
  EAEC <- ('aggR' %chin% vector.of.present.genes & any(c('aggA', 'aafA', 'agg3A', 'agg4A', 'agg5A') %chin% vector.of.present.genes)) | 'C719-09_CS22' %chin% vector.of.present.genes
  EIEC <- 'ipaH' %chin% vector.of.present.genes
  AEEC <- 'eae' %chin% vector.of.present.genes & !(any(c('stx1', 'stx2') %chin% vector.of.present.genes))
  STEC <- any(c('stx1', 'stx2') %in% vector.of.present.genes)
  ExPECjj <- (sum(c('pap', 'iutA', 'kpsMII') %chin% vector.of.present.genes) + sfa.foc.matches + afa.dra.matches) > 1
  UPEC <- sum(c('chuA', 'fyuA', 'vat', 'yfcV') %chin% vector.of.present.genes) > 2
  EPEC <- 'bfpA' %chin% vector.of.present.genes & 'eae' %chin% vector.of.present.genes
  APEC <- all(c('iutA','hlyF', 'iss', 'iroN', 'ompT') %chin% vector.of.present.genes)
  DAEC <- any(c('afa','Afa', 'Dra','bma', 'Daa') %chin% vector.of.present.genes)
  match.vector <- c('ETEC', 'EAEC', 'EIEC', 'AEEC', 'STEC', 'ExPECjj', 'UPEC', 'EPEC', 'APEC', 'DAEC')
  match.vector %>% lapply(function(x){
    if (get(x)){ # evaluate name of variable
      init.match <<- c(init.match, x)
    }
  })
  #if (length(init.match) < 1){
  #  init.match <- 'unassigned'
  #}
  return(init.match)
}

# generate the pathotype vector and append it to the metadata file
pathotypes_per_assembly <- super_sets_per_sample %>% lapply(assign_pathotype)

pathotype_strings <- pathotypes_per_assembly %>% lapply(function(x){ifelse(length(x) < 1, 'unassigned', paste0(x, collapse=','))}) %>% unlist()
pathotype_by_names <- setNames(pathotype_strings, assembly_gene_matrix$assembly_barcode)
meta_final %<>% mutate(pathotype_from_assembly=unname(pathotype_by_names[`Assembly barcode`]))
print(paste0('writing metadata file with appended pathotypes to ', PROJECT_DIR))
write_tsv(meta_final, paste0(PROJECT_DIR, 'joined_meta_filtered_no_missing_with_quast_results_and_pathotypes.tsv'))

