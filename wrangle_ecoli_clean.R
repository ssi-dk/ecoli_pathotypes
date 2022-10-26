library(tidyverse)
library(magrittr)
library(jsonlite)

# load tsvs and join
args = commandArgs(trailingOnly=T)
PROJECT_DIR = args[1] # directory containing files and subfolders, must end with a /

ecoli_meta <- read_tsv(paste0(PROJECT_DIR, 'meta_with_assemblystats.tsv'), guess_max=100000, col_types = cols(.default='?', 'Collection Time' = 'c'))
meta_cgmlst <- read_tsv(paste0(PROJECT_DIR, 'meta_with_cgmlst.tsv'), guess_max=100000, col_types = cols(.default='?', 'Collection Time' = 'c')) %>% 
  rename('cgMLST_Differences'='Differences', 'cgMLST_ST'='ST')
meta_wgmlst <- read_tsv(paste0(PROJECT_DIR, 'meta_with_wgmlst.tsv'), guess_max=100000, col_types = cols(.default='?', 'Collection Time' = 'c')) %>% 
  rename('wgMLST_Differences'='Differences', 'wgMLST_ST'='ST')
meta_phylotypes <- read_tsv(paste0(PROJECT_DIR, 'meta_with_phylotypes.tsv'), guess_max=100000, col_types = cols(.default='?', 'Collection Time' = 'c')) %>% 
  rename('phylotypeST'='ST')
meta_achtman <- read_tsv(paste0(PROJECT_DIR, 'meta_with_achtman.tsv'), guess_max=100000, col_types = cols(.default='?', 'Collection Time' = 'c'))


core_cols <- intersect(ecoli_meta %>% names, meta_cgmlst %>% names)

joined_meta <- ecoli_meta %>% 
  left_join(meta_cgmlst, by=core_cols) %>% 
  left_join(meta_wgmlst, by=core_cols) %>% 
  left_join(meta_phylotypes, by=core_cols) %>% 
  left_join(meta_achtman, by=core_cols)

## filter rows
# want a source
any_source <- (!is.na(joined_meta$`Source Details`) | !is.na(joined_meta$`Source Niche`)) | !is.na(joined_meta$`Source Details`)
all_source <- (!is.na(joined_meta$`Source Details`) & !is.na(joined_meta$`Source Niche`)) & !is.na(joined_meta$`Source Details`)

joined_meta_filt <- joined_meta %>% filter(!is.na(`Collection Year`) & 
                                             !is.na(Country) &
                                             any_source)

unique_cols <- names(joined_meta_filt)[names(joined_meta_filt) %>% lapply(function(x){length(joined_meta_filt[[x]] %>% unique())==nrow(joined_meta_filt)}) %>% unlist()]
# write metadata tsv
write_tsv(joined_meta_filt, file = paste0(PROJECT_DIR, 'ecoli_joined_meta_filtered.tsv'))


