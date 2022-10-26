# get the final tsv
library(tidyverse)
library(magrittr)

args = commandArgs(trailingOnly=T)
PROJECT_DIR = args[1] # directory containing files and subfolders, must end with a /

# file_sizes file generated in bash with the script generate_file_size_table.sh
# note that before filtering on unavailable fastas you should make sure that you have re-run
# samples that failed in the case that they were unavailable due to enterobase instability
# rather than other reasons

fasta_sizes <- read.table(paste0(PROJECT_DIR, '/fasta_file_sizes.txt'))
range(fasta_sizes$V1)
hist(fasta_sizes$V1)

joined_meta_filt <- read_tsv(paste0(PROJECT_DIR, 'ecoli_joined_meta_filtered.tsv'), guess_max=100000)

sum(fasta_sizes$V2 %>% gsub('.fasta', '', .) %in% joined_meta_filt$`Assembly barcode`)

# filter away the samples whose fastas were not available for download
meta_without_missing_fastas <- joined_meta_filt %>% filter(`Assembly barcode` %in% gsub('.fasta', '', fasta_sizes$V2))


# load concatenated quast output tsv and compare with metadata from enterobase
# quast files were generated with the command python basic_fasta_stats.py -i PROJECT_DIR/fasta_files -o PROJECT_DIR/quast_outputs
# and then concatenated using the collect_quast_tsvs.py script
quast_tsv <- read_tsv(paste0(PROJECT_DIR, 'concatenated_quast.tsv')) %>% rename('N50_quast'='N50') # orig tsv contains a N50 col

meta_without_missing_fastas_with_quast <- meta_without_missing_fastas %>% left_join(quast_tsv, by=c('Assembly barcode'='Assembly'))

#all.equal(meta_without_missing_fastas_with_quast$`Total length`, meta_without_missing_fastas_with_quast$`Total length (>= 0 bp)`)

# check a few values, these will not be identical as we don't know how the original statistics
# were generated, but it will give us something to briefly look at
(meta_without_missing_fastas_with_quast$N50/meta_without_missing_fastas_with_quast$N50_quast) %>% hist()
(meta_without_missing_fastas_with_quast$N50/meta_without_missing_fastas_with_quast$N50_quast) %>% summary()

(meta_without_missing_fastas_with_quast$Length/meta_without_missing_fastas_with_quast$`Total length (>= 0 bp)`) %>% hist()
(meta_without_missing_fastas_with_quast$Length/meta_without_missing_fastas_with_quast$`Total length (>= 0 bp)`) %>% summary()

(meta_without_missing_fastas_with_quast$Length/meta_without_missing_fastas_with_quast$`Total length (>= 0 bp)`) %>% hist()
(meta_without_missing_fastas_with_quast$Length/meta_without_missing_fastas_with_quast$`Total length (>= 0 bp)`) %>% summary()

write_tsv(meta_without_missing_fastas_with_quast, paste0(PROJECT_DIR, 'joined_meta_filtered_no_missing_with_quast_results.tsv'))
