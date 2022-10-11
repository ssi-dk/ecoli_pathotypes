library(tidyverse)
library(magrittr)
library(jsonlite)

# temporary hack for running both mounted and on srver
if(file.exists("/srv/data/")) {
  BASE_PATH <- "/srv/data/MPV/"
} else {
  BASE_PATH <- "/Volumes/data/MPV/"
}
PROJECT_DIR <- paste0(BASE_PATH, 'projects/ecoli_pathotypes/')
joined_meta_filt <- read_tsv(paste0(BASE_PATH, 'LEBC/ecoli_joined_meta_filtered.tsv'))

# we have used the API to download fasta links, load them and filter and check that they're all there
fasta_tsv <- fromJSON(paste0(BASE_PATH, 'LEBC/enterobase_meta_with_links.json')) %>% bind_rows()

#length(intersect(fasta_tsv$assembly_barcode, joined_meta_filt$`Assembly barcode`))/nrow(joined_meta_filt)

#joined_meta_filt_with_fasta <- joined_meta_filt %>% left_join(fasta_tsv %>% select(assembly_barcode, download_fasta_link), by=c('Assembly barcode'='assembly_barcode'))

fasta_in_desired_seqs <- fasta_tsv %>% filter(assembly_barcode %in% joined_meta_filt$`Assembly barcode`)
# generate a subset download
KEY_PATH <-  paste0(PROJECT_DIR, 'ecoli_key.tsv')
if (file.exists(KEY_PATH)){
  ecoli_key <- read_tsv(KEY_PATH)
} else {
  ecoli_key <- NA
}

starting_row <- ifelse(any(is.na(ecoli_key)), 1, nrow(ecoli_key)+1)
ending_row <- min(starting_row+1000-1, nrow(fasta_in_desired_seqs)) # enterobase specified increments of 1000
curr_fasta_subset <- fasta_in_desired_seqs[starting_row:ending_row, ] %>% select(assembly_barcode, download_fasta_link, strain_barcode)
curr_fasta_file_name <- paste0('curr_fasta_subset_', starting_row, '_', ending_row, '.tsv')
write_tsv(curr_fasta_subset, paste0(PROJECT_DIR, curr_fasta_file_name))
if (any(is.na(ecoli_key))){
  ecoli_key <- curr_fasta_subset
} else {
  write_tsv(ecoli_key, paste0(PROJECT_DIR, 'ecoli_key_bk.tsv'))
  ecoli_key <- bind_rows(ecoli_key, curr_fasta_subset)
}
write_tsv(ecoli_key, paste0(PROJECT_DIR, 'ecoli_key.tsv'))
DL_SCRIPT_PATH <- paste0(PROJECT_DIR, 'download_fastas.py')
INPUT_TSV <- paste0(PROJECT_DIR, curr_fasta_file_name)
OUTPUT_FOLDER <- paste0(PROJECT_DIR, 'fasta_files')
job_time <- '4:0:0'
python_cmd <- paste0('python -u ', DL_SCRIPT_PATH, ' -i ', INPUT_TSV, ' -o ', OUTPUT_FOLDER)
sbatch_cmd <- sprintf('sbatch -D %s -c 2 --mem=2G -J ecoli_download_%s_%s -p daytime --time %s --wrap=\'%s\'', PROJECT_DIR, starting_row, ending_row, job_time, python_cmd)
print('sbatching with..')
print(sbatch_cmd)
system(sbatch_cmd)
