library(tidyverse)
library(magrittr)
library(Biostrings)
library(data.table)
library(cowplot)
library(RColorBrewer)
library(UpSetR)
library(scales)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-p", "--project_dir", help="directory containing all the relevant files and folders")

PROJECT_DIR <- parser$project_dir
kma_res <- read_tsv(paste0(PROJECT_DIR, 'concatenated_genes_kma_results.tsv'))
kma_res %<>% mutate('assembly_barcode'=gsub('.fasta.res', '', fasta_name))
# also filter based on template identity and template coverage
kma_res %<>% filter(Template_Identity>=90, Template_Coverage>=90)

assembly_gene_matrix <- read_tsv(paste0(PROJECT_DIR, 'assembly_gene_matrix.tsv'))

meta_final <- read_tsv(paste0(PROJECT_DIR, 'joined_meta_filtered_no_missing_with_quast_results_and_pathotypes.tsv'), guess_max = 10000) #%>% 
  #filter(!is.na(ST))


# do some interactive shit using kma results, matrix and metadata with annotated pathotypes
meta_final %<>% mutate(one_pathotype=pathotype_from_assembly %>% str_split(',') %>% lapply(function(x){x[[1]]}) %>% unlist())
orig_pathotypes_formatted <- meta_final$Pathovar %>% sapply(function(x){
  y <- x %>% gsub('Shigella.*', 'EIEC', .) %>% gsub('E. coli - ', '', .) %>% gsub('/', ',', .) %>% 
    gsub('ND|-', 'unassigned', .) %>% gsub('EHEC', 'STEC', .)
  y <- ifelse(is.na(y), 'unassigned', y)
  return(y)
})
pathotype_table <- meta_final$pathotype_from_assembly %>% table()
pathotype_orig_table <- orig_pathotypes_formatted %>% table()

meta_final %<>% mutate(st_no_na=ifelse(is.na(meta_final$ST), 'missing', meta_final$ST))
st_table <- table(meta_final$st_no_na)
meta_final %<>% mutate(st_no_na = ifelse(st_table[st_no_na] < 300, 'less_than_300', st_no_na))
meta_final %>% ggplot(aes(x=pathotype_from_assembly))+geom_histogram(stat='count') + theme(axis.text.x = element_text(angle = -90, vjust = 1, hjust=1))
meta_final %>% ggplot(aes(x=st_no_na))+geom_histogram(stat='count') + theme(axis.text.x = element_text(angle = -90, vjust = 1, hjust=1)) # there's just too many STs to visualize

sum(orig_pathotypes_formatted==meta_final$pathotype_from_assembly)/nrow(meta_final) # only comparing the strings
concordance_check <- function(x,y, sep=','){
  #browser()
  z <- x %>% str_split(',') %>% unlist()
  u <- y %>% str_split(',') %>% unlist()
  n.intersect <- length(intersect(z, u))
  return(n.intersect)
}
pathotype_concordance <- mapply(concordance_check, meta_final$pathotype_from_assembly, orig_pathotypes_formatted)
print(paste0('amount of assemblies for which there is at least one common label between the original pathovar and the assigned pathotype: ', sum(pathotype_concordance > 0)))
# without unassigned
no_unassigned_at_all <- meta_final$pathotype_from_assembly!='unassigned' & orig_pathotypes_formatted != 'unassigned'
pathotype_concordance_no_unassigned <- mapply(concordance_check, 
                                              meta_final$pathotype_from_assembly[no_unassigned_at_all], 
                                              orig_pathotypes_formatted[no_unassigned_at_all])
print(paste0('amount of assemblies for which there is at least one common label between the original pathovar and the assigned pathotype when filtering out all unassigned from both vectors: ', 
             sum(pathotype_concordance_no_unassigned > 0)/sum(no_unassigned_at_all)))

n50_len <- meta_final %>% ggplot(aes(x=N50, y=`Total length`, col=one_pathotype))+ geom_point()
contigs_len <- meta_final %>% ggplot(aes(x=`# contigs`, y=`Total length`, col=one_pathotype))+ geom_point()
len_vp <- meta_final %>% ggplot(aes(y=`Total length`, x=one_pathotype))+geom_violin()
len_release <- meta_final %>% ggplot(aes(y=`Total length`, x=`Release Date`))+geom_point()

contigs_len
# filter the metadata
meta_filt <- meta_final %>% filter((`# contigs` < 800 & `Total length` > 4e+06) & `Total length` < 65e+05)
nrow(meta_final) - nrow(meta_filt)
meta_filt %>% ggplot(aes(x=`# contigs`, y=`Total length`, col=one_pathotype))+ geom_point()
pathotype_table_filtered <- meta_filt$pathotype_from_assembly %>% table()
has_been_filtered <- meta_final %>% filter(`# contigs` >= 800 | `Total length` <= 4e+06)
was_filtered_table <- has_been_filtered$pathotype_from_assembly %>% table()
was_filtered_table %>% names %>% sapply(function(x){was_filtered_table[x]/pathotype_table[x]})

# generate a new pathotype column where we collapse upec and expec
collapse_expec_upec <- function(pathotype.string){
  subbed.string <- pathotype.string %>% gsub('UPEC', 'ExPECjj', .)
  substrings <- subbed.string %>% str_split(',') %>% unlist() %>% unique()
  return(paste(substrings, collapse=','))
}
combined_expec_upec <- pathotype_table_filtered %>% names %>% sapply(collapse_expec_upec)
meta_filt %<>% mutate(pathotype_combined = combined_expec_upec[pathotype_from_assembly])

### plot pathotype frequency for the most common STs
sorted_st_table <- meta_filt %>% pluck('ST') %>% table() %>% sort(decreasing = T)
top_sts <- sorted_st_table[1:50] %>% names()

# sort by n genomes per pathotypes mby
create_table_dfs <- function(st, merge.below=0.05){
  meta.st <- meta_filt %>% filter(ST==st)
  pathotype.table <- meta.st$pathotype_combined %>% table()
  pathotype.fractions <- pathotype.table/sum(pathotype.table)
  pathotype.few.sample <- pathotype.table[pathotype.fractions < merge.below] %>% names()
  meta.st %<>% mutate(pathotype.clean = ifelse(pathotype_combined %in% pathotype.few.sample, paste0('less_than_', as.character(merge.below)), pathotype_combined))
  pathotype.table.clean <- meta.st$pathotype.clean %>% table()
  distribution <- pathotype.table.clean %>% as.numeric()
  fraction <- distribution/sum(distribution)
  pathotypes <- pathotype.table.clean %>% names()
  total <- sum(distribution)
  st.df <- data.frame(distribution=distribution, fraction=fraction, pathotypes=as.factor(pathotypes), ST=st, total=total)
  return(st.df)
}


hand_picked_palettes <- c('Reds', 'Purples', 'Greys', 'Blues', 'Greens', 'RdPu')
#qual_col_pals <- brewer.pal.info[rownames(brewer.pal.info) %in% hand_picked_palettes, ]
qual_col_pals = brewer.pal.info[brewer.pal.info$category %in% c('qual'),]
col_vector = unlist(mapply(brewer.pal, 6, rownames(qual_col_pals)))
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
n_cols <- length(col_vector)
pie(rep(1,n_cols), col=col_vector) # can show us which to excluse
#col_vector_filtered <- col_vector[-c(1, 2, 3, 10, 18, 37, 46, 45, 47, 48, 29, 30, 38, 31, 20, 19, 12, 21)]
col_vector_filtered <- col_vector[-c(11, 16, 14, 2, 3, 4, 5, 6, 7, 9, 18, 27, 16, 24, 23, 15, 19, 41, 22, 38, 21, 20, 30, 1)]
pie(rep(1,length(col_vector_filtered)), col=col_vector_filtered) # can show us which to excluse

table_dfs <- top_sts %>% lapply(create_table_dfs) %>% bind_rows() %>% 
  mutate(ST=factor(ST, levels = top_sts)) %>% 
  mutate(st_n = factor(paste0('ST=', ST, ', n=', total),
                       levels=paste0('ST=', ST, ', n=', total) %>% unique()))


table_dfs %>% ggplot(aes(x='', y=distribution, color=pathotypes, fill=pathotypes)) + 
  geom_bar(stat = 'identity', width = 1, position = position_fill()) + 
  coord_polar(theta="y", start=0, clip='off') +
  facet_wrap(~st_n,nrow=5) +
  theme_void() + 
  scale_fill_manual(values=col_vector_filtered) + 
  scale_color_manual(values=col_vector_filtered) #+ 
  #theme(legend.position="top", legend.direction='horizontal')

table_dfs %>% ggplot(aes(x=pathotypes, y=fraction)) + geom_bar(stat='identity')  +facet_wrap(~ST,nrow=10) + theme(axis.text.x = element_text(angle = -90, vjust = 1, hjust=1))
table_dfs %>% ggplot(aes(x=pathotypes, y=distribution)) + geom_point() + facet_wrap(~ST,nrow=10) + theme(axis.text.x = element_text(angle = -90, vjust = 1, hjust=1))
table_dfs %>% ggplot(aes(x=st_n, y=fraction, fill=pathotypes)) + geom_bar(stat='identity', position='stack') + 
  theme(axis.text.x = element_text(angle = -90, vjust = 1, hjust=1))+
  scale_fill_manual(values=col_vector_filtered) #+ 
  #geom_text(stat='count', aes(label=after_stat(count)), vjust=-1)

### creating a graph, this will give us an idea of connectedness and if applying a venn diagram is feasible
library(igraph)
create_adjacency_matrix <- function(some.table){
  #browser()
  unique.values <- some.table %>% names %>% lapply(str_split, ',') %>% unlist() %>% unique
  n.items <- length(unique.values)
  init.matrix <- matrix(0, nrow = n.items, ncol = n.items)
  rownames(init.matrix) <- unique.values
  colnames(init.matrix) <- unique.values
  for (i in 1:length(some.table)) {
    verts <- names(some.table)[i]
    verts <- str_split(verts, ',') %>% unlist()
    if (length(verts) == 2){
      init.matrix[verts[1], verts[2]] <- 1
      init.matrix[verts[2], verts[1]] <- 1
    }
  }
  return(init.matrix)
}

some_graph <- graph_from_adjacency_matrix(create_adjacency_matrix(pathotype_table_filtered), mode='undirected')
plot(some_graph)

edges_per_pathotype <- unique_values %>% sapply(function(x){str_detect(names(pathotype_table_filtered), x)} %>% sum())
###

### the pathotypes are fully connected, we must use an upset plot instead of basic venn
pathotype_for_upsetr <- setNames(pathotype_table_filtered %>% as.numeric(), 
                                 names(pathotype_table_filtered) %>% gsub(',', '&', .))
n_intersects <- pathotype_for_upsetr %>% names %>% str_detect('&') %>% sum()

unique_values <- pathotype_table_filtered %>% names %>% lapply(str_split, ',') %>% unlist() %>% unique
pathotypes_total_sorted <- setNames(unique_values %>% lapply(function(x){meta_filt$pathotype_from_assembly %>% str_detect(x) %>% sum()}) %>% unlist(),
                                    unique_values) %>% sort(decreasing=T)
upset(fromExpression(pathotype_for_upsetr), 
      nintersects = n_intersects, 
      nsets = length(unique_values), 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1, 
      point.size = 2.8, 
      line.size = 1
)
###

### make a new genome x gene matrix but based only on the gene bases i.e. not variants, this code is from generate_gene_matrix.R
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

create_gene_matrix <- function(assemblies){
  kma.res <- kma_res %>% mutate(gene_bases = ifelse(gene_bases=='Afa', 'afa', gene_bases)) %>% filter(!is.na(gene_bases))
  #browser()
  kma.res.subset <- kma.res %>% filter(assembly_barcode %in% assemblies) # using all genes for more granularity
  unique.genes <- kma.res %>% pluck('gene_bases') %>% unique() # need to have gene bases initialized remember this
  kma.res.by.assembly <- kma.res.subset %$% split(., assembly_barcode)
  gene.matrix <- kma.res.by.assembly %>% sapply(function(x){
    binary.vec <- as.integer(unique.genes %in% x$gene_bases)
  }) %>% t()
  colnames(gene.matrix) <- unique.genes
  return(gene.matrix)
}

# try to visualzie all the points for marc
full_matrix <- meta_filt$`Assembly barcode` %>% create_gene_matrix()
unique_full_matrix <- full_matrix %>% unique()
distance_matrix <- dist(unique_full_matrix, method='manhattan')
full_genomes_mds <- cmdscale(distance_matrix) %>% as_tibble(rownames='assembly_barcode')
full_genomes_mds %<>% left_join(meta_filt %>% select(`Assembly barcode`, pathotype_combined, ST),
                                by = c('assembly_barcode'='Assembly barcode')) %>% 
  rename(c('V1'='mds1', 'V2'='mds2'))
full_genomes_mds %>% ggplot(aes(x=mds1, y=mds2, col=pathotype_combined))+geom_point()
# assign black to unassigned
hexes <- hue_pal()(length(unique(full_genomes_mds$pathotype_combined)))
hexes[length(hexes)] <- '#000000'
full_genomes_mds %>% ggplot(aes(x=mds1, y=mds2, col=pathotype_combined))+geom_point()+
  scale_color_manual(values=hexes)
show_col(hexes)
###



### make entropy plot for the sts
table_df_full_pathotypes <- top_sts %>% lapply(create_table_dfs, merge.below=0) %>% bind_rows() %>% 
  mutate(ST=factor(ST, levels = top_sts)) %>% 
  mutate(st_n = factor(paste(ST, total, sep=', n='),
                       levels=paste(ST, total, sep=', n=') %>% unique()))

# we will just treat unassigned as a pathotype for the probability
shannon_entropy <- function(probabilities){
  entropy <- -sum(probabilities*log(probabilities))
  return(entropy)
}

st_entropies <- table_df_full_pathotypes %$% split(., ST) %>% lapply(function(x){shannon_entropy(x$fraction)}) %>% unlist()
table_df_full_pathotypes %<>% mutate(entropy=st_entropies[ST])
table_df_full_pathotypes %>% ggplot(aes(x=ST, y=entropy))+geom_point()


# generate mds plot with all (core or pan?) genes

gene_absence_tsv <- read_csv(paste0(PROJECT_DIR, 'roary_stuff/roary_output/gene_presence_absence.csv'), guess_max = 100000)
sample_columns <- (gene_absence_tsv %>% names)[15:ncol(gene_absence_tsv)]
gene_matrix_from_cols <- function(assembly.barcode, roary.df){
  return(as.integer(!(is.na(roary.df[[assembly.barcode]]))))
}
roary_gene_mat <- sample_columns %>% sapply(gene_matrix_from_cols, gene_absence_tsv) %>% t()
colnames(roary_gene_mat) <- gene_absence_tsv$Gene
st_718_dist <- dist(roary_gene_mat, method='manhattan')

st_718_mds <- cmdscale(st_718_dist) %>% as_tibble(rownames='assembly_barcode')
st_718_mds %<>% left_join(meta_filt %>% select(`Assembly barcode`, pathotype_combined), 
                           by = c('assembly_barcode'='Assembly barcode')) %>% 
  rename(c('V1'='mds1', 'V2'='mds2'))
st_718_mds %>% ggplot(aes(x=mds1, y=mds2, col=pathotype_combined))+geom_point()
full_genomes_mds %>% filter(assembly_barcode %in% st_718_mds$assembly_barcode) %>% ggplot(aes(x=mds1, y=mds2, col=pathotype_combined))+geom_point()

# heatmap
# get genes with highest variance
gene_vars <- roary_gene_mat %>% colnames() %>% sapply(function(x){var(roary_gene_mat[, x])})
most_var <- gene_vars %>% sort(decreasing=T) %>% head(100) %>% names()
most_var_mat <- roary_gene_mat[, most_var]
pathotype_frequencies <- st_718_mds %$% split(., pathotype_combined) %>% lapply(function(x){
  submat <- most_var_mat[x$assembly_barcode, ]
  return(colMeans(submat))
}) %>% do.call(rbind, .)
test_heatmap <- heatmap(most_var_mat)
pheatmap(pathotype_frequencies)

