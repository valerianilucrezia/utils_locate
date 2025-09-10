library(ProCESS)
library(dplyr)
library(optparse)

setwd('/orfeo/cephfs/scratch/area/lvaleriani/utils_locate/simulations_rRACES/out')
source('../utils_races.R')

option_list = list(
  make_option(c("-s", "--simulation"), type="character", default="sim_1", 
              help="simulation ID", metavar="character"),
  
  make_option(c("-p", "--path"), type="character", 
              default="/orfeo/cephfs/scratch/area/lvaleriani/utils_locate/simulations_rRACES/out", 
              help="simulation path", metavar="character"),
  
  make_option(c("-t", "--type"), type="character", default="clonal", 
              help="type of simulation", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

res_dir <- paste0(opt$path, '/', opt$type, '/', opt$simulation, '/')

phylo_forest <- load_phylogenetic_forest(paste0(res_dir, "phylo_forest.sff"))
basic_seq <- BasicIlluminaSequencer(4e-3)

## simulate seq ####
coverage = c(30, 50, 70)
purity = c(0.9, 0.6, 0.3)

cov= coverage[[1]]
pur = purity[[1]]

cna_seg <- readRDS(paste0(res_dir, 'cna_event.RDS')) %>% 
         mutate(seg_id = paste(chr, begin, end, major, minor, sep = ':'),
                CN = paste(major, minor, sep = ':'))

for (cov in coverage){
  for (pur in purity){
    cli::cli_text(cli::col_yellow('Coverage = {cov} & Purity = {pur}'))
    
    data_out <- paste0(res_dir, 'cov_', cov, '_p_', pur)
    dir.create(data_out, recursive = T, showWarnings = F)

    seq_results <- simulate_seq(phylo_forest, coverage = cov, write_SAM = F, purity = pur, sequencer = basic_seq, with_normal_sample = T)
    saveRDS(object = seq_results$mutations, file = paste0(data_out, '/raw_seq_res.RDS'))
    saveRDS(object = seq_results$parameters, file = paste0(data_out, '/params_seq_res.RDS'))
    
    seq <- seq_results$mutations
    seq_cna <- seq %>% 
      left_join(cna_seg) %>% 
      filter(chr_pos >= begin, chr_pos <= end)
    seq_long <- seq_to_long_cna(seq_cna)
    saveRDS(object = seq_long, file = paste0(data_out, '/seq_res.RDS'))
  }
}