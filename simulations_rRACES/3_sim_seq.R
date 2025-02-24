library(rRACES)
library(dplyr)
library(optparse)

source('/Users/lucreziavaleriani/Documents/GitHub/utils_locate/simulations_rRACES/utils_races.R')

option_list = list(
  make_option(c("-s", "--simulation"), type="character", default="sim_1", 
              help="simulation ID", metavar="character"),
  
  make_option(c("-p", "--path"), type="character", default="/Users/lucreziavaleriani/Documents/GitHub/utils_locate/simulations_rRACES/results", 
              help="simulation path", metavar="character"),
  
  make_option(c("-t", "--type"), type="character", default="clonal", 
              help="type of simulation", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

res_dir <- paste0(opt$path, '/', opt$type, '/', opt$simulation, '/')

phylo_forest <- load_phylogenetic_forest(paste0(res_dir, "phylo_forest.sff"))
basic_seq <- BasicIlluminaSequencer(4e-3)
cna_seg <- readRDS(paste0(res_dir, 'cna_event.RDS')) %>% 
         mutate(cna_id = paste(chr, begin, end, major, minor, sep = ':'),
         cna = paste(major, minor, sep = ':'))


## simulate seq ####
coverage = c(50, 70, 100)
purity = c(0.9, 0.6, 0.3)

for (cov in coverage){
  for (pur in purity){
    cli::cli_text(cli::col_yellow('Coverage = {cov} & Purity = {pur}'))
    
    data_out <- paste0(res_dir, 'cov_', cov, '_p_', pur)
    dir.create(data_out, recursive = T)

    seq_results <- simulate_seq(phylo_forest, coverage = cov, write_SAM = FALSE, purity = pur, sequencer = basic_seq)
    saveRDS(object = seq_results$mutations, file = paste0(data_out, '/raw_seq_res.RDS'))
    saveRDS(object = seq_results$parameters, file = paste0(data_out, '/params_seq_res.RDS'))
    
    seq <- seq_results$mutations %>% mutate(cna_id = NA, cna = NA)
    for (cna in seq(1,nrow(cna_seg))){
      seg <- cna_seg[cna,]
      seq <- seq %>% mutate(cna_id = ifelse(chr_pos %in% seq(seg$begin, seg$end), seg$cna_id, cna_id),
                                            cna = ifelse(chr_pos %in% seq(seg$begin, seg$end), seg$cna, cna),
                                            M = ifelse(chr_pos %in% seq(seg$begin, seg$end), seg$major, cna),
                                            m = ifelse(chr_pos %in% seq(seg$begin, seg$end), seg$minor, cna))
    }
    
    seq_long <- seq_to_long_cna(seq)
    saveRDS(object = seq_long, file = paste0(data_out, '/seq_res.RDS'))
  }
}
