library(rRACES)
library(dplyr)

library("optparse")

option_list = list(
  make_option(c("-s", "--simulation"), type="character", default="sim_3", 
              help="simulation ID", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


dir_orfeo <- "/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/sim_data_races/simulate_data/"
setwd(dir_orfeo)

name_sim <- opt$simulation
type <- 'clonal'
res_dir <- paste0("/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/sim_data_races/data_races/",type, '/', name_sim, '/')

phylo_forest <- load_phylogenetic_forest(paste0(res_dir, "phylo_forest.sff"))
basic_seq <- BasicIlluminaSequencer(4e-3)
cna_seg <- readRDS(paste0(res_dir, 'cna_event.RDS')) %>% 
  mutate(cna_id = paste(chr, begin, end, major, minor, sep = ':'),
         cna = paste(major, minor, sep = ':'))


source('./utils_races.R')
## simulate seq ####
coverage = c(50, 70, 100)
purity = c(0.9, 0.6, 0.3)

for (cov in coverage){
  for (pur in purity){
    cli::cli_text(cli::col_yellow('Coverage = {cov} & Purity = {pur}'))
    
    data_out <- paste0(res_dir, 'cov_', cov, '_p_', pur)
    dir.create(data_out, recursive = T)

    seq_results <- simulate_seq(phylo_forest, coverage = cov, write_SAM = FALSE, purity = pur, sequencer = basic_seq)
    saveRDS(seq_results, file = paste0(data_out, '/raw_seq_res.RDS'))
    
    seq_results <- seq_results %>% mutate(cna_id = NA,
                                          cna = NA)
    for (cna in seq(1,nrow(cna_seg))){
      seg <- cna_seg[cna,]
      seq_results <- seq_results %>% mutate(cna_id = ifelse(chr_pos %in% seq(seg$begin, seg$end), seg$cna_id, cna_id),
                                            cna = ifelse(chr_pos %in% seq(seg$begin, seg$end), seg$cna, cna))
    }
    
    seq_results <- seq_to_long_cna(seq_results)
    saveRDS(seq_results, file = paste0(data_out, '/seq_res.RDS'))
  }
}
