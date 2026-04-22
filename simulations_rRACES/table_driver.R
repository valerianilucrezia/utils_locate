library(rRACES)
library(dplyr)
library(optparse)
library(ggplot2)
library(patchwork)
setwd('~/Documents/GitHub/utils_locate/simulations_rRACES/')

res_dir <- '~/Documents/GitHub/utils_locate/simulations_rRACES/results/sub-clonal/sim_1/'
seq_results <- readRDS('~/Documents/GitHub/utils_locate/simulations_rRACES/results/sub-clonal/sim_1/cov_50_p_0.9/raw_seq_res.RDS')
phylo_forest <- rRACES::load_phylogenetic_forest('/Users/lucreziavaleriani/Documents/GitHub/utils_locate/simulations_rRACES/results/sub-clonal/sim_1/phylo_forest.sff')

driver <- phylo_forest$get_driver_mutations()
info <- phylo_forest$get_species_info()

get_event <- function(driver){
  event <- driver %>% 
    filter(type == 'CNA') %>% 
    select(-ref, -alt) %>% 
    group_by(mutant, chr, start, end, type) %>% 
    summarize(A0 = sum(CNA_type == 'A' & src_allele == '0'),
              A1 = sum(CNA_type == 'A' & src_allele == '1'),
              D0 = sum(CNA_type == 'D' & allele == '0'),
              D1 = sum(CNA_type == 'D' & allele == '1')) %>% 
    mutate(major = 1 + A0 - D0,
           minor = 1 + A1 - D1) %>% 
    mutate(CN = paste(major, minor, sep = ':'),
           length = end - start) %>% 
    select(mutant, type, chr,  start, end, length, major, minor, CN) 
  
  return(event)
}

cna_event <- get_event(driver)
mut_event <- driver %>% filter(type != 'CNA') %>% select(-CNA_type, -order, -allele, -src_allele)


patchwork::wrap_table(cna_event) + patchwork::wrap_table(mut_event)  + patchwork::wrap_table(info) + plot_layout(nrow = 3)
