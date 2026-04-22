library(rRACES)
library(dplyr)
library(optparse)
library(ggplot2)
library(patchwork)
setwd('~/Documents/GitHub/utils_locate/simulations_rRACES/')
source('~/Documents/GitHub/utils_locate/simulations_rRACES/utils_plot.R')
source('~/Documents/GitHub/utils_locate/simulations_rRACES/utils_races.R')

res_dir <- '~/Documents/GitHub/utils_locate/simulations_rRACES/results/sub-clonal/sim_1/'
seq_results <- readRDS('~/Documents/GitHub/utils_locate/simulations_rRACES/results/sub-clonal/sim_1/cov_50_p_0.9/raw_seq_res.RDS')
phylo_forest <- rRACES::load_phylogenetic_forest('/Users/lucreziavaleriani/Documents/GitHub/utils_locate/simulations_rRACES/results/sub-clonal/sim_1/phylo_forest.sff')


seq <- seq_to_long(seq_results)
cna_seg <- readRDS(paste0(res_dir, 'cna_event.RDS'))


muts_cn <- list()
for (sample in names(cna_seg)){
  tmp_cna <- cna_seg[[sample]] %>% rename(CN_type = type)
  
  for (cna in seq(1,nrow(tmp_cna))){
    seg <- tmp_cna[cna,]
    
    muts_cn[[sample]] <- seq %>% 
      filter(sample_name == sample) %>% 
      inner_join(tmp_cna, by = c("chr"), relationship = "many-to-many") %>%  
      filter(from >= begin & to <= end) %>% 
      mutate(seg_id = paste0("chr", chr, ":", begin, ":", end))
  }
}

offset <- 0.05
muts_cn <- muts_cn %>% bind_rows()
cn <- muts_cn %>% select(chr, sample_name, begin, end, CN_type, CN, seg_id, ratio, major, minor) %>% unique() %>%  
  mutate(CN_type = ifelse(CN_type == 'clonal', 'clonal', 'sub-clonal')) %>% 
  mutate(chr = paste0('chr', chr))
cn <- absolute_to_relative_coordinates(cn)

genome_subclonal <- blank_genome(chromosomes = c('chr22')) + 
  geom_rect(data =  cn, aes(xmin = begin, xmax = end, ymin = -Inf, ymax = Inf, fill = factor(CN)), alpha = 0.3) +
  geom_segment(data =  cn, aes(x = begin, xend = end, y = major+offset, yend = major+offset), col = 'red', size = 2) +
  geom_segment(data =  cn, aes(x = begin, xend = end, y = minor-offset, yend = minor-offset), col = 'steelblue', size = 2) +
  scale_fill_manual(values = karyo_colors) + 
  ylab('CN') + 
  xlab('position') + 
  ggplot2::guides(fill = ggplot2::guide_legend('CN', override.aes = list(alpha = 1))) + 
  ylim(-0.2, 3.2) + 
  ggh4x::facet_nested(sample_name + CN_type + ratio ~.)

forest <- rRACES::load_phylogenetic_forest('results/sub-clonal/sim_1/phylo_forest.sff')
cna_seg <- lapply(names(cna_seg), FUN = function(sample){
  cna_seg[[sample]] %>% mutate(sample = sample)
})
cna_segs <- cna_seg %>% bind_rows() %>% select(-ratio) %>% 
  mutate(chr = paste0('chr', chr))   %>% 
  mutate(seg_id =paste(chr, begin, end, sep = ':'))

cell_cna <- forest$get_sampled_cell_CNAs() %>% 
  mutate(chr = paste0('chr', chr))   %>% 
  mutate(seg_id = paste(chr, begin, end, sep = ':')) %>% 
  select(-allele, -src.allele, -class, -type) 

cell_muts <- forest$get_sampled_cell_mutations() %>% 
  select(-allele, -class, -cause) %>% 
  mutate(chr = paste0('chr', chr))

sampled_cells <- forest$get_nodes() %>%
  filter(!is.na(.data$sample)) %>% 
  select(cell_id, mutant, sample)

cell_muts <- cell_muts %>% left_join(sampled_cells)
cell_cna <- cell_cna %>% left_join(sampled_cells) 

cell_cna <- left_join(cell_cna, cna_segs, relationship = "many-to-many") %>% unique() 


muts_cn %>% 
  filter(classes != 'germinal') %>% 
  filter(CN_type != 'clonal') %>% 
  #mutate(id = paste(CN, ratio, sep = ':'))  %>% 
  ggplot() +
  geom_histogram(aes(x = VAF, fill = CN), binwidth = 0.02, position = 'identity') +
  scale_fill_manual(values = karyo_colors) + 
  xlim(-0.01, 1) +
  ggh4x::facet_nested(sample_name + CN_type + ratio ~ chr + seg_id, scales = 'free_y', independent = 'y') +
  CNAqc:::my_ggplot_theme(cex = cex) 


muts_cn %>% 
  filter(classes != 'germinal') %>% 
  filter(CN_type == 'clonal') %>% 
  #mutate(id = paste(CN, ratio, sep = ':'))  %>% 
  ggplot() +
  geom_histogram(aes(x = VAF, fill = CN), binwidth = 0.01, position = 'identity') +
  scale_fill_manual(values = karyo_colors) + 
  xlim(-0.01, 1) +
  ggh4x::facet_nested(sample_name + CN_type + ratio ~ chr, scales = 'free_y', independent = 'y') +
  CNAqc:::my_ggplot_theme(cex = cex) 

  
