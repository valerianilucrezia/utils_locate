library(ProCESS)
library(dplyr)
library(optparse)
library(ggplot2)
library(patchwork)
setwd('/orfeo/cephfs/scratch/area/lvaleriani/utils_locate/simulations_rRACES/out')

source('../utils_races.R')
source('../plot_races.R')

for (s in 1:60){
  sim_n = paste0('sim_',s)
  cna_events <- readRDS(paste0('/orfeo/cephfs/scratch/area/lvaleriani/utils_locate/simulations_rRACES/out/clonal/',sim_n,'/cna_event.RDS'))
  tmp <- cna_events %>% mutate(len = end-begin, CN = paste(major, minor, sep = ':'))
  if (nrow(tmp < 4)){
    if (sum(tmp$CN %in% c('3:2', '3:3', '3:0', '3:1') %>% sum())  == 0){
      print(sim_n)
    }
  }
}
sim_n = "sim_43"
cna_events <- readRDS(paste0('/orfeo/cephfs/scratch/area/lvaleriani/utils_locate/simulations_rRACES/out/clonal/',sim_n,'/cna_event.RDS'))
frag <- plot_allelic_fragmentation(events = cna_events) 

seq_results <- readRDS('clonal/sim_43/cov_70_p_0.9/raw_seq_res.RDS')
cna_seg <- cna_events %>% mutate(seg_id = paste(chr, begin, end, major, minor, sep = ':'), CN = paste(major, minor, sep = ':'))
seq_cna <- seq_results %>% 
  left_join(cna_seg) %>% 
  filter(chr_pos >= begin, chr_pos <= end)
data <- seq_to_long_cna(seq_cna)

snv <- get_somatic_data(data, sample = 'Sample')
snp <- get_germline_data(data, sample = 'Sample') 
bps <- cna_seg %>% tidyr::pivot_longer(c(begin, end)) %>% pull(value)

filt_snp <- snp %>%
  filter(VAF.normal > 0.35, VAF.normal < 0.75) %>%
  filter(VAF.tumour > 0.1, VAF.tumour < 0.99) %>%
  filter(ref.tumour  %in% c('A', 'C', 'T', 'G')) %>%
  filter(ALT.tumour  %in% c('A', 'C', 'T', 'G')) %>%
  filter(ALT.tumour != ref.tumour)
filt_snp <- filt_snp %>% sample_n(size = 1e4)

cna_new  = cna_events %>%
  mutate(begin = ifelse(begin == 1, min(filt_snp$from.tumour), begin),
         end = ifelse(end == 51304566, max(filt_snp$from.tumour), end))
frag_new <- plot_allelic_fragmentation(events = cna_new)+ 
  xlab('')  + 
  theme(legend.position = 'none', axis.text.x = element_blank())
  

plt_1 <- frag_new +
  plot_gw(germline = filt_snp, somatic = snv, bp = NULL) + 
  plot_layout(design = 'A\nB\nB\nB')

plt <- wrap_plots(ggplot(), plt_1, plot_hist_BAF(filt_snp),
           plot_hist_VAF(snv),
           design = 'AAAA\nBBCD\nBBCD',
           guides = 'collect')
#plt
ggsave(filename = '/orfeo/cephfs/scratch/area/lvaleriani/tesi/prova.png', plot = plt, width = 10, height = 9, units = 'in', dpi = 400)
ggsave(filename = '/orfeo/cephfs/scratch/area/lvaleriani/tesi/prova.pdf', plot = plt, width = 10, height = 9, units = 'in')
