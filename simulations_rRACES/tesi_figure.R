library(ProCESS)
library(dplyr)
library(optparse)
library(ggplot2)
library(patchwork)
setwd('/orfeo/cephfs/scratch/area/lvaleriani/utils_locate/simulations_rRACES/out')

source('../utils_races.R')
source('../plot_races.R')

sim_n = 'sim_21'
cna_events <- readRDS(paste0('/orfeo/cephfs/scratch/area/lvaleriani/utils_locate/simulations_rRACES/out/clonal/',sim_n,'/cna_event.RDS'))
frag <- plot_allelic_fragmentation(events = cna_events)

seq_results <- readRDS('clonal/sim_21/cov_70_p_0.9/raw_seq_res.RDS')
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


plt_1 <- frag + xlab('chromosome') + theme(legend.position = 'none') +
  plot_gw(germline = filt_snp %>% sample_n(size = 1e4), somatic = snv, bp = NULL) + 
  plot_layout(design = 'A\nB\nB\nB')

plt <- wrap_plots(sim_plot, plt_1, plot_hist_BAF(filt_snp),
           plot_hist_VAF(snv),
           design = 'AAAA\nBBCD\nBBCD',
           guides = 'collect')
#plt
ggsave(filename = 'tesi/prova.png', plot = plt, width = 10, height = 9, units = 'in', dpi = 400)
ggsave(filename = 'tesi/prova.pdf', plot = plt, width = 10, height = 9, units = 'in')
