library(dplyr)
library(ggplot2)
library(patchwork)
library("optparse")


option_list = list(make_option(c("-s", "--simulation"), type="character", default="sim_1", 
              help="simulation ID", metavar="character")); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

dir_orfeo <- "/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/sim_data_races/simulate_data/"
setwd(dir_orfeo)

source('./utils_smooth.R')
source('./plot_races.R')
source('./utils_races.R')

name_sim <- opt$simulation

path <- paste0("/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/sim_data_races/data_races/", name_sim, '/')
combinations <- list.dirs(path, full.names = F)

for (comb in combinations){
  if (grepl('cov', comb)){
      cli::cli_text('Plotting {name_sim}-{comb}')

      res_dir <- paste0(path, comb)
      event <- readRDS(paste0(path, '/cna_event.RDS'))
      bps <- event %>% tidyr::pivot_longer(c(start, end)) %>% pull(value)
      

      data <- readRDS(paste0(res_dir, '/seq_res.RDS'))
      snv <- get_somatic_data(data, sample = 'Sample2')
      snp <- get_germline_data(data, sample = 'Sample2') 
      
      filt_snp <- snp %>% 
        filter(VAF.normal > 0.2, VAF.normal < 0.8) %>% 
        filter(VAF.tumour > 0.1, VAF.tumour < 0.99) %>% 
        filter(ref.tumour  %in% c('A', 'C', 'T', 'G')) %>% 
        filter(ALT.tumour  %in% c('A', 'C', 'T', 'G')) %>% 
        filter(ALT.tumour != ref.tumour)
        
      vaf_smooth <- mirro_and_smoothing(sv = snv, 
                                  sp = filt_snp, 
                                  save = FALSE)
      saveRDS(file = paste0(res_dir, '/mirr_smooth_snv.RDS'), object = vaf_smooth)
      write.csv(vaf_smooth, file = paste0(res_dir, '/mirr_smooth_snv.csv'), quote = F, row.names = F)

      
      # plt <- patchwork::wrap_plots(plot_gw(filt_snp, snv, bps), plt_smooth_vaf_median(vaf_smooth, bps))
      # ggsave(filename = paste0(res_dir, '/data.png'), plot = plt)
    }
}

