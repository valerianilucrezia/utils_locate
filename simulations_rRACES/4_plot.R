library(dplyr)
library(ggplot2)
library(patchwork)
library(ProCESS)
library("optparse")

setwd('/orfeo/cephfs/scratch/area/lvaleriani/utils_locate/simulations_rRACES/out/')
source('../utils_smooth.R')
source('../plot_races.R')
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
type = opt$type
simulation = opt$simulation

dir <- paste0('/orfeo/cephfs/scratch/area/lvaleriani/utils_locate/simulations_rRACES/out/', type, '/', simulation, '/')

combinations <- list.dirs(dir, full.names = F)
if (type == 'clonal'){
  for (comb in combinations){
    if (grepl('cov', comb)){
        cli::cli_text('Processing {simulation}-{comb}')
  
        res_dir <- paste0(dir, comb)
        # event <- readRDS(paste0(dir, '/cna_event.RDS'))
        # bps <- event %>% tidyr::pivot_longer(c(begin, end)) %>% pull(value) %>% unique()
        # 
        # data <- readRDS(paste0(res_dir, '/seq_res.RDS'))
        # snv <- get_somatic_data(data, sample = 'Sample')
        # snp <- get_germline_data(data, sample = 'Sample') 
        # 
        # filt_snp <- snp %>% 
        #   filter(VAF.normal > 0.35, VAF.normal < 0.75) %>%
        #   filter(VAF.tumour > 0.1, VAF.tumour < 0.99) %>%
        #   filter(ref.tumour  %in% c('A', 'C', 'T', 'G')) %>%
        #   filter(ALT.tumour  %in% c('A', 'C', 'T', 'G')) %>%
        #   filter(ALT.tumour != ref.tumour)
        #   
        # vaf_smooth <- mirro_and_smoothing(sv = snv, 
        #                             sp = filt_snp, 
        #                             save = FALSE)
      
        vaf_smooth <- readRDS(paste0(res_dir, '/mirr_smooth_snv.RDS'))
        bps <- vaf_smooth %>% 
          mutate(new_pos = 1:n()) %>% 
          group_by(seg_id,CN) %>% 
          summarise(start = min(new_pos), end = max(new_pos))  %>% 
          ungroup() %>% 
          select(-seg_id) %>% 
          arrange(start)
        write.csv(bps, file = paste0(res_dir, '/mirr_smooth_bps.csv'), quote = F, row.names = F)
        
        # wrap_plots(plot_gw(germline = vaf_smooth %>% rename(VAF.tumour = median_baf, DR = median_dr, CN.tumour = CN, from.tumour = pos), 
        #         somatic = vaf_smooth %>% rename(VAF = vaf, from = pos) %>% filter(VAF > 0.1), bp = NULL),
        # plot_gw(germline = vaf_smooth %>% rename(VAF.tumour = mean_baf, DR = mean_dr, CN.tumour = CN, from.tumour = pos), 
        #           somatic = vaf_smooth %>% rename(VAF = vaf, from = pos) %>% filter(VAF > 0.1), bp = NULL))
        # 
        # saveRDS(file = paste0(res_dir, '/mirr_smooth_snv.RDS'), object = vaf_smooth)
        # write.csv(vaf_smooth, file = paste0(res_dir, '/mirr_smooth_snv.csv'), quote = F, row.names = F)
        
        # plt <- patchwork::wrap_plots(plot_gw(germline = filt_snp, somatic = snv, bp = bps), 
        #                              plot_hist_BAF(filt_snp), 
        #                              plot_hist_VAF(snv), 
        #                              design = 'AA\nCD', 
        #                              guides = 'collect'
        #                             )
        #ggsave(filename = paste0(res_dir, '/data.png'), plot = plt, width = 10, height = 8, units = 'in', dpi = 300)
      }
  }
}

