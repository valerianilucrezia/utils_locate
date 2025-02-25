library(dplyr)
library(ggplot2)
library(patchwork)
library("optparse")


option_list = list(make_option(c("-s", "--simulation"), type="character", default="sim_1", 
              help="simulation ID", metavar="character"),
              make_option(c("-t", "--type"), type="character", default="clonal", 
                          help="type of simulation", metavar="character")
              ); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

dir <- paste0('/Users/lucreziavaleriani/Documents/GitHub/utils_locate/simulations_rRACES/results/', opt$type, '/', opt$simulation, '/')

source('/Users/lucreziavaleriani/Documents/GitHub/utils_locate/simulations_rRACES/utils_smooth.R')
source('/Users/lucreziavaleriani/Documents/GitHub/utils_locate/simulations_rRACES/plot_races.R')
source('/Users/lucreziavaleriani/Documents/GitHub/utils_locate/simulations_rRACES/utils_races.R')


combinations <- list.dirs(dir, full.names = F)

if (opt$type == 'clonal'){
  for (comb in combinations){
    if (grepl('cov', comb)){
        cli::cli_text('Plotting {opt$simulation}-{comb}')
  
        res_dir <- paste0(dir, comb)
        event <- readRDS(paste0(dir, '/cna_event.RDS'))
        bps <- event %>% tidyr::pivot_longer(c(begin, end)) %>% pull(value)
        
        data <- readRDS(paste0(res_dir, '/seq_res.RDS'))
        snv <- get_somatic_data(data, sample = 'Sample')
        snp <- get_germline_data(data, sample = 'Sample') 
        
        filt_snp <- snp %>% 
          filter(VAF.normal > 0.2, VAF.normal < 0.8) %>% 
          filter(VAF.tumour > 0.1, VAF.tumour < 0.99) %>% 
          filter(ref.tumour  %in% c('A', 'C', 'T', 'G')) %>% 
          filter(ALT.tumour  %in% c('A', 'C', 'T', 'G')) %>% 
          filter(ALT.tumour != ref.tumour)
          
        vaf_smooth <- mirro_and_smoothing(sv = snv, 
                                    sp = filt_snp, 
                                    save = FALSE)
        #saveRDS(file = paste0(res_dir, '/mirr_smooth_snv.RDS'), object = vaf_smooth)
        #write.csv(vaf_smooth, file = paste0(res_dir, '/mirr_smooth_snv.csv'), quote = F, row.names = F)
  
        
        plt <- patchwork::wrap_plots(plot_gw(filt_snp, snv, bps), 
                                     #plt_smooth_vaf_median(vaf_smooth, bps),
                                     plot_hist_BAF(filt_snp), 
                                     plot_hist_VAF(snv), 
                                     design = 'AA\nCD', 
                                     guides = 'collect'
                                    )
        # ggsave(filename = paste0(res_dir, '/data.png'), plot = plt)
      }
  }
} else{
  
}

