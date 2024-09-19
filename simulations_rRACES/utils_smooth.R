# smooting ####
data_smoothing <- function(sp, res_path, bin_size = 50000, save = TRUE){
  smooth_snp <- sp %>%
    rename(pos = from.tumour) %>%
    mutate(group = ceiling(pos/ bin_size)) %>% 
    group_by(group, cna_id.tumour) %>% 
    summarize(mean_BAF = mean(VAF.tumour),
              median_BAF = median(VAF.tumour),
              mean_DR = mean(DR),
              median_DR = median(DR), 
              nSNP = n_distinct(pos), 
              minSNP = min(pos),
              maxSNP = max(pos),
    )  
  
  
  if (save == TRUE){
    saveRDS(smooth_snp, paste0(res_path, 'smooth_snp.RDS') )
    return(smooth_snp)
  } else if( save == FALSE){
    return(smooth_snp)
  }
}


vaf_smoothing <- function(sv, sp, res_path, save = TRUE, vaf_th = 0.15, wd = 20000){
  snv <- sv %>% filter(VAF >= vaf_th) %>% rename(pos = from) %>% rename(vaf = VAF)
  snp <- sp %>%  rename(pos = from.tumour)
  
  smooth_vaf <- tibble() 
  for (i in seq(1, nrow(snv))){
    tmp <- snv[i,]
    mean_baf <- snp %>% filter(pos > tmp$pos - wd, pos < tmp$pos + wd) %>% pull(VAF.tumour) %>% mean()
    mean_dr <- snp %>% filter(pos > tmp$pos - wd, pos < tmp$pos + wd) %>% pull(DR) %>% mean()
    median_baf <- snp %>% filter(pos > tmp$pos - wd, pos < tmp$pos + wd) %>% pull(VAF.tumour) %>% median()
    median_dr <- snp %>% filter(pos > tmp$pos - wd, pos < tmp$pos + wd) %>% pull(DR) %>% median()
    tmp$mean_dr = mean_dr
    tmp$mean_baf = mean_baf
    tmp$median_dr = median_dr
    tmp$median_baf = median_baf
    smooth_vaf =  bind_rows(smooth_vaf, tmp)
  }
  
  smooth_vaf <- smooth_vaf %>% arrange(desc(pos))
  smooth_vaf <- smooth_vaf %>% filter(!is.na(mean_baf))

  
  if (save == TRUE){
    saveRDS(smooth_vaf, paste0(res_path, 'smooth_vaf.RDS'))
    return(smooth_vaf)
  } else if (save == FALSE){
    return(smooth_vaf)
  }
}

mirro_and_smoothing <- function(sv, sp, res_path, save = TRUE, vaf_th = 0.15, wd = 20000){
  snv <- sv %>% rename(pos = from) %>% rename(vaf = VAF)
  snp <- sp %>%  rename(pos = from.tumour) %>% 
    mutate(new_baf = ifelse(VAF.tumour > 0.6, 1-VAF.tumour, VAF.tumour))
  
  smooth_vaf <- tibble() 
  for (i in seq(1, nrow(snv))){
    tmp <- snv[i,]
    mean_baf <- snp %>% filter(pos > tmp$pos - wd, pos < tmp$pos + wd) %>% pull(new_baf) %>% mean()
    mean_dr <- snp %>% filter(pos > tmp$pos - wd, pos < tmp$pos + wd) %>% pull(DR) %>% mean()
    median_baf <- snp %>% filter(pos > tmp$pos - wd, pos < tmp$pos + wd) %>% pull(new_baf) %>% median()
    median_dr <- snp %>% filter(pos > tmp$pos - wd, pos < tmp$pos + wd) %>% pull(DR) %>% median()
    mean_dp <- snp %>% filter(pos > tmp$pos - wd, pos < tmp$pos + wd) %>% pull(DP.tumour) %>% mean()
    median_dp <- snp %>% filter(pos > tmp$pos - wd, pos < tmp$pos + wd) %>% pull(DP.tumour) %>% median()
    tmp$mean_dr = mean_dr
    tmp$mean_baf = mean_baf
    tmp$median_dr = median_dr
    tmp$median_baf = median_baf
    tmp$mean_dp = mean_dp
    tmp$median_dp = median_dp
    smooth_vaf =  bind_rows(smooth_vaf, tmp)
  }
  
  smooth_vaf <- smooth_vaf %>% arrange(desc(pos))
  smooth_vaf <- smooth_vaf %>% filter(!is.na(mean_baf))
  
  
  if (save == TRUE){
    saveRDS(smooth_vaf, paste0(res_path, 'mirr_smooth_vaf.RDS'))
    write.csv(smooth_vaf, paste0(res_path, 'mirr_smooth_vaf.csv'), quote = F, row.names = F)
    return(smooth_vaf)
  } else if (save == FALSE){
    return(smooth_vaf)
  }
}


