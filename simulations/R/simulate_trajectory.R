source('Dropbox/CNA/1. Subclonal CN caller/scripts/all_scripts/evolution.R')

karyotypes <- c('1:0', '1:1', '2:0', '2:1', '2:2')
l_ccf <- seq(0.1, 0.9, by = 0.1)
l_purity <- seq(0.05, 1, by = 0.05)

res_dir <- 'Desktop/MultiSegmenter/sub_clonal_peaks/'

karyotypes_1 <- c("1:0", "1:1", "2:0", "2:1", "2:2")
karyotypes_2 <- c("1:0", "1:1", "2:0", "2:1", "2:2")

k1_used <- c()
total_list_of_processes <- list()
for (k1 in karyotypes_1) {
  
  for (k2 in karyotypes_2) {
    if (!(k2 %in% k1_used)) {
      if (sum(karyotype_distances(k1, k2)) <= 2) {
        
        list_of_processes <- get_evolutionary_processes('1:1', k1, k2)
        model_id <- paste(k1, k2, sep='-')
        total_list_of_processes[[model_id]] <- list_of_processes
        
      }
    }
  }
  k1_used <- c(k1_used, k1)
}

names(total_list_of_processes)
obtain_models <- function(total_list_of_processes){
  a <- list()
  
  for (n in names(total_list_of_processes)) {
    dir.create(paste0(res_dir, n, '/'))
    print(n)
    c <- list()

    processes <- total_list_of_processes[[n]]
    karyotypes <- strsplit(n[1], "-")[[1]]
    k1 <- karyotypes[1]
    k2 <- karyotypes[2]
    
    for (ccf in l_ccf){
      for (pr in l_purity){
        
        lab <- paste(ccf, pr ,sep='-')
        
        dd <- expectations_subclonal(list_of_processes = processes, CCF_1 = ccf, purity = pr)
        d <- dd %>% select(peak, model_id, genotype_1, genotype_2, model)
        write.table(d, file = paste0(res_dir, n, '/', lab, '_peaks.tsv'), quote = FALSE, sep = '\t', col.names = TRUE)
        
        c[lab] <- list(d)
      }
    }
    a[n] <- list(c)
  }
  return(a)
}
res <- obtain_models(total_list_of_processes)
saveRDS(res, paste0(res_dir, 'all_peaks.RDS'))
