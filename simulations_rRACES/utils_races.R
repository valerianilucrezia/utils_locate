seq_to_long_cna <- function(seq_results) {
  # Extract sample names from column names
  sample_names <- strsplit(colnames(seq_results)[grepl(".VAF", colnames(seq_results), fixed = TRUE)], ".VAF") %>% unlist()
  
  seq_df <- lapply(sample_names, function(sn) {
    cc <- c("chr", "chr_pos", "ref", "alt", "causes", "classes", 'cna_id', colnames(seq_results)[grepl(paste0(sn, "."), colnames(seq_results), fixed = TRUE)])
    seq_results[, cc] %>%
      `colnames<-`(c("chr", "chr_pos", "ref", "alt", "causes", "classes", "cna_id", "occurences", "coverage", "VAF")) %>%
      dplyr::mutate(sample_name = sn)
  }) %>% do.call("bind_rows", .)
  
  seq_df %>%
    dplyr::rename(chr = chr, from = chr_pos, DP = coverage, NV = occurences, ALT = alt) %>%
    dplyr::mutate(to = from)
}


get_somatic_data <- function(seq_res, sample, sample_n = nrow(seq_res)){
  somatic <- seq_res %>% 
    filter(classes != 'germinal') %>% 
    filter(sample_name == sample) %>% 
    filter(VAF >= 0.05) %>% 
    slice_sample(n = sample_n)
  return(somatic)
}


get_germline_data <- function(seq_res, sample){
  germline <- seq_res %>% 
    filter(classes == 'germinal') 
  
  normal_data <- germline %>% dplyr::filter(sample_name == "normal_sample") %>% mutate(mut_id =  paste(chr, from, to, sep = ":"))
  tumour_data <- germline %>% dplyr::filter(sample_name == sample)  %>% mutate(mut_id =  paste(chr, from, to, sep = ":"))
  
  germline <- tumour_data %>%
    dplyr::left_join(normal_data, suffix = c(".tumour", ".normal"), by = "mut_id") %>%
    dplyr::mutate(DR = DP.tumour / DP.normal) %>% 
    filter(ref.normal != ALT.normal) %>% 
    filter(ref.normal %in% c('A', 'C', 'T', 'G')) %>% 
    filter(ALT.normal %in% c('A', 'C', 'T', 'G')) 
  
  return(germline)
}

get_bp <- function(cna_id){
  bps <- cna_id %>% tidyr::pivot_longer(cols = c(start,end)) %>% pull(value)
  return(bps)
}

