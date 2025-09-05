seq_to_long_cna <- function(seq_results) {
  # Extract sample names from column names
  sample_names <- strsplit(colnames(seq_results)[grepl(".VAF", colnames(seq_results), fixed = TRUE)], ".VAF") %>% unlist()
  
  seq_df <- lapply(sample_names, function(sn) {
    cc <- c("chr", "chr_pos", "ref", "alt", "causes", "classes", "major", "minor", "ratio", "seg_id", "CN", colnames(seq_results)[grepl(paste0(sn, "."), colnames(seq_results), fixed = TRUE)])
    seq_results[, cc] %>%
      `colnames<-`(c("chr", "chr_pos", "ref", "alt", "causes", "classes", "major", "minor", "ratio", "seg_id", "CN", "occurences", "coverage", "VAF")) %>%
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


get_sample_info <- function(sim, forest){
  info = sim$get_samples_info()
  nodes = forest$get_nodes()
  clones = nodes %>% 
    dplyr::filter(!is.na(sample)) %>% 
    dplyr::group_by(sample, mutant) %>% 
    dplyr::pull(mutant) %>% 
    unique()
  clones_of_origin = nodes %>%
    dplyr::filter(!is.na(sample)) %>% 
    dplyr::group_by(sample, mutant) %>% 
    # dplyr::mutate(mutant = gsub(" ", "_", mutant)) %>% 
    dplyr::count(mutant) %>% 
    tidyr::pivot_wider(values_from = n, names_from = mutant, values_fill = 0) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(sample_type = sum(c_across(clones) == 0)) %>% 
    dplyr::mutate(sample_type = ifelse(sample_type == (length(clones)-1), "Monoclonal", "Polyclonal"))
  
  info = dplyr::full_join(info, clones_of_origin, by = join_by("name" == "sample"))
  
  info = info %>% 
    dplyr::group_by(time) %>% 
    dplyr::group_split() 
  info = lapply(1:length(info), function(x) {
    oo = info[[x]] %>% 
      dplyr::mutate(t = paste0("t~", x, "~"))
  }) %>% 
    bind_rows() %>% select(-id, -xmin, -xmax, -ymin, -ymax, -tumour_cells_in_bbox, -t) %>%
    mutate(across(starts_with("Clone"), ~ . / tumour_cells, .names = "CFF_{.col}"))
}
