library(rRACES)
library(dplyr)
event_list <- list('1:0' = 1,
                   '2:1' = 1,
                   '2:2' = 2,
                   '2:0' = 2,
                   '3:1' = 2,
                   '3:0' = 3,
                   '3:2' = 3,
                   '3:3' = 4)


get_cna_event = function(nE, Sizes, Start, End, Events){
  Len = End - Start
  
  full_len = 0 
  while (full_len < Len){
    lenghts = sample(Sizes, size = nE, replace = T)
    full_len = sum(lenghts)
    if (full_len < Len){
      break
    }
  }
  
  # define start pos of segments
  missing_pos = Len - full_len - 150
  normal_seg = c(2,3,4)
  nNorm = sample(normal_seg, 1)
  lenNorm = round(missing_pos/nNorm)
  
  all_len = c(rep(lenNorm, nNorm), lenghts)
  all_len = sample(all_len)

  
  tmp = Start
  start_pos = c()
  end_pos = c()
  for (seg in seq(1, length(all_len))){
    start_pos = c(start_pos, tmp)
    e_pos = tmp + all_len[seg]
    end_pos = c(end_pos, e_pos)
    
    tmp = tmp + all_len[seg] + 1
  }

  all_event = dplyr::tibble(start = start_pos, 
                            end = end_pos, 
                            len = all_len)
  
  idx_event = which(all_len != lenNorm)
  pos_event = start_pos[idx_event]
  len_event = all_len[idx_event]
  
  current_ev = sample(Events, size = nE, replace = F)
  
  event_df = dplyr::tibble(events = current_ev, 
                           length = len_event, 
                           start = pos_event) %>%
              mutate(end = start + length) %>% 
              select(events, start, end , length) %>% 
              mutate(nAdd = event_list[current_ev] %>% unlist())
  
  
  CNAs = list()
  
  nAllele = 2
  for (event in seq(1, nE)){
    tmp_event = event_df[event,]
    CN = tmp_event$events
    size = tmp_event$length
    
    if (CN == '1:0') {
      CNAs = append(CNAs, CNA(type = "D", "22", chr_pos = pos_event[event], len = size))
      
    } else if (CN == '2:1'){
      CNAs = append(CNAs, CNA(type = "A", "22", chr_pos = pos_event[event], len = size, src_allele = 0, allele = nAllele))
      nAllele = nAllele + 1 
  
    } else if (CN == '2:2'){
      CNAs = append(CNAs, CNA(type = "A", "22", chr_pos = pos_event[event], len = size, src_allele = 0, allele = nAllele))
      nAllele = nAllele + 1 
      
      CNAs = append(CNAs,CNA(type = "A", "22", chr_pos = pos_event[event], len = size, src_allele = 1, allele = nAllele))
      nAllele = nAllele + 1 
  
    } else if(CN == '2:0') {
      CNAs = append(CNAs, CNA(type = "A", "22", chr_pos = pos_event[event], len = size, src_allele = 0, allele = nAllele))
      nAllele = nAllele + 1 
      
      CNAs = append(CNAs, CNA(type = "D", "22", chr_pos = pos_event[event], len = size,  allele = 1))
      
    } else if(CN == '3:1'){
      CNAs = append(CNAs, CNA(type = "A", "22", chr_pos = pos_event[event], len = size, src_allele = 0, allele = nAllele))
      nAllele = nAllele + 1 
      CNAs = append(CNAs, CNA(type = "A", "22", chr_pos = pos_event[event], len = size,  src_allele = 0, allele = nAllele))
      nAllele = nAllele + 1 
      
    } else if(CN == '3:2'){
      print(nAllele)
      
      CNAs = append(CNAs, CNA(type = "A", "22", chr_pos = pos_event[event], len = size, src_allele = 0, allele = nAllele))
      nAllele = nAllele + 1 
      CNAs = append(CNAs, CNA(type = "A", "22", chr_pos = pos_event[event], len = size,  src_allele = 0, allele = nAllele))
      nAllele = nAllele + 1 
      CNAs = append(CNAs, CNA(type = "A", "22", chr_pos = pos_event[event], len = size,  src_allele = 1, allele = nAllele))
      nAllele = nAllele + 1 

    } else if(CN == '3:0'){
      CNAs = append(CNAs, CNA(type = "A", "22", chr_pos = pos_event[event], len = size, src_allele = 0, allele = nAllele))
      nAllele = nAllele + 1 
      CNAs = append(CNAs, CNA(type = "A", "22", chr_pos = pos_event[event], len = size,  src_allele = 0, allele = nAllele))
      nAllele = nAllele + 1 
      CNAs = append(CNAs, CNA(type = "D", "22", chr_pos = pos_event[event], len = size,  allele = 1))
      
    } else if(CN == '3:3'){
      CNAs = append(CNAs, CNA(type = "A", "22", chr_pos = pos_event[event], len = size, src_allele = 0, allele = nAllele))
      nAllele = nAllele + 1 
      CNAs = append(CNAs, CNA(type = "A", "22", chr_pos = pos_event[event], len = size,  src_allele = 0, allele = nAllele))
      nAllele = nAllele + 1 
      CNAs = append(CNAs, CNA(type = "A", "22", chr_pos = pos_event[event], len = size,  src_allele = 1, allele = nAllele))
      nAllele = nAllele + 1 
      CNAs = append(CNAs, CNA(type = "A", "22", chr_pos = pos_event[event], len = size,  src_allele = 1, allele = nAllele))
      nAllele = nAllele + 1 
    } 
    
  }
  return(list(CNAs, event_df, all_event))
}


