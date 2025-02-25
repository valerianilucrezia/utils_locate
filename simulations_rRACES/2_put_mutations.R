library(rRACES)
library(dplyr)
library(optparse)
library(ggplot2)
setwd('/Users/lucreziavaleriani/Documents/GitHub/utils_locate/simulations_rRACES/results/')

# Source auxiliary file for creating random CN events
source('/Users/lucreziavaleriani/Documents/GitHub/utils_locate/simulations_rRACES/put_cn_random.R')
source('/Users/lucreziavaleriani/Documents/GitHub/utils_locate/simulations_rRACES/plot_races.R')

option_list = list(
  make_option(c("-o", "--outdir"), type="character", default="/Users/lucreziavaleriani/Documents/GitHub/utils_locate/simulations_rRACES/results", 
              help="simulation path", metavar="character"),
  
  make_option(c("-t", "--type"), type="character", default="sub-clonal", 
              help="type of simulation", metavar="character")
  
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# Define output directory
base_dir <- paste0(opt$outdir, '/', opt$type, '/')

# Set parameters for simulating CN events, for the moment only in chr22
Nsims = seq(1, 30)                                  # how many simulations you want
Nevents = 1:6                    # possible number of events in each simulation
Sizes = c(1e7, 9e6, 6e6, 3e6, 1e6)                  # possible lengths of events
Start = 16000000                                    # where chr22 start
End = 51304566                                      # where chr22 end
Events = c('2:1', '2:2', '2:0', '1:0',
           '3:0', '3:1', '3:2', '3:3')              # types of events that we are able to simulate

mr_SNV = 1e-8                                    # set mutation rate for SNV

for (sim in Nsims){
  cli::cli_text(cli::col_blue('Simulation {sim}'))
  
  name_sim <- paste0('sim_', sim)
  res_dir <- paste0(base_dir, name_sim, '/')
  dir.create(res_dir, recursive = T, showWarnings = F)
  
  
  if (opt$type == 'clonal'){
    m_engine <- MutationEngine(setup_code = "demo")
    
    Nev <- sample(Nevents, 1)  
    event <- get_cna_event(Nev, Sizes, Start, End, Events)
    saveRDS(object = event, file = paste0(res_dir, 'original_cna_event.RDS'))
    
    m_engine$add_mutant(mutant_name = "Clone 1",
                        passenger_rates = c(SNV = mr_SNV,
                                            CNA = 0),
                        drivers = event[[1]])
    
    m_engine$add_mutant(mutant_name = "Clone 2",
                        passenger_rates = c(SNV = mr_SNV, 
                                            CNA = 0),
                        drivers = list("NF2 R221*"))
    
    m_engine$add_exposure(c(SBS5 = 0.3, SBS3 = 0.7))
    
    forest <- load_samples_forest(paste0(base_dir, "samples_forest.sff"))
    
    phylo_forest <- m_engine$place_mutations(forest, 
                                             num_of_preneoplatic_SNVs = 800,
                                             num_of_preneoplatic_indels = 0)
    # store CNAs events
    events <- phylo_forest$get_bulk_allelic_fragmentation('Sample')
    saveRDS(object = events, file = paste0(res_dir, 'cna_event.RDS'))
    
    plot_events <- plot_allelic_fragmentation(events)
    ggsave(filename = paste0(res_dir, "cna_events.png"), plot = plot_events, dpi = 300)
    
    phylo_forest$save(paste0(res_dir, "phylo_forest.sff"))
    
  } else if (opt$type == 'sub-clonal'){
    m_engine <- MutationEngine(setup_code = "demo")
    cli::cli_text(cli::col_blue('Simulation {sim}'))
    
    Nevents <- 3:8
    Nev <- sample(Nevents, 1)  
    
    events_all <- get_cna_event(Nev, Sizes, Start, End, Events)
    saveRDS(object = events_all, file = paste0(res_dir, 'original_cna_event.RDS'))
    
    N1 <- sample(Nev-1, 1)
    N2 <- Nev - N1
    
    i <- 1
    list <- lapply(events_all[[2]]$nAdd, FUN = function(n){
      end <- n+i-1
      c <- events_all[[1]][i:end]
      i <<- i + n 
      return(c)
    })
    

    event_1 <- events_all[[2]][1:N1,]
    list_1 <- list[1:N1] %>% unlist() %>% unname()
    
    event_2 <- events_all[[2]][N1+1:N2,]
    list_2 <- list[N1+1:N2] %>% unlist() %>% unname()
    
    m_engine$add_mutant(mutant_name = "Clone 1",
                        passenger_rates = c(SNV = mr_SNV,
                                            CNA = 0),
                        drivers = list_1)
    
    m_engine$add_mutant(mutant_name = "Clone 2",
                        passenger_rates = c(SNV = mr_SNV, 
                                            CNA = 0),
                        drivers = list("NF2 R221*"))
    
    m_engine$add_mutant(mutant_name = "Clone 3",
                        passenger_rates = c(SNV = 9e-8, 
                                            CNA = 0),
                        drivers = list_2)
    
    m_engine$add_exposure(c(SBS5 = 0.3, SBS3 = 0.7))
    
    forest <- load_samples_forest(paste0(base_dir, "samples_forest.sff"))
    phylo_forest <- m_engine$place_mutations(forest, 
                                             num_of_preneoplatic_SNVs = 800,
                                             num_of_preneoplatic_indels = 0)
    
    events_1 <- phylo_forest$get_bulk_allelic_fragmentation('Sample_Clonal') %>%
      mutate(type = ifelse(ratio < 0.9, 'sub-clonal', 'clonal'),
             CN = paste(major, minor, sep = ':')) 
    
    events_2 <- phylo_forest$get_bulk_allelic_fragmentation('Sample_Subclonal') %>% 
      mutate(type = ifelse(ratio < 0.9, 'sub-clonal', 'clonal'),
             CN = paste(major, minor, sep = ':'))  %>% 
      mutate(type = ifelse(type == 'sub-clonal' & CN == '1:1', 'sub-clonal_1', type),
             type = ifelse(type == 'sub-clonal' & CN != '1:1', 'sub-clonal_2', type))
    
    
    events <- list('Sample_Clonal' = events_1, 
                   'Sample_Subclonal' = events_2)
    saveRDS(object = events, file = paste0(res_dir, 'cna_event.RDS'))
    
    plot_events <- plot_allelic_fragmentation(events_1, name = 'Sample_Clonal')  + plot_allelic_fragmentation_subclonal(events_2, name = 'Sample_Subclonal') 
    ggsave(filename = paste0(res_dir, "cna_events.png"), plot = plot_events, dpi = 300)
    
    phylo_forest$save(paste0(res_dir, "phylo_forest.sff"))
  }
}

