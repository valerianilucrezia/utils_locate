library(rRACES)
library(dplyr)

dir_orfeo <- "/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/sim_data_races/simulate_data/"
setwd(dir_orfeo)
source('./put_mut_random.R')

type <- 'clonal'
base_dir <- paste0('/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/sim_data_races/data_races/', type, '/')

Nsims = seq(30, 30)
Nevents = c(1, 2, 3, 4, 5, 6)
Sizes = c(1e7, 9e6, 6e6, 3e6, 1e6) #1e6, 9e5)
Start = 16000000
End = 51304566
Events = c('2:1', '2:2', '2:0', '1:0', 
           '3:0', '3:1', '3:2', '3:3')


for (sim in Nsims){
  cli::cli_text(cli::col_blue('Simulation {sim}'))
  
  Nev <- sample(Nevents, 1)  
  name_sim <- paste0('sim_', sim)
  res_dir <- paste0(base_dir, name_sim, '/')
  dir.create(res_dir, recursive = T, showWarnings = F)
  
  m_engine <- MutationEngine(setup_code = "demo")
  
  if (type == 'clonal'){
    event <- get_cna_event(Nev, Sizes, Start, End, Events)
    m_engine$add_mutant(mutant_name = "Clone 1",
                        passenger_rates = c(SNV = 9e-8,
                                            CNA = 0),
                        drivers = event[[1]])
    
    m_engine$add_mutant(mutant_name = "Clone 2",
                        passenger_rates = c(SNV = 9e-8, 
                                            CNA = 0),
                        drivers = list("NF2 R221*"))
    
    m_engine$add_exposure(c(SBS5 = 0.3, SBS3 = 0.7))
    
    forest <- load_samples_forest(paste0(base_dir, "samples_forest.sff"))
    
    phylo_forest <- m_engine$place_mutations(forest, 
                                             num_of_preneoplatic_SNVs = 800,
                                             num_of_preneoplatic_indels = 0)
    events <- phylo_forest$get_bulk_allelic_fragmentation('Sample')
    
    saveRDS(events, file = paste0(res_dir, 'cna_event.RDS'))
    phylo_forest$save(paste0(res_dir, "phylo_forest.sff"))
    
  } else if (type == 'sub-clonal'){
    event_1 <- get_cna_event(Nev, Sizes, Start, End, Events)
    m_engine$add_mutant(mutant_name = "Clone 1",
                        passenger_rates = c(SNV = 9e-8,
                                            CNA = 0),
                        drivers = event_1[[1]])
    
    event_2 <- get_cna_event(Nev, Sizes, Start, End, Events)
    m_engine$add_mutant(mutant_name = "Clone 2",
                        passenger_rates = c(SNV = 9e-8, 
                                            CNA = 0),
                        drivers = event_2[[1]])
    
    m_engine$add_mutant(mutant_name = "Clone 3",
                        passenger_rates = c(SNV = 9e-8, 
                                            CNA = 0),
                        drivers = list("NF2 R221*"))
    
    m_engine$add_exposure(c(SBS5 = 0.3, SBS3 = 0.7))
    
    forest <- load_samples_forest(paste0(base_dir, "samples_forest.sff"))
    phylo_forest <- m_engine$place_mutations(forest, 
                                             num_of_preneoplatic_SNVs = 800,
                                             num_of_preneoplatic_indels = 0)
    events <- phylo_forest$get_bulk_allelic_fragmentation('Sample')
    
    
    saveRDS(events, file = paste0(res_dir, 'cna_event.RDS'))
    phylo_forest$save(paste0(res_dir, "phylo_forest.sff"))
  }
}

