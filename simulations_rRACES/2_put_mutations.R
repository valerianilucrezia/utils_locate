library(rRACES)
library(dplyr)
library(optparse)

# Source auxiliary file for creating random CN events
source('./put_cn_random.R')

option_list = list(
  make_option(c("-o", "--outdir"), type="character", default="", 
              help="simulation path", metavar="character"),
  
  make_option(c("-t", "--type"), type="character", default="clonal", 
              help="type of simulation", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Define output directory
base_dir <- paste0(opt$outdir, opt$type, '/')

# Set parameters for simulating CN events, for the moment only in chr22
Nsims = seq(1, 30)                                  # how many simulations you want
Nevents = c(1, 2, 3, 4, 5, 6)                       # possible number of events in each simulation
Sizes = c(1e7, 9e6, 6e6, 3e6, 1e6)                  # possible lengths of events
Start = 16000000                                    # where chr22 start
End = 51304566                                      # where chr22 end
Events = c('2:1', '2:2', '2:0', '1:0',
           '3:0', '3:1', '3:2', '3:3')              # types of events that we are able to simulate

mr_SNV = 9e-8                                       # set mutation rate for SNV

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

