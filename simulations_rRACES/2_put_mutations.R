library(ProCESS)
library(dplyr)
library(optparse)
library(ggplot2)
setwd('/orfeo/cephfs/scratch/area/lvaleriani/utils_locate/simulations_rRACES/out')

source('../put_cn_random.R')
source('../plot_races.R')

set.seed(1234)

option_list = list(
  make_option(c("-o", "--outdir"), type="character", 
              default="/orfeo/cephfs/scratch/area/lvaleriani/utils_locate/simulations_rRACES/out", 
              help="simulation path", metavar="character"),
  
  make_option(c("-t", "--type"), type="character", default="clonal", 
              help="type of simulation", metavar="character")
  
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

base_dir <- file.path(opt$outdir, opt$type)

Nsims = seq(1,60)                                 
Nevents = 1:6                                  
Sizes = c(2e7,1e7,9e6,7e6,5e6,3e6,2e6,1e6,9e5,7e5,5e5)   
Start = 16000000                                    
End = 51304566    
Events = c('2:1', '2:2', '2:0', '1:0','3:0', '3:1', '3:2', '3:3')            
mr_SNV = 8e-8                                  

sim=1
for (sim in Nsims){
  cli::cli_text(cli::col_blue('Simulation {sim}'))
  
  name_sim <- paste0('sim_', sim)
  res_dir <- file.path(base_dir, name_sim)
  dir.create(res_dir, recursive = T, showWarnings = F)
  
  m_engine <-  MutationEngine(setup_code = "demo", 
                              context_sampling = 20)
  
  Nev <- sample(Nevents, 1)  
  event <- get_cna_event(nE = Nev, Sizes = Sizes, Start = Start, End = End, Events = Events, Chr = "22")
  saveRDS(object = event, file = paste0(res_dir, '/original_cna_event.RDS'))
  
  m_engine$add_mutant(mutant_name = "Clone 1",
                      passenger_rates = c(SNV = mr_SNV,
                                          CNA = 0),
                      drivers = event[[1]])
  
  m_engine$add_mutant(mutant_name = "Clone 2",
                      passenger_rates = c(SNV = mr_SNV, 
                                          CNA = 0 ),
                      drivers = list("NF2 R262*"))
  
  m_engine$add_exposure(c(SBS1 = 0.3, SBS5 = 0.7, ID1 = 1))
  
  forest <- load_sample_forest(paste0(base_dir, "/samples_forest.sff"))
  
  phylo_forest <- m_engine$place_mutations(forest, 
                                           num_of_preneoplatic_SNVs = 800,
                                           num_of_preneoplatic_indels = 200)
  
  # store CNAs events
  events <- phylo_forest$get_bulk_allelic_fragmentation('Sample') %>% filter(chr == 1)
  saveRDS(object = events, file = paste0(res_dir, '/cna_event.RDS'))
  
  plot_events <- plot_allelic_fragmentation(events = events)
  ggsave(filename = paste0(res_dir, "/cna_events.png"), plot = plot_events, dpi = 300)
  
  phylo_forest$save(paste0(res_dir, "/phylo_forest.sff"))
}

