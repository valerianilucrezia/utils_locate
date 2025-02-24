library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)
library(optparse)


option_list = list(
  make_option(c("-o", "--outdir"), type="character", default="/Users/lucreziavaleriani/Documents/GitHub/utils_locate/simulations_rRACES/results", 
              help="simulation path", metavar="character"),
  
  make_option(c("-t", "--type"), type="character", default="clonal", 
              help="type of simulation", metavar="character")
  
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# Shoose which type of tumor you want to simulate between clonal and sub-clonal
# For now it works ONLY clonal
if ( opt$type == 'clonal'){ 
  
  res_dir <- paste0(opt$outdir, '/', opt$type, '/')
  dir.create(res_dir, recursive = T, showWarnings = F)
  
  setwd(res_dir)
  
  # Simulate tissue ####
  sim <- SpatialSimulation(name =  opt$type, seed = 12345, save_snapshots = T)
  sim$history_delta <- 1
  sim$death_activation_level <- 50
  
  sim$add_mutant(name = "Clone 1", growth_rates = 0.3, death_rates = 0.01)
  sim$place_cell("Clone 1", 500, 500)
  sim$run_up_to_size("Clone 1", 5000)
  
  sim$add_mutant(name = "Clone 2", growth_rates = 1.5, death_rates = 0.01)
  sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
  sim$run_up_to_size(species = 'Clone 2', 5000)
  
  sim$update_rates('Clone 1', c(growth = 0.05, death = 0.05))
  sim$run_up_to_size(species = 'Clone 2', 15000)
  
  sim$update_rates('Clone 1', c(growth = 0.01, death = 0.1))
  sim$run_up_to_size(species = 'Clone 2', 35000)
  
  sim$update_rates('Clone 1', c(growth = 0.001, death = 1))
  sim$run_up_to_size(species = 'Clone 2', 60000)
  
  bbox <- sim$search_sample(c("Clone 2" = 1000), 50, 50)
  sim$sample_cells("Sample", bbox$lower_corner, bbox$upper_corner)
  
  muller <- plot_muller(sim)
  tissue <- plot_tissue(sim)
  sim_plot <-  muller + tissue
  sim_plot
  sim
  ggsave(filename = paste0(res_dir, 'sim.png'), plot = sim_plot)
  
  forest <- sim$get_samples_forest()
  plot_forest(forest)
  
  forest$save(paste0(res_dir, "samples_forest.sff"))

} else if (type == 'sub-clonal'){
  res_dir <- paste0(opt$outdir, opt$type, '/')
  dir.create(res_dir, recursive = T, showWarnings = F)
  
  # Simulate tissue ####
  sim <- SpatialSimulation(name = name_sim, seed = 12345)
  sim$history_delta <- 1
  sim$death_activation_level <- 50
  
  sim$add_mutant(name = "Clone 1", growth_rates = 0.3, death_rates = 0.01)
  sim$place_cell("Clone 1", 100, 100)
  sim$run_up_to_size("Clone 1", 500)
  
  sim$add_mutant(name = "Clone 2", growth_rates = 1, death_rates = 0.01)
  sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
  sim$run_up_to_size(species = 'Clone 2', 1000)
  
  sim$update_rates('Clone 1', c(growth = 0.01, death = 1.2))
  sim$run_up_to_size(species = 'Clone 2', 3500)
  
  sim$add_mutant(name = "Clone 3", growth_rates = 3, death_rates = 0.01)
  sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")
  sim$run_up_to_size(species = 'Clone 3', 4000) 
  
  bbox <- sim$search_sample(c("Clone 2" = 150, "Clone 3" = 150), 20, 20)
  sim$sample_cells("Sample", bbox$lower_corner, bbox$upper_corner)
  
  muller <- plot_muller(sim)
  tissue <- plot_tissue(sim)
  sim_plot <-  muller + tissue
  ggsave(paste0(res_dir, 'sim.png'), sim_plot)
  
  forest <- sim$get_samples_forest()
  forest$save(paste0(res_dir, "samples_forest.sff"))
  }
