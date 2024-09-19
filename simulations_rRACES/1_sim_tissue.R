library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)

dir_orfeo <- "/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/sim_data_races/simulate_data/"
setwd(dir_orfeo)

type <- 'clonal'

if (type == 'clonal'){ 
  name_sim <- type
  res_dir <- paste0("/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/sim_data_races/data_races/", name_sim, '/')
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
  
  bbox <- sim$search_sample(c("Clone 2" = 200), 20, 20)
  sim$sample_cells("Sample", bbox$lower_corner, bbox$upper_corner)
  
  muller <- plot_muller(sim)
  tissue <- plot_tissue(sim)
  sim_plot <-  muller + tissue
  ggsave(paste0(res_dir, 'sim.png'), sim_plot)
  
  forest <- sim$get_samples_forest()
  forest$save(paste0(res_dir, "samples_forest.sff"))

} else if (type == 'sub-clonal'){
  name_sim <- type
  res_dir <- paste0("/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/sim_data_races/data_races/", name_sim, '/')
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
