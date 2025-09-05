library(ProCESS)
library(dplyr)
library(ggplot2)
library(patchwork)
library(optparse)


option_list = list(
  make_option(c("-o", "--outdir"), type="character", default="/orfeo/cephfs/scratch/area/lvaleriani/utils_locate/simulations_rRACES/out", 
              help="simulation path", metavar="character"),
  
  make_option(c("-t", "--type"), type="character", default="clonal", 
              help="type of simulation", metavar="character")
  
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# Shoose which type of tumor you want to simulate between clonal and sub-clonal
# For now it works ONLY clonal
if ( opt$type == 'clonal'){ 
  
  res_dir <- file.path(opt$outdir, opt$type)
  dir.create(res_dir, recursive = T, showWarnings = F)
  
  setwd(res_dir)
  
  # Simulate tissue ####
  sim <- TissueSimulation(name = opt$type, seed = 12345, save_snapshots = T)
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
  forest <- sim$get_sample_forest()
  forest$save(paste0(res_dir, "/samples_forest.sff"))
  
  sim_plot <-  muller + tissue + plot_forest(forest) %>% annotate_forest(forest) & theme_minimal() & theme(legend.position = 'none')
  ggsave(filename = file.path(res_dir, '/sim.png'), plot = sim_plot, width = 8, height = 3.5, units = 'in', dpi = 200)
  ggsave(filename = file.path(res_dir, '/sim.pdf'), plot = sim_plot, width = 8, height = 3.5, units = 'in', dpi = 400)

} else if (type == 'sub-clonal'){
  res_dir <- paste0(opt$outdir, '/', opt$type, '/')
  dir.create(res_dir, recursive = T, showWarnings = F)
  
  # Simulate tissue ####
  sim <- SpatialSimulation(name = opt$type, seed = 12345)
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
  sim$sample_cells("Sample_Clonal", bbox$lower_corner, bbox$upper_corner)
  
  sim$add_mutant(name = "Clone 3", growth_rates = 3, death_rates = 0.01)
  sim$update_rates('Clone 2', c(growth = .5, death = 0.1))
  sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")
  sim$run_up_to_size(species = 'Clone 3', 40000) 
  
  bbox <- sim$search_sample(c("Clone 2" = 500, "Clone 3" = 400), 50, 50)
  sim$sample_cells("Sample_Subclonal", bbox$lower_corner, bbox$upper_corner)
  
  muller <- plot_muller(sim)
  tissue1 <- plot_tissue(sim, at_sample = 'Sample_Clonal') 
  tissue2 <- plot_tissue(sim, at_sample = 'Sample_Subclonal') 
  sim_plot <-  muller + tissue1 + tissue2
  ggsave(paste0(res_dir, 'sim.png'), sim_plot)
  
  forest <- sim$get_samples_forest()
  #plot_forest(forest)
  forest$save(paste0(res_dir, "samples_forest.sff"))
  
  info <- get_sample_info(sim, forest)
  info
  saveRDS(object = info, file = paste0(res_dir, "sample_info.rds"))
  }
