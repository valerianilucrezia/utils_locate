library(ggplot2)
library(patchwork)
library(tidyverse)

col_seg <- list('steelblue', 'hotpink4', 'turquoise4', 'chocolate3', 'palegreen4')
path <- 'Desktop/simulations/subclonal_w_tail/'
s <- 1

for (s in seq(0, 19)){
  print(s)
  snp <-  read.table(paste0(path,'/sim_', s, '/', s, '_snp_sim.tsv'), sep =  '\t', header = T) %>% 
    arrange(pos) %>% 
    mutate(id = paste0(CN_1, ' ', CN_2, ' ', ccf, ' ', purity))
  
  snv <-  read.table(paste0(path, '/sim_', s, '/', s, '_snv_sim.tsv'), sep =  '\t', header = T)  %>% 
    arrange(pos) %>% 
    mutate(id = paste0(CN_1, ' ', CN_2, ' ', ccf, ' ', purity, ' ', model))
  
  seg <- read.table(paste0(path, '/sim_', s, '/', s, '_segs.tsv'), sep =  '\t', header = T) 
 
  plot <- snp %>% 
          ggplot() +
          geom_point(aes(x = pos, y = baf, color = id), size = 0.2) +
          scale_color_manual(values = col_seg) + 
          theme_bw() + 
          theme(legend.position = 'none') +
          ylim(0,1) + 
    
          snp %>% 
            ggplot() +
            geom_point(aes(x = pos, y = dr, color = id), size = 0.2) +
            scale_color_manual(values = col_seg) + 
            theme_bw() + 
            theme(legend.position = 'none') +
            ylim(0,3) +
          
          snv %>% 
            ggplot() +
            geom_point(aes(x = pos, y = vaf, color = id), size = 0.4) +
            scale_color_manual(values = col_seg) + 
            theme_bw() + 
            theme(legend.position = 'none') +
            ylim(0,1) + 
          
          snv %>% 
            ggplot() +
            geom_histogram(aes(vaf, fill = id), alpha = 0.4, position = 'identity', binwidth = 0.015) +
            scale_fill_manual(values = col_seg) +
            theme_bw() + 
            xlim(0,1) +
  
          snv %>% 
            ggplot() +
            geom_histogram(aes(vaf, fill = id), alpha = 0.4, position = 'identity', binwidth = 0.015) +
            scale_fill_manual(values = col_seg) +
            theme_bw() + 
            xlim(0,1) +
            facet_wrap(id~., scales = 'free') +
            theme_bw() + 
            theme(legend.position = 'none') +
    plot_layout(design = 'ABC
                          DEE', guides = 'collect') 
  
  result <- paste0(path, 'sim_', s, '/')
  dir.create(result, recursive = TRUE)
  ggsave(filename = paste0(result, 'sim_', s, '_data.png'), 
         plot = plot, 
         dpi = 400, 
         width = 15, 
         height = 10, 
         units = 'in')
}

