col_seg  <- c('lightsteelblue', 'rosybrown2', 'wheat2', 'snow3','darkseagreen')

source('~/Documents/GitHub/utils_locate/simulations_rRACES/utils_plot.R')
plot_allelic_fragmentation <- function(events){
  data <- absolute_to_relative_coordinates(events %>% mutate(chr = paste0('chr', chr)))
  plt <- blank_genome(chromosomes = unique(data$chr)) +
    geom_segment(data = data, aes(x = begin, xend = end, y = major+0.03), col = '#01796F', linewidth = 1.5) +
    geom_segment(data = data, aes(x = begin, xend = end, y = minor-0.03), col = 'goldenrod', linewidth = 1.5) +
    ylab('CN')
  return(plt)
}


plot_hist_BAF <- function(germline, cut = 0.99){
  plt <- germline %>% 
    filter(VAF.tumour < cut)  %>% 
    ggplot() + 
    geom_histogram(aes(x = VAF.tumour, fill = cna_id.tumour), binwidth = 0.01) +
    xlab('BAF') +
    xlim(-0.01, 1.01) +
    facet_grid(cna_id.tumour~sample_name.tumour , scales = 'free') + 
    theme_bw() +  theme(legend.position = 'None') 
  return(plt)
}

plot_hist_VAF <- function(somatic, cut = 0.05){
  plt <- somatic %>% 
    filter(VAF > cut) %>% 
    ggplot() + 
    geom_histogram(aes(x = VAF, fill = cna_id), binwidth = 0.01) +
    xlim(-0.01, 1.01) +
    facet_grid(cna_id ~ sample_name, scales = 'free') + 
    theme_bw() +  theme(legend.position = 'None')   
  
  return(plt)
}


plot_gw <- function(germline, somatic, bp, info = ''){
  sample <- somatic %>% pull(sample_name) %>% unique()
  plt <-     ggplot() +
    geom_point(aes(x = germline$from.tumour, y = germline$DR, color = germline$cna_id.tumour), alpha = 0.3, size = 0.5) + 
    geom_vline(aes(xintercept = bp), color = 'gray') +
    scale_color_manual(values = col_seg) + 
    ylab('DR') + 
    xlab('pos') + 
    theme_bw() +  theme(legend.position = 'None') +
    
    ggplot() +
    geom_point(aes(x = germline$from.tumour, y = germline$VAF.tumour, color = germline$cna_id.tumour), alpha = 0.3, size = 0.5) + 
    geom_vline(aes(xintercept = bp), color = 'gray') +
    scale_color_manual(values = col_seg) + 
    xlab('pos') + 
    ylim(0,1) +
    ylab('BAF') + 
    theme_bw() +  theme(legend.position = 'None')  +
    
    
    ggplot() +
    geom_point(aes(x = somatic$from, y = somatic$VAF, color = somatic$cna_id), alpha = 0.8, size = 0.5) + 
    geom_vline(aes(xintercept = bp), color = 'gray') +
    scale_color_manual(values = col_seg) + 
    ylim(0,1) +
    xlab('pos') + 
    ylab('VAF') + 
    theme_bw() +  theme(legend.position = 'None') +
    plot_layout(nrow = 3) + plot_annotation(title = info) 
  return(plt)
}

plt_smooth_vaf_median <- function(vaf_smooth, bp = c(), info = ''){
  plt <- ggplot() +
    geom_point(aes(x = vaf_smooth$pos, y = vaf_smooth$median_dr, color = vaf_smooth$cna_id), size = 1) + 
    geom_vline(aes(xintercept = bp), color = 'gray') +
    scale_color_manual(values = col_seg) + 
    ylab('median DR') + 
    xlab('pos') + 
    theme_bw() +  theme(legend.position = 'None') +
    
    ggplot() +
    geom_point(aes(x = vaf_smooth$pos, y = vaf_smooth$median_baf, color = vaf_smooth$cna_id), size = 1) + 
    geom_vline(aes(xintercept = bp), color = 'gray') +
    scale_color_manual(values = col_seg) + 
    xlab('pos') + 
    ylim(0,1) +
    ylab('median BAF') + 
    theme_bw() +  theme(legend.position = 'None')  +
    
    
    ggplot() +
    geom_point(aes(x = vaf_smooth$pos, y = vaf_smooth$vaf, color = vaf_smooth$cna_id), size = 1) + 
    geom_vline(aes(xintercept = bp), color = 'gray') +
    scale_color_manual(values = col_seg) + 
    ylim(0,1) +
    xlab('pos') + 
    ylab('VAF') + 
    theme_bw() +  theme(legend.position = 'None') +
    plot_layout(nrow = 3) + plot_annotation(title = info) 
  
  return(plt)
}


plt_smooth_vaf_mean <- function(vaf_smooth, bp = c(), info = ''){
  plt <- ggplot() +
    geom_point(aes(x = vaf_smooth$pos, y = vaf_smooth$mean_dr, color = vaf_smooth$cna_id), size = 1) + 
    geom_vline(aes(xintercept = bp), color = 'gray') +
    scale_color_manual(values = col_seg) + 
    ylab('mean DR') + 
    xlab('pos') + 
    theme_bw() +  theme(legend.position = 'None') +
    
    ggplot() +
    geom_point(aes(x = vaf_smooth$pos, y = vaf_smooth$mean_baf, color = vaf_smooth$cna_id), size = 1) + 
    geom_vline(aes(xintercept = bp), color = 'gray') +
    scale_color_manual(values = col_seg) + 
    xlab('pos') + 
    ylim(0,1) +
    ylab('mean BAF') + 
    theme_bw() +  theme(legend.position = 'None')  +
    
    
    ggplot() +
    geom_point(aes(x = vaf_smooth$pos, y = vaf_smooth$vaf, color = vaf_smooth$cna_id), size = 1) + 
    geom_vline(aes(xintercept = bp), color = 'gray') +
    scale_color_manual(values = col_seg) + 
    ylim(0,1) +
    xlab('pos') + 
    ylab('VAF') + 
    theme_bw() +  theme(legend.position = 'None') +
    plot_layout(nrow = 3) + plot_annotation(title = info) 
  
  return(plt)
}
