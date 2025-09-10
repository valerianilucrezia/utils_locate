col_seg  <- list('1:1'='indianred3', 
                 '2:1' = 'deepskyblue4', 
                 '2:0'='olivedrab',
                 '1:0'='sienna2',
                 '2:2'='goldenrod', 
                 '3:1'='palevioletred3',
                 '3:0'='plum4',
                 '3:2'='darkslategray4',
                 '3:3'='tan3')

source('../utils_plot.R')
plot_allelic_fragmentation <- function(events, name = ''){
  data <- absolute_to_relative_coordinates(events %>% mutate(chr = paste0('chr', chr)))
  data <- data %>% mutate(CN = paste(major, minor, sep = ':'))
  
  plt <- blank_genome(chromosomes = unique(data$chr)) +
    geom_rect(data = data, aes(xmin =begin, xmax = end, ymin = -Inf, ymax = Inf, fill = CN), alpha = .4) +
    geom_segment(data = data, aes(x = begin, xend = end, y = major+0.03), col = 'deepskyblue4', linewidth = 1.5) +
    geom_segment(data = data, aes(x = begin, xend = end, y = minor-0.03), col = 'indianred3', linewidth = 1.5) +
    ggtitle(name) + 
    scale_fill_manual(values = col_seg) + 
    ylab('CN')
  return(plt)
}

plot_allelic_fragmentation_subclonal <- function(events, name = ''){
  data <- absolute_to_relative_coordinates(events %>% mutate(chr = paste0('chr', chr)))
  plt <- blank_genome(chromosomes = unique(data$chr)) +
    geom_segment(data = data, aes(x = begin, xend = end, y = major+0.03), col = '#01796F', linewidth = 1.5) +
    geom_segment(data = data, aes(x = begin, xend = end, y = minor-0.03), col = 'goldenrod', linewidth = 1.5) +
    ylab('CN') +
    ggtitle(name) + 
    facet_wrap(type~., nrow = 3)
  return(plt)
}



plot_hist_BAF <- function(germline, cut = 0.99){
  plt <- germline %>% 
    filter(VAF.tumour < cut)  %>% 
    ggplot() + 
    geom_histogram(aes(x = VAF.tumour, fill = CN.tumour), binwidth = 0.01) +
    scale_fill_manual('CN', values = col_seg) + 
    xlab('BAF') +
    xlim(-0.01, 1.01) +
    facet_grid(CN.tumour~sample_name.tumour , scales = 'free') + 
    theme_bw()
  return(plt)
}

plot_hist_VAF <- function(somatic, cut = 0){
  plt <- somatic %>% 
    filter(VAF >= cut) %>% 
    ggplot() + 
    geom_histogram(aes(x = VAF, fill = CN), binwidth = 0.01) +
    scale_fill_manual('CN', values = col_seg) + 
    xlim(-0.01, 1.01) +
    facet_grid(CN ~ sample_name, scales = 'free') + 
    theme_bw() 
  
  return(plt)
}


plot_gw <- function(germline, somatic, bp, info = ''){
  sample <- somatic %>% pull(sample_name) %>% unique()
  plt <- ggplot() +
    geom_point(aes(x = germline$from.tumour, y = germline$DR, color = germline$CN.tumour), alpha = 0.5, size = 0.1) + 
    geom_vline(aes(xintercept = bp), color = 'gray') +
    scale_color_manual('CN', values = col_seg) + 
    ylab('DR') + 
    xlab('pos') + 
    theme_bw() +
    scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) + 
    ylim(0,3) +
    
    ggplot() +
    geom_point(aes(x = germline$from.tumour, y = germline$VAF.tumour, color = germline$CN.tumour), alpha = 0.5, size = 0.1) + 
    geom_vline(aes(xintercept = bp), color = 'gray') +
    scale_color_manual('CN', values = col_seg) + 
    xlab('pos') + 
    ylim(0,1) +
    ylab('BAF') + 
    scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) + 
    theme_bw() +
    
    
    ggplot() +
    geom_point(aes(x = somatic$from, y = somatic$VAF, color = somatic$CN), alpha = 0.5, size = 0.1) + 
    geom_vline(aes(xintercept = bp), color = 'gray') +
    scale_color_manual('CN', values = col_seg) + 
    ylim(0,1) +
    xlab('pos') + 
    ylab('VAF') + 
    theme_bw()  +
    scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) + 
    plot_layout(nrow = 3, guides = 'collect') + plot_annotation(title = info) & theme(legend.position = 'none')
  return(plt)
}

plt_smooth_vaf_median <- function(vaf_smooth, bp = c(), info = ''){
  plt <- ggplot() +
    geom_point(aes(x = vaf_smooth$pos, y = vaf_smooth$median_dr, color = vaf_smooth$cna_id), size = .7) + 
    geom_vline(aes(xintercept = bp), color = 'gray') +
    scale_color_manual(values = col_seg) + 
    ylab('median DR') + 
    xlab('pos') + 
    theme_bw() +  theme(legend.position = 'None') +
    
    ggplot() +
    geom_point(aes(x = vaf_smooth$pos, y = vaf_smooth$median_baf, color = vaf_smooth$cna_id), size = .7) + 
    geom_vline(aes(xintercept = bp), color = 'gray') +
    scale_color_manual(values = col_seg) + 
    xlab('pos') + 
    ylim(0,1) +
    ylab('median BAF') + 
    theme_bw() +  theme(legend.position = 'None')  +
    
    
    ggplot() +
    geom_point(aes(x = vaf_smooth$pos, y = vaf_smooth$vaf, color = vaf_smooth$cna_id), size = .7) + 
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
