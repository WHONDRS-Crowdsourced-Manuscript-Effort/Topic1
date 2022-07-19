library(tidyverse) # data wrangling
library(patchwork) # plotting with layouts

# Set directory to source file directory (for ease of use in Rstudio. remove if using from command line)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### FUNCTIONS
plot_vk <- function(cross_table, sed_or_water, threshold_method, pathways){
  #' plot Van-Krevelen diagram for Core-Satellite species
  #' cross_table: Input dataset.
  #' sed_or_water: Pick either 'sed' or 'water.
  #' threshold_method: Pick either 'emergent', 'pca' or 'rf'.
  #' pathways: Display direction of biogeochemical pathways or not.
  
  p <- cross_table %>%
    mutate(cs_flag = !!as.name(paste0('cs.flag.',threshold_method,'_', sed_or_water))) %>%
    filter(cs_flag %in% c('Core','Satellite')) %>%
    ggplot() +
    geom_point(data = . %>% filter(cs_flag == 'Satellite'),
               aes(x = OtoC_ratio, y = HtoC_ratio, col = cs_flag), size = 4, alpha = 0.45) +
    geom_point(data = . %>% filter(cs_flag == 'Core'),
               aes(x = OtoC_ratio, y = HtoC_ratio, col = cs_flag), size = 4, alpha = 0.45) +
    scale_color_discrete('CS Flag') +
    theme(legend.position = 'top') + xlab('O/C ratio') + ylab('H/C ratio') +
    ggtitle(paste0('VK plot for ', toupper(sed_or_water), ' with ', toupper(threshold_method), ' thresholds'))
  
  if(pathways == 1){
      p <- p + annotate("text", x = -0.025, y = 1.1, label = "Oxidation") +
      geom_segment(aes(x = 0, y = 2, xend = -0.05, yend = 2.1),
                   arrow = arrow(length = unit(0.25, "cm"))) +
      annotate("text", x = -0.025, y = 2.2, label = "Methylation") +
      geom_segment(aes(x = 0.5, y = 2.7, xend = 0.5, yend = 2.5),
                   arrow = arrow(length = unit(0.25, "cm"))) +
      annotate("text", x = 0.575, y = 2.6, label = "Hydrogen") +
      geom_segment(aes(x = -0.05, y = 0.05, xend = 0, yend = 0.15),
                   arrow = arrow(length = unit(0.25, "cm"))) +
      annotate("text", x = -0.025, y = 0.25, label = "Hydration") +
      geom_segment(aes(x = -0.05, y = 1, xend = 0, yend = 1),
                     arrow = arrow(length = unit(0.25, "cm"))) 
  }
  return(p)
}

plot_kmd <- function(cross_table, sed_or_water, threshold_method) {
  # ' plot Kendrick mass defect analysis for Core-Satellite species
  #' cross_table: Input dataset.
  #' sed_or_water: Pick either 'sed' or 'water.
  #' threshold_method: Pick either 'emergent', 'pca' or 'rf'.

  p <- cross_table %>%
    mutate(cs_flag = !!as.name(paste0('cs.flag.',threshold_method,'_', sed_or_water))) %>%
    filter(cs_flag %in% c('Core','Satellite')) %>%
    ggplot() +
    geom_point(data = . %>% filter(cs_flag == 'Satellite'),
               aes(x = kmass.CH2, y = kdefect.CH2, col = cs_flag), size = 3, alpha = 0.45) +
    geom_point(data = . %>% filter(cs_flag == 'Core'),
               aes(x = kmass.CH2, y = kdefect.CH2, col = cs_flag), size = 3, alpha = 0.45) +
    # geom_point(aes(x = kmass.CH2, y = kdefect.CH2, color = cs_flag, alpha = cs_flag)) +
    # scale_alpha_manual('CS Flag',values = c(1,0.15)) +
    scale_color_discrete('CS Flag') +
    theme(legend.position = 'top') + xlab('Kendrick mass (CH2)') + ylab('Kendrick mass defect') +
    ggtitle(paste0('KMD analysis for ', toupper(sed_or_water), ' with ', toupper(threshold_method), ' thresholds')) +
    geom_segment(aes(x = 400, y = 0.75, xend = 300, yend = 0.625),
                 arrow = arrow(length = unit(0.5, "cm"))) +
    annotate("text", x = 350, y = 0.775, label = "2H replaced by O") 
  return(p)
}

plot_dens <- function(cross_table, threshold_method, sed_or_water, axis_name){
  # ' plot densities to use along VK diagram for Core-Satellite species
  #' cross_table: Input dataset.
  #' sed_or_water: Pick either 'sed' or 'water.
  #' threshold_method: Pick either 'emergent', 'pca' or 'rf'.
  #' axis: pick 'HC' or 'OC' OR 'Kmass' or 'Kdefect' to use with kendrick mass analysis
  if(axis_name =='HC'){
    p <- cross_table %>%
      mutate(cs_flag = !!as.name(paste0('cs.flag.',threshold_method,'_', sed_or_water))) %>%
      filter(cs_flag %in% c('Core','Satellite')) %>%
      ggplot() +
      geom_density(aes(HtoC_ratio, fill = cs_flag), alpha = .5) +
      scale_fill_discrete('CS Flag') +
      theme(legend.position = 'none') +
      xlab('H/C ratio') + ylab('Density') +
      # ggtitle(paste0('Densities ', toupper(sed_or_water), ' with ', toupper(threshold_method), ' thresholds')) + 
      coord_flip()
  } else if(axis_name == 'OC'){
    p <- cross_table %>%
      mutate(cs_flag = !!as.name(paste0('cs.flag.',threshold_method,'_', sed_or_water))) %>%
      filter(cs_flag %in% c('Core','Satellite')) %>%
      ggplot() +
      geom_density(aes(OtoC_ratio, fill = cs_flag), alpha = .5) +
      scale_fill_discrete('CS Flag') +
      theme(legend.position = 'none') +
      xlab('O/C ratio') + ylab('Density') 
      # ggtitle(paste0('Densities ', toupper(sed_or_water), ' with ', toupper(threshold_method), ' thresholds')) 
  } else if(axis_name == 'Kmass'){
    p <- cross_table %>%
      mutate(cs_flag = !!as.name(paste0('cs.flag.',threshold_method,'_', sed_or_water))) %>%
      filter(cs_flag %in% c('Core','Satellite')) %>%
      ggplot() +
      geom_density(aes(kmass.CH2, fill = cs_flag), alpha = .5) +
      scale_fill_discrete('CS Flag') +
      theme(legend.position = 'none') +
      xlab('Kendrick mass') + ylab('Density') 
    # ggtitle(paste0('Densities ', toupper(sed_or_water), ' with ', toupper(threshold_method), ' thresholds')) 
  } else if(axis_name == 'Kdefect'){
    p <- cross_table %>%
      mutate(cs_flag = !!as.name(paste0('cs.flag.',threshold_method,'_', sed_or_water))) %>%
      filter(cs_flag %in% c('Core','Satellite')) %>%
      ggplot() +
      geom_density(aes(kdefect.CH2, fill = cs_flag), alpha = .5) +
      scale_fill_discrete('CS Flag') +
      theme(legend.position = 'none') +
      xlab('Kendrick defect') + ylab('Density') +
      coord_flip()
    # ggtitle(paste0('Densities ', toupper(sed_or_water), ' with ', toupper(threshold_method), ' thresholds')) 
  }
  return(p)
}

plot_with_layout <- function(p1, p2, p3){
  #' plot scatter and density plots together
  p1 + plot_spacer() + p2 + p3 + 
    plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
}

### PLOTTING
# read in the data
ct1 <- read.csv('./Data/FTICR_crosstable_rep.merged1_all_em.thres_2022-03-18.csv') # extract molecular information
ct2 <- read.csv('./Data/FTICR_crosstable_rep.merged2_all_em.thres_2022-03-18.csv') # extract molecular information

# VK
pdf('vk_diagram_emergent_merged2.pdf')
sample_type = 'water'
p1 <- plot_vk(cross_table = ct2, sed_or_water = sample_type, threshold_method = 'emergent', pathways = 1) 
p1d1 <- plot_dens(cross_table = ct2, threshold_method = 'emergent', sed_or_water = sample_type, axis_name = 'OC')
p1d2 <- plot_dens(cross_table = ct2, threshold_method = 'emergent', sed_or_water = sample_type, axis_name = 'HC')

plot_with_layout(p1d1, p1, p1d2)

# KMD
p2 <- plot_kmd(cross_table = ct2, sed_or_water = sample_type, threshold_method = 'emergent')
plot(p2)

dev.off()


  