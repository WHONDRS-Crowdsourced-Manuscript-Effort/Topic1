##############################################
# WHONDRS CROWDSOURCED TOPIC 1: CORE SATELLITE 

# This script was written to determine the best threshold value to divide core and satellite species

# Approach is to find threshold value that maximizes explanation of variance by core satellite variable
# (css.flag) in unsupervised clustering methods (i.e. PCA and unsupervised random forest)

# Requires csv files:
  # molecular properties of each molecular formula
  # presence absence data of each molecular formula

# Author: Kadir Bice

# Load necessary packages
library(tidyverse)    # data wrangling
library(cluster)      # clustering algorithms
library(factoextra)   # clustering algorithms & visualization
library(ggfortify)    # for PCA and visualization
library(randomForest) # random forest

# ################## #
#### Load in data ####
# ################## #

# Set directory to source file directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# get the list of input files
input_files <- 
  list.files(pattern = '.csv', path = c('./Data/rep.merged1'), full.names = T) %>%
  matrix(ncol = 2) %>%
  `colnames<-`(c('commat', 'crosstable')) %>%
  data.frame() %>%
  bind_rows(matrix(data = list.files(pattern = '.csv', path = c('./Data/rep.merged2'), full.names = T), ncol = 2) %>%
              `colnames<-`(c('commat', 'crosstable')) %>% data.frame()) %>%
  mutate(merge = substr(crosstable, start = 12, stop = 18), rarity = rep(c('all', 'rar1', 'rar2'),2))

# loop over input files with different merging and rarity cutoffs
for(file_index in 2:nrow(input_files)){
  # Load files
  data <- as.data.frame(t(read.csv(input_files$commat[file_index], row.names = 1 ))) # extract presence information 
  mol <- read.csv(input_files$crosstable[file_index]) # extract molecular information
  mol <- mol[order(mol$MolForm),]
  rownames(mol) <- mol$MolForm
  
  # checking mol and data files
  if(!identical(row.names(data), row.names(mol))){
    stop("data file and mol file do not match")
  }
  
  # ################## #
  #### Preprocessing ####
  # ################## #
  
  # Setting data to presence/absence (already done in Masumi's new data)
  data[data > 0] =  1 # Given that intensities are not necessarily related to concentration, we set values to presence/absence
  
  # loop over type of the sample (SED or SW)
  for(sample_type in c('SED', 'SW')) {
    # filter Sediment OR Surface Water samples
    data_filtered <- data %>%
      select(contains(sample_type))
    
    # Calculate frequencies of each species
    freq_css <- data.frame(frequencies = (rowSums(data_filtered)/ncol(data_filtered))*100) # Finding frequencies (number of occurences in samples)
    
    # # just to see if frequencies looks alright
    # freq_css %>%
    #   count(frequencies) %>%
    #   ggplot() +
    #   geom_histogram(aes(x = frequencies))
    
    # pick molecular properties to include in the test
    pick_mol <- as.data.frame(apply(mol[,c('NeutralMass', 'Error_ppm', 'Candidates',
                                           'AI', 'AI_Mod', 'DBE', 'DBE_O', 'DBE_AI', 'GFE', 'kmass.CH2',
                                           'kdefect.CH2', 'NOSC', 'OtoC_ratio', 'HtoC_ratio', 'NtoC_ratio',
                                           'PtoC_ratio', 'delGcoxPerCmol', 'delGd0', 'delGd',
                                           'n.mf')], 2, as.numeric))
    
    # filter frequencies = 0 and scale
    num_mol <- pick_mol %>%
      bind_cols(freq_css) %>%
      filter(frequencies > 0) %>%
      mutate_at(vars(-frequencies), scale)
    
    ## output for python implementation
    # write.csv(format(num_mol, digits = 17), paste0('nummol_all_', sample_type, '_', input_files$merge[file_index], '.csv'))
    # readnmol = read.csv('nummol.csv')
    ##
    
    # ################## #
    #### Loop over thresholds ####
    # ################## #
    
    # initialization of storage
    freq.mol_css_store <- list() # list to store each matrix
    pca_results_store <- c() # list to store each PCA
    rf_results_store <- c() # list to store each RF
    counts_store <- c() # list to store counts per each class
    css_threshold <- c(seq(from = 2, to = 50, by = 1)) # series of threshold values
    
    # loop over thresholds
    for(i in 1:length(css_threshold)){ 
      # determine Core vs Satellite compounds
      freq.mol_css <- num_mol %>%
        mutate(.before = NeutralMass, css.flag = factor(ifelse(frequencies >= css_threshold[i], 0, 1))) %>% # if a compound has frequency higher than threshold make it Core (0), otherwise Satellite (1)
        mutate(css.flag = as.numeric(as.character(css.flag))) %>%
        select(-frequencies)
      
      # store counts
      counts_store <- rbind(counts_store,freq.mol_css %>%
                              count(css.flag) %>%
                              mutate(threshold = css_threshold[i]))
      
      # store matrix
      freq.mol_css_store[[i]] <- freq.mol_css
      
      ############ 
      # # plot CDF
      # cdf_ai = ecdf(freq.mol_css$GFE) # get empirical CDF on a variable
      # 
      # freq.mol_css %>%
      #   mutate(css.flag = as.factor(css.flag)) %>%
      #   cbind(cdf = cdf_ai(freq.mol_css$GFE)) %>%
      #   # distinct(OtoC_ratio, HtoC_ratio, .keep_all = T) %>%
      #   # distinct(MolForm, .keep_all = T) %>%
      #   ggplot() +
      #   geom_point(aes(x = GFE, y = cdf, color = css.flag, alpha = css.flag, size = css.flag), stat = 'unique') + #, shape = time_point)) +
      #   # scale_color_manual(values=c(Core = "#333BFF", Satellite = "#CC6600")) +
      #   scale_alpha_manual(values=c('0' = 1, '1' = 0.25)) +
      #   scale_size_manual(values=c('0' = 2, '1' = 1))
      
      # # # kmeans cluster and plot
      # cluster_kmeans <- kmeans(freq.mol_css, centers = 2)
      # fviz_cluster(cluster_kmeans, data = freq.mol_css, geom = 'point' )
      # 
      # # get_dist(freq.mol_css, method = 'pearson')
      
      # # # Function to compute total within-cluster sum of square (to test what number of clusters is best)
      # wssplot <- function(data, nc, seed=1234){
      #   wss <- (nrow(data)-1)*sum(apply(data,2,var))
      #   for (j in 2:nc){
      #     set.seed(seed)
      #     wss[j] <- sum(kmeans(data, centers=j)$withinss)}
      #   plot(1:nc, wss, type="b", xlab="Number of Clusters",
      #        ylab="Within groups sum of squares")
      # }
      # # plotting values for each cluster starting from 1 to 15
      # wssplot(freq.mol_css, nc = 15)
      # fviz_nbclust(freq.mol_css, kmeans, method = "silhouette")
      
      ##############################
      ###### run and plot PCA ######
      ##############################
      
      # tried ONE-HOT encoding instead of integer encoding but does not give that much of meaning
      # css.flag2 = as.numeric(freq.mol_css$css.flag == 0)
    
      pca_1 <- prcomp(freq.mol_css)
      autoplot(pca_1, data = freq.mol_css, colour = 'css.flag', loadings = T, loadings.label = T)

      # Scree plot to see percentage variance explained by each principal component
      fviz_eig(pca_1, addlabels = TRUE, ylim = c(0, 50))

      # Individual contributions of variables for each principal component
      get_pca_var(pca_1)$contrib

      # plot variable contributions
      fviz_pca_var(pca_1,
                   col.var = "contrib", # Color by contributions to the PC
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   repel = TRUE   # Avoid text overlapping
      ) + labs(title = paste("PCA variables for threshold = ", css_threshold[i]))
      # ggsave(paste('threshold', css_threshold[i],'_PCA.png'),bg = '#FFFFFF')

      # get variance explained of each component
      eigs <- pca_1$sdev^2
      # calculate total (variance) weighted contribution of css flag for first 3 principal components
      cont_pca_css <- sum(get_pca_var(pca_1)$contrib[1,1:3] %*% (eigs[1:3] / sum(eigs)))

      print(paste('contributions of css flag to first 3 dimensions with threshold = ', css_threshold[i]))
      print(cont_pca_css)

      # store pca contribution
      pca_results_store <- rbind(pca_results_store, cbind(contribution = cont_pca_css, threshold = css_threshold[i]))
      
      ##############################
      ###### run unsupervised RF ###
      ##############################
      # (gives RAM error in the large dataset)
      rf1 <- randomForest(freq.mol_css)

      df_rf1 <- data.frame(vars = rownames(rf1$importance), importance = as.vector(rf1$importance))

      print(paste('contributions of css flag in random forest with threshold = ', css_threshold[i]))
      cont_rf_css <- df_rf1 %>%
        filter(df_rf1 == 'css.flag') %>%
        select(importance) %>%
        as.numeric()
      print(cont_rf_css)
      rf_results_store <- rbind(rf_results_store, cbind(contribution = cont_rf_css, threshold = css_threshold[i]))

      # ggplot(df_rf1, aes(y = vars, x = importance)) + geom_bar(stat="identity", fill = 4)
      rm(rf1)
      ######
      
    }
    ##################
    # plotting results
    ##################
    
    # pca_results_store %>%
    #   ggplot() + geom_line(aes(x = threshold, y = contribution), size = 1.2) +
    #   geom_vline(xintercept = pca_results_store[which.max(pca_results_store[,1]),2], col = 2) +
    #   annotate(geom = 'text', x = 20, y = 0, label = paste('max at threshold = ', pca_results_store[which.max(pca_results_store[,1]),2]), col = 2) +
    #   ggtitle(paste0(input_files$merge[file_index], '_', input_files$rarity[file_index],'_pca_',sample_type))
    # ggsave(paste0(input_files$merge[file_index], '_', input_files$rarity[file_index],'_pca_',sample_type, '.png'))
    rf_results_store %>%
      ggplot() + geom_line(aes(x = threshold, y = contribution), size = 1.2) +
      geom_vline(xintercept = rf_results_store[which.max(rf_results_store[,1]),2], col = 2) +
      annotate(geom = 'text', x = 20, y = 0, label = paste('max at threshold = ', rf_results_store[which.max(rf_results_store[,1]),2]), col = 2) +
      ggtitle(paste0(input_files$merge[file_index], '_', input_files$rarity[file_index],'_rf_',sample_type))
    ggsave(paste0(input_files$merge[file_index], '_', input_files$rarity[file_index],'_rf_',sample_type, '.png'))
    # counts_store %>%
    #   ggplot() + geom_line(aes(x = threshold, y = n, col = css.flag), size = 1.2)
    
    # store if needed
    #write.csv(rf_results_store, file = paste0(input_files$merge[file_index], '_', input_files$rarity[file_index],'_rf_',sample_type, '.csv'))
  }
}