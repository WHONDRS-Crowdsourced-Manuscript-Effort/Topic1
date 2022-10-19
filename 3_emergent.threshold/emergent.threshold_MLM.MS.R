# Apply emergent threshold

# Original script from Michaela de Melo (in emthres-exploration_MLM.R)
# Re-organized for parallel processing by Masumi Stadler

#--------------------------------------------------------------------------------

# Data overview:
# 12 datasets
# 6 for molecular formulae, 6 for peaks
# Within each:
# 2 replicate merging thresholds = rep.merged1, rep.merged2
# 3 rarity cutoffs = none, rar1, rar2

# Number of peaks/MF retained:
# rep.merged1 > rep.merged2
# none > rar1 > rar2
# (Summary in Github repo or in 1_data.cleaning/2_rarity_cutoff_MS.R)

# Read in data through plyr to allow parallel processing:
# Original script from ll. 237-onwards in emthres-exploration_MLM.R

### Packages -------------------------------------------------------------------------------
pckgs <- list("plyr","tidyverse","data.table","doMC")

### Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
# install missing packages with following line. Please uncomment:
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)

### Load packages
invisible(lapply(pckgs, library, character.only = T))
rm(pckgs, miss.pckgs)

### Parallel set-up -----------------------------------------------------------------------
# register number of cores for parallel computing with 'apply' family
detectCores() # 24 (private PC: 4), we do not have 24 cores. It's 12 cores, 24 threads.
registerDoMC()
getDoParWorkers() # 12 (private PC: 2)

# prepare for parallel processing (with 'parallel') for 'vegan' functions
numCores <- detectCores()
cl <- makeCluster(numCores, type = "FORK") # using forking

# Set seed for session and reproducibility of permutations
# (NMDS were done, but not part of main workflow anymore)
set.seed(3)

# Read in molecular formulae -----------------------------------------------------------------
# Read in community matrices
mf.files <- list.files(path = './1_data.cleaning/output/', pattern = "commat_rep.merged")
#only extract new ones
mf.files <- mf.files[str_detect(mf.files, pattern = "2022-01")]
length(mf.files) # should be 6, yes

# Read in cross table
mf.cross <- list.files(path = './1_data.cleaning/output/', pattern = "crosstable_rep.merged")
length(mf.cross)

# Check if order corresponds the same sequence
data.frame(mf.files, mf.cross)

# Read in meta data
meta <- read.csv("./1_data.cleaning/output/FTICR_meta_eachriver_2022-01-19.csv",
                 sep = ',', stringsAsFactors = F)

mf.ls <- vector("list", 6)
# Run loop to read in data
for(i in 1:6){
  data_commat <- read.csv(paste0("./1_data.cleaning/output/", mf.files[i]),
           sep = ",", stringsAsFactors = F)
  # data_cross <- read.csv(paste0("./1_data.cleaning/output/", mf.cross[i]),
  #                        sep = ",", stringsAsFactors = F)
  
  # assign to list bin
  mf.ls[[i]] <- data_commat
}

# Some sanity check
lapply(mf.ls, dim)
# Add names to list
names(mf.ls) <- c("mf_rep.merged1_all", "mf_rep.merged1_rar1", "mf_rep.merged1_rar2",
                  "mf_rep.merged2_all", "mf_rep.merged2_rar1", "mf_rep.merged2_rar2")

# Read in peaks -----------------------------------------------------------------
# Read in community matrices
mz.files <- list.files(path = './1_data.cleaning/output/peaks/', pattern = "commat_rep.merged")
length(mz.files) # should be 6, yes

# There is no cross table for peaks, as there are some peaks with no molecular formula

mz.ls <- vector("list", 6)
# Run loop to read in data
for(i in 1:6){
  data_commat <- read.csv(paste0("./1_data.cleaning/output/peaks/", mz.files[i]),
                          sep = ",", stringsAsFactors = F)
  
  # assign to list bin
  mz.ls[[i]] <- data_commat
}

# Some sanity check
lapply(mz.ls, dim)
# Add names to list
names(mz.ls) <- c("mz_rep.merged1_all", "mz_rep.merged1_rar1", "mz_rep.merged1_rar2",
                  "mz_rep.merged2_all", "mz_rep.merged2_rar1", "mz_rep.merged2_rar2")

# Run re-formatting on all data ----------------------------------------------------------
all.ls <- c(mf.ls, mz.ls)

proc.ls <- llply(all.ls, function(x){
  setDT(x)
  # Separate sediment and surface water by ID and remove ID column
  sed <- x[str_detect(ID, "SED"),][,-1]
  water <- x[str_detect(ID, "SW"),][,-1]
  
  # Process surface water
  
  # how often were MF found across water samples?
  sum_water <- data.frame(colSums(water)) # sum of presence-absence = number of sites present
  # Make into percentage
  df_water <- tibble(table(sum_water)) # transform to counts of MF that were found in 1, 2, 3 sites and so on...
  #table(sum_water) #for visualization of results
  df_water <- df_water[-1,] #remove First cell, which is number of formulae = 0. Means we found 9,575 formulae were absent in water samples (compared to the merged table with sediment samples)
  df_water$occupancy <- seq(1:nrow(df_water)) #95 is the number of water samples
  colnames(df_water)[1] <- "frequency"
  
  rm(sum_water)
  
  # Process sediment
  # Sequence of code is same as water samples.
  sum_sed <- data.frame(colSums(sed))
  df_sed <- tibble(table(sum_sed))
  #table(sum_sed) #for visualization of results
  df_sed <- df_sed[-1,] #remove First cell, which is number of formulae = 0. Means we found 9,575 formulae were absent in water samples (compared to the merged table with sediment samples)
  df_sed$occupancy <- seq(1:nrow(df_sed)) #93 is the number of water samples
  colnames(df_sed)[1] <- "frequency"
  
  return(list(df_water, df_sed))
  
  }, .parallel = T)

# Plot frequency - occupancy curves -----------------------------------------------------------
for(i in 1:length(proc.ls)){
  
  p <- ggplot(proc.ls[[i]][[1]], 
              aes(x=occupancy, y=log(frequency))) + 
    geom_point() + 
    labs(x="Occupancy", y = "Frequency (log)",title=paste("Water", names(proc.ls)[i], sep = " "))
  ggsave(paste0("3_emergent.threshold/prelim_figures/FreqXoccupancy_water_", 
                names(proc.ls)[i],".png"),
         p, dpi=250, height = 10, width = 10, unit = "cm")
  
  p <- ggplot(proc.ls[[i]][[2]], 
              aes(x=occupancy, y=log(frequency))) + 
    geom_point() + 
    labs(x="Occupancy", y = "Frequency (log)",title=paste("Sediment", names(proc.ls)[i], sep = " "))
  ggsave(paste0("3_emergent.threshold/prelim_figures/FreqXoccupancy_sed_", 
                names(proc.ls)[i],".png"),
         p, dpi=250, height = 10, width = 10, unit = "cm")
}

# Remove unnecessary objects
# rm(sed, water, sum_sed, sum_water, water, x, p, data_commat, df_sed, df_water)

# Emergent threshold --------------------------------------------------------------------------------
#-- Start: Addition by MS
# Find points of maximum deceleration to identify emergent thresholds -----------------------
library(ggpubr)

# Define functions to find points of maximum acceleration/deceleration ----------------------
# Functions published in:
# Stadler M, del Giorgio PA. ISME J 2021. 
# Originally modified from Tommy's answer at: https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima

# Identify positions of local maxima
# this function is very sensitive some quality control has to be done afterwards
localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  # Use .Machine$integer.max instead of Inf for integer
  y <- diff(c(-Inf, x)) > 0L
  # returns logical vector with which numbers in 1st derivative are positive
  rle(y)$lengths
  # returns how often the same value (i.e. TRUE) repeats in a sequence without being interrupted
  # and when there is a switch from T to F
  y <- cumsum(rle(y)$lengths)
  # returns vector with cumulating the sum of TRUE and FALSE observations = returns location of switch from T to F
  y <- y[seq.int(1L, length(y), 2L)]
  # samples every second location (i.e. switch to TRUE)
  
  y
  # return locations of local maxima
}

localMinima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  # Use .Machine$integer.max instead of Inf for integer
  y <- diff(c(Inf, x)) < 0L
  # returns logical vector with which numbers in 1st derivative are positive
  rle(y)$lengths
  # returns how often the same value (i.e. TRUE) repeats in a sequence without being interrupted
  # and when there is a switch from T to F
  y <- cumsum(rle(y)$lengths)
  # returns vector with cumulating the sum of TRUE and FALSE observations = returns location of switch from T to F
  y <- y[seq.int(1L, length(y), 2L)]
  # samples every second location (i.e. switch to TRUE)
  
  y
  # return locations of local minima
}

# Merge all data sets into a list for parallel processing --------------------
# Flatten list into a one-level list
flattened.ls <- purrr::flatten(proc.ls)
# add names
names(flattened.ls) <- c("mf_rep.merged1_all_water", "mf_rep.merged1_all_sed",
                         "mf_rep.merged1_rar1_water", "mf_rep.merged1_rar1_sed",
                         "mf_rep.merged1_rar2_water", "mf_rep.merged1_rar2_sed",
                         "mf_rep.merged2_all_water", "mf_rep.merged2_all_sed",
                         "mf_rep.merged2_rar1_water","mf_rep.merged2_rar1_sed", 
                         "mf_rep.merged2_rar2_water","mf_rep.merged2_rar2_sed",
                         "mz_rep.merged1_all_water","mz_rep.merged1_all_sed",
                         "mz_rep.merged1_rar1_water","mz_rep.merged1_rar1_sed",
                         "mz_rep.merged1_rar2_water","mz_rep.merged1_rar2_sed",
                         "mz_rep.merged2_all_water","mz_rep.merged2_all_sed",
                         "mz_rep.merged2_rar1_water","mz_rep.merged2_rar1_sed",
                         "mz_rep.merged2_rar2_water","mz_rep.merged2_rar2_sed")

# Example data set to write function
#x<- flattened.ls[[1]]

# Add percentage instead of occupancy to keep consistent with other thresholds
flattened.ls <- lapply(flattened.ls, function(x){
  x$percentage <- round(x$occupancy * 100 / nrow(x), 0)
  return(x)
})

# Extract only MF and MZ rep.merged1 and rep.merged2 - all for Christof to compare the effect
# of smoothing

export <- names(flattened.ls)[str_detect(names(flattened.ls), "all")]
export.ls <- flattened.ls[names(flattened.ls) %in% export]

export.df <- ldply(export.ls, bind_rows)
export.df <- export.df %>% separate(.id, into = c("dataset","rep.merged","rarity.cutoff","sample.type"),
                       remove = F, sep = "_")
export.df$dataset <- factor(export.df$dataset, levels = c("mf","mz"), labels = c("molecular.formulae",
                                                                                 "peaks"))

# write.table(export.df, "./3_emergent.threshold/output/frequency_occupancy_mfmz.csv",
#             sep = ",", dec = ".", row.names = F)

# Some sanity check
# Are percentage values uniquely rounded?
# lapply(flattened.ls, function(x){
#   any(duplicated(x$percentage))
# }) # all FALSE --> yes.

# We're using smooth spline here instead of defining a function and getting a derivative from the function
# as described in https://cran.r-project.org/web/packages/inflection/vignettes/inflectionMissionImpossible.html
# as it is hard to come up with a function describing the observed pattern
# not exponential, not logarithmic etc

em.list <- llply(flattened.ls, function(x){
  # take log
  x$log.freq <- log(x$frequency)
  # make a smooth curve
  spl <- smooth.spline(x$percentage, x$log.freq, spar = 0.5)
  # predict to get fit
  pred <- predict(spl)
  # get second derivative
  sec <- predict(spl, deriv = 2)
  
  # get number of peaks to classify core
  min.tail.peaks <- c(length(localMinima(sec$y)), length(localMinima(sec$y)) -1)
  max.tail.peaks <- c(length(localMaxima(sec$y)), length(localMaxima(sec$y)) -1)
  
  options(scipen = 999) # avoid scientific annotations
  
  # Core classification
  # Sediments have a max peak after min peak in contrast to
  # Surface water which has a min peak after max peak
  # Classifciation is different between these two sample types
  if(localMaxima(sec$y)[max.tail.peaks[1]] > localMinima(sec$y)[min.tail.peaks[1]]){
    # get one plot with raw numbers (non-log)
    raw <- ggplot() +
      theme_pubr() +
      annotate(xmax = x$percentage[localMinima(sec$y)[2]],
               xmin = min(x$percentage),
               ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "tomato") +
      annotate(xmax = max(x$percentage),
               xmin = x$percentage[localMinima(sec$y)[2]],
               ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "aquamarine") +
      annotate(xmax = x$percentage[localMaxima(sec$y)[max.tail.peaks[2]]],
               xmin = max(x$percentage),
               ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "tomato") +
      geom_line(aes(x = x$percentage, y = x$frequency)) +
      geom_point(aes(x = x$percentage[localMaxima(sec$y)[1:2]],
                     y = x$frequency[localMaxima(sec$y)[1:2]]), colour = "tomato") +
      geom_point(aes(x = x$percentage[localMaxima(sec$y)[max.tail.peaks[2]:max.tail.peaks[1]]],
                     y = x$frequency[localMaxima(sec$y)[max.tail.peaks[2]:max.tail.peaks[1]]]),
                 colour = "tomato") +
      geom_point(aes(x = x$percentage[localMinima(sec$y)[1:2]],
                     y = x$frequency[localMinima(sec$y)[1:2]]), colour = "royalblue") +
      geom_point(aes(x = x$percentage[localMinima(sec$y)[min.tail.peaks[2]:min.tail.peaks[1]]],
                     y = x$frequency[localMinima(sec$y)[min.tail.peaks[2]:min.tail.peaks[1]]]),
                 colour = "royalblue") +
      labs(x = "", y = "Frequency") +
      theme(axis.title = element_text(size = 12))
    
    # get another plot with data used to identify the points of maximum acceleration and minimum deceleration
    logged <- ggplot() +
      theme_pubr() +
      annotate(xmax = pred$x[localMinima(sec$y)[2]],
               xmin = min(x$percentage),
               ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "tomato") +
      annotate(xmax = max(pred$x),
               xmin = pred$x[localMinima(sec$y)[2]],
               ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "aquamarine") +
      annotate(xmax = pred$x[localMaxima(sec$y)[max.tail.peaks[2]]],
               xmin = max(x$percentage),
               ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "tomato") +
      geom_line(aes(x = pred$x, y = pred$y)) +
      geom_point(aes(x = x$percentage, y = x$log.freq), alpha = 0.5, colour = "black") +
      geom_point(aes(x = pred$x[localMaxima(sec$y)[1:2]],
                     y = pred$y[localMaxima(sec$y)[1:2]]), colour = "tomato") +
      geom_point(aes(x = pred$x[localMinima(sec$y)[1:2]],
                     y = pred$y[localMinima(sec$y)[1:2]]), colour = "royalblue") +
      geom_point(aes(x = pred$x[localMaxima(sec$y)[max.tail.peaks[2]:max.tail.peaks[1]]],
                     y = pred$y[localMaxima(sec$y)[max.tail.peaks[2]:max.tail.peaks[1]]]),
                 colour = "tomato") +
      geom_point(aes(x = pred$x[localMinima(sec$y)[min.tail.peaks[2]:min.tail.peaks[1]]],
                     y = pred$y[localMinima(sec$y)[min.tail.peaks[2]:min.tail.peaks[1]]]),
                 colour = "royalblue") +
      labs(x = "", y = "Frequency (log)") +
      theme(axis.title = element_text(size = 12))
    
    # Show second derivative used to find the points
    deriv <- ggplot() +
      theme_pubr() +
      geom_line(aes(x = sec$x, y = sec$y * 1000)) +
      geom_point(aes(x = sec$x[localMaxima(sec$y)[1:2]],
                     y = sec$y[localMaxima(sec$y)[1:2]]* 1000), colour = "tomato") +
      geom_point(aes(x = sec$x[localMinima(sec$y)[1:2]],
                     y = sec$y[localMinima(sec$y)[1:2]]* 1000), colour = "royalblue") +
      geom_point(aes(x = sec$x[localMaxima(sec$y)[max.tail.peaks[2]:max.tail.peaks[1]]],
                     y = sec$y[localMaxima(sec$y)[max.tail.peaks[2]:max.tail.peaks[1]]]* 1000),
                 colour = "tomato") +
      geom_point(aes(x = sec$x[localMinima(sec$y)[min.tail.peaks[2]:min.tail.peaks[1]]],
                     y = sec$y[localMinima(sec$y)[min.tail.peaks[2]:min.tail.peaks[1]]]* 1000),
                 colour = "royalblue") +
      labs(x = "", y = expression(paste("Acceleration x10"^3, " (2"^"nd", " derivative)"))) +
      theme(axis.title = element_text(size = 12))
    
    p <- ggarrange(raw, logged, deriv, ncol = 3, labels = c("d","e","f"))
    # add x axis title to be in the middle of two panels
    
    # Extract identified threshold
    thres.df <- data.frame(cs.flag = c("Satellite", "In-between", "Core"),
                           occup.thres.min = c(min(x$percentage), 
                                               localMinima(sec$y)[2]+1,
                                               localMaxima(sec$y)[max.tail.peaks[2]]),
                           occup.thres.max = c(localMinima(sec$y)[2],
                                               localMaxima(sec$y)[max.tail.peaks[2]]-1,
                                               max(x$percentage)
                           ))
    
  } else {
    # get one plot with raw numbers (non-log)
    raw <- ggplot() +
      theme_pubr() +
      annotate(xmax = x$percentage[localMinima(sec$y)[2]],
               xmin = min(x$percentage),
               ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "tomato") +
      annotate(xmax = max(x$percentage),
               xmin = x$percentage[localMinima(sec$y)[2]],
               ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "aquamarine") +
      annotate(xmax = x$percentage[localMinima(sec$y)[min.tail.peaks[2]]],
               xmin = max(x$percentage),
               ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "tomato") +
      geom_line(aes(x = x$percentage, y = x$frequency)) +
      geom_point(aes(x = x$percentage[localMaxima(sec$y)[1:2]],
                     y = x$frequency[localMaxima(sec$y)[1:2]]), colour = "tomato") +
      geom_point(aes(x = x$percentage[localMaxima(sec$y)[max.tail.peaks[2]:max.tail.peaks[1]]],
                     y = x$frequency[localMaxima(sec$y)[max.tail.peaks[2]:max.tail.peaks[1]]]),
                 colour = "tomato") +
      geom_point(aes(x = x$percentage[localMinima(sec$y)[1:2]],
                     y = x$frequency[localMinima(sec$y)[1:2]]), colour = "royalblue") +
      geom_point(aes(x = x$percentage[localMinima(sec$y)[min.tail.peaks[2]:min.tail.peaks[1]]],
                     y = x$frequency[localMinima(sec$y)[min.tail.peaks[2]:min.tail.peaks[1]]]),
                 colour = "royalblue") +
      labs(x = "", y = "Frequency") +
      theme(axis.title = element_text(size = 12))
    
    # get another plot with data used to identify the points of maximum acceleration and minimum deceleration
    logged <- ggplot() +
      theme_pubr() +
      annotate(xmax = pred$x[localMinima(sec$y)[2]],
               xmin = min(x$percentage),
               ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "tomato") +
      annotate(xmax = max(pred$x),
               xmin = pred$x[localMinima(sec$y)[2]],
               ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "aquamarine") +
      annotate(xmax = pred$x[localMinima(sec$y)[min.tail.peaks[2]]],
               xmin = max(x$percentage),
               ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "tomato") +
      geom_line(aes(x = pred$x, y = pred$y)) +
      geom_point(aes(x = x$percentage, y = x$log.freq), alpha = 0.5, colour = "black") +
      geom_point(aes(x = pred$x[localMaxima(sec$y)[1:2]],
                     y = pred$y[localMaxima(sec$y)[1:2]]), colour = "tomato") +
      geom_point(aes(x = pred$x[localMinima(sec$y)[1:2]],
                     y = pred$y[localMinima(sec$y)[1:2]]), colour = "royalblue") +
      geom_point(aes(x = pred$x[localMaxima(sec$y)[max.tail.peaks[2]:max.tail.peaks[1]]],
                     y = pred$y[localMaxima(sec$y)[max.tail.peaks[2]:max.tail.peaks[1]]]),
                 colour = "tomato") +
      geom_point(aes(x = pred$x[localMinima(sec$y)[min.tail.peaks[2]:min.tail.peaks[1]]],
                     y = pred$y[localMinima(sec$y)[min.tail.peaks[2]:min.tail.peaks[1]]]),
                 colour = "royalblue") +
      labs(x = "", y = "Frequency (log)") +
      theme(axis.title = element_text(size = 12))
    
    # Show second derivative used to find the points
    deriv <- ggplot() +
      theme_pubr() +
      geom_line(aes(x = sec$x, y = sec$y * 1000)) +
      geom_point(aes(x = sec$x[localMaxima(sec$y)[1:2]],
                     y = sec$y[localMaxima(sec$y)[1:2]]* 1000), colour = "tomato") +
      geom_point(aes(x = sec$x[localMinima(sec$y)[1:2]],
                     y = sec$y[localMinima(sec$y)[1:2]]* 1000), colour = "royalblue") +
      geom_point(aes(x = sec$x[localMaxima(sec$y)[max.tail.peaks[2]:max.tail.peaks[1]]],
                     y = sec$y[localMaxima(sec$y)[max.tail.peaks[2]:max.tail.peaks[1]]]* 1000),
                 colour = "tomato") +
      geom_point(aes(x = sec$x[localMinima(sec$y)[min.tail.peaks[2]:min.tail.peaks[1]]],
                     y = sec$y[localMinima(sec$y)[min.tail.peaks[2]:min.tail.peaks[1]]]* 1000),
                 colour = "royalblue") +
      labs(x = "", y = expression(paste("Acceleration x10"^3, " (2"^"nd", " derivative)"))) +
      theme(axis.title = element_text(size = 12))
    
    p <- ggarrange(raw, logged, deriv, ncol = 3, labels =  c("a","b","c"))
    # add x axis title to be in the middle of two panels
    
    # Extract identified threshold
    thres.df <- data.frame(cs.flag = c("Satellite", "In-between", "Core"),
                           occup.thres.min = c(min(x$percentage), 
                                               localMinima(sec$y)[2]+1,
                                               localMinima(sec$y)[min.tail.peaks[2]]),
                           occup.thres.max = c(localMinima(sec$y)[2],
                                               localMinima(sec$y)[min.tail.peaks[2]]-1,
                                               max(x$percentage)
                           ))
    
  }
  
  out <- list(thres.df, p)
  return(out)
})

# Plot surface water vs sediment merge rep 1 and all for supplements
p <- ggarrange(annotate_figure(em.list[[1]][[2]], top = "Surface water"),
          annotate_figure(em.list[[2]][[2]], top = "Sediment", bottom = "Occupancy (%)"), nrow = 2)

ggsave(paste0("./3_emergent.threshold/em.thres_curves_sed.surf.png"), p, width = 25, height = 15, unit = "cm", dpi = 250)

# Extract smoothing
fitted.df <- ldply(flattened.ls, function(x){
  x$log.freq <- log(x$frequency)
  # make a smooth curve
  spl <- smooth.spline(x$percentage, x$log.freq, spar = 0.5)
  # predict to get fit
  x$loess.5 <- predict(spl)$y
  # get second derivative
  x$loess.5.secderiv <- predict(spl, deriv = 2)$y
  return(x)
})

# exclude, rarity cutoffs
fitted.df <- separate(fitted.df, col = '.id', into = c("dataset", "replicate.merging", "rarity.cutoff", "sample.type"),
                     remove = F, sep = "_")
fitted.df$dataset <- factor(fitted.df$dataset, levels = c("mf","mz"), labels = c("Molecular formulae",
                                                                               "Peaks"))
fitted.df <- fitted.df[fitted.df$rarity.cutoff == "all",]

# Extract results --------------------------------------------------------------------------
# Identified thresholds
thres.df <- ldply(em.list, "[[", 1)
thres.df <- separate(thres.df, col = '.id', into = c("dataset", "replicate.merging", "rarity.cutoff", "sample.type"),
         remove = F, sep = "_")
thres.df$dataset <- factor(thres.df$dataset, levels = c("mf","mz"), labels = c("Molecular formulae",
                                                                       "Peaks"))

# Take only replicate merging, remove rarity cutoffs
thres.df <- thres.df[thres.df$rarity.cutoff == "all",]

# Save as table
# write.table(thres.df, paste0("./4_gather.thresholds/keys/emergent_tresholds_key_",Sys.Date(), ".csv"),
#             sep = ",",
#             row.names = F)

# Add to fitted data frame
# pivot wider, threshold key
setDT(thres.df); setDT(fitted.df)
thres.df[, thres := paste(occup.thres.min, occup.thres.max, sep = "-")]
thres.wide <- dcast(thres.df, formula = .id ~ cs.flag, value.var = "thres")
fitted.df <- fitted.df[thres.wide, , on = .(.id)]
fitted.df[, occupancy := NULL]

# Save fitted data
# write.table(fitted.df, "./3_emergent.threshold/output/em.thres_loess_spar0.5_fitted.csv",
#             sep = ",", dec = ".", row.names = F)

# Knit a table for Github
# thres.df <- read.csv("./3_emergent.threshold/output/emergent_tresholds_key_2022-03-01.csv",
#                      sep = ",", stringsAsFactors = F)

# Select necessary columns
# thres.df <- thres.df %>% dplyr::select(-.id)
# 
# # Re-label
# colnames(thres.df) <- c("Dataset", "Replicate merging", "Rarity cutoff", "Sample type", "C-S flag",
#                           "Threshold start (Occupancy)", "Threshold end (Occupancy)", "Threshold Range")
# # Re-level
# thres.df$`Sample type` <- factor(thres.df$`Sample type`, levels = c('water',
#                                                                     'sed'),
#                                  labels = c('Water','Sediment'))
# knitr::kable(thres.df)
# 
# setDT(thres.df)
# thres.df[, max(`Threshold start (Occupancy)`), by = .(Dataset, `Sample type`, `C-S flag`)]

# Save generated plots
title <- c("Molecular formulae: Surface water - rep.merged1 - all",
           "Molecular formulae: Sediment - rep.merged1 - all",
           "Molecular formulae: Surface water - rep.merged1 - rar1",
           "Molecular formulae: Sediment - rep.merged1 - rar1",
           "Molecular formulae: Surface water - rep.merged1 - rar2",
           "Molecular formulae: Sediment - rep.merged1 - rar2",
           "Molecular formulae: Surface water - rep.merged2 - all",
           "Molecular formulae: Sediment - rep.merged2 - all",
           "Molecular formulae: Surface water - rep.merged2 - rar1",
           "Molecular formulae: Sediment - rep.merged2 - rar1",
           "Molecular formulae: Surface water - rep.merged2 - rar2",
           "Molecular formulae: Sediment - rep.merged2 - rar2",
           "Molecular formulae: Surface water - rep.merged1 - all",
           "Peaks: Sediment - rep.merged1 - all",
           "Peaks: Surface water - rep.merged1 - rar1",
           "Peaks: Sediment - rep.merged1 - rar1",
           "Peaks: Surface water - rep.merged1 - rar2",
           "Peaks: Sediment - rep.merged1 - rar2",
           "Peaks: Surface water - rep.merged2 - all",
           "Peaks: Sediment - rep.merged2 - all",
           "Peaks: Surface water - rep.merged2 - rar1",
           "Peaks: Sediment - rep.merged2 - rar1",
           "Peaks: Surface water - rep.merged2 - rar2",
           "Peaks: Sediment - rep.merged2 - rar2")

for(i in 1:length(em.list)){
  p <- annotate_figure(em.list[[i]][[2]], top = text_grob(title[i]),
                       bottom = "Occupancy (%)")
  ggsave(paste0("./3_emergent.threshold/prelim_figures/em.thres_", names(em.list)[i],
                ".png"), p, width = 25, height = 9, unit = "cm", dpi = 250)
}


# Apply threshold -----------------------------------------------------------------------------------------------------
# We're done with getting the threshold now, we must apply the threshold to the dataset and
# export in a way that other users can use it

# What we will do is, calculate the percentage occupancy for each peak/molecular formulae in dataset
# and then we apply our thresholds
# Output will be saved in the cross table, so it's easy to use by other users

# Calcualte percentage occupancy for each peak/molecular formulae


occup.ls <- llply(all.ls, function(x){
  setDT(x)
  # Separate sediment and surface water by ID and remove ID column
  sed <- x[str_detect(ID, "SED"),][,-1]
  water <- x[str_detect(ID, "SW"),][,-1]
  
  # Process surface water
  
  # how often were MF found across water samples?
  sum_water <- data.frame(occupancy = colSums(water)) # sum of presence-absence = number of sites present
  setDT(sum_water, keep.rownames = "mf.mz")
  
  # Calculate percentage
  sum_water[, perc.occup := round(occupancy * 100 / nrow(water), 2)] # divide by n
  
  # Process sediment
  # Sequence of code is same as water samples.
  sum_sed <- data.frame(occupancy = colSums(sed))
  setDT(sum_sed, keep.rownames = "mf.mz")
  
  # Calculate percentage
  sum_sed[, perc.occup := round(occupancy * 100 / nrow(sed), 2)] # divide by n
  
  # combine water and sediment together
  # sum_water[sum_sed, c("occupancy.sed",
  #                      "perc.occup.sed") := list(occupancy.sed, perc.occup.sed), on = .(mf.mz)]
  # sanity check
  # ncol(x)-1 ==  length(unique(sum_water$mf.mz))
  
  return(list(sum_water, sum_sed))
  
}, .parallel = T)

# Flatten list into a one-level list
occup.ls <- purrr::flatten(occup.ls)
# add names
names(occup.ls) <- c("mf_rep.merged1_all_water", "mf_rep.merged1_all_sed",
                         "mf_rep.merged1_rar1_water", "mf_rep.merged1_rar1_sed",
                         "mf_rep.merged1_rar2_water", "mf_rep.merged1_rar2_sed",
                         "mf_rep.merged2_all_water", "mf_rep.merged2_all_sed",
                         "mf_rep.merged2_rar1_water","mf_rep.merged2_rar1_sed", 
                         "mf_rep.merged2_rar2_water","mf_rep.merged2_rar2_sed",
                         "mz_rep.merged1_all_water","mz_rep.merged1_all_sed",
                         "mz_rep.merged1_rar1_water","mz_rep.merged1_rar1_sed",
                         "mz_rep.merged1_rar2_water","mz_rep.merged1_rar2_sed",
                         "mz_rep.merged2_all_water","mz_rep.merged2_all_sed",
                         "mz_rep.merged2_rar1_water","mz_rep.merged2_rar1_sed",
                         "mz_rep.merged2_rar2_water","mz_rep.merged2_rar2_sed")

# Apply threshold values
occup.df <- bind_rows(occup.ls, .id = ".id")

# Make another wide threshold df
thres.wide <- dcast(thres.df, formula = .id ~ cs.flag, value.var = c("occup.thres.min","occup.thres.max"))
# merge
occup.df <- occup.df[thres.wide, , on = .(.id)]

# Identify cs group for each mf/mz
occup.df[round(perc.occup,0) <= occup.thres.max_Satellite, cs.flag := "Satellite"]
occup.df[round(perc.occup,0) >= occup.thres.min_Core, cs.flag := "Core"]
occup.df[round(perc.occup,0) > occup.thres.max_Satellite & round(perc.occup,0) < occup.thres.min_Core, cs.flag :=  "In-between"]

# sanity check
occup.df[is.na(cs.flag),]

# Clean df
occup.df <- occup.df %>% dplyr::select(.id:perc.occup, cs.flag)

# Read to crosstable
mf.ls <- vector("list", 6)
# Run loop to read in data
for(i in 1:6){
  data_cross <- read.csv(paste0("./1_data.cleaning/output/", mf.cross[i]),
                          sep = ",", stringsAsFactors = F)
  
  # assign to list bin
  mf.ls[[i]] <- data_cross
}

# Some sanity check
lapply(mf.ls, dim)
# Add names to list
names(mf.ls) <- c("mf_rep.merged1_all", "mf_rep.merged1_rar1", "mf_rep.merged1_rar2",
                  "mf_rep.merged2_all", "mf_rep.merged2_rar1", "mf_rep.merged2_rar2")

# bind into one dataframe
cross <- bind_rows(mf.ls, .id = ".id")

# For molecular formulae, rearrange occup.df and merge
# Separate sed and water

occup.df <- separate(occup.df, col = '.id', into = c("dataset", "replicate.merging", "rarity.cutoff", "sample.type"),
                     remove = F, sep = "_")

out <- occup.df[rarity.cutoff == "all",] %>% setDT()
# make ID same as cross table
out[, .id := paste(dataset, replicate.merging, rarity.cutoff, sep = "_")]

out.mf <- out[dataset == "mf",]
out.mz <- out[dataset == "mz",]

# change column names
colnames(out.mf)[c(6,9)] <- c("MolForm", "cs.flag.emergent")
out.mf <- dcast(out.mf, formula = .id + MolForm ~ sample.type, value.var = c("occupancy", "perc.occup", "cs.flag.emergent"))

# merge with cross table
setDT(cross)
cross <- cross[out.mf, , on = .(MolForm, .id)]

# Separate into different data frames by ID for saving
cross.ls <- split(cross, by = ".id")

# Do same for peaks
# change column names
colnames(out.mz)[c(6,9)] <- c("Mass", "cs.flag.emergent")
# Change mass into numeric
out.mz[, Mass := str_remove(Mass, "X")]
out.mz <- dcast(out.mz, formula = .id + Mass ~ sample.type, value.var = c("occupancy", "perc.occup", "cs.flag.emergent"))

# Separate into different data frames by ID for saving
cross.mz.ls <- split(out.mz, by = ".id")

# Saving paths
# create list with names
save.vec <- c(paste0("./4_gather.thresholds/FTICR_crosstable_rep.merged1_all_em.thres_",Sys.Date(),".csv"),
              paste0("./4_gather.thresholds/FTICR_crosstable_rep.merged2_all_em.thres_",Sys.Date(),".csv"),
              paste0("./4_gather.thresholds/FTICR_peaks_crosstable_rep.merged1_all_em.thres_",Sys.Date(),".csv"),
              paste0("./4_gather.thresholds/FTICR_peaks_crosstable_rep.merged2_all_em.thres_",Sys.Date(),".csv"))

out.ls <- c(cross.ls, cross.mz.ls)
# Save community matrix
mapply(write.table, out.ls, save.vec, MoreArgs = list(sep = ",", dec = ".", row.names = F))


# Compare smoothing methods -----------------------------------------------------------------------------------------

smooth.methods <- read.csv("./3_emergent.threshold/output/em.thres_loess_spar0.5_fitted_added_2nd.csv",
                           sep = ",", dec = ".", stringsAsFactors = F)


# Calculate thresholds

# Plot smoothing
# make long format for ggplot

plot.df <- smooth.methods %>% dplyr::select(.id:loess.5, moving.avg.1:knot55.1) %>% setDT()

# widen
plot.df <- melt(plot.df[, c(".id", "frequency") := list(NULL, NULL)],
     id.vars = c("dataset", "replicate.merging", "rarity.cutoff", "sample.type", "percentage", "log.freq"),
     variable.name = "method", value.name = "smooth.value")

# MF - Sediment
mf <- plot.df[dataset == "Molecular formulae",]
ggplot() +
  theme_bw()+
  geom_point(data = mf[sample.type == "sed",], aes(x = percentage, y = log.freq),
             size = 1, alpha = 0.5) +
  geom_line(data = mf[sample.type == "sed",], aes(x = percentage, y = smooth.value, colour = method), size =2) +
  scale_colour_viridis_d(option = "magma", end = 0.9) +
  facet_wrap(vars(replicate.merging, method), scale = "free_y") +
  
  labs(x = "Occupancy (%)", y = "Frequency", title = "MF - Sediment")

# MF - water
ggplot() +
  theme_bw()+
  geom_point(data = mf[sample.type == "water",], aes(x = percentage, y = log.freq),
             size = 1, alpha = 0.5) +
  geom_line(data = mf[sample.type == "water",], aes(x = percentage, y = smooth.value, colour = method)) +
  scale_colour_viridis_d(option = "magma", end = 0.9) +
  facet_wrap(vars(replicate.merging, method), scale = "free_y") +
  labs(x = "Occupancy (%)", y = "Frequency", title = "MF - Surface water")

# Peaks - Sediment
mz <- plot.df[dataset == "Peaks",]
ggplot() +
  theme_bw()+
  geom_point(data = mz[sample.type == "sed",], aes(x = percentage, y = log.freq),
             size = 1, alpha = 0.5) +
  geom_line(data = mz[sample.type == "sed",], aes(x = percentage, y = smooth.value, colour = method)) +
  scale_colour_viridis_d(option = "magma", end = 0.9) +
  facet_wrap(vars(replicate.merging, method), scale = "free_y") +
  labs(x = "Occupancy (%)", y = "Frequency", title = "Peaks - Sediment")

# Peaks - Surface water
ggplot() +
  theme_bw()+
  geom_point(data = mz[sample.type == "water",], aes(x = percentage, y = log.freq),
             size = 1, alpha = 0.5) +
  geom_line(data = mz[sample.type == "water",], aes(x = percentage, y = smooth.value, colour = method)) +
  scale_colour_viridis_d(option = "magma", end = 0.9) +
  facet_wrap(vars(replicate.merging, method), scale = "free_y") +
  labs(x = "Occupancy (%)", y = "Frequency", title = "Peaks - Surface water")

# Get thresholds ------------------------------------------------------------------------------------------------------------
# Run function by group

thres.df <- smooth.methods %>% dplyr::select(.id:loess.5.secderiv, moving.avg.1:knot55.1, moving.avg.2:knot55.2) %>% setDT()

# into long format
thres.df <- melt(thres.df[, c(".id", "frequency") := list(NULL, NULL)],
                id.vars = c("dataset", "replicate.merging", "rarity.cutoff", "sample.type", "percentage", "log.freq"),
                variable.name = "method", value.name = "smooth.value")

# extract patterns
thres.df[, method.uni := sub("\\..*", "", method)]
thres.df[method == "loess.5", method.uni := "loess.5"]
thres.df[method == "loess.5.secderiv", method.uni := "loess.5"]
thres.df[, deriv := ifelse(str_detect(method, ".2"), "deriv", "pred")]
thres.df[method == "loess.5.secderiv", deriv := "deriv"]

# widen
# dcast(thres.df, dataset + replicate.merging + rarity.cutoff + sample.type + percentage + log.freq ~ method.uni + deriv,
#       value.var = "smooth.value")

#x <- thres.df[dataset == "Peaks" & sample.type == "sed" & method.uni == "knot25",]

all.thres <- ddply(thres.df, .(dataset, replicate.merging, sample.type, method.uni), function(x){
  setDT(x)
  pred <- x[deriv == "pred",]$smooth.value
  deriv <- x[deriv == "deriv",]$smooth.value 
  
  # get number of peaks to classify core
  min.tail.peaks <- c(length(localMinima(deriv)), length(localMinima(deriv)) -1)
  max.tail.peaks <- c(length(localMaxima(deriv)), length(localMaxima(deriv)) -1)
  
  if(localMaxima(deriv)[max.tail.peaks[1]] > localMinima(deriv)[min.tail.peaks[1]]){
    # Extract identified threshold
    out <- data.frame(cs.flag = c("Satellite", "In-between", "Core"),
                           occup.thres.min = c(min(x$percentage), 
                                               localMinima(deriv)[2]+1,
                                               localMaxima(deriv)[max.tail.peaks[2]]),
                           occup.thres.max = c(localMinima(deriv)[2],
                                               localMaxima(deriv)[max.tail.peaks[2]]-1,
                                               max(x$percentage)
                           ))
  } else {
    # Extract identified threshold
    out <- data.frame(cs.flag = c("Satellite", "In-between", "Core"),
                           occup.thres.min = c(min(x$percentage), 
                                               localMinima(deriv)[2]+1,
                                               localMinima(deriv)[min.tail.peaks[2]]),
                           occup.thres.max = c(localMinima(deriv)[2],
                                               localMinima(deriv)[min.tail.peaks[2]]-1,
                                               max(x$percentage)
                           ))
  }
  return(out)

})

# Verifying smoothing method --------------------------------------------------------------------
setDT(all.thres)
all.thres[cs.flag == "Core", thres := occup.thres.min]
all.thres[cs.flag == "Satellite", thres := occup.thres.max]
all.thres[thres > 100, thres := NA]

(p <- ggplot(all.thres[cs.flag != "In-between",], aes(x = cs.flag, y = thres)) +
  theme_bw() +
  geom_boxplot() +
  scale_y_continuous(limits = c(0,100)) +
  facet_grid(.~dataset) +
  labs(x = "", y = "Occupancy (%)"))

ggsave("./3_emergent.threshold/threshold.values_boxplot.png", p,
       width = 8, height = 5, units = "cm", dpi = 300)


# Identify outliers
library(rstatix)

outliers <- rbind(all.thres %>%
  group_by(cs.flag, dataset) %>%
  identify_outliers(occup.thres.min) %>%
  select(dataset, replicate.merging, sample.type, method.uni) %>%
  unique(),
  all.thres %>%
  group_by(cs.flag, dataset) %>%
  identify_outliers(occup.thres.max) %>%
  select(dataset, replicate.merging, sample.type, method.uni) %>%
  unique()) %>% unique()

View(all.thres %>%
  group_by(cs.flag, dataset) %>%
  identify_outliers(occup.thres.min))

# how many of the smoothing methods are outliers?
nlevels(factor(all.thres$method.uni)) *2 *2 *2
# 104

nrow(outliers) # 22 out of 104 were identified as outliers, that's not too bad

# add a .1 to method to be able to merge
outliers <- outliers %>% mutate(method.uni = paste0(method.uni, ".1"))
# make an ID
outliers <- outliers %>% mutate(ID = paste(replicate.merging, sample.type, method.uni, sep = "_"))

# seems like three smoothing methods are the only outliers
mf <- mf %>% mutate(ID = paste(replicate.merging, sample.type, method, sep = "_"))
mz <- mz %>% mutate(ID = paste(replicate.merging, sample.type, method, sep = "_"))

# select only those that are outliers
sub.mf <- mf %>% filter(ID %in% outliers[outliers$dataset == "Molecular formulae",] $ID)
sub.mz <- mz %>% filter(ID %in% outliers[outliers$dataset == "Peaks",] $ID)

# plot smooth curves of outliers
outly.df <- rbind(sub.mf, sub.mz)

ggplot() +
  theme_bw()+
  geom_point(data = outly.df, aes(x = percentage, y = log.freq),
             size = 1, alpha = 0.5) +
  geom_line(data = outly.df, aes(x = percentage, y = smooth.value, colour = method)) +
  scale_colour_viridis_d(option = "magma", end = 0.9) +
  facet_wrap(vars(dataset, sample.type, replicate.merging, method), scale = "free_y") +
  labs(x = "Occupancy (%)", y = "Frequency", title = "Peaks - Sediment")


# plot derivatives to see what is wrong
knot5 <- thres.df[deriv == "deriv" & method.uni == "knot5",]
knot15 <- thres.df[deriv == "deriv" & method.uni == "knot15",]
knot25 <- thres.df[deriv == "deriv" & method.uni == "knot25",]
moving <- thres.df[deriv == "deriv" & method.uni == "moving",]
loess5 <- thres.df[deriv == "deriv" & method.uni == "loess.5",]

ggplot(knot5, aes(x = percentage, y = smooth.value)) +
  facet_wrap(dataset ~ sample.type) +
  geom_line(aes(colour = replicate.merging)) +
  labs(title ="2nd derivative of knot5")

ggplot(knot15, aes(x = percentage, y = smooth.value)) +
  facet_wrap(dataset ~ sample.type) +
  geom_line(aes(colour = replicate.merging)) +
  labs(title ="2nd derivative of knot15")

ggplot(knot25, aes(x = percentage, y = smooth.value)) +
  facet_wrap(dataset ~ sample.type) +
  geom_line(aes(colour = replicate.merging)) +
  labs(title ="2nd derivative of knot25")

ggplot(moving, aes(x = percentage, y = smooth.value)) +
  facet_wrap(dataset ~ sample.type) +
  geom_line(aes(colour = replicate.merging))  +
  labs(title ="2nd derivative of moving")

ggplot(loess5, aes(x = percentage, y = smooth.value)) +
  facet_wrap(dataset ~ sample.type) +
  geom_line(aes(colour = replicate.merging))+
  labs(title ="2nd derivative of loess (0.5)")

# check actual points of acceleration
dlply(knot15, .(dataset, sample.type, replicate.merging), function(x){
  localMinima(x$smooth.value)
})

# so what is the issue in knot methods?
# it seems that knot5 has very few acceleration points and knot25 has too many. Hence the threshold identifying method does not work
# for knot15 it seems that the satellite decline doesn't have a hump shape and hence, it doesn't work.


# Identify overlap groups -------------------------------------------------------------------------------
# Change Satellites that are 0 but classified as satellite

cross <- c(list.files("./4_gather.thresholds/", pattern = "2022-03-07.csv"),
           list.files("./4_gather.thresholds/", pattern = "peaks"))
files <- list()

for(i in 1:length(cross)){
  files[[i]] <- read.csv(paste0("./4_gather.thresholds/", cross[i]),
           sep = ",", dec = ".", stringsAsFactors = F) %>% setDT()
}

cor.files <- lapply(files, function(x){
  cols <- colnames(x)[str_detect(colnames(x), pattern = "cs.flag")]
  seds <- cols[str_detect(cols, pattern = "sed")]
  wats <- cols[str_detect(cols, pattern = "water")]
  x<-x[perc.occup_sed == 0, (seds) := as.list(rep(NA, length(seds)))]
  x<-x[perc.occup_water == 0, (wats) := as.list(rep(NA, length(wats)))]
  return(x)
})

# sanity check
cor.files[[1]][perc.occup_water == 0,]

# write files
for(i in 1:length(cor.files)){
  write.table(cor.files[[i]],
              paste0("./4_gather.thresholds/", str_replace(cross[i], pattern = "2022-03-07|2022-03-01", 
                                                           replacement = paste(Sys.Date()))),
              sep = ',', dec = '.', row.names = F)
}

# Find overlap -------------------------------------------------------------------------------------------------

cross <- c(list.files("./4_gather.thresholds/", pattern = "2022-03-18.csv"))
files <- list()

for(i in 1:length(cross)){
  files[[i]] <- read.csv(paste0("./4_gather.thresholds/", cross[i]),
                         sep = ",", dec = ".", stringsAsFactors = F) %>% setDT()
}

# add new columns with overlap IDs
cor.files <- lapply(files, function(x){
  
  # Add a column of habitat overlap
  x[occupancy_sed > 0 & occupancy_water > 0, habitat.overlap := "y"]
  x[is.na(habitat.overlap), habitat.overlap := "n"]
  
  cols <- colnames(x)[str_detect(colnames(x), pattern = "cs.flag")]
  #em <- cols[str_detect(cols, pattern = "emergent")]
  pca <- cols[str_detect(cols, pattern = "pca")]
  
  # Emergent is in both MF and peaks dataset
  # make new columns
  x[cs.flag.emergent_sed == "Core" & cs.flag.emergent_water == "Core", "cs.flag.emergent_overlap" := "Global core"]
  x[cs.flag.emergent_sed == "Satellite" & cs.flag.emergent_water == "Satellite", "cs.flag.emergent_overlap" := "Global satellite"]
  x[cs.flag.emergent_sed == "In-between" & cs.flag.emergent_water == "In-between", "cs.flag.emergent_overlap" := "Global in-between"]
  
  # Add another column less detailed
  x[, cs.flag.emergent_general.overlap := cs.flag.emergent_overlap]
  
  # Identify those that shift
  x[cs.flag.emergent_sed == "Core" & cs.flag.emergent_water == "Satellite", "cs.flag.emergent_overlap" := "Sed Core - Water Sat"]
  x[cs.flag.emergent_sed == "Core" & cs.flag.emergent_water == "In-between", "cs.flag.emergent_overlap" := "Sed Core - Water Inbetween"]
  x[cs.flag.emergent_sed == "Satellite" & cs.flag.emergent_water == "Core", "cs.flag.emergent_overlap" := "Sed Sat - Water Core"]
  x[cs.flag.emergent_sed == "Satellite" & cs.flag.emergent_water == "In-between", "cs.flag.emergent_overlap" := "Sed Sat - Water Inbetween"]
  x[cs.flag.emergent_sed == "In-between" & cs.flag.emergent_water == "Core", "cs.flag.emergent_overlap" := "Sed Inbetween - Water Core"]
  x[cs.flag.emergent_sed == "In-between" & cs.flag.emergent_water == "Satellite", "cs.flag.emergent_overlap" := "Sed Inbetween - Water Sat"]
  
  x[cs.flag.emergent_sed != cs.flag.emergent_water,
    cs.flag.emergent_general.overlap := "Shifter"]
  
  # sanity check
  # x[habitat.overlap == "y" & is.na(cs.flag.emergent_overlap),] # should be NA
  # x[habitat.overlap == "y" & is.na(cs.flag.emergent_overlap),] # should be NA
  
  # For MF dataset do...
  if(length(pca) != 0L){
    # Do for other thresholds
    x[cs.flag.pca_sed == "Core" & cs.flag.pca_water == "Core", cs.flag.pca_general.overlap := "Global core"]
    x[cs.flag.pca_sed == "Satellite" & cs.flag.pca_water == "Satellite", cs.flag.pca_general.overlap := "Global satellite"]
    x[cs.flag.pca_sed != cs.flag.pca_water, cs.flag.pca_general.overlap := "Shifter"]
    
    x[cs.flag.rf_sed == "Core" & cs.flag.rf_water == "Core", cs.flag.rf_general.overlap := "Global core"]
    x[cs.flag.rf_sed == "Satellite" & cs.flag.rf_water == "Satellite", cs.flag.rf_general.overlap := "Global satellite"]
    x[cs.flag.rf_sed != cs.flag.rf_water, cs.flag.rf_general.overlap := "Shifter"]
    
    # Sanity check
    # x[habitat.overlap == "y" & is.na(cs.flag.pca_general.overlap),] # should be NA
    # x[habitat.overlap == "y" & is.na(cs.flag.rf_general.overlap),] # should be NA
  }
  
  return(x)
})

# write files
for(i in 1:length(cor.files)){
  write.table(cor.files[[i]],
              paste0("./4_gather.thresholds/", str_replace(cross[i], pattern = "2022-03-18", 
                                                           replacement = paste(Sys.Date()))),
              sep = ',', dec = '.', row.names = F)
}

# Check formula assignment quality ----------------------------------------------------------------------------------------------------------
# Have commonly removed MFs been removed?

# After Herzsprung et al 2014 Anal Bioanal Chem 406:30

cross <- c(list.files("./4_gather.thresholds/", pattern = "2022-03-23.csv"))
cross <- cross[!str_detect(cross, "peaks")]
files <- list()

for(i in 1:length(cross)){
  files[[i]] <- read.csv(paste0("./4_gather.thresholds/", cross[i]),
                         sep = ",", dec = ".", stringsAsFactors = F) %>% setDT()
}

# What will be removed?
any(files[[1]][OtoC_ratio < 0,]$cs.flag.emergent_sed == "Core")
any(files[[1]][OtoC_ratio < 0,]$cs.flag.emergent_water == "Core")

any(files[[1]][HtoC_ratio > 2,]$cs.flag.emergent_sed == "Core")
any(files[[1]][HtoC_ratio > 2,]$cs.flag.emergent_water == "Core")

any(files[[1]][OtoC_ratio > 1,]$cs.flag.emergent_sed == "Core")
any(files[[1]][OtoC_ratio > 1,]$cs.flag.emergent_water == "Core")

any(files[[1]][DBE_O > 10,]$cs.flag.emergent_sed == "Core") #36
any(files[[1]][DBE_O > 10,]$cs.flag.emergent_water == "Core")

# Apply filters
cor.files <- lapply(files, function(x){
  # apply MF elemental ranges
  # after T. Riedel, T. Dittmar, A method detection limit for the analysis of natural
  # organic matter via Fourier transform ion cyclotron resonance mass spectrometry. Anal. Chem. 86, 8376–8382 (2014).

  x <- x[C >= 1 & C <= 130,]
  x <- x[H >= 1 & H <= 200,]
  x <- x[O >= 1 & O <= 50,]
  x <- x[N >= 0 & N <= 4,]
  x <- x[S >= 0 & S <= 2,]
  x <- x[P >= 0 & P <= 1,]
  
  # remove certain MF that are unreliable
  # Hawkes JA, D’Andrilli J, Agar JN, Barrow MP, Berg SM, Catalán N, et al. An international laboratory comparison of
  # dissolved organic matter composition by high resolution mass spectrometry: Are we getting the same answer? Limnol Oceanogr Methods 2020; 18: 235–258. 
  x <- x[OtoC_ratio > 0 & OtoC_ratio <= 1.2,]
  x <- x[HtoC_ratio >= 0.3 & HtoC_ratio <= 2.2,]
  
  # Herzsprung P, Hertkorn N, Von Tumpling W, Harir M, Friese K, Schmitt-Kopplin P. 
  # Understanding molecular formula assignment of Fourier transform ion cyclotron resonance
  # mass spectrometry data of natural organic matter from a chemical point of view. Anal Bioanal Chem 2014; 406: 7977–7987. 
  x <- x[DBE_O <= 10,]
  
  # Add a few more variables
  x[, molweight := ((12.0107*C) + (1.00784*H) + (14.0067*N) + (15.999*O) + (32.065*S) + (30.974*P))]
  
  x[HtoC_ratio >= 1.5 & HtoC_ratio <= 2 & OtoC_ratio > 0.9 & N > 0, Group := "Peptides"]
  x[OtoC_ratio >= 0.9 & AI_Mod < 0.5, Group := "Sugars"]
  x[HtoC_ratio >= 2 & OtoC_ratio >= 0.9, Group := "Saturated fatty acids"]
  x[HtoC_ratio >= 1.5 & HtoC_ratio <= 2 & N == 0, Group := "Unsaturated aliphatics"]
  x[HtoC_ratio < 1.5 & OtoC_ratio < 0.9 & AI_Mod < 0.5, Group := "Highly unsaturated compounds"]
  x[AI_Mod >= 0.5 & C < 12, Group := "Phenols"]
  x[AI_Mod >= 0.5 & C >= 12, Group := "Polyphenols"]
  x[AI_Mod > 0.66, Group := "Combusted polycondensed aromatics"]
  x[,IOS := F]
  x[HtoC_ratio <= 1.3 & HtoC_ratio >= 1.04 & OtoC_ratio <= 0.62 & OtoC_ratio >= 0.42 & molweight <= 388 & molweight >= 332, IOS := T]
  x[HtoC_ratio <= 1.3 & HtoC_ratio >= 1.04 & OtoC_ratio <= 0.62 & OtoC_ratio >= 0.42 & molweight <= 548 & molweight >= 446, IOS := T]
  
  return(x)
})

# check if many were removed
lapply(cor.files, nrow)
lapply(files, nrow)

ggplot(files[[2]], aes(x = OtoC_ratio, HtoC_ratio)) +
  geom_point()

# Read in processed data ------------------------------------------------------------------

# This code was for when we had multiple thresholds
# cross <- c(list.files("./4_gather.thresholds/", pattern = "2022-05-05.csv"))
# cross <- cross[!str_detect(cross, "peaks")]
# files <- list()
# 
# for(i in 1:length(cross)){
#   files[[i]] <- read.csv(paste0("./4_gather.thresholds/", cross[i]),
#                          sep = ",", dec = ".", stringsAsFactors = F) %>% setDT()
# }

cross <- read.csv("./4_gather.thresholds/FTICR_crosstable_rep.merged1_all_em.thres_2022-05-05.csv",
                  sep = ",", dec = ".", stringsAsFactors = F) %>% setDT()

cross[, .(n = .N), by = .(cs.flag.emergent_sed)]
commat <- read.csv("./4_gather.thresholds/FTICR_commat_rep.merged1_2022-07-19.csv",
                   sep = ",", dec = ".", stringsAsFactors = F)

# only keep MF in cross table
commat <- commat[, which(colnames(commat) %in% cross$MolForm)]

cross[]

# write.table(commat, paste0("./4_gather.thresholds/FTICR_commat_rep.merged1_", Sys.Date(),".csv"),
#             sep = ",", dec = ".", row.names = F)

ios.df <- data.frame(HtoC_ratio = c(1.3,1.04), OtoC_ratio = c(0.62,0.42))


files[[1]][, cs.flag.emergent_overlap.mod := factor(cs.flag.emergent_overlap,
                                                        levels = c("Global core","Global in-between",
                                                                   "Global satellite",
                                                                   "Sed Core - Water Sat", "Sed Core - Water Inbetween", 
                                                                   "Sed Sat - Water Core", "Sed Inbetween - Water Core",
                                                                   "Sed Sat - Water Inbetween", "Sed Inbetween - Water Sat"),
                                                        labels = c("Global core","Global in-between",
                                                                   "Global satellite",
                                                                   "Sed Core - Water Sat/In-bet.", "Sed Core - Water Sat/In-bet.", 
                                                                   "Water Core - Sed Sat/In-bet.", "Water Core - Sed Sat/In-bet.", 
                                                                   "Sed Sat - Water Inbetween", "Water Sat - Sed Inbetween"))]

comp_groups <- data.frame(group = c("Saturated fatty acids",
                                    "Peptides and\nunsaturated\naliphatics",
                                    "Highly unsaturated\ncompounds",
                                    "Vascular plant-\nderived polyphenols\nand phenols"),
                          HtoC_ratio = c(2, 1.7, 1.25, 0.5),
                          cs.flag.emergent_overlap = "Global core") %>% setDT()

comp_groups[, cs.flag.emergent_overlap := factor(cs.flag.emergent_overlap,
                                                 levels = c("Global core","Global in-between",
                                                            "Global satellite",
                                                            "Sed Core - Water Sat", "Sed Core - Water Inbetween", 
                                                            "Sed Sat - Water Core", "Sed Inbetween - Water Core",
                                                            "Sed Sat - Water Inbetween", "Sed Inbetween - Water Sat"),
                                                 labels = c("Global core","Global in-between",
                                                            "Global satellite",
                                                            "Sed Core - Water Sat/In-bet.", "Sed Core - Water Sat/In-bet.", 
                                                            "Water Core - Sed Sat/In-bet.", "Water Core - Sed Sat/In-bet.", 
                                                            "Sed Sat - Water Inbetween", "Water Sat - Sed Inbetween"))]

cross[, cs.flag.emergent_overlap := factor(cs.flag.emergent_overlap,
                                           levels = c("Global core","Global in-between",
                                                      "Global satellite",
                                                      "Sed Core - Water Sat", "Sed Core - Water Inbetween", 
                                                      "Sed Sat - Water Core", "Sed Inbetween - Water Core",
                                                      "Sed Sat - Water Inbetween", "Sed Inbetween - Water Sat"),
                                           labels = c("Global core","Global in-between",
                                                      "Global satellite",
                                                      "Sed Core - Water Sat/In-bet.", "Sed Core - Water Sat/In-bet.", 
                                                      "Water Core - Sed Sat/In-bet.", "Water Core - Sed Sat/In-bet.", 
                                                      "Sed Sat - Water Inbetween", "Water Sat - Sed Inbetween"))]


# check if data looks fine.
# Do some Van Krevelen
library(grid)
p <- ggplot(cross[!is.na(cs.flag.emergent_overlap),], aes(x = OtoC_ratio, HtoC_ratio)) +
  theme_bw() +
  geom_point(aes(colour = cs.flag.emergent_overlap), alpha = 0.7) +
  scale_colour_manual(values = c("#999999", "#661100", "#0072B2",
                                 "#FFB000", "#785EF0", "#FE6100", "#648FFF"),
                      name = "Overlapping categories") +
  facet_wrap(.~cs.flag.emergent_overlap, ncol =2) +
  geom_hline(yintercept = 1.5, linetype = "dashed") +
  geom_abline(intercept = 1.1, slope = -0.3, linetype = "dashed") +
  #stat_ellipse(data = ios.df, aes(x = OtoC_ratio, y = HtoC_ratio), type = "norm",
  #              colour = "tomato", linetype = "dashed") +
  labs(x = "O/C", y = "H/C") +
  lims(x = c(0,1.25), y = c(0,2)) +
  annotation_custom(grob=circleGrob(r=unit(1,"npc"),
                                    gp = gpar(col = 'black', lty = 3, fill = 'transparent')),
                    xmin=ios.df$OtoC_ratio[2], xmax=ios.df$OtoC_ratio[1],
                    ymin=ios.df$HtoC_ratio[2], ymax=ios.df$HtoC_ratio[1])
  

(p <- p + geom_text(data = comp_groups, aes(x = 0.78, y = HtoC_ratio, label = group), size = 3, 
              hjust = 0, lineheight = 0.7))

# save
# ggsave("./Supplementary/VK_by_overlap.groups.png", p, dpi = 300,
#        height = 20, width = 28, units = "cm")
ggsave("./Supplementary/VK_by_overlap.groups.png", p, dpi = 300,
       height = 25, width = 23, units = "cm")

# Save filtered dataset
# write files
for(i in 1:length(cor.files)){
  write.table(cor.files[[i]],
              paste0("./4_gather.thresholds/", str_replace(cross[i], pattern = "2022-03-23", 
                                                           replacement = paste(Sys.Date()))),
              sep = ',', dec = '.', row.names = F)
}

#-- End: Addition by MS
