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
meta <- read.csv("./1_data.cleaning/output/FTICR_meta_eachriver_2021-11-03.csv",
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
x<- flattened.ls[[1]]

# We're using smooth spline here instead of defining a function and getting a derivative from the function
# as described in https://cran.r-project.org/web/packages/inflection/vignettes/inflectionMissionImpossible.html
# as it is hard to come up with a function describing the observed pattern
# not exponential, not logarithmic etc

em.list <- llply(flattened.ls, function(x){
  # take log
  x$log.freq <- log(x$frequency)
  # make a smooth curve
  spl <- smooth.spline(x$occupancy, x$log.freq, spar = 0.5)
  # predict to get fit
  pred <- predict(spl)
  # get second derivative
  sec <- predict(spl, deriv = 2)
  
  # get number of peaks to classify core
  min.tail.peaks <- c(length(localMinima(sec$y)), length(localMinima(sec$y)) -1)
  max.tail.peaks <- c(length(localMaxima(sec$y)), length(localMaxima(sec$y)) -1)
  
  options(scipen = 999) # avoid scientific annotations
  # get one plot with raw numbers (non-log)
  raw <- ggplot() +
    theme_pubr() +
    annotate(xmax = x$occupancy[localMinima(sec$y)[2]],
             xmin = min(x$occupancy),
             ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "tomato") +
    annotate(xmax = max(x$occupancy),
             xmin = x$occupancy[localMinima(sec$y)[2]],
             ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "aquamarine") +
    annotate(xmax = x$occupancy[localMaxima(sec$y)[max.tail.peaks[1]]],
             xmin = max(x$occupancy),
             ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "tomato") +
    geom_line(aes(x = x$occupancy, y = x$frequency)) +
    geom_point(aes(x = x$occupancy[localMaxima(sec$y)[1:2]],
                   y = x$frequency[localMaxima(sec$y)[1:2]]), colour = "tomato") +
    geom_point(aes(x = x$occupancy[localMaxima(sec$y)[max.tail.peaks[2]:max.tail.peaks[1]]],
                   y = x$frequency[localMaxima(sec$y)[max.tail.peaks[2]:max.tail.peaks[1]]]),
               colour = "tomato") +
    geom_point(aes(x = x$occupancy[localMinima(sec$y)[1:2]],
                   y = x$frequency[localMinima(sec$y)[1:2]]), colour = "royalblue") +
    geom_point(aes(x = x$occupancy[localMinima(sec$y)[min.tail.peaks[2]:min.tail.peaks[1]]],
                   y = x$frequency[localMinima(sec$y)[min.tail.peaks[2]:min.tail.peaks[1]]]),
               colour = "royalblue") +
    labs(x = "", y = "Frequency") +
    theme(axis.title = element_text(size = 9))
  
  # get another plot with data used to identify the points of maximum acceleration and minimum deceleration
  logged <- ggplot() +
    theme_pubr() +
    annotate(xmax = pred$x[localMinima(sec$y)[2]],
             xmin = min(x$occupancy),
             ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "tomato") +
    annotate(xmax = max(pred$x),
             xmin = pred$x[localMinima(sec$y)[2]],
             ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "aquamarine") +
    annotate(xmax = pred$x[localMaxima(sec$y)[max.tail.peaks[1]]],
             xmin = max(x$occupancy),
             ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "tomato") +
    geom_line(aes(x = pred$x, y = pred$y)) +
    geom_point(aes(x = x$occupancy, y = x$log.freq), alpha = 0.5, colour = "black") +
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
    theme(axis.title = element_text(size = 9))
  
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
    theme(axis.title = element_text(size = 9))
  
  p <- ggarrange(raw, logged, deriv, ncol = 3, labels = "auto")
  # add x axis title to be in the middle of two panels
  
  # Extract identified threshold
  thres.df <- data.frame(cs.flag = c("Satellite", "In-between", "Core"),
                         occup.thres.min = c(min(x$occupancy), 
                                             localMinima(sec$y)[2]+1,
                                             localMaxima(sec$y)[max.tail.peaks[1]]),
                         occup.thres.max = c(localMinima(sec$y)[2],
                                             localMaxima(sec$y)[max.tail.peaks[1]]-1,
                                             max(x$occupancy)
                                             ))
  
  out <- list(thres.df, p)
  return(out)
})

# Extract results --------------------------------------------------------------------------
# Identified thresholds
thres.df <- ldply(em.list, "[[", 1)
thres.df <- separate(thres.df, col = '.id', into = c("dataset", "replicate.merging", "rarity.cutoff", "sample.type"),
         remove = F, sep = "_")
thres.df$dataset <- factor(thres.df$dataset, levels = c("mf","mz"), labels = c("Molecular formulae",
                                                                       "Peaks"))

# Save as table
write.table(thres.df, paste0("./3_emergent.threshold/output/emergent_tresholds_",Sys.Date(), ".csv"),
            sep = ",",
            row.names = F)

thres.df <- read.csv("./3_emergent.threshold/output/emergent_tresholds_2022-01-24.csv",
                     sep = ",", stringsAsFactors = F)
# Select necessary columns
thres.df <- thres.df %>% dplyr::select(-.id)

# Re-label
colnames(thres.df) <- c("Dataset", "Replicate merging", "Rarity cutoff", "Sample type", "C-S flag",
                          "Threshold start (Occupancy)", "Threshold end (Occupancy)")
# Re-level
thres.df$`Sample type` <- factor(thres.df$`Sample type`, levels = c('water',
                                                                    'sed'),
                                 labels = c('Water','Sediment'))
knitr::kable(thres.df)

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
                       bottom = "Occupancy")
  ggsave(paste0("./3_emergent.threshold/prelim_figures/em.thres_", names(em.list)[i],
                ".png"), p, width = 25, height = 9, unit = "cm", dpi = 250)
}

#-- End: Addition by MS
