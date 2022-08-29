rm(list=ls());graphics.off()
# Loading libraries
library (vegan); library(tidyverse)
library("dplyr")                                    
library("plyr")                                     
library("readr")                                    

options(digits=10)# to ensure we bring in all the decimals
# Setting working directories
home.dir = "C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/Topic1/Fig4_cross.envCS_VK/Transformations/"

all.peaks.path = "All_peaks/Transformations per Peak/S19S_all_peaks/"
mf.path = "only_MF_assigned/Transformations per Peak/S19S_MF_assigned/"

input.path = "C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/Topic1/4_gather.thresholds/"

# Read in data  
file = list.files(path = input.path,pattern = "all_em.thres_2022-05-05.csv", full.names = T) # File to use "FTICR_crosstable_rep.merged1_all_em.thres_2022-05-05.csv"
data = read.csv(file)

# Transformations
# Read in transformations per peak and calculate the total number of transformations in each peak

trans.all.peaks <- list.files(path = paste0(home.dir,all.peaks.path), pattern = "*.csv",full.names = TRUE) %>%
  lapply(read_csv) %>% 
  reduce(full_join, by = "peak")%>%
  select(-starts_with("sam")) %>%
  dplyr::mutate(Total_all_peaks = rowSums(dplyr::across(where(is.numeric)&starts_with("nu")), na.rm = TRUE)) %>%
  select(peak,Total_all_peaks)

trans.all.sed <- list.files(path = paste0(home.dir,all.peaks.path), pattern = "Sed_Field_ICR",full.names = TRUE) %>%
  lapply(read_csv) %>% 
  reduce(full_join, by = "peak")%>%
  select(-starts_with("sam")) %>%
  dplyr::mutate(Total_all_sed = rowSums(dplyr::across(where(is.numeric)&starts_with("nu")), na.rm = TRUE)) %>%
  select(peak,Total_all_sed)

trans.all.water <- list.files(path = paste0(home.dir,all.peaks.path), pattern = "*.csv",full.names = TRUE)[!grepl(pattern = "Sed_Field_ICR",list.files(path = paste0(home.dir,all.peaks.path), pattern = "*.csv",full.names = TRUE))] %>%
  lapply(read_csv) %>% 
  reduce(full_join, by = "peak")%>%
  select(-starts_with("sam")) %>%
  dplyr::mutate(Total_all_water= rowSums(dplyr::across(where(is.numeric)&starts_with("nu")), na.rm = TRUE)) %>%
  select(peak,Total_all_water)

####### Repeating steps for Trans with MF assigned #####
trans.mf.peaks <- list.files(path = paste0(home.dir,mf.path), pattern = "*.csv",full.names = TRUE) %>%
  lapply(read_csv) %>% 
  reduce(full_join, by = "peak")%>%
  select(-starts_with("sam")) %>%
  dplyr::mutate(Total_mf_peaks = rowSums(dplyr::across(where(is.numeric)&starts_with("nu")), na.rm = TRUE)) %>%
  select(peak,Total_mf_peaks)


trans.mf.sed <- list.files(path = paste0(home.dir,mf.path), pattern = "Sed_Field_ICR",full.names = TRUE) %>%
  lapply(read_csv) %>% 
  reduce(full_join, by = "peak")%>%
  select(-starts_with("sam")) %>%
  dplyr::mutate(Total_mf_sed = rowSums(dplyr::across(where(is.numeric)&starts_with("nu")), na.rm = TRUE)) %>%
  select(peak,Total_mf_sed)

trans.mf.water <- list.files(path = paste0(home.dir,mf.path), pattern = "*.csv",full.names = TRUE)[!grepl(pattern = "Sed_Field_ICR",list.files(path = paste0(home.dir,all.peaks.path), pattern = "*.csv",full.names = TRUE))] %>%
  lapply(read_csv) %>% 
  reduce(full_join, by = "peak")%>%
  select(-starts_with("sam")) %>%
  dplyr::mutate(Total_mf_water= rowSums(dplyr::across(where(is.numeric)&starts_with("nu")), na.rm = TRUE)) %>%
  select(peak,Total_mf_water)

# Merge all of the datasets by peaks
trans.all.final = join_all(list(trans.all.peaks,trans.all.sed,trans.all.water), by = "peak")

trans.mf.final = join_all(list(trans.mf.peaks,trans.mf.sed,trans.mf.water), by = "peak")

trans.all.final$peak = formatC(trans.all.final$peak,digits = 7, format = "f")

trans.final = merge(trans.all.final,trans.mf.final)
names(trans.final)[1] = "Mass"

write.csv(trans.final,paste0(home.dir,"Total_number_of_transformations_per_peak.csv"),row.names = F)
# In order to merge with the rest of the dataset the masses and the digits across them need to be consistent
# Changing digits on the fly to be able to make a figur. We would need to cut to 4 digits but this is not accurate enough for a final figure

trans.final$Mass = formatC(trans.final$Mass,digits = 3, format = "f")
data$Mass = formatC(data$Mass,digits = 10, format = "f")

data.final = merge(data,trans.final, by = "Mass")
