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

trans.final$Mass = formatC(trans.final$Mass,digits = 6, format = "f")

key = read.csv(paste0(input.path,"/keys/match_all.masses_to_MF.csv"))

data.int = merge(key,trans.final, by = "Mass")

# Subset the data for only unique MF 
data.int$info = duplicated(data.int$MolForm)

data.int.uni = subset(data.int,data.int$info == FALSE)

# Subset the data for only unique MF and Calculate mean values 

dat = as.data.frame(matrix(NA, nrow = nrow(data.int),ncol = ncol(data.int)))
colnames(dat) = colnames(data.int)
 for (i in 1:nrow(data.int.uni)){
   temp = subset(data.int,data.int$MolForm==data.int.uni$MolForm[i])
 if (nrow(temp)>1){
   dat$Mass[i] = temp$Mass[1]
   dat$MolForm[i] = temp$MolForm[1]
   dat$Total_all_peaks[i] = round(sum(temp$Total_all_peaks),0)
   dat$Total_all_sed [i] = round(sum(temp$Total_all_sed),0)
   dat$Total_all_water [i] = round(sum(temp$Total_all_water),0)
   dat$Total_mf_peaks[i] = round(sum(temp$Total_mf_peaks),0)
   dat$Total_mf_sed [i] = round(sum(temp$Total_mf_sed),0)
   dat$Total_mf_water [i] = round(sum(temp$Total_mf_water),0)
 } else{
   dat$Mass[i] = temp$Mass[1]
   dat$MolForm[i] = temp$MolForm[1]
   dat$Total_all_peaks[i] = temp$Total_all_peaks[1]
   dat$Total_all_sed [i] = temp$Total_all_sed[1]
   dat$Total_all_water [i] = temp$Total_all_water[1]
   dat$Total_mf_peaks[i] = temp$Total_mf_peaks[1]
   dat$Total_mf_sed [i] = temp$Total_mf_sed[1]
   dat$Total_mf_water [i] = temp$Total_mf_water[1]
 }
  
   }

dat= dat[!is.na(dat$MolForm),]
data.final = merge(data,dat, by = "MolForm")

data.final2 = data.final %>% select(Mass.y,Mass.x,MolForm,cs.flag.emergent_overlap,Total_all_peaks,Total_all_sed,Total_all_water,Total_mf_peaks,Total_mf_sed,Total_mf_water)
  
write.csv(data.final2,paste0(home.dir,"Total_number_of_transformations_per_threshold.csv"),row.names = F)
