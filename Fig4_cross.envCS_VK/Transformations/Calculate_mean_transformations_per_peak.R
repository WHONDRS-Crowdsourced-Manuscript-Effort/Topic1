rm(list=ls());graphics.off()
# Loading libraries
library (vegan); library(tidyverse)
library("dplyr")                                    
library("plyr")                                     
library("readr")                                    

options(digits=10)# to ensure we bring in all the decimals
# Setting working directories
home.dir = "C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/Topic1/Fig4_cross.envCS_VK/Transformations/"


mf.path = "MF_assigned_merged_data/Transformations per Peak/S19S_MF_assigned_merged_data/"

input.path = "C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/Topic1/4_gather.thresholds/"

# Read in data  
file = list.files(path = input.path,pattern = "all_em.thres_2022-05-05.csv", full.names = T) # File to use "FTICR_crosstable_rep.merged1_all_em.thres_2022-05-05.csv"
data = read.csv(file)

# Transformations
# Read in transformations per peak and calculate the mean,median and std number of transformations in each peak

trans.mf.peaks <- list.files(path = paste0(home.dir,mf.path), pattern = "*.csv",full.names = TRUE) %>%
  lapply(read_csv) %>% 
  reduce(full_join, by = "peak")%>%
  select(-starts_with("sam"))

stats.mf.peaks = as.data.frame(matrix(NA,ncol = 4, nrow = nrow(trans.mf.peaks)))
colnames(stats.mf.peaks) = c("peak","Mean","Median","SD")

for (i in 1:nrow(trans.mf.peaks)){
  stats.mf.peaks$peak[i] = trans.mf.peaks$peak[i]
  stats.mf.peaks$Mean[i] = mean(na.omit(t(trans.mf.peaks[i,grep("nu",colnames(trans.mf.peaks))])))
  stats.mf.peaks$Median[i] = median(na.omit(t(trans.mf.peaks[i,grep("nu",colnames(trans.mf.peaks))])))
  stats.mf.peaks$SD[i] = sd(na.omit(t(trans.mf.peaks[i,grep("nu",colnames(trans.mf.peaks))])))
}

# Only sediment transformations
trans.mf.sed <- list.files(path = paste0(home.dir,mf.path), pattern = "SED",full.names = TRUE) %>%
  lapply(read_csv) %>% 
  reduce(full_join, by = "peak")%>%
  select(-starts_with("sam")) 

stats.sed.peaks = as.data.frame(matrix(NA,ncol = 4, nrow = nrow(trans.mf.sed)))
colnames(stats.sed.peaks) = c("peak","SED_Mean","SED_Median","SED_SD")

for (i in 1:nrow(trans.mf.sed)){
  stats.sed.peaks$peak[i] = trans.mf.sed$peak[i]
  stats.sed.peaks$SED_Mean[i] = mean(na.omit(t(trans.mf.sed[i,grep("nu",colnames(trans.mf.sed))])))
  stats.sed.peaks$SED_Median[i] = median(na.omit(t(trans.mf.sed[i,grep("nu",colnames(trans.mf.sed))])))
  stats.sed.peaks$SED_SD[i] = sd(na.omit(t(trans.mf.sed[i,grep("nu",colnames(trans.mf.sed))])))
}

# Finally water
trans.mf.water <- list.files(path = paste0(home.dir,mf.path), pattern = "SW",full.names = TRUE) %>%
  lapply(read_csv) %>% 
  reduce(full_join, by = "peak")%>%
  select(-starts_with("sam")) 


stats.water.peaks = as.data.frame(matrix(NA,ncol = 4, nrow = nrow(trans.mf.water)))
colnames(stats.water.peaks) = c("peak","SW_Mean","SW_Median","SW_SD")

for (i in 1:nrow(trans.mf.water)){
  stats.water.peaks$peak[i] = trans.mf.water$peak[i]
  stats.water.peaks$SW_Mean[i] = mean(na.omit(t(trans.mf.water[i,grep("nu",colnames(trans.mf.water))])))
  stats.water.peaks$SW_Median[i] = median(na.omit(t(trans.mf.water[i,grep("nu",colnames(trans.mf.water))])))
  stats.water.peaks$SW_SD[i] = sd(na.omit(t(trans.mf.water[i,grep("nu",colnames(trans.mf.water))])))
}

# Merge all of the datasets by peaks
trans.final = join_all(list(stats.mf.peaks,stats.sed.peaks,stats.water.peaks), by = "peak")


names(trans.final)[1] = "Mass"
trans.final2 = trans.final
trans.final2[2:10] = round(trans.final2[2:10],2)
write.csv(trans.final2,paste0(home.dir,"Mean_number_of_transformations_per_peak.csv"),row.names = F)

trans.final$Mass = formatC(trans.final$Mass,digits = 6, format = "f")

key = read.csv(paste0(input.path,"/keys/match_all.masses_to_MF.csv"))

data.int = merge(key,trans.final, by = "Mass")

# Subset the data for duplicated MF 
data.int$info = duplicated(data.int$MolForm)

data.int.dup = subset(data.int,data.int$info == TRUE)
# Unique MFs 
data.int.uni = subset(data.int,data.int$info == FALSE)
# Then we will merge the avg transformations per peak with the molecular formula and classifications were using for the paper and we will check for duplicate MFs with different masses, If we encounter these we will sum them such that specific MF can have both  avg transformations representes

# Subset the data for only unique MF and Calculate mean values 

dat = as.data.frame(matrix(NA, nrow = nrow(data.int.dup),ncol = ncol(data.int.dup)))
colnames(dat) = colnames(data.int.dup)
 for (i in 1:nrow(data.int.dup)){
   dat$Mass[i] = "Double"
   dat$MolForm[i] = data.int.dup$MolForm[i]
   dat$Mean[i] = sum(na.omit(data.int.dup$Mean[i]),na.omit(data.int.uni$Mean[grep(data.int.dup$MolForm[i],data.int.uni$MolForm)]))
   dat$Median[i] = sum(na.omit(data.int.dup$Median[i]),na.omit(data.int.uni$Median[grep(data.int.dup$MolForm[i],data.int.uni$MolForm)]))
  
   dat$SED_Mean[i] = sum(na.omit(data.int.dup$SED_Mean[i]),na.omit(data.int.uni$SED_Mean[grep(data.int.dup$MolForm[i],data.int.uni$MolForm)]))
   dat$SED_Median[i] = sum(na.omit(data.int.dup$SED_Median[i]),na.omit(data.int.uni$SED_Median[grep(data.int.dup$MolForm[i],data.int.uni$MolForm)]) )
   
   dat$SW_Mean[i] = sum(na.omit(data.int.dup$SW_Mean[i]),na.omit(data.int.uni$SW_Mean[grep(data.int.dup$MolForm[i],data.int.uni$MolForm)]))
   dat$SW_Median[i] = sum(na.omit(data.int.dup$SW_Median[i]),na.omit(data.int.uni$SW_Median[grep(data.int.dup$MolForm[i],data.int.uni$MolForm)]))
   
   dat$SD[i] = sum(na.omit(data.int.dup$SD[i]),na.omit(data.int.uni$SD[grep(data.int.dup$MolForm[i],data.int.uni$MolForm)]))
   dat$SED_SD[i] = sum(na.omit(data.int.dup$SED_SD[i]),na.omit(data.int.uni$SED_SD[grep(data.int.dup$MolForm[i],data.int.uni$MolForm)]))
   dat$SW_SD[i] = sum(na.omit(data.int.dup$SW_SD[i]),na.omit(data.int.uni$SW_SD[grep(data.int.dup$MolForm[i],data.int.uni$MolForm)]))
 }


data.final = rbind(data.int.uni,dat)
key2 = read.csv(paste0(input.path,"/FTICR_crosstable_rep.merged1_all_em.thres_2022-05-05.csv"))
key2 = key2 %>% dplyr::select(MolForm,cs.flag.emergent_overlap)

data.final2 = merge(key2,data.final, by = "MolForm")
data.final2[3:11] = round(data.final2[3:11],2)

write.csv(data.final2,paste0(home.dir,"Total_number_of_transformations_per_threshold.csv"),row.names = F)
