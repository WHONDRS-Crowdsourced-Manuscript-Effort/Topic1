#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#SCRIPT EXPLORATORY ANALYSIS Topic 1: DOM CORE-SAT PAPER (WHONDRS crowdsource)
#AUTHOR: Michaela de Melo
#CONTACT:michaelaldemelo@gmail.com

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###########################################
   ########## FTICR data ##############
#Upload Dataset Masumi organized
setwd("/Volumes/GoogleDrive/My Drive/WHONDRS/Topic 1")
data_meta=read.csv("FT-ICR_meta_all_2021-08-19.csv", row.names = 1)
data_commat=read.csv("FTICR_commat_2021-08-19.csv", row.names = 1)
data_cross=read.csv("FTICR_crosstable_2021-08-19.csv", row.names = 2) #37,528 formula 


## Libraries ##
library(vegan)
library(ggplot2)
library(dplyr)
library(nord)

## Run NMDS Analysis ##
#Organize input to metaMDS
tab_merge=merge(data.frame(data_commat), data.frame(data_meta), by.x="row.names", by.y="row.names") #merge formulae with metadata
rownames(tab_merge)=tab_merge[,1]
dom = tab_merge[,2:37529] 
map = tab_merge[,37530:37598] 
rownames(map)==rownames(dom) #checking if merge works

#Run it!
ord = metaMDS(dom, dist="bray", k=2, trymax=200, autotransform=FALSE, noshare=FALSE, wascores=FALSE)
ord$stress #0.085
nmds_score = scores(ord)
mapa_dist = merge(map, nmds_score, by.x="row.names", by.y="row.names")

## Plot NMDS results ##

#By sample type
ggplot(mapa_dist, aes(NMDS1, NMDS2)) +
  geom_point(aes(fill = sample.type, shape=sample.type), size=4) + 
  scale_shape_manual(values=c(21, 22, 23, 24, 25))+ 
  scale_fill_nord("aurora")+
  theme_bw()+
  theme(aspect.ratio=1, legend.title=element_blank(), 
        legend.text=element_text(size=16), axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 16))+
  annotate("text", x = 0.3, y = 3.8, label = "Stress=0.08", size=6)
setwd("/Volumes/GoogleDrive/My Drive/WHONDRS/Topic 1")
ggsave("NMDS_DOM_sed_water.png",  width = 7, height = 7, dpi = 400) 

#We can see an outlier
ggplot(mapa_dist, aes(NMDS1, NMDS2, label=row.names(mapa_dist))) +
  geom_point( size=1) + geom_text()
#It seems to be the row 193
#remove it
tab_merge2=tab_merge[-193,]
rownames(tab_merge2) = seq(length=nrow(tab_merge2))
dom2 = tab_merge2[,2:37529] 
map2 = tab_merge2[,37530:37598] 
rownames(map2)==rownames(dom2) #checking if merge works

ord2 = metaMDS(dom2, dist="bray", k=2, trymax=200, autotransform=FALSE, noshare=FALSE, wascores=FALSE)
ord2$stress #0.085
nmds_score2 = scores(ord2)
mapa_dist2 = merge(map2, nmds_score2, by.x="row.names", by.y="row.names")

#Plot without OUTLIER
ggplot(mapa_dist2, aes(NMDS1, NMDS2)) +
  geom_point(aes(fill = sample.type, shape=sample.type), size=4) + 
  scale_shape_manual(values=c(21, 22, 23, 24, 25))+ 
  scale_fill_nord("aurora")+
  theme_bw()+
  theme(aspect.ratio=1, legend.title=element_blank(), 
        legend.text=element_text(size=16), axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 16))+
  annotate("text", x = 0.3, y = 1.3, label = "Stress=0.08", size=6)
setwd("/Volumes/GoogleDrive/My Drive/WHONDRS/Topic 1")
ggsave("NMDS_DOM_sed_water_without_outlier.png",  width = 7, height = 7, dpi = 400) 

# There is a big difference between DOM composition from water and sediment, so better separate them 

########################################################################
          ###########   WATER DOM  #############
water=tab_merge2 [tab_merge2$sample.type %in% c("SW"), ] #select only water samples 

#NMDS water samples: organizing data and running metaMDS
water2=water[,-1] #remove first column
rownames(water2) = seq(length=nrow(water2)) #reset rows order
dom_water = water2[,2:37528] 
map_water = water2[,37529:37597] 
rownames(map_water)==rownames(dom_water) #checking if merging worked
ord_water = metaMDS(dom_water, dist="bray", k=2, trymax=200, autotransform=FALSE, noshare=FALSE, wascores=FALSE)
ord_water$stress #0.15
nmds_score_water = scores(ord_water)
mapa_dist_water = merge(map_water, nmds_score_water, by.x="row.names", by.y="row.names")

#By stream order
ggplot(mapa_dist_water, aes(NMDS1, NMDS2)) +
  geom_point(aes(color =as.character(Stream_Order)), size=4) + 
  scale_shape_manual(values=c(21, 22, 23, 24, 25))+ 
  scale_color_nord("baie_mouton")+ 
  theme_bw()+
  theme(aspect.ratio=1, legend.title=element_blank(), 
        legend.text=element_text(size=16), axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 16))
  #annotate("text", x = 0.3, y = 3.8, label = "Stress=0.15", size=6)
setwd("/Volumes/GoogleDrive/My Drive/WHONDRS/Topic 1")
ggsave("NMDS_Water_streamorder.png",  width = 7, height = 7, dpi = 400) 
#It seems river with smaller orders spread more in the 2-D space, but >5 order cluster more

#Latitude--> why "Latitude_dec.deg" is empty?? 
ggplot(mapa_dist_water, aes(NMDS1, NMDS2, colour =Gauge_Latitude_dec.deg)) +
  geom_point( size=4) + 
  scale_colour_gradientn(colours = terrain.colors(5))+
  theme_bw()+
  theme(aspect.ratio=1, legend.title=element_text(size=18), 
        legend.text=element_text(size=16), axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 16))
setwd("/Volumes/GoogleDrive/My Drive/WHONDRS/Topic 1")
ggsave("NMDS_Water_Latitude.png",  width = 7, height = 7, dpi = 400) 


#Longitude--> why "Latitude_dec.deg" is empty?? 
ggplot(mapa_dist_water, aes(NMDS1, NMDS2, colour =Gauge_Longitude_dec.deg)) +
  geom_point( size=4) + 
  scale_colour_gradientn(colours = terrain.colors(5))+
  theme_bw()+
  theme(aspect.ratio=1, legend.title=element_text(size=18), 
        legend.text=element_text(size=16), axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 16))


#Contamination
ggplot(mapa_dist_water, aes(NMDS1, NMDS2, colour =Contamination.Source.Upstream)) +
  geom_point( size=4) + 
  scale_colour_nord("baie_mouton")+
  theme_bw()+
  theme(aspect.ratio=1, legend.title=element_text(size=18), 
        legend.text=element_text(size=16), axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 16))


#Temperature 
ggplot(mapa_dist_water, aes(NMDS1, NMDS2, colour =as.numeric(Temp_degC))) +
  geom_point( size=4) + 
  scale_colour_gradientn(colours = terrain.colors(8))+
  theme_bw()+
  theme(aspect.ratio=1, legend.title=element_text(size=18), 
        legend.text=element_text(size=16), axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 16))
ggsave("NMDS_Water_Temperature.png",  width = 7, height = 7, dpi = 400)




###########   SEDIMENT DOM  #############
sed=tab_merge2 [tab_merge2$sample.type %in% c("SED"), ]

#NMDS sediment samples: organizing and running it
sed=sed[,-1]
dom_sed = sed[,2:37528] 
map_sed = sed[,37529:37597] 
rownames(map_sed)==rownames(dom_sed) #checking if merging worked
ord_sed = metaMDS(dom_sed, dist="bray", k=2, trymax=200, autotransform=FALSE, noshare=FALSE, wascores=FALSE)
ord_sed$stress #0.17
nmds_score_sed = scores(ord_sed)
mapa_dist_sed = merge(map_sed, nmds_score_sed, by.x="row.names", by.y="row.names")


#By stream order
ggplot(mapa_dist_sed, aes(NMDS1, NMDS2)) +
  geom_point(aes(color =as.character(Stream_Order)), size=4) + 
  scale_shape_manual(values=c(21, 22, 23, 24, 25))+ 
  scale_color_nord("baie_mouton")+
  theme_bw()+
  theme(aspect.ratio=1, legend.title=element_blank(), 
        legend.text=element_text(size=16), axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 16))
#annotate("text", x = 0.3, y = 3.8, label = "Stress=0.15", size=6)
setwd("/Volumes/GoogleDrive/My Drive/WHONDRS/Topic 1")
ggsave("NMDS_Sediment_streamorder.png",  width = 7, height = 7, dpi = 400) 
#It seems river with smaller orders spread more in the 2-D space, but >5 order cluster more

#Latitude
ggplot(mapa_dist_sed, aes(NMDS1, NMDS2, colour =Latitude_dec.deg)) +
  geom_point( size=4) + 
  scale_colour_gradientn(colours = terrain.colors(5))+
  theme_bw()+
  theme(aspect.ratio=1, legend.title=element_text(size=18), 
        legend.text=element_text(size=16), axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 16))
setwd("/Volumes/GoogleDrive/My Drive/WHONDRS/Topic 1")
ggsave("NMDS_Sediment_Latitude.png",  width = 7, height = 7, dpi = 400) 
#NMDS2 divided low  (green) to high  (white-ish) latitude

#Longitude
ggplot(mapa_dist_sed, aes(NMDS1, NMDS2, colour =Longitude_dec.deg)) +
  geom_point( size=4) + 
  scale_colour_gradientn(colours = terrain.colors(5))+
  theme_bw()+
  theme(aspect.ratio=1, legend.title=element_text(size=18), 
        legend.text=element_text(size=16), axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 16))
setwd("/Volumes/GoogleDrive/My Drive/WHONDRS/Topic 1")
ggsave("NMDS_Sediment_Longitude.png",  width = 7, height = 7, dpi = 400) 

#Contamination
ggplot(mapa_dist_sed, aes(NMDS1, NMDS2, colour =Contamination.Source.Upstream)) +
  geom_point( size=4) + 
  scale_colour_nord("baie_mouton")+
  theme_bw()+
  theme(aspect.ratio=1, legend.title=element_text(size=18), 
        legend.text=element_text(size=16), axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 16))


#exploring data
ggplot(mapa_dist_sed, aes(NMDS1, NMDS2, colour =Macrophyte.Coverage)) +
  geom_point( size=4) + 
  scale_colour_nord("baie_mouton")+
  theme_bw()+
  theme(aspect.ratio=1, legend.title=element_text(size=18), 
        legend.text=element_text(size=6), axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 16))


ggplot(mapa_dist_sed, aes(NMDS1, NMDS2, colour =as.numeric(Water.Column.Height_cm))) +
  geom_point( size=4) + 
  scale_colour_gradientn(colours = terrain.colors(5))+
  theme_bw()+
  theme(aspect.ratio=1, legend.title=element_text(size=18), 
        legend.text=element_text(size=16), axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 16))

ggplot(mapa_dist_sed, aes(NMDS1, NMDS2, colour =location.id)) +
  geom_point( size=4) + 
  scale_colour_nord("baie_mouton")+
  theme_bw()+
  theme(aspect.ratio=1, legend.title=element_text(size=18), 
        legend.text=element_text(size=6), axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 16))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#####################################################
      ##  CORE AND SAT ###
####################################################@
#Separation using presence and absence
setwd("/Volumes/GoogleDrive/My Drive/WHONDRS/Topic 1")
data_commat=read.csv("FTICR_commat_2021-08-19.csv", row.names = 1)
data_cross=read.csv("FTICR_crosstable_2021-08-19.csv", row.names = 2) #37,528 formula 
data_meta=read.csv("FT-ICR_meta_all_2021-08-19.csv", row.names = 1)


#Merge with metadata
tab_core_sat=merge(data.frame(data_commat), data.frame(data_meta), by.x="row.names", by.y="row.names")
rownames(tab_core_sat)=tab_core_sat[,1]

#2 tables: water and sediment
water=tab_core_sat [tab_core_sat$sample.type %in% c("SW"), ] #265 obs means 265 water samples
sed=tab_core_sat [tab_core_sat$sample.type %in% c("SED"), ] #239 obs means 239 water samples

#remove metadata
water2=water[,-1]
dom_water = water2[,2:37527] 

sed2=sed[,-1]
dom_sed = sed2[,2:37528] 


#Check in how many samples molecules were found
sum_water=data.frame(colSums(dom_water))
sat_water=filter(sum_water, colSums.dom_water.==1) #7,782 obs (formulae) were present in only a water site (sparce distribution)
#we can define different thresholds here for core
x=filter(sum_water, colSums.dom_water.>0) #Just to see how many formulae there is in water samples = 23,568

sum_sed=data.frame(colSums(dom_sed))
sat_sed=filter(sum_sed, colSums.dom_sed.==1) #9,476 obs (formulae) from  were present in only a water site (sparce distribution)
y=filter(sum_sed, colSums.dom_sed.>0) #Just to see how many formulae there is in water samples = 22,830


#Occupancy-Frequency WATER#
df_water=data_frame(table(sum_water))
table(sum_water) #for visualization of results
df_water2=df_water[-1,]#remove First cell, which is number of formulae = 0. Means we found 7,464 formulae were absent in water samples (compared to the merged table with sediment samples)
df_water2$occupancy=seq(1:265) #265 is the number of water samples
colnames(df_water2)[1]="frequency"

#Plots
ggplot(df_water2, aes(x=occupancy, y=log(frequency)))+ geom_point()
barplot(frequency ~ occupancy, data = df_water2)
barplot(log(frequency) ~ occupancy, data = df_water2, xlab="Occupancy (number of sites)", ylab="Log Frequency (number of DOM formulae)", main="Surface water samples")
#Saved as "Water_occup_freq_bars.png
#MOStest(df_water2$frequency, df_water2$occupancy)


#Occupancy-Frequency  SEDIMENT
df_sed=data_frame(table(sum_sed))
table(sum_sed) #for visualization of results
df_sed2=df_sed[-1,] #remove First cell, which is number of formulae = 0. Means we found 9,575 formulae were absent in water samples (compared to the merged table with sediment samples)
df_sed2$occupancy=seq(1:239) #239 is the number of sediment samples
colnames(df_sed2)[1]="frequency"

#Plots
ggplot(df_sed2, aes(x=occupancy, y=log(frequency)))+ geom_point()
barplot(frequency ~ occupancy, data = df_sed2)
barplot(log(frequency) ~ occupancy, data = df_sed2, xlab="Occupancy (number of sites)", ylab="Log Frequency (number of DOM formulae)", main="Sediment samples")
#Saved as "Sediment_occup_freq_bars.png




#######################################################################
#Date: 10 January 2022

#Updated data

#Add data

#Path:
#setwd("C:/Users/micha/OneDrive/Documentos/GitHub/whondrs-coresat-topic1/Dataset/")

#%%%%%#
# ALL #
#%%%%%#
# data_commat=read.csv("FTICR_commat_rep.merged_all_2021-11-03.csv", row.names = 1)
# data_cross=read.csv("FTICR_cross.table_rep.merged_all_2021-11-03.csv", row.names = 2) #37,528 formula 
# data_meta=read.csv("FTICR_meta_all_2021-09-29.csv")

# Adding relative paths (-MS)
data_commat=read.csv("./1_data.cleaning/output/FTICR_commat_rep.merged_all_2021-11-03.csv", row.names = 1)
data_cross=read.csv("./1_data.cleaning/output/FTICR_cross.table_rep.merged_all_2021-11-03.csv", row.names = 2) #37,528 formula 
data_meta=read.csv("./1_data.cleaning/output/FTICR_meta_all_2021-09-29.csv")

#Script to find inflection based on second derivative:
#https://cran.r-project.org/web/packages/inflection/vignettes/inflectionMissionImpossible.html

#2 tables: water and sediment
rownames(data_commat)
sed=data_commat [1:93, ] #93 samples
rownames(sed)
water=data_commat [94:188, ]
rownames(water) #95 samples


#Occupancy-Frequency  WATER
# how often were MF found across water samples?
sum_water=data.frame(colSums(water)) # sum of presence-absence = number of sites present
df_water=data_frame(table(sum_water)) # transform to counts of MF that were found in 1, 2, 3 sites and so on...
table(sum_water) #for visualization of results
df_water2=df_water[-1,] #remove First cell, which is number of formulae = 0. Means we found 9,575 formulae were absent in water samples (compared to the merged table with sediment samples)
df_water2$occupancy=seq(1:95) #95 is the number of water samples
colnames(df_water2)[1]="frequency"

#SEDIMENT
# how often were MF found across sediment samples?
# Sequence of code is same as water samples.
sum_sed=data.frame(colSums(sed))
df_sed=data_frame(table(sum_sed))
table(sum_sed) #for visualization of results
df_sed2=df_sed[-1,] #remove First cell, which is number of formulae = 0. Means we found 9,575 formulae were absent in water samples (compared to the merged table with sediment samples)
df_sed2$occupancy=seq(1:93) #93 is the number of water samples
colnames(df_sed2)[1]="frequency"


#Plots
ggplot(df_water2, aes(x=occupancy, y=log(frequency)))+ geom_point() + labs(x="occupancy", title="Water, whole dataset - all")
#setwd("C:/Users/micha/OneDrive/Documentos/GitHub/whondrs-coresat-topic1/Analyses/Michaela/")
ggsave("3_emergent.threshold/prelim_figures/FreqXoccupancy_Water_All.png", dpi=250)


ggplot(df_sed2, aes(x=occupancy, y=log(frequency)))+ geom_point() + labs(x="occupancy", title="Sediment, whole dataset - all")
ggsave("3_emergent.threshold/prelim_figures/FreqXoccupancy_Sed_All.png", dpi=250)


# Re-do the whole process with the other two rarity cutoffs
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55#
# Rar1 - removes 36.6% of MF#
#%%%%%%%%%%%%%%%%%%%%%%%%%%#

#Path:
#setwd("C:/Users/micha/OneDrive/Documentos/GitHub/whondrs-coresat-topic1/Dataset/")

data_commat_rar1=read.csv("./1_data.cleaning/output/FTICR_commat_rep.merged_rar1_2021-11-03.csv", row.names = 1)
data_cross_rar1=read.csv("./1_data.cleaning/output/FTICR_cross.table_rep.merged_rar1_2021-11-03.csv", row.names = 2) #37,528 formula 
data_meta=read.csv("./1_data.cleaning/output/FTICR_meta_all_2021-09-29.csv")


#2 tables: water and sediment
rownames(data_commat_rar1)
sed_rar1=data_commat_rar1 [1:93, ] #93 samples
rownames(sed_rar1)
water_rar1=data_commat_rar1 [94:188, ]
rownames(water_rar1) #95 samples



#Occupancy-Frequency  WATER
sum_water_rar1=data.frame(colSums(water_rar1))
df_water_rar1=data_frame(table(sum_water_rar1))
table(sum_water_rar1) #for visualization of results
df_water2_rar1=df_water_rar1[-1,] #remove First cell, which is number of formulae = 0. Means we found 9,575 formulae were absent in water samples (compared to the merged table with sediment samples)
df_water2_rar1$occupancy=seq(1:95) #95 is the number of water samples
colnames(df_water2_rar1)[1]="frequency"

#SEDIMENT
sum_sed_rar1=data.frame(colSums(sed_rar1))
df_sed_rar1=data_frame(table(sum_sed_rar1))
table(sum_sed_rar1) #for visualization of results
df_sed2_rar1=df_sed_rar1[-1,] #remove First cell, which is number of formulae = 0. Means we found 9,575 formulae were absent in water samples (compared to the merged table with sediment samples)
df_sed2_rar1$occupancy=seq(1:93) #93 is the number of water samples
colnames(df_sed2_rar1)[1]="frequency"


#Plots
ggplot(df_water2_rar1, aes(x=occupancy, y=log(frequency)))+ geom_point() + labs(x="occupancy", title="Water, Remove present 1 site")
#setwd("C:/Users/micha/OneDrive/Documentos/GitHub/whondrs-coresat-topic1/Analyses/Michaela/")
ggsave("./3_emergent.threshold/prelim_figures/FreqXoccupancy_Water_Rar1.png", dpi=250)

ggplot(df_sed2_rar1, aes(x=occupancy, y=log(frequency)))+ geom_point() + labs(x="occupancy", title="Sediment, Remove present 1 site")
ggsave("./3_emergent.threshold/prelim_figures/FreqXoccupancy_Sed_Rar1.png", dpi=250)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55#
# Rar2 - removes 48.6% of MF#
#%%%%%%%%%%%%%%%%%%%%%%%%%%#

#Path:
#setwd("C:/Users/micha/OneDrive/Documentos/GitHub/whondrs-coresat-topic1/Dataset/")
data_commat_rar2=read.csv("./1_data.cleaning/output/FTICR_commat_rep.merged_rar2_2021-11-03.csv", row.names = 1)
data_cross_rar2=read.csv("./1_data.cleaning/output/FTICR_cross.table_rep.merged_rar2_2021-11-03.csv", row.names = 2) #37,528 formula 
data_meta=read.csv("./1_data.cleaning/output/FTICR_meta_all_2021-09-29.csv")


#2 tables: water and sediment
rownames(data_commat_rar2)
sed_rar2=data_commat_rar2 [1:93, ] #93 samples
rownames(sed_rar2)
water_rar2=data_commat_rar2 [94:188, ]
rownames(water_rar2) #95 samples


#Occupancy-Frequency  WATER
sum_water_rar2=data.frame(colSums(water_rar2))
df_water_rar2=data_frame(table(sum_water_rar2))
table(sum_water_rar2) #for visualization of results
df_water2_rar2=df_water_rar2[-1,] #remove First cell, which is number of formulae = 0. Means we found 9,575 formulae were absent in water samples (compared to the merged table with sediment samples)
df_water2_rar2$occupancy=seq(1:95) #95 is the number of water samples
colnames(df_water2_rar2)[1]="frequency"

#SEDIMENT
sum_sed_rar2=data.frame(colSums(sed_rar2))
df_sed_rar2=data_frame(table(sum_sed_rar2))
table(sum_sed_rar2) #for visualization of results
df_sed2_rar2=df_sed_rar2[-1,] #remove First cell, which is number of formulae = 0. Means we found 9,575 formulae were absent in water samples (compared to the merged table with sediment samples)
df_sed2_rar2$occupancy=seq(1:93) #93 is the number of water samples
colnames(df_sed2_rar2)[1]="frequency"


#Plots
ggplot(df_water2_rar2, aes(x=occupancy, y=log(frequency)))+ geom_point() + labs(x="occupancy", title="Water, Remove present 2 sites")
#setwd("C:/Users/micha/OneDrive/Documentos/GitHub/whondrs-coresat-topic1/Analyses/Michaela/")
ggsave("./3_emergent.threshold/prelim_figures/FreqXoccupancy_Water_Rar2.png", dpi=250)


ggplot(df_sed2_rar2, aes(x=occupancy, y=log(frequency)))+ geom_point() + labs(x="occupancy", title="Sediment, Remove present 2 sites")
ggsave("./3_emergent.threshold/prelim_figures/FreqXoccupancy_Sed_Rar2.png", dpi=250)


#-- Start: Addition by MS
# Find points of maximum deceleration to identify emergent thresholds -----------------------
library(data.table)
library(plyr)
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
data.list <- list(df_water2, df_sed2,
                  df_water2_rar1, df_sed2_rar1,
                  df_water2_rar2, df_sed2_rar2)
names(data.list) <- c("SW_all", "SED_all",
                      "SW_rar1", "SED_rar1",
                      "SW_rar2", "SED_rar2")

# Example data set to write function
#x<- df_water2_rar1

# We're using smooth spline here instead of defining a function and getting a derivative from the function
# as described in https://cran.r-project.org/web/packages/inflection/vignettes/inflectionMissionImpossible.html
# as it is hard to come up with a function describing the observed pattern
# not exponential, not logarithmic etc

em.list <- llply(data.list, function(x){
  # take log
  x$log.freq <- log(x$frequency)
  # make a smooth curve
  spl <- smooth.spline(x$occupancy, x$log.freq, spar = 0.5)
  # predict to get fit
  pred <- predict(spl)
  # get second derivative
  sec <- predict(spl, deriv = 2) 
  
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
    geom_line(aes(x = x$occupancy, y = x$frequency)) +
    geom_point(aes(x = x$occupancy[localMaxima(sec$y)[1:2]],
                   y = x$frequency[localMaxima(sec$y)[1:2]]), colour = "tomato") +
    geom_point(aes(x = x$occupancy[localMinima(sec$y)[1:2]],
                   y = x$frequency[localMinima(sec$y)[1:2]]), colour = "royalblue") +
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
    geom_line(aes(x = pred$x, y = pred$y)) +
    geom_point(aes(x = x$occupancy, y = x$log.freq), alpha = 0.5, colour = "black") +
    geom_point(aes(x = pred$x[localMaxima(sec$y)[1:2]],
                   y = pred$y[localMaxima(sec$y)[1:2]]), colour = "tomato") +
    geom_point(aes(x = pred$x[localMinima(sec$y)[1:2]],
                   y = pred$y[localMinima(sec$y)[1:2]]), colour = "royalblue") +
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
    labs(x = "", y = expression(paste("Acceleration x10"^3, " (2"^"nd", " derivative)"))) +
    theme(axis.title = element_text(size = 9))
  
  p <- ggarrange(raw, logged, deriv, ncol = 3, labels = "auto")
  # add x axis title to be in the middle of two panels
  
  # Extract identified threshold
  thres.df <- data.frame(cs.flag = c("Satellite", "Core"),
                    occup.thres.min = c(min(x$occupancy), localMinima(sec$y)[2]+1),
                    occup.thres.max = c(localMinima(sec$y)[2], max(x$occupancy)))
  
  out <- list(thres.df, p)
  return(out)
})

# Extract results --------------------------------------------------------------------------
# Identified thresholds
thres.df <- ldply(em.list, "[[", 1)
# Save as table
write.table(thres.df, "./3_emergent.threshold/output/emergent_tresholds.csv", sep = ",",
            row.names = F)

# Save generated plots
title <- c("Surface water - all","Sediment - all",
           "Surface water - rar1","Sediment - rar1",
           "Surface water - rar2","Sediment - rar2")

for(i in 1:length(em.list)){
  p <- annotate_figure(em.list[[i]][[2]], top = text_grob(title[i]),
                  bottom = "Occupancy")
  ggsave(paste0("./3_emergent.threshold/prelim_figures/em.thres_", names(em.list)[i],
         ".png"), p, width = 25, height = 9, unit = "cm", dpi = 250)
}

#-- End: Addition by MS
