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
library(nord)
library(dplyr)

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