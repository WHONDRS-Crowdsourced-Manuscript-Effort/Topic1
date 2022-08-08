library (vegan)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotrix)

#Dataset obtained from rep.merged1 downloaded from 

setwd("C:/Users/erick/Desktop/WHONDRS/")                
molecules = read.csv("FTICR_commat_rep.merged1_2022-07-19.csv", header = F) #used dataset
mol = read.csv("FTICR_commat_rep.merged1_2022-01-18.csv", header = F, row.names = 1) #just to get rownames
rownames(molecules) = row.names(mol) #just to get rownames
molecules = t(molecules)

#Cross table was previously edited to take of molecules set as NA in the 
#"cs.flag.emergent_general.overlap" collumn, and delete unused collumns,
#as the R software in my PC was unable to lead with the complete dataset.
crosstab = read.csv("FTICR_crosstable_rep.merged1_all_em.thres_2022-05-05noNA.csv", header = T) 

mergedtab <- merge (crosstab,molecules,by.x = "MolForm",by.y = "ID",0,all = F)
data <- filter(mergedtab, cs.flag.emergent_sed == "Satellite" | cs.flag.emergent_water == "Satellite") #change here to "Core", and wherever else there is "Satellite" or "sat" to the analysis of Core molcules  
mergeddata <- data[,8:195]
row.names(mergeddata) <- data$MolForm

write.csv(mergeddata, file="mergeddata_sat.csv") 
#Save and upload the table without any editing. It was the only way I found to the software recognize the data
mergeddata2 = read.csv("mergeddata_sat.csv", header = T, row.names = 1)
dist_jac = vegdist(mergeddata2, method="jaccard", binary = TRUE)
beta_dist_site=betadisper(dist_jac, data$cs.flag.emergent_overlap)

#Significance tests
permutest(beta_dist_site, control=permControl (nperm=1000))
anova(beta_dist_site)
TukeyHSD(beta_dist_site)
plot(TukeyHSD(beta_dist_site), las = 1)

#extract axis' importance
valor = beta_dist_site$eig
valor2 = (beta_dist_site$eig/sum(beta_dist_site$eig))*100
head(valor2)

#extract axis' values
PCoA_score = beta_dist_site$vectors
PCoA_score = PCoA_score[,1:3]
PCoA_score[,3] <- row.names(PCoA_score)
colnames(PCoA_score) [3] = "MolForm"

map_dist = merge(PCoA_score, data,by.x = "MolForm",by.y = "MolForm")
#colnames(map_dist)[1] <-"MolForm"
write.csv(map_dist, file="pcoa_data_sat.csv")
#Save and upload the table without any editing. It was the only way I found to the software recognize the data
map_dist2 = read.csv("pcoa_data_sat.csv", header = T, row.names = 1)
attach(map_dist)

#draw the main chart
pdf("main_chart_sat.pdf", width = 9,height = 5)
ggplot(map_dist2, aes(PCoA1, PCoA2, color = cs.flag.emergent_overlap)) + 
  geom_jitter(size=1.5) +  
  #scale_shape_manual(values=c(1,2,14,15,3,4,16,17,10)) + 
  scale_color_manual(name="Sediment - Water overlap", 
                     values=c("red2","green4", "navy","dodgerblue2","magenta2","chocolate3","purple"))+ 
  #stat_ellipse(aes(group=cs.flag.emergent_overlap,color=cs.flag.emergent_overlap),size=0.75)+
  #annotate("text", x = -2, y = -2, label = "Stress = 0.139")+
  #annotate("text", x = -2, y = -2.2, label = "PERMANOVA = 0.001")+
  theme_classic()+
  xlab("PCoA 1 (14.15%)")+
  ylab("PCoA 2 (7.45%)")
dev.off()

#draw the pie charts
piedonut=data %>% group_by(cs.flag.emergent_overlap,Class) %>% summarize(n=n())
write.csv(piedonut,"pie_sat.csv")

piedonut=read.csv("pie_sat.csv", header = T, row.names = 1)
#1.Organize data (1/5)
piedonut_filtered <- filter(piedonut, cs.flag.emergent_overlap == "Global satellite")
overlap_data <-aggregate(piedonut_filtered$n,by = list(overlap = piedonut_filtered$cs.flag.emergent_overlap),FUN = sum)
class_data <- piedonut_filtered[order(piedonut_filtered$cs.flag.emergent_overlap), ]
class_labels <- paste(class_data$Class, " (", class_data$n,")",  sep = "")

#2.Define colors
class_colors <- c("lightcoral","springgreen", "seashell1","lightgoldenrod3","aquamarine2","grey60","deeppink4","olivedrab4","midnightblue") #Change color to match with classes

#3.Draw the chart
pdf("sat_pie_Global-satellite.pdf", width = 6,height = 6)
center_x <- 0.5
center_y <- 0.5
plot.new()
class_chart <-floating.pie(xpos = center_x,ypos = center_y,x = class_data$n,radius = 0.35,border = "white",col = class_colors)
pie.labels(x = center_x,y = center_y,angles = class_chart,labels = class_labels,radius = 0.36,bg = NULL,cex = 1,font = 1.5,col = "gray40")
overlap_chart <-draw.circle(0.5,0.5,radius = 0.18,border = "white",col = "white")
dev.off()

#1.Organize data (2/5)
piedonut_filtered <- filter(piedonut, cs.flag.emergent_overlap == "Sed Core - Water Sat")
overlap_data <-aggregate(piedonut_filtered$n,by = list(overlap = piedonut_filtered$cs.flag.emergent_overlap),FUN = sum)
class_data <- piedonut_filtered[order(piedonut_filtered$cs.flag.emergent_overlap), ]
class_labels <- paste(class_data$Class, " (", class_data$n,")",  sep = "")

#2.Define colors
class_colors <- c("lightcoral","springgreen", "seashell1","lightgoldenrod3","aquamarine2","grey60","deeppink4","olivedrab4","midnightblue") #Change color to match with classes

#3.Draw the chart
pdf("sat_pie_core-sat.pdf", width = 6,height = 6)
center_x <- 0.5
center_y <- 0.5
plot.new()
class_chart <-floating.pie(xpos = center_x,ypos = center_y,x = class_data$n,radius = 0.35,border = "white",col = class_colors)
pie.labels(x = center_x,y = center_y,angles = class_chart,labels = class_labels,radius = 0.36,bg = NULL,cex = 1,font = 1.5,col = "gray40")
overlap_chart <-draw.circle(0.5,0.5,radius = 0.18,border = "white",col = "white")
dev.off()

#1.Organize data (3/5)
piedonut_filtered <- filter(piedonut, cs.flag.emergent_overlap == "Sed Inbetween - Water Sat")
overlap_data <-aggregate(piedonut_filtered$n,by = list(overlap = piedonut_filtered$cs.flag.emergent_overlap),FUN = sum)
class_data <- piedonut_filtered[order(piedonut_filtered$cs.flag.emergent_overlap), ]
class_labels <- paste(class_data$Class, " (", class_data$n,")",  sep = "")

#2.Define colors
class_colors <- c("lightcoral","springgreen", "seashell1","lightgoldenrod3","aquamarine2","grey60","deeppink4","olivedrab4","midnightblue") #Change color to match with classes

#3.Draw the chart
pdf("sat_pie_between-sat.pdf", width = 6,height = 6)
center_x <- 0.5
center_y <- 0.5
plot.new()
class_chart <-floating.pie(xpos = center_x,ypos = center_y,x = class_data$n,radius = 0.35,border = "white",col = class_colors)
pie.labels(x = center_x,y = center_y,angles = class_chart,labels = class_labels,radius = 0.36,bg = NULL,cex = 1,font = 1.5,col = "gray40")
overlap_chart <-draw.circle(0.5,0.5,radius = 0.18,border = "white",col = "white")
dev.off()

#1.Organize data (4/5)
piedonut_filtered <- filter(piedonut, cs.flag.emergent_overlap == "Sed Sat - Water Core")
overlap_data <-aggregate(piedonut_filtered$n,by = list(overlap = piedonut_filtered$cs.flag.emergent_overlap),FUN = sum)
class_data <- piedonut_filtered[order(piedonut_filtered$cs.flag.emergent_overlap), ]
class_labels <- paste(class_data$Class, " (", class_data$n,")",  sep = "")

#2.Define colors
class_colors <- c("lightcoral","springgreen", "seashell1","lightgoldenrod3","aquamarine2","grey60","deeppink4","olivedrab4","midnightblue") #Change color to match with classes

#3.Draw the chart
pdf("sat_pie_sat-core.pdf", width = 6,height = 6)
center_x <- 0.5
center_y <- 0.5
plot.new()
class_chart <-floating.pie(xpos = center_x,ypos = center_y,x = class_data$n,radius = 0.35,border = "white",col = class_colors)
pie.labels(x = center_x,y = center_y,angles = class_chart,labels = class_labels,radius = 0.36,bg = NULL,cex = 1,font = 1.5,col = "gray40")
overlap_chart <-draw.circle(0.5,0.5,radius = 0.18,border = "white",col = "white")
dev.off()

#1.Organize data (5/5)
piedonut_filtered <- filter(piedonut, cs.flag.emergent_overlap == "Sed Sat - Water Inbetween")
overlap_data <-aggregate(piedonut_filtered$n,by = list(overlap = piedonut_filtered$cs.flag.emergent_overlap),FUN = sum)
class_data <- piedonut_filtered[order(piedonut_filtered$cs.flag.emergent_overlap), ]
class_labels <- paste(class_data$Class, " (", class_data$n,")",  sep = "")

#2.Define colors
class_colors <- c("lightcoral","springgreen", "seashell1","lightgoldenrod3","aquamarine2","grey60","deeppink4","olivedrab4","midnightblue") #Change color to match with classes

#3.Draw the chart
pdf("sat_pie_Sat-between.pdf", width = 6,height = 6)
center_x <- 0.5
center_y <- 0.5
plot.new()
class_chart <-floating.pie(xpos = center_x,ypos = center_y,x = class_data$n,radius = 0.35,border = "white",col = class_colors)
pie.labels(x = center_x,y = center_y,angles = class_chart,labels = class_labels,radius = 0.36,bg = NULL,cex = 1,font = 1.5,col = "gray40")
overlap_chart <-draw.circle(0.5,0.5,radius = 0.18,border = "white",col = "white")
dev.off()
