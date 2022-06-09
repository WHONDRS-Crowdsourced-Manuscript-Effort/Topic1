library (vegan)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotrix)
#library(ggforce)
################################
###Upload and data processing###
################################
molecules = read.csv("/../FTICR_commat_rep.merged1_2022-01-18.csv", row.names = 1, header = F)
molecules = t(molecules)
crosstab = read.csv("/../FTICR_crosstable_rep.merged1_all_em.thres_2022-03-18.csv", header = T, row.names = 1) %>% 
  unite(col = "cs.flag.emergent_overlap", cs.flag.emergent_sed:cs.flag.emergent_water ,sep = "-", remove = F)

mergedtab <- merge (crosstab,molecules,by.x = "MolForm",by.y = "ID",0,all = F)
data <- filter(mergedtab, cs.flag.emergent_sed == "Core" | cs.flag.emergent_water == "Core")
mergeddata <- data[,52:239]
row.names(mergeddata) <- data$MolForm
write.csv(mergeddata, file="/../mergeddata.csv") #It was the only way I found to force R to understand the dataset properly
mergeddata2 = read.csv("/../mergeddata.csv", header = T, row.names = 1)
dist_jac = vegdist(mergeddata2, method="jaccard", binary = TRUE)

#############################
###Betadispersion analysis###
#############################
beta_dist_site=betadisper(dist_jac, data$cs.flag.emergent_overlap)
permutest(beta_dist_site, control=permControl (nperm=1000))
anova(beta_dist_site)
valor = beta_dist_site$eig
(beta_dist_site$eig/sum(beta_dist_site$eig))*100

PCoA_score = scores(beta_dist_site$vectors)
PCoA_score = PCoA_score[,1:3]
PCoA_score[,3] <- row.names(PCoA_score)
map_dist = merge(PCoA_score, data,by.x = "PCoA3",by.y = "MolForm")
colnames(map_dist)[1] <-"MolForm"
write.csv(map_dist, file="/../pcoa_data.csv")
map_dist2 = read.csv("/../pcoa_data.csv", header = T, row.names = 1)
attach(map_dist)

#######################
###Main figure plot###
#######################
ggplot(map_dist2, aes(PCoA1, PCoA2, color = cs.flag.emergent_overlap)) + 
  geom_jitter(size=1.5) +  
  #scale_shape_manual(values=c(1,2,14,15,3,4,16,17,10)) + 
  scale_color_manual(name="Sediment - Water overlap",labels=c("Core - Core","Core - in-between","Core - NA", "Core - Satellite","in-between - Core", "NA - Core", "Satellite - Core"), 
                     values=c("red2","green4", "navy","dodgerblue2","magenta2","chocolate3","purple"))+ 
  #stat_ellipse(aes(group=cs.flag.emergent_overlap,color=cs.flag.emergent_overlap),size=0.75)+
  #annotate("text", x = -2, y = -2, label = "Stress = 0.139")+
  #annotate("text", x = -2, y = -2.2, label = "PERMANOVA = 0.001")+
  theme_classic()+
  xlab("PCoA 1 (78.96%)")+
  ylab("PCoA 2 (6.73%)")

#####################
###PieDonut charts###
#####################
###Considering one figure for each overlap group

piedonut=data %>% group_by(cs.flag.emergent_overlap,Class) %>% summarize(n=n())
#1.Organize data
piedonut_filtered <- filter(piedonut, cs.flag.emergent_overlap == "Satellite-Core")
##count the number of times the same group was in the data
overlap_data <-
  aggregate(piedonut_filtered$n,by = list(overlap = piedonut_filtered$cs.flag.emergent_overlap),FUN = sum)
##order version data by browser so it will line up with browser pie chart
class_data <- piedonut_filtered[order(piedonut_filtered$cs.flag.emergent_overlap), ]

##format labels to display version and % market share
class_labels <- paste(class_data$Class, " (", class_data$n,")",  sep = "")


#2.Define colors
##adjust these as desired (currently colors all versions the same as browser)
class_colors <- c("lightcoral","springgreen", "seashell1","lightgoldenrod3","aquamarine2","grey60","deeppink4","olivedrab4","midnightblue")

#3.Draw the chart
##coordinates for the center of the chart
center_x <- 0.5
center_y <- 0.5
plot.new()
class_chart <-floating.pie(xpos = center_x,ypos = center_y,
                           x = class_data$n,radius = 0.35,border = "white",col = class_colors)
pie.labels(x = center_x,y = center_y,angles = class_chart,labels = class_labels,
           radius = 0.36,bg = NULL,cex = 1,font = 1.5,col = "gray40")
overlap_chart <-draw.circle(0.5,0.5,
                             radius = 0.18,border = "white",col = "white")

##The final figure was composed using an external editor
