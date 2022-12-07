rm(list=ls());graphics.off()
# Loading libraries
library (vegan); library(tidyverse)
library("dplyr") ;                     library('FSA')              
library("plyr")
library("readr") 
library("reshape2")
library(ggpubr)
options(digits=10)# to ensure we bring in all the decimals
# Setting working directories
home.dir = "C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/Topic1/Fig4_cross.envCS_VK/Transformations/"
setwd(home.dir)

# Read in transformations with emergent thresholds 
data = read.csv("Total_number_of_transformations_per_threshold.csv")

# dat = data %>% dplyr::select(cs.flag.emergent_overlap,Total_all_peaks,Total_all_sed,Total_all_water,Total_mf_peaks,Total_mf_sed,Total_mf_water)

# dat = data %>% dplyr::select(MolForm,cs.flag.emergent_overlap,Total_mf_peaks)

dat = data
dat = dat[!is.na(dat$cs.flag.emergent_overlap),]
df.mf = melt(dat)

df.mf$cs.flag.emergent_overlap= factor(df.mf$cs.flag.emergent_overlap,
                                       levels = c("Global core","Global in-between",
                                                  "Global satellite",
                                                  "Sed Core - Water Sat", "Sed Core - Water Inbetween",
                                                  "Sed Sat - Water Core", "Sed Inbetween - Water Core",
                                                  "Sed Sat - Water Inbetween", "Sed Inbetween - Water Sat"), 
                                       labels = c("Global core","Global in-between",
                                                  "Global satellite",
                                                  "Sed Core - Water Sat/In-bet.", "Sed Core - Water Sat/In-bet.",
                                                  "Water Core - Sed Sat/In-bet.", "Water Core - Sed Sat/In-bet.",
                                                  "Sed Sat - Water Inbetween", "Water Sat - Sed Inbetween"))

class = levels(unique(df.mf$cs.flag.emergent_overlap))
df.mf2 = subset(df.mf,df.mf$variable == "Mean")

ggplot(df.mf2, aes(x = cs.flag.emergent_overlap,fill = cs.flag.emergent_overlap , y = (value))) +
coord_cartesian() + geom_boxplot() +
stat_summary(fun.y = mean, geom = "point", shape = 18, size = 2, colour = "yellow",  position = position_dodge(width=0.75), alpha = 0.8)+ labs(x=expression()) +
  labs(y = expression(Mean~Transformations))+
  theme_bw()+theme(legend.position="")+theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()  ) +
  theme(legend.title = element_blank(),  legend.background = element_rect(fill = 'NA'), legend.text = element_text(size=8))+
  #theme(legend.position=c(0.8,0.9))+
  scale_fill_manual(values = c("#999999", "#661100", "#0072B2", "#FFB000",  "#FE6100", "#648FFF", "#DC267F" ))+
  theme(axis.text.x=element_text(colour = c("black","black")))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x=element_text(size=8))+
  theme(axis.text.y=element_text(size=8))+
  theme(axis.title.x =element_text(size=8))+
  theme(axis.title =element_text(size=8))+
  theme(axis.title.y =element_text(size=8))+annotate("text", x = 2, y = 120, label = "abc")+
   annotate("text", x = 7, y = 120, label = "bce")+annotate("text", x = 4, y = 120, label = "abd")+annotate("text", x = 6, y = 120, label = "b")+annotate("text", x = 5, y = 120, label = "d")+annotate("text", x = 3, y = 120, label = "e")

ggsave(file="Mean_Transformations_stats2.pdf")
ggsave(file="Mean_Transformations_stats2.png")

## Statistics
#Kruskal-Wallis: Non-parametric test
kruskal.test(value ~ cs.flag.emergent_overlap, data = df.mf2) #p-value<0.001 p = 2.2204e-16


#If the results of a Kruskal-Wallis test are statistically significant, then it's appropriate to conduct Dunn's Test to determine exactly which groups are different.
test = dunnTest(value ~ cs.flag.emergent_overlap, data = df.mf2,   method="holm") 
test2= test[["res"]]


write.csv(test2,"DunnTest_Mean_Transf.csv",row.names = F)


