rm(list=ls());graphics.off()
# Loading libraries
library (vegan); library(tidyverse)
library("dplyr")                                    
library("plyr")                                     
library("readr") 
library("reshape2")

options(digits=10)# to ensure we bring in all the decimals
# Setting working directories
home.dir = "C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/Topic1/Fig4_cross.envCS_VK/Transformations/"
setwd(home.dir)

# Read in transfomations with emergent tresholds 
data = read.csv("Total_number_of_transformations_per_threshold.csv")

dat = data %>% dplyr::select(cs.flag.emergent_overlap,Total_all_peaks,Total_all_sed,Total_all_water,Total_mf_peaks,Total_mf_sed,Total_mf_water)

dat = dat[!is.na(dat$cs.flag.emergent_overlap),]
df = melt(dat)

df.all = df[grep("all",df$variable),]
df.mf = df[grep("mf",df$variable),]

ggplot(df.all, aes(x = cs.flag.emergent_overlap , fill = variable, y = value)) +
coord_cartesian() + geom_boxplot() +
stat_summary(fun.y = mean, geom = "point", shape = 18, size = 2, colour = "yellow",  position = position_dodge(width=0.75), alpha = 0.8)+
labs(x=expression()) +
labs(y = expression(Number~of~Transformations))+ 
  # stat_compare_means(method="t-test",hide.ns = TRUE,method.args = list(alternative = "less"),label="p.signif", label.y=0)+
  theme_bw()+theme(legend.position="top")+theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()  ) +
  theme(legend.title = element_blank(),  legend.background = element_rect(fill = 'NA'), legend.text = element_text(size=8))+
  #theme(legend.position=c(0.8,0.9))+
  scale_fill_brewer(palette="Set2")+ 
  theme(axis.text.x=element_text(colour = c("black","black")))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x=element_text(size=8))+
  theme(axis.text.y=element_text(size=8))+
  theme(axis.title.x =element_text(size=8))+
  theme(axis.title =element_text(size=8))+
  theme(axis.title.y =element_text(size=8))
ggsave(file="Transformation_all_thresholds.pdf")

ggplot(df.mf, aes(x = cs.flag.emergent_overlap , fill = variable, y = value)) +
  coord_cartesian() + geom_boxplot() +
  stat_summary(fun.y = mean, geom = "point", shape = 18, size = 2, colour = "yellow",  position = position_dodge(width=0.75), alpha = 0.8)+
  labs(x=expression()) +
  labs(y = expression(Number~of~Transformations))+ 
  # stat_compare_means(method="t-test",hide.ns = TRUE,method.args = list(alternative = "less"),label="p.signif", label.y=0)+
  theme_bw()+theme(legend.position="top")+theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()  ) +
  theme(legend.title = element_blank(),  legend.background = element_rect(fill = 'NA'), legend.text = element_text(size=8))+
  #theme(legend.position=c(0.8,0.9))+
  scale_fill_brewer(palette="Set2")+ 
  theme(axis.text.x=element_text(colour = c("black","black")))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x=element_text(size=8))+
  theme(axis.text.y=element_text(size=8))+
  theme(axis.title.x =element_text(size=8))+
  theme(axis.title =element_text(size=8))+
  theme(axis.title.y =element_text(size=8))
ggsave(file="Transformation_mf_thresholds.pdf")


ggplot(df.all, aes(fill = cs.flag.emergent_overlap , x = variable, y = value)) +
  coord_cartesian() + geom_boxplot() +
  stat_summary(fun.y = mean, geom = "point", shape = 18, size = 2, colour = "yellow",  position = position_dodge(width=0.75), alpha = 0.8)+
  labs(x=expression()) +
  labs(y = expression(Number~of~Transformations))+ 
  # stat_compare_means(method="t-test",hide.ns = TRUE,method.args = list(alternative = "less"),label="p.signif", label.y=0)+
  theme_bw()+theme(legend.position="top")+theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()  ) +
  theme(legend.title = element_blank(),  legend.background = element_rect(fill = 'NA'), legend.text = element_text(size=8))+
  #theme(legend.position=c(0.8,0.9))+
  scale_fill_brewer(palette="Set2")+ 
  theme(axis.text.x=element_text(colour = c("black","black")))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x=element_text(size=8))+
  theme(axis.text.y=element_text(size=8))+
  theme(axis.title.x =element_text(size=8))+
  theme(axis.title =element_text(size=8))+
  theme(axis.title.y =element_text(size=8))
ggsave(file="Transformation_all_thresholds_v2.pdf")

ggplot(df.mf, aes(fill = cs.flag.emergent_overlap , x = variable, y = value)) +
  coord_cartesian() + geom_boxplot() +
  stat_summary(fun.y = mean, geom = "point", shape = 18, size = 2, colour = "yellow",  position = position_dodge(width=0.75), alpha = 0.8)+
  labs(x=expression()) +
  labs(y = expression(Number~of~Transformations))+ 
  # stat_compare_means(method="t-test",hide.ns = TRUE,method.args = list(alternative = "less"),label="p.signif", label.y=0)+
  theme_bw()+theme(legend.position="top")+theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()  ) +
  theme(legend.title = element_blank(),  legend.background = element_rect(fill = 'NA'), legend.text = element_text(size=8))+
  #theme(legend.position=c(0.8,0.9))+
  scale_fill_brewer(palette="Set2")+ 
  theme(axis.text.x=element_text(colour = c("black","black")))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x=element_text(size=8))+
  theme(axis.text.y=element_text(size=8))+
  theme(axis.title.x =element_text(size=8))+
  theme(axis.title =element_text(size=8))+
  theme(axis.title.y =element_text(size=8))
ggsave(file="Transformation_mf_thresholds_v2.pdf")