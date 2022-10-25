#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#SCRIPT EXPLORATORY ANALYSIS Topic 1: DOM CORE-SAT PAPER (WHONDRS crowdsource)
#Supplementary figures: boxplots and density plots
#AUTHOR: Michaela de Melo
#CONTACT:michaelaldemelo@gmail.com

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##############################
# Analysis of Core-Sat flags #
##############################
#Load libraries
library(car)
library(ggplot2)
library(vegan)
library(FSA)
library(viridis)
library(ggpubr)
library(rstatix)
library(stats)
library(vegan)
library(ggpattern)
library("gridExtra")
#Path
setwd('C:/Users/micha/OneDrive/Documentos/GitHub/Topic1/4_gather.thresholds')
setwd("/Users/newuser/Downloads/Topic1-main/4_gather.thresholds")

#Input table
cross= read.csv('FTICR_crosstable_rep.merged1_all_em.thres_2022-05-05.csv')

#Checking the tables
colnames(cross)
unique(cross$cs.flag.emergent_sed) #Checking categories
unique(cross$cs.flag.emergent_water)

cross$cs.flag.emergent_sed= factor(cross$cs.flag.emergent_sed,
                                       levels = c("Core","In-between","Satellite" ), 
                                       labels = c("Core-Sed","In-Sed", "Sat-Sed"))

cross$cs.flag.emergent_water= factor(cross$cs.flag.emergent_water,
                                   levels = c("Core","In-between","Satellite" ), 
                                   labels = c("Core-Water","In-Water", "Sat-Water"))

#I wanna plot water and sediment boxplots in the same figure, so I renamed the flags and now I wil create a new table

sed=cross[,c (1:41,46 )]
colnames(sed)[42]="cs_flag"
sed$cat="sediment"
water=cross[,c (1:41,47 )]
colnames(water)[42]= "cs_flag"
water$cat="water"

bind=rbind(sed, water)

na_bind=na.omit(bind)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
           #Boxplots and violin Plots
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
 
########################
#1) NOSC
########################

res.kruskal <- na_bind %>% kruskal.test(cs_flag)
res.kruskal
stat.test <-na_bind %>% dunn_test( NOSC~cs_flag, p.adjust.method = "holm") 
stat.test
stat.test <- stat.test %>% add_xy_position(x = "cs_flag")

gnosc=ggboxplot(na_bind, x = "cs_flag" , y = "NOSC", fill = "cs_flag") +
  #stat_pvalue_manual(stat.test, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.text.x = element_blank())+
  scale_fill_manual(values=c("#CB6778", "#9785E0", "#8CCEED","#CB6778", "#9785E0", "#8CCEED" ))+
  labs(x=" ", y="NOSC")+ annotate("text", x=4, y=2.7, label="a")+annotate("text", x=5, y=2.7, label="a")


########################
#2) DBE
########################
stat.test <-na_bind %>% dunn_test( DBE~cs_flag, p.adjust.method = "holm") 
stat.test
gdbe=ggboxplot(na_bind, x = "cs_flag" , y = "DBE", fill = "cs_flag") +
  #stat_pvalue_manual(stat.test, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.text.x = element_blank())+
  scale_fill_manual(values=c("#CB6778", "#9785E0", "#8CCEED","#CB6778", "#9785E0", "#8CCEED" ))+
  labs(x=" ", y="DBE")+ annotate("text", x=4, y=29, label="a")+
  annotate("text", x=5, y=29, label="a")+annotate("text", x=6, y=29, label="a")


########################
#3) GFE
########################
stat.test <-na_bind %>% dunn_test( GFE~cs_flag, p.adjust.method = "holm") 
stat.test
ggfe=ggboxplot(na_bind, x = "cs_flag" , y = "GFE", fill = "cs_flag") +
  #stat_pvalue_manual(stat.test, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.text.x = element_text(angle=90), axis.ticks.x = element_blank())+ 
  scale_fill_manual(values=c("#CB6778", "#9785E0", "#8CCEED","#CB6778", "#9785E0", "#8CCEED" ))+
  labs(x=" ", y="GFE")+ annotate("text", x=4, y=125, label="a")+
  annotate("text", x=5, y=125, label="a")

#############
#4) AI_mod
############
stat.test <-na_bind %>% dunn_test( AI_Mod~cs_flag, p.adjust.method = "holm") 
stat.test
gaimod=ggboxplot(na_bind, x = "cs_flag" , y = "AI_Mod", fill = "cs_flag") +
 #stat_pvalue_manual(stat.test, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.text.x = element_text(angle=90), axis.ticks.x = element_blank())+ 
  scale_fill_manual(values=c("#CB6778", "#9785E0", "#8CCEED","#CB6778", "#9785E0", "#8CCEED" ))+
 annotate("text", x=3, y=2.5, label="a")+ annotate("text", x=6, y=2.5, label="a")+
annotate("text", x=4, y=2.5, label="b")+ annotate("text", x=5, y=2.5, label="b")+
  labs(x=" ", y="AI modified")

########################
#5) Mass
########################
stat.test <-na_bind %>% dunn_test( Mass~cs_flag, p.adjust.method = "holm") 
stat.test
gmass=ggboxplot(na_bind, x = "cs_flag" , y = "Mass", fill = "cs_flag") +
  #stat_pvalue_manual(stat.test, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank())+ 
  scale_fill_manual(values=c("#CB6778", "#9785E0", "#8CCEED","#CB6778", "#9785E0", "#8CCEED" ))+
  labs(x=" ", y="Mass")


##-> Arrange and save it!

g=arrangeGrob(gmass,  gnosc, gdbe, gaimod, ggfe, ncol=3, heights=c(3,4))

ggsave("Panel_DOMproperties.png", dpi=300, height = 6.5,  width = 8, g)
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#    Density plots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dens_mass=ggplot(data=na_bind, aes(x =Mass, fill = cs_flag)) +geom_density(alpha=0.5)+ facet_grid(~cat)+
  theme_bw() +theme(legend.position = c(0.35, 0.78), legend.title=element_blank(), legend.text=element_text(size=5), legend.key.size = unit(0.3, 'cm'))+
  scale_fill_manual(values=c("#CB6778", "#9785E0", "#8CCEED","#CB6778", "#9785E0", "#8CCEED" ))+
  labs(x="Mass", y="Density")

dens_gfe=ggplot(data=na_bind, aes(x =GFE, fill = cs_flag)) +geom_density(alpha=0.5)+ facet_grid(~cat)+
  theme_bw() +theme(legend.position = "none")+
  scale_fill_manual(values=c("#CB6778", "#9785E0", "#8CCEED","#CB6778", "#9785E0", "#8CCEED" ))+
  labs(x="GFE", y="Density", fill=" ")

dens_nosc=ggplot(data=na_bind, aes(x =NOSC, fill = cs_flag)) +geom_density(alpha=0.5)+ facet_grid(~cat)+
  theme_bw() +theme(legend.position = "none")+
  scale_fill_manual(values=c("#CB6778", "#9785E0", "#8CCEED","#CB6778", "#9785E0", "#8CCEED" ))+
  labs(x="NOSC", y="Density", fill=" ")

dens_dbe=ggplot(data=na_bind, aes(x =DBE, fill = cs_flag)) +geom_density(alpha=0.5)+ facet_grid(~cat)+
  theme_bw() +theme(legend.position = "none")+
  scale_fill_manual(values=c("#CB6778", "#9785E0", "#8CCEED","#CB6778", "#9785E0", "#8CCEED" ))+
  labs(x="DBE", y="Density", fill=" ")

## Arrange and save it!
dens=arrangeGrob(dens_mass, dens_gfe, dens_nosc, dens_dbe, nrow=2)

ggsave("Panel_Density_Plots.png", dpi=300, height = 6,  width = 8, dens)
