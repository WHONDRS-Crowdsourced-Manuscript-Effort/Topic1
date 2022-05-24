#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#SCRIPT EXPLORATORY ANALYSIS Topic 1: DOM CORE-SAT PAPER (WHONDRS crowdsource)
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
library(gridExtra)
library(rstatix)
library(dplyr)


#-->  INPUT DATA #

#Path
setwd('C:/Users/micha/OneDrive/Área de Trabalho/Topic1/4_gather.thresholds')


#Input tables
cross_merge1= read.csv('FTICR_crosstable_rep.merged1_all_em.thres_2022-05-05.csv')
names(cross_merge1)

cross_merge2= read.csv('FTICR_crosstable_rep.merged2_all_em.thres_2022-05-05.csv')
names(cross_merge2)




#Remove rows with Flags=NA
cross_merge1=subset(cross_merge1, !is.na(cs.flag.emergent_sed))
cross_merge1=subset(cross_merge1, !is.na(cs.flag.emergent_water))
cross_merge2=subset(cross_merge1, !is.na(cs.flag.emergent_sed))
cross_merge2=subset(cross_merge1, !is.na(cs.flag.emergent_water))


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
           #Boxplots and violin Plots
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
 
########################
#1) Merged 1 - SEDIMENTS
########################
#Plot
  geom_boxplot()

# To compare if metrics are statistically different between groups, we need to run a test and check the assumptions
#Normality and Variances to decide for a parametric or non-parametric test

#Testing the homegeneity of Variances with Levene's Test in 2+x samples
## Levene's test
#Null hypothesis: group variances are equal
leveneTest(NOSC ~ cs.flag.emergent_sed, data = cross_merge1)
leveneTest(DBE ~ cs.flag.emergent_sed, data = cross_merge1)
leveneTest(AI_Mod ~ cs.flag.emergent_sed, data = cross_merge1)
leveneTest(OtoC_ratio ~ cs.flag.emergent_sed, data = cross_merge1)
leveneTest(HtoC_ratio ~ cs.flag.emergent_sed, data = cross_merge1)
leveneTest(Mass ~ cs.flag.emergent_sed, data = cross_merge1)



#Non homogeneous --> need to use a non-parametric test

#A Kruskal-Wallis test is used to determine whether or not there is a statistically significant difference between the medians of three or more independent groups. 
#It is considered to be the non-parametric equivalent of the One-Way ANOVA.
kruskal.test(NOSC ~ cs.flag.emergent_sed, data = cross_merge1)

#If the results of a Kruskal-Wallis test are statistically significant, then it's appropriate to conduct Dunn's Test to determine exactly which groups are different.
dunnTest(NOSC ~ cs.flag.emergent_sed, data = cross_merge1,
         method="holm")
#All the 3 groups are statistically significantly different from each other 

##Boxplot
p=ggplot(data=cross_merge1, aes(x=cs.flag.emergent_sed, y=NOSC, fill=cs.flag.emergent_sed))+   
  geom_boxplot() + theme_bw() +theme(legend.position="none")+ scale_fill_viridis(discrete = TRUE)
p+stat_compare_means(method = "kruskal.test")


########################################################################################
#I can plot everything together in a simpler command

##NOSC - SED MERGE1
res.kruskal <- cross_merge1 %>% kruskal_test( NOSC~cs.flag.emergent_sed)
res.kruskal
stat.test <-cross_merge1 %>% dunn_test( NOSC~cs.flag.emergent_sed, p.adjust.method = "holm") 
stat.test
stat.test <- stat.test %>% add_xy_position(x = "cs.flag.emergent_sed")

cross_merge1$cs.flag.emergent_sed= factor(cross_merge1$cs.flag.emergent_sed, levels= c("Core", "In-between", "Satellite", "NA"))
plot1=ggboxplot(cross_merge1, x = "cs.flag.emergent_sed" , y = "NOSC", fill = "cs.flag.emergent_sed") +
  stat_pvalue_manual(stat.test, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.title.x=element_blank())+ scale_fill_viridis(discrete = TRUE)
#setwd("C:/Users/micha/OneDrive/Documentos/GitHub/Topic1/7_molecular.traits/Boxplots/")
#ggsave("Boxplot_NOSC_sediment_merged1.png", dpi=300, width = 6, height = 4)

########################################################################################
## DBE - SED MERGE1
#I can plot everything together
res.kruskal2 <- cross_merge1 %>% kruskal_test( DBE~cs.flag.emergent_sed)
res.kruskal2
stat.test2 <-cross_merge1 %>% dunn_test( DBE~cs.flag.emergent_sed, p.adjust.method = "holm") 
stat.test2
stat.test2 <- stat.test2 %>% add_xy_position(x = "cs.flag.emergent_sed")


plot2=ggboxplot(cross_merge1, x = "cs.flag.emergent_sed" , y = "DBE", fill = "cs.flag.emergent_sed") +
  stat_pvalue_manual(stat.test2, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.title.x=element_blank())+ scale_fill_viridis(discrete = TRUE)

#ggsave("Boxplot_DBE_sediment_merged1.png", dpi=300, width = 6, height = 4)

######################################################################################
## Ai Mod - SED MERGE1
#I can plot everything together
res.kruskal3 <- cross_merge1 %>% kruskal_test( AI_Mod~cs.flag.emergent_sed)
res.kruskal3
stat.test3 <-cross_merge1 %>% dunn_test( AI_Mod~cs.flag.emergent_sed, p.adjust.method = "holm") 
stat.test3
stat.test3 <- stat.test3 %>% add_xy_position(x = "cs.flag.emergent_sed")

plot3=ggboxplot(cross_merge1, x = "cs.flag.emergent_sed" , y = "AI_Mod", fill = "cs.flag.emergent_sed") +
  stat_pvalue_manual(stat.test3, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.title.x=element_blank())+ scale_fill_viridis(discrete = TRUE)

#ggsave("Boxplot_AI-Mod_sediment_merged1.png", dpi=300, width = 6, height = 4)


######################################################################################
## O/C_ratio - SED MERGE1
#I can plot everything together
res.kruskal3a <- cross_merge1 %>% kruskal_test( OtoC_ratio~cs.flag.emergent_sed)
res.kruskal3a
stat.test3a <-cross_merge1 %>% dunn_test( OtoC_ratio~cs.flag.emergent_sed, p.adjust.method = "holm") 
stat.test3a
stat.test3a <- stat.test3a %>% add_xy_position(x = "cs.flag.emergent_sed")

plot3a=ggboxplot(cross_merge1, x = "cs.flag.emergent_sed" , y = "OtoC_ratio", fill = "cs.flag.emergent_sed") +
  stat_pvalue_manual(stat.test3a, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.title.x=element_blank())+ scale_fill_viridis(discrete = TRUE)

#ggsave("Boxplot_OtoC_ratio_sediment_merged1.png", dpi=300, width = 6, height = 4)

######################################################################################
## H/C_ratio - SED MERGE1
#I can plot everything together
res.kruskal3b <- cross_merge1 %>% kruskal_test( HtoC_ratio~cs.flag.emergent_sed)
res.kruskal3b
stat.test3b <-cross_merge1 %>% dunn_test( HtoC_ratio~cs.flag.emergent_sed, p.adjust.method = "holm") 
stat.test3b
stat.test3b <- stat.test3b %>% add_xy_position(x = "cs.flag.emergent_sed")

plot3b=ggboxplot(cross_merge1, x = "cs.flag.emergent_sed" , y = "HtoC_ratio", fill = "cs.flag.emergent_sed") +
  stat_pvalue_manual(stat.test3b, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.title.x=element_blank())+ scale_fill_viridis(discrete = TRUE)

#ggsave("Boxplot_HtoC_ratio_sediment_merged1.png", dpi=300, width = 6, height = 4)

######################################################################################
## Mass
#I can plot everything together
res.kruskal3c <- cross_merge1 %>% kruskal_test( Mass~cs.flag.emergent_sed)
res.kruskal3c
stat.test3c <-cross_merge1 %>% dunn_test( Mass~cs.flag.emergent_sed, p.adjust.method = "holm") 
stat.test3c
stat.test3c <- stat.test3c %>% add_xy_position(x = "cs.flag.emergent_sed")

plot3c=ggboxplot(cross_merge1, x = "cs.flag.emergent_sed" , y = "Mass", fill = "cs.flag.emergent_sed") +
  stat_pvalue_manual(stat.test3c, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.title.x=element_blank())+ scale_fill_viridis(discrete = TRUE)

#ggsave("Boxplot_Mass_sediment_merged1.png", dpi=300, width = 6, height = 4)



#Arrange plots together

setwd("C:/Users/micha/OneDrive/Área de Trabalho/Topic1/7_molecular.traits/Boxplots/")
#Arrange plots together
grob_sed1=gridExtra::grid.arrange(plot1, plot2, plot3, plot3c, nrow=2, top = "Sediment")
ggsave("Boxplots_arranged_sediment_merged1.png", dpi=300, width = 6, height = 4, grob_sed1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
########################
#2) Merged 1 - WATER
########################
## Levene's test
leveneTest(NOSC ~ cs.flag.emergent_water, data = cross_merge1)
leveneTest(DBE ~ cs.flag.emergent_water, data = cross_merge1)
leveneTest(AI_Mod ~ cs.flag.emergent_water, data = cross_merge1)
leveneTest(OtoC_ratio ~ cs.flag.emergent_water, data = cross_merge1)
leveneTest(HtoC_ratio ~ cs.flag.emergent_water, data = cross_merge1)
leveneTest(Mass ~ cs.flag.emergent_water, data = cross_merge1)



##NOSC - WATER MERGE1
res.kruskal4 <- cross_merge1 %>% kruskal_test( NOSC~cs.flag.emergent_water)
res.kruskal4
stat.test4 <-cross_merge1 %>% dunn_test( NOSC~cs.flag.emergent_water, p.adjust.method = "holm") 
stat.test4
stat.test4 <- stat.test4 %>% add_xy_position(x = "cs.flag.emergent_water")

cross_merge1$cs.flag.emergent_water= factor(cross_merge1$cs.flag.emergent_water, levels= c("Core", "In-between", "Satellite", "NA"))
plot4=ggboxplot(cross_merge1, x = "cs.flag.emergent_water" , y = "NOSC", fill = "cs.flag.emergent_water") +
  stat_pvalue_manual(stat.test4, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.title.x = element_blank())+ scale_fill_viridis(discrete = TRUE)
#setwd("C:/Users/micha/OneDrive/Documentos/GitHub/Topic1/7_molecular.traits/Boxplots/")
#ggsave("Boxplot_NOSC_water_merged1.png", dpi=300, width = 6, height = 4)

########################################################################################
## DBE - WATER MERGE1
#I can plot everything together
res.kruskal5 <- cross_merge1 %>% kruskal_test( DBE~cs.flag.emergent_water)
res.kruskal5
stat.test5 <-cross_merge1 %>% dunn_test( DBE~cs.flag.emergent_water, p.adjust.method = "holm") 
stat.test5
stat.test5 <- stat.test5 %>% add_xy_position(x = "cs.flag.emergent_water")

plot5=ggboxplot(cross_merge1, x = "cs.flag.emergent_water" , y = "DBE", fill = "cs.flag.emergent_water") +
  stat_pvalue_manual(stat.test5, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.title.x = element_blank())+ scale_fill_viridis(discrete = TRUE)

#ggsave("Boxplot_DBE_water_merged1.png", dpi=300, width = 6, height = 4)

######################################################################################
## Ai Mod - WATER MERGE1
#I can plot everything together
res.kruskal6 <- cross_merge1 %>% kruskal_test( AI_Mod~cs.flag.emergent_water)
res.kruskal6
stat.test6 <-cross_merge1 %>% dunn_test( AI_Mod~cs.flag.emergent_water, p.adjust.method = "holm") 
stat.test6
stat.test6 <- stat.test6 %>% add_xy_position(x = "cs.flag.emergent_water")

plot6=ggboxplot(cross_merge1, x = "cs.flag.emergent_water" , y = "AI_Mod", fill = "cs.flag.emergent_water") +
  stat_pvalue_manual(stat.test6, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.title.x = element_blank())+ scale_fill_viridis(discrete = TRUE)
#setwd("C:/Users/micha/OneDrive/Documentos/GitHub/Topic1/7_molecular.traits/Boxplots/")
#ggsave("Boxplot_AI-Mod_water_merged1.png", dpi=300, width = 6, height = 4)

######################################################################################
## O/C ratio - WATER MERGE1
#I can plot everything together
res.kruskal6a <- cross_merge1 %>% kruskal_test( OtoC_ratio~cs.flag.emergent_water)
res.kruskal6a
stat.test6a <-cross_merge1 %>% dunn_test( OtoC_ratio~cs.flag.emergent_water, p.adjust.method = "holm") 
stat.test6a
stat.test6a <- stat.test6a %>% add_xy_position(x = "cs.flag.emergent_water")

ggboxplot(cross_merge1, x = "cs.flag.emergent_water" , y = "OtoC_ratio", fill = "cs.flag.emergent_water") +
  stat_pvalue_manual(stat.test6a, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none")+ scale_fill_viridis(discrete = TRUE)
#setwd("C:/Users/micha/OneDrive/Documentos/GitHub/Topic1/7_molecular.traits/Boxplots/")
#ggsave("Boxplot_OtoC_ratio_water_merged1.png", dpi=300, width = 6, height = 4)


######################################################################################
## H/C ratio - WATER MERGE1
#I can plot everything together
res.kruskal6b <- cross_merge1 %>% kruskal_test( HtoC_ratio~cs.flag.emergent_water)
res.kruskal6b
stat.test6b <-cross_merge1 %>% dunn_test( HtoC_ratio~cs.flag.emergent_water, p.adjust.method = "holm") 
stat.test6b
stat.test6b <- stat.test6b %>% add_xy_position(x = "cs.flag.emergent_water")

ggboxplot(cross_merge1, x = "cs.flag.emergent_water" , y = "HtoC_ratio", fill = "cs.flag.emergent_water") +
  stat_pvalue_manual(stat.test6b, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none")+ scale_fill_viridis(discrete = TRUE)
#setwd("C:/Users/micha/OneDrive/Documentos/GitHub/Topic1/7_molecular.traits/Boxplots/")
#ggsave("Boxplot_HtoC_ratio_water_merged1.png", dpi=300, width = 6, height = 4)



######################################################################################
## Mass - WATER MERGE1
#I can plot everything together
res.kruskal6c <- cross_merge1 %>% kruskal_test( Mass~cs.flag.emergent_water)
res.kruskal6c
stat.test6c <-cross_merge1 %>% dunn_test( Mass~cs.flag.emergent_water, p.adjust.method = "holm") 
stat.test6c
stat.test6c <- stat.test6c %>% add_xy_position(x = "cs.flag.emergent_water")

plot6c=ggboxplot(cross_merge1, x = "cs.flag.emergent_water" , y = "Mass", fill = "cs.flag.emergent_water") +
  stat_pvalue_manual(stat.test6c, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.title.x = element_blank())+ scale_fill_viridis(discrete = TRUE)
#setwd("C:/Users/micha/OneDrive/Documentos/GitHub/Topic1/7_molecular.traits/Boxplots/")
#ggsave("Boxplot_Mass_water_merged1.png", dpi=300, width = 6, height = 4)

setwd("C:/Users/micha/OneDrive/Área de Trabalho/Topic1/7_molecular.traits/Boxplots/")
#Arrange plots together
grob_water1=gridExtra::grid.arrange(plot4, plot5, plot6, plot6c, nrow=2, top = "Water")
ggsave("Boxplots_arranged_water_merged1.png", dpi=300, width = 6, height = 4, grob_water1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
########################
#2) Merged 2 - SEDIMENT
########################
leveneTest(NOSC ~ cs.flag.emergent_water, data = cross_merge2)
leveneTest(DBE ~ cs.flag.emergent_water, data = cross_merge2)
leveneTest(AI_Mod ~ cs.flag.emergent_water, data = cross_merge2)
leveneTest(OtoC_ratio ~ cs.flag.emergent_water, data = cross_merge2)
leveneTest(HtoC_ratio ~ cs.flag.emergent_water, data = cross_merge2)
leveneTest(Mass ~ cs.flag.emergent_water, data = cross_merge2)


##NOSC - SED MERGE2
res.kruskal8 <- cross_merge2 %>% kruskal_test( NOSC~cs.flag.emergent_sed)
res.kruskal8
stat.test8 <-cross_merge2 %>% dunn_test( NOSC~cs.flag.emergent_sed, p.adjust.method = "holm") 
stat.test8
stat.test8<- stat.test8 %>% add_xy_position(x = "cs.flag.emergent_sed")

cross_merge2$cs.flag.emergent_sed= factor(cross_merge2$cs.flag.emergent_sed, levels= c("Core", "In-between", "Satellite", "NA")) #Order flags in the boxplot
ggboxplot(cross_merge2, x = "cs.flag.emergent_sed" , y = "NOSC", fill = "cs.flag.emergent_sed") +
  stat_pvalue_manual(stat.test8, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none")+ scale_fill_viridis(discrete = TRUE)
setwd("C:/Users/micha/OneDrive/Documentos/GitHub/Topic1/7_molecular.traits/Boxplots/")
ggsave("Boxplot_NOSC_sediment_merged2.png", dpi=300, width = 6, height = 4)

########################################################################################
## DBE - SED MERGE2
#I can plot everything together
res.kruskal9 <- cross_merge2 %>% kruskal_test( DBE~cs.flag.emergent_sed)
res.kruskal9
stat.test9 <-cross_merge2 %>% dunn_test( DBE~cs.flag.emergent_sed, p.adjust.method = "holm") 
stat.test9
stat.test9 <- stat.test9 %>% add_xy_position(x = "cs.flag.emergent_sed")

ggboxplot(cross_merge2, x = "cs.flag.emergent_sed" , y = "DBE", fill = "cs.flag.emergent_sed") +
  stat_pvalue_manual(stat.test9, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none")+ scale_fill_viridis(discrete = TRUE)

ggsave("Boxplot_DBE_sediment_merged2.png", dpi=300, width = 6, height = 4)

######################################################################################
## Ai Mod - SED MERGE2
#I can plot everything together
res.kruskal10 <- cross_merge2 %>% kruskal_test( AI_Mod~cs.flag.emergent_sed)
res.kruskal10
stat.test10 <-cross_merge2 %>% dunn_test( AI_Mod~cs.flag.emergent_sed, p.adjust.method = "holm") 
stat.test10
stat.test10 <- stat.test10 %>% add_xy_position(x = "cs.flag.emergent_sed")

ggboxplot(cross_merge2, x = "cs.flag.emergent_sed" , y = "AI_Mod", fill = "cs.flag.emergent_sed") +
  stat_pvalue_manual(stat.test10, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none")+ scale_fill_viridis(discrete = TRUE)

ggsave("Boxplot_AI-Mod_sediment_merged2.png", dpi=300, width = 6, height = 4)


######################################################################################
## O/C ratio - SED MERGE2
#I can plot everything together
res.kruskal10a <- cross_merge2 %>% kruskal_test( Mass~cs.flag.emergent_sed)
res.kruskal10a
stat.test10a <-cross_merge2 %>% dunn_test( Mass~cs.flag.emergent_sed, p.adjust.method = "holm") 
stat.test10a
stat.test10a <- stat.test10a %>% add_xy_position(x = "cs.flag.emergent_sed")

ggboxplot(cross_merge2, x = "cs.flag.emergent_sed" , y = "Mass", fill = "cs.flag.emergent_sed") +
  stat_pvalue_manual(stat.test10a, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none")+ scale_fill_viridis(discrete = TRUE)

ggsave("Boxplot_Mass_sediment_merged2.png", dpi=300, width = 6, height = 4)


######################################################################################
## H/C ratio - SED MERGE2
#I can plot everything together
res.kruskal10b <- cross_merge2 %>% kruskal_test( HtoC_ratio~cs.flag.emergent_sed)
res.kruskal10b
stat.test10b <-cross_merge2 %>% dunn_test( HtoC_ratio~cs.flag.emergent_sed, p.adjust.method = "holm") 
stat.test10b
stat.test10b <- stat.test10b %>% add_xy_position(x = "cs.flag.emergent_sed")

ggboxplot(cross_merge2, x = "cs.flag.emergent_sed" , y = "HtoC_ratio", fill = "cs.flag.emergent_sed") +
  stat_pvalue_manual(stat.test10b, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none")+ scale_fill_viridis(discrete = TRUE)

ggsave("Boxplot_HtoC_ratio_sediment_merged2.png", dpi=300, width = 6, height = 4)


######################################################################################
## Mass - SED MERGE2
#I can plot everything together
res.kruskal10c <- cross_merge2 %>% kruskal_test( OtoC_ratio~cs.flag.emergent_sed)
res.kruskal10c
stat.test10c <-cross_merge2 %>% dunn_test( OtoC_ratio~cs.flag.emergent_sed, p.adjust.method = "holm") 
stat.test10c
stat.test10c <- stat.test10c %>% add_xy_position(x = "cs.flag.emergent_sed")

ggboxplot(cross_merge2, x = "cs.flag.emergent_sed" , y = "OtoC_ratio", fill = "cs.flag.emergent_sed") +
  stat_pvalue_manual(stat.test10c, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none")+ scale_fill_viridis(discrete = TRUE)

ggsave("Boxplot_OtoC_ratio_sediment_merged2.png", dpi=300, width = 6, height = 4)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

########################
#2) Merged 2 - WATER
########################
## Levene's test
leveneTest(NOSC~ cs.flag.emergent_water, data = cross_merge2)
leveneTest(DBE ~ cs.flag.emergent_water, data = cross_merge2)
leveneTest(AI_Mod ~ cs.flag.emergent_water, data = cross_merge2)
leveneTest(OtoC_ratio ~ cs.flag.emergent_water, data = cross_merge2)
leveneTest(HtoC_ratio ~ cs.flag.emergent_water, data = cross_merge2)
leveneTest(Mass ~ cs.flag.emergent_water, data = cross_merge2)


##NOSC - WATER MERGE1
res.kruskal11 <- cross_merge2 %>% kruskal_test( NOSC~cs.flag.emergent_water)
res.kruskal11
stat.test11 <-cross_merge2 %>% dunn_test( NOSC~cs.flag.emergent_water, p.adjust.method = "holm") 
stat.test11
stat.test11 <- stat.test11 %>% add_xy_position(x = "cs.flag.emergent_water")

cross_merge2$cs.flag.emergent_water= factor(cross_merge2$cs.flag.emergent_water, levels= c("Core", "In-between", "Satellite", "NA")) #Order flags in the boxplot
ggboxplot(cross_merge2, x = "cs.flag.emergent_water" , y = "NOSC", fill = "cs.flag.emergent_water") +
  stat_pvalue_manual(stat.test11, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none")+ scale_fill_viridis(discrete = TRUE)
setwd("C:/Users/micha/OneDrive/Documentos/GitHub/Topic1/7_molecular.traits/Boxplots/")
ggsave("Boxplot_NOSC_water_merged2.png", dpi=300, width = 6, height = 4)

########################################################################################
## DBE - WATER MERGE2
#I can plot everything together
res.kruskal12 <- cross_merge2 %>% kruskal_test( DBE~cs.flag.emergent_water)
res.kruskal12
stat.test12 <-cross_merge2 %>% dunn_test( DBE~cs.flag.emergent_water, p.adjust.method = "holm") 
stat.test12
stat.test12 <- stat.test12 %>% add_xy_position(x = "cs.flag.emergent_water")

ggboxplot(cross_merge2, x = "cs.flag.emergent_water" , y = "DBE", fill = "cs.flag.emergent_water") +
  stat_pvalue_manual(stat.test12, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none")+ scale_fill_viridis(discrete = TRUE)

setwd("C:/Users/micha/OneDrive/Documentos/GitHub/Topic1/7_molecular.traits/Boxplots/")
ggsave("Boxplot_DBE_water_merged2.png", dpi=300, width = 6, height = 4)

######################################################################################
## Ai Mod - WATER MERGE2
#I can plot everything together
res.kruskal13 <- cross_merge2 %>% kruskal_test( AI_Mod~cs.flag.emergent_water)
res.kruskal13
stat.test13 <-cross_merge2 %>% dunn_test( AI_Mod~cs.flag.emergent_water, p.adjust.method = "holm") 
stat.test13
stat.test13 <- stat.test13 %>% add_xy_position(x = "cs.flag.emergent_water")

ggboxplot(cross_merge2, x = "cs.flag.emergent_water" , y = "AI_Mod", fill = "cs.flag.emergent_water") +
  stat_pvalue_manual(stat.test13, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none")+ scale_fill_viridis(discrete = TRUE)

ggsave("Boxplot_AI-Mod_water_merged2.png", dpi=300, width = 6, height = 4)


######################################################################################
## O/C ratio - WATER MERGE1
#I can plot everything together
res.kruskal13a <- cross_merge2 %>% kruskal_test( OtoC_ratio~cs.flag.emergent_water)
res.kruskal13a
stat.test13a <-cross_merge2 %>% dunn_test( OtoC_ratio~cs.flag.emergent_water, p.adjust.method = "holm") 
stat.test13a
stat.test13a <- stat.test13a %>% add_xy_position(x = "cs.flag.emergent_water")

ggboxplot(cross_merge2, x = "cs.flag.emergent_water" , y = "OtoC_ratio", fill = "cs.flag.emergent_water") +
  stat_pvalue_manual(stat.test13a, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none")+ scale_fill_viridis(discrete = TRUE)

ggsave("Boxplot_OtoC_ratio_water_merged2.png", dpi=300, width = 6, height = 4)


######################################################################################
## H/C ratio - WATER MERGE1
#I can plot everything together
res.kruskal13b <- cross_merge2 %>% kruskal_test( HtoC_ratio~cs.flag.emergent_water)
res.kruskal13b
stat.test13b <-cross_merge2 %>% dunn_test( HtoC_ratio~cs.flag.emergent_water, p.adjust.method = "holm") 
stat.test13b
stat.test13b <- stat.test13b %>% add_xy_position(x = "cs.flag.emergent_water")

ggboxplot(cross_merge2, x = "cs.flag.emergent_water" , y = "HtoC_ratio", fill = "cs.flag.emergent_water") +
  stat_pvalue_manual(stat.test13b, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none")+ scale_fill_viridis(discrete = TRUE)

ggsave("Boxplot_HtoC_ratio_water_merged2.png", dpi=300, width = 6, height = 4)

######################################################################################
## Mass - WATER MERGE1
#I can plot everything together
res.kruskal13c <- cross_merge2 %>% kruskal_test( Mass~cs.flag.emergent_water)
res.kruskal13c
stat.test13c <-cross_merge2 %>% dunn_test( Mass~cs.flag.emergent_water, p.adjust.method = "holm") 
stat.test13c
stat.test13c <- stat.test13c %>% add_xy_position(x = "cs.flag.emergent_water")

ggboxplot(cross_merge2, x = "cs.flag.emergent_water" , y = "Mass", fill = "cs.flag.emergent_water") +
  stat_pvalue_manual(stat.test13c, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none")+ scale_fill_viridis(discrete = TRUE)

ggsave("Boxplot_Mass_water_merged2.png", dpi=300, width = 6, height = 4)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#Date: 2022-05-12

#Andrew's idea, need to work on that

# normalise relative to contribution in pool
cross_merge1$rel_occ_sed <- cross_merge1$perc.occup_sed/sum(cross_merge1$perc.occup_sed)
cross_merge1$rel_occ_water <- cross_merge1$perc.occup_wat/sum(cross_merge1$perc.occup_wat)

m1 <- lm(NOSC~cs.flag.emergent_sed+scale(rel_occ_sed)+scale(rel_occ_water),data=cross_merge1)
summary(m1)

# Andrew is now adding some extra code here
# the goal is to estimate group means that account for differences in occupancy 
# the present analysis assumes that each molecule contributes equally
# but some are found 95% of the time vs 1% of the time, so we want to weight by that contribution to the pool
library(emmeans)
emmeans(m1,'cs.flag.emergent_sed')
# we would then want to repeat these analyses across the different traits 




    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#     DENSITY PLOTS AND LINEAR REGRESSIONS
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


##-> Sediment Merge 1
#cs.flag.emergent_sed

#DBE
dens_sed1=ggplot(cross_merge1, aes(x=DBE, fill=cs.flag.emergent_sed))+geom_density(alpha=0.5)

sed1=dens_sed1+scale_fill_manual(values=c("#4682B4" ,"#B4464B", "#B4AF46"))+theme_bw()+  
  xlab(bquote(DBE))+ ylab(bquote(Density))+ #EDIT!
  theme(aspect.ratio=1, strip.text.x = element_text(size = 12),legend.title=element_blank(), legend.text=element_text(size=14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), axis.title.x = element_text(size=16), axis.text.x = element_text(size = 14))+
  theme(legend.position = c(0.75, 0.8),legend.box="horizontal", legend.background = element_rect(fill = "white"), 
        legend.text=element_text(size=9))

#AImod
dens_sed2=ggplot(cross_merge1, aes(x=AI_Mod, fill=cs.flag.emergent_sed))+geom_density(alpha=0.5)

sed2=dens_sed2+scale_fill_manual(values=c("#4682B4" ,"#B4464B", "#B4AF46"))+theme_bw()+  
  xlab(bquote(AI[mod]))+ ylab(bquote(Density))+ #EDIT!
  theme(aspect.ratio=1, strip.text.x = element_text(size = 12),legend.title=element_blank(), legend.text=element_text(size=14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), axis.title.x = element_text(size=16), axis.text.x = element_text(size = 14))+
  theme(legend.position = c(0.8, 0.85),legend.box="horizontal", legend.background = element_rect(fill = "white"), 
        legend.text=element_text(size=9))+theme(legend.position = "none")

#NOSC
dens_sed3=ggplot(cross_merge1, aes(x=NOSC, fill=cs.flag.emergent_sed))+geom_density(alpha=0.5)

sed3=dens_sed3+scale_fill_manual(values=c("#4682B4" ,"#B4464B", "#B4AF46"))+theme_bw()+  
  xlab(bquote(NOSC))+ ylab(bquote(Density))+ #EDIT!
  theme(aspect.ratio=1, strip.text.x = element_text(size = 12),legend.title=element_blank(), legend.text=element_text(size=14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), axis.title.x = element_text(size=16), axis.text.x = element_text(size = 14))+
  theme(legend.position = c(0.8, 0.85),legend.box="horizontal", legend.background = element_rect(fill = "white"), 
        legend.text=element_text(size=9))+theme(legend.position = "none")

#Mass
dens_sed4=ggplot(cross_merge1, aes(x=Mass, fill=cs.flag.emergent_sed))+geom_density(alpha=0.5)

sed4=dens_sed4+scale_fill_manual(values=c("#4682B4" ,"#B4464B", "#B4AF46"))+theme_bw()+  
  xlab(bquote(Mass))+ ylab(bquote(Density))+ #EDIT!
  theme(aspect.ratio=1, strip.text.x = element_text(size = 12),legend.title=element_blank(), legend.text=element_text(size=14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), axis.title.x = element_text(size=16), axis.text.x = element_text(size = 14))+
  theme(legend.position = c(0.8, 0.85),legend.box="horizontal", legend.background = element_rect(fill = "white"), 
        legend.text=element_text(size=9))+theme(legend.position = "none")


#Arrange Plots same figure
grob_sed1=gridExtra::grid.arrange(sed1, sed2, sed3, sed4, nrow=2, top = "Emergent Merge 1 Sediment")
setwd("C:/Users/micha/OneDrive/Documentos/GitHub/Topic1/7_molecular.traits/Boxplots")
ggsave("Panel_Density_sediment_merge1.png", dpi=300, height = 7, width=8, grob_sed1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##-> Water Merge 1
#cs.flag.emergent_sed

#DBE
dens_wat1=ggplot(cross_merge1, aes(x=DBE, fill=cs.flag.emergent_water))+geom_density(alpha=0.5)

wat1=dens_wat1+scale_fill_manual(values=c("#4682B4" ,"#B4464B", "#B4AF46"))+theme_bw()+  
  xlab(bquote(DBE))+ ylab(bquote(Density))+ #EDIT!
  theme(aspect.ratio=1, strip.text.x = element_text(size = 12),legend.title=element_blank(), legend.text=element_text(size=14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), axis.title.x = element_text(size=16), axis.text.x = element_text(size = 14))+
  theme(legend.position = c(0.75, 0.8),legend.box="horizontal", legend.background = element_rect(fill = "white"), 
        legend.text=element_text(size=9))

#AImod
dens_wat2=ggplot(cross_merge1, aes(x=AI_Mod, fill=cs.flag.emergent_water))+geom_density(alpha=0.5)

wat2=dens_wat2+scale_fill_manual(values=c("#4682B4" ,"#B4464B", "#B4AF46"))+theme_bw()+  
  xlab(bquote(AI[mod]))+ ylab(bquote(Density))+ #EDIT!
  theme(aspect.ratio=1, strip.text.x = element_text(size = 12),legend.title=element_blank(), legend.text=element_text(size=14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), axis.title.x = element_text(size=16), axis.text.x = element_text(size = 14))+
  theme(legend.position = c(0.8, 0.85),legend.box="horizontal", legend.background = element_rect(fill = "white"), 
        legend.text=element_text(size=9))+theme(legend.position = "none")

#NOSC
dens_wat3=ggplot(cross_merge1, aes(x=NOSC, fill=cs.flag.emergent_water))+geom_density(alpha=0.5)

wat3=dens_wat3+scale_fill_manual(values=c("#4682B4" ,"#B4464B", "#B4AF46"))+theme_bw()+  
  xlab(bquote(NOSC))+ ylab(bquote(Density))+ #EDIT!
  theme(aspect.ratio=1, strip.text.x = element_text(size = 12),legend.title=element_blank(), legend.text=element_text(size=14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), axis.title.x = element_text(size=16), axis.text.x = element_text(size = 14))+
  theme(legend.position = c(0.8, 0.85),legend.box="horizontal", legend.background = element_rect(fill = "white"), 
        legend.text=element_text(size=9))+theme(legend.position = "none")

#Mass
dens_wat4=ggplot(cross_merge1, aes(x=Mass, fill=cs.flag.emergent_water))+geom_density(alpha=0.5)

wat4=dens_wat4+scale_fill_manual(values=c("#4682B4" ,"#B4464B", "#B4AF46"))+theme_bw()+  
  xlab(bquote(Mass))+ ylab(bquote(Density))+ #EDIT!
  theme(aspect.ratio=1, strip.text.x = element_text(size = 12),legend.title=element_blank(), legend.text=element_text(size=14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), axis.title.x = element_text(size=16), axis.text.x = element_text(size = 14))+
  theme(legend.position = c(0.8, 0.85),legend.box="horizontal", legend.background = element_rect(fill = "white"), 
        legend.text=element_text(size=9))+theme(legend.position = "none")


#Arrange Plots same figure
grob_wat1=gridExtra::grid.arrange(wat1, wat2, wat3, wat4, nrow=2, top = "Emergent Merge 1 Water")
setwd("C:/Users/micha/OneDrive/Documentos/GitHub/Topic1/7_molecular.traits/Boxplots")
ggsave("Panel_Density_water_merge1.png", dpi=300, height = 7, width=8, grob_wat1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Class Histograms Subdivided


g_sed=ggplot(cross_merge1, aes(x = cs.flag.emergent_sed, fill = Class)) +theme_classic()+
  geom_bar(position = "fill") +labs(title="Sediment", X=" ", y="Relative Contribution (%)")+
  theme(aspect.ratio=1, strip.text.x = element_text(size = 12), axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title.x = element_blank(), axis.text.x = element_text(size = 12))+
    scale_fill_viridis(discrete=TRUE)+theme(legend.position = "none")

g_water=ggplot(cross_merge1, aes(x = cs.flag.emergent_water, fill = Class)) +theme_classic()+
  geom_bar(position = "fill") +labs(title="Water", X=" ", y="Relative Contribution (%)")+
  theme(aspect.ratio=1, strip.text.x = element_text(size = 12),legend.title=element_text(12), legend.text=element_text(size=10), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank() , axis.title.x = element_blank(), axis.text.x = element_text(size = 12))+
  scale_fill_viridis(discrete=TRUE)

g=arrangeGrob(g_sed, g_water, nrow=1, widths=c(2.5,3))
setwd("C:/Users/micha/OneDrive/Documentos/GitHub/Topic1/7_molecular.traits/Boxplots")
ggsave("Panel_Histogram_Class_merge1.png", dpi=300, height = 5, width=8, g)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# addition from AJ Tanentzap
# what proportion of core molecules are identical between sediment and water versus satellite?
sum(cross_merge1$cs.flag.emergent_sed == 'Core' & cross_merge1$cs.flag.emergent_water == 'Core') / sum(cross_merge1$cs.flag.emergent_sed == 'Core' | cross_merge1$cs.flag.emergent_water == 'Core')
sum(cross_merge1$cs.flag.emergent_sed == 'Satellite' & cross_merge1$cs.flag.emergent_water == 'Satellite') / sum(cross_merge1$cs.flag.emergent_sed == 'Satellite' | cross_merge1$cs.flag.emergent_water == 'Satellite')
# is the percentage of identical core compounds higher when we look at those that are everywhere
sum((cross_merge1$cs.flag.emergent_sed == 'Core' & cross_merge1$perc.occup_sed == 100) & (cross_merge1$cs.flag.emergent_water == 'Core' & cross_merge1$perc.occup_water == 100)) /
	sum((cross_merge1$cs.flag.emergent_sed == 'Core' & cross_merge1$perc.occup_sed == 100) | (cross_merge1$cs.flag.emergent_water == 'Core' & cross_merge1$perc.occup_water == 100))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#DO NOT RUN



#Plot showing number of samples per boxplot

#Add number samples and mean in each boxplot
#stat_box_data <- function(y, upper_limit = max(cross_merge1$DBE) * 1.15) {
#  return( 
#    data.frame(
#      y = 0.95 * upper_limit,
#      label = paste('count =', length(y), '\n',
#                    'mean =', round(mean(y), 1), '\n')
#    )
#  )
#}

#plot_DBE= ggplot(data=cross_merge1, aes(x=cs.flag.emergent_sed, y=DBE))+   
#  geom_boxplot() +stat_summary(
#  fun.data = stat_box_data, 
#  geom = "text", 
#  hjust = 0.5,
#  vjust = 0.9
#)+  geom_jitter(color="black", size=0.4, alpha=0.9)

#It is difficult to compare because the number of formula among groups is very different

#Violin plot
#ggplot(data=cross_merge1, aes(x=cs.flag.emergent_sed, y=DBE))+   
#  geom_violin() +stat_summary(
#    fun.data = stat_box_data, 
#    geom = "text", 
#    hjust = 0.5,
#    vjust = 0.9
#  )+  geom_jitter(color="black", size=0.4, alpha=0.9)
#ggsave("Violin_DBE_sediment_merged1.png", dpi=300)
