#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   Figure 5: Boxplots Indices 
#(NOSC, DBE, Gibbs, Aromaticity, maybe putative transformations) 
#in relation to cross-environment CS-classification  #

#Author: Michaela de Melo
#Project: Core-Sat DOM-WHONDRS
#Date: July 2022
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### Load libraries
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



#Path --> https://github.com/WHONDRS-Crowdsourced-Manuscript-Effort/Topic1/tree/main/4_gather.thresholds
setwd("")


# Tables used to generate figures
cross=read.csv("FTICR_crosstable_rep.merged1_all_em.thres_2022-05-05.csv")
cross=subset(cross, !is.na(cs.flag.emergent_overlap)) #remove NAs

### Check categories
factor(cross$cs.flag.emergent_overlap) 


## No statistics for now, just exploratory
cross$cs.flag.emergent_overlap= factor(cross$cs.flag.emergent_overlap,
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


plot1=ggboxplot(cross, x = "cs.flag.emergent_overlap" , y = "NOSC", fill = "cs.flag.emergent_overlap") +
  #stat_pvalue_manual(stat.test, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.text.x=element_text(angle=45), axis.title.x=element_blank())+ 
  scale_fill_viridis(discrete = TRUE)

plot2=ggboxplot(cross, x = "cs.flag.emergent_overlap" , y = "DBE", fill = "cs.flag.emergent_overlap") +
  #stat_pvalue_manual(stat.test, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.text.x=element_text(angle=45), axis.title.x=element_blank())+ 
  scale_fill_viridis(discrete = TRUE)

plot3=ggboxplot(cross, x = "cs.flag.emergent_overlap" , y = "AI_Mod", fill = "cs.flag.emergent_overlap") +
  #stat_pvalue_manual(stat.test, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.text.x=element_text(angle=45), axis.title.x=element_blank())+ 
  scale_fill_viridis(discrete = TRUE)

plot4=ggboxplot(cross, x = "cs.flag.emergent_overlap" , y = "Mass", fill = "cs.flag.emergent_overlap") +
  #stat_pvalue_manual(stat.test, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.text.x=element_text(angle=45), axis.title.x=element_blank())+ 
  scale_fill_viridis(discrete = TRUE)


##Arrange and save it!
g=arrangeGrob(plot1, plot2, plot3, plot4,  ncol=2)
ggsave("Fig.5_Panel_boxplots.png", dpi=300, height = 10, width=12, g)

