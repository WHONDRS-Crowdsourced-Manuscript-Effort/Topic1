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
  scale_fill_manual(values = c("#999999", "#661100", "#0072B2", "#FFB000",  "#FE6100", "#648FFF", "#DC267F" ))


plot2=ggboxplot(cross, x = "cs.flag.emergent_overlap" , y = "DBE", fill = "cs.flag.emergent_overlap") +
  #stat_pvalue_manual(stat.test, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.text.x=element_text(angle=45), axis.title.x=element_blank())+ 
  scale_fill_manual(values = c("#999999", "#661100", "#0072B2", "#FFB000",  "#FE6100", "#648FFF", "#DC267F" ))


plot3=ggboxplot(cross, x = "cs.flag.emergent_overlap" , y = "AI_Mod", fill = "cs.flag.emergent_overlap") +
  #stat_pvalue_manual(stat.test, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.text.x=element_text(angle=45), axis.title.x=element_blank())+ 
  scale_fill_manual(values = c("#999999", "#661100", "#0072B2", "#FFB000",  "#FE6100", "#648FFF", "#DC267F" ))


plot4=ggboxplot(cross, x = "cs.flag.emergent_overlap" , y = "Mass", fill = "cs.flag.emergent_overlap") +
  #stat_pvalue_manual(stat.test, hide.ns = FALSE) +
  theme_bw() +theme(legend.position="none", axis.text.x=element_text(angle=45), axis.title.x=element_blank())+ 
  scale_fill_manual(values = c("#999999", "#661100", "#0072B2", "#FFB000",  "#FE6100", "#648FFF", "#DC267F" ))



##Arrange and save it!
g=arrangeGrob(plot1, plot2, plot3, plot4,  ncol=2)
ggsave("Fig.5_Panel_boxplots.png", dpi=300, height = 10, width=12, g)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##Statistics
# Length of data
levels(cross$cs.flag.emergent_overlap)
aggregate( NOSC~cs.flag.emergent_overlap, data =cross, FUN=length)
aggregate( DBE~cs.flag.emergent_overlap, data =cross, FUN=length)
aggregate( AI_Mod~cs.flag.emergent_overlap, data =cross, FUN=length)
aggregate( Mass~cs.flag.emergent_overlap, data =cross, FUN=length)

# Checking assumptions 
##Levene test: quality of variances
leveneTest(NOSC ~ cs.flag.emergent_overlap, data = cross) #p<0.001

#Shapiro wilk: normality test
#Test by category
global_core=filter(cross, cs.flag.emergent_overlap=="Global core")
shapiro.test(global_core$NOSC)
shapiro.test(global_core$DBE)
shapiro.test(global_core$AI_Mod)
shapiro.test(global_core$Mass)
#Not normal, I will use non-parametric test

#Kruskal-Wallis: Non-parametric test
kruskal.test(NOSC ~ cs.flag.emergent_overlap, data = cross) #p-value<0.001
kruskal.test(DBE ~ cs.flag.emergent_overlap, data = cross) #p-value<0.001
kruskal.test(AI_Mod ~ cs.flag.emergent_overlap, data = cross) #p-value<0.001
kruskal.test(Mass ~ cs.flag.emergent_overlap, data = cross) #p-value<0.001

#If the results of a Kruskal-Wallis test are statistically significant, then it's appropriate to conduct Dunn's Test to determine exactly which groups are different.
dunnTest(NOSC ~ cs.flag.emergent_overlap, data = cross,   method="holm") #all are different

dunnTest(DBE ~ cs.flag.emergent_overlap, data = cross,   method="holm") 

#                    Comparison                                Z       P.unadj         P.adj
# Sed Sat - Water Inbetween - Water Core - Sed Sat/In-bet.  -0.2852473  7.754547e-01  7.754547e-01

dunnTest(AI_Mod ~ cs.flag.emergent_overlap, data = cross,   method="holm") 
# Global core - Global in-between  -1.123171  2.613650e-01  2.613650e-01
#Global core - Global satellite   1.482623  1.381747e-01  2.763494e-01

dunnTest(Mass ~ cs.flag.emergent_overlap, data = cross,   method="holm") #all are different
#Global core - Sed Core - Water Sat/In-bet.  -1.9796630  4.774141e-02  1.432242e-01
#Global in-between - Sed Core - Water Sat/In-bet.   1.4842431  1.377444e-01  2.754889e-01
#Global satellite - Sed Sat - Water Inbetween   0.5817314  5.607476e-01  5.607476e-01

#Plots with statistics
plot1a=plot1+annotate("text", x = 1:7, y = 3, label =c("n=666", "n=543", "n=2645", "n=524", "n=1367", "n=1151", "n=881"))
plot2a=plot2+annotate("text", x = 5, y = 20, label = "a")+annotate("text", x = 6, y = 21, label = "a")
plot3a=plot3+annotate("text", x = 1, y = 1.1, label = "ab")+annotate("text", x = 2, y = 1.1, label = "a")+
annotate("text", x = 3, y = 1.1, label = "b")
plot4a=plot4+annotate("text", x = 1, y = 905, label = "a")+annotate("text", x = 4, y = 905, label = "ab")+
  annotate("text", x = 2, y = 905, label = "b")+ 
  annotate("text", x = 3, y = 905, label = "c")+ annotate("text", x = 6, y = 905, label = "c")


##Arrange and save it!
ga=arrangeGrob(plot1a, plot2a, plot3a, plot4a,  ncol=2)
ggsave("Fig.5_Panel_boxplots_statistics.png", dpi=300, height = 10, width=12, ga)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# END

