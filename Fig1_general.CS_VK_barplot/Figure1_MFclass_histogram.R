#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   Figure 1: Class Histograms Subdivided  #

#Author: Michaela de Melo
#Project: Core-Sat DOM-WHONDRS
#Date: July 2022
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

### Load libraries
library(ggplot2)
library(vegan)
library(viridis)
library(gridExtra)
library(pals)
library(dplyr)


### Import data

#Path --> https://github.com/WHONDRS-Crowdsourced-Manuscript-Effort/Topic1/tree/main/4_gather.thresholds
setwd("")


# Tables used to generate figures

cross_merge1= read.csv('FTICR_crosstable_rep.merged1_all_em.thres_2022-05-05.csv')
names(cross_merge1)


### Remove rows with Flags=NA
cross_merge1=subset(cross_merge1, !is.na(cs.flag.emergent_sed))
cross_merge1=subset(cross_merge1, !is.na(cs.flag.emergent_water))

#Recode to spell out the whole name of molecular classes
cross_merge1$Class=recode(cross_merge1$Class, "ConHC" = "Condensed aromatics", "UnsatHC" = "Unsaturated hydrocarbons", "Carb" = "Carbohydrates")

# Order classes 
cross_merge1$Class<- factor(cross_merge1$Class, levels = c("AminoSugar", "Carbohydrates", "Condensed aromatics", "Lignin", "Lipid", "Protein", "Tannin", "Unsaturated hydrocarbons", "Other"))

### PLot figures
g_sed=ggplot(cross_merge1, aes(x = cs.flag.emergent_sed, fill = Class)) +theme_classic()+
  geom_bar(position = "fill") +labs(title="Sediment", X=" ", y="Relative Contribution ")+
  theme(aspect.ratio=1, strip.text.x = element_text(size = 12), axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title.x = element_blank(), axis.text.x = element_text(size = 12))+
  scale_fill_manual(values = c("#88CCEE", "#CC6677", "#DDCC77","#44AA99", "#117733", "#332288", "#AA4499", "#999933", "#661100"))+
  theme(legend.position = "none")

g_water=ggplot(cross_merge1, aes(x = cs.flag.emergent_water, fill = Class)) + theme_classic()+
  geom_bar(position = "fill") +labs(title="Water", X=" ", y="Relative Contribution (%)") +
  theme(aspect.ratio=1, strip.text.x = element_text(size = 12),  axis.text.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(size = 12),  axis.ticks.y = element_blank(), axis.title.y = element_blank())+
  scale_fill_manual(values = c("#88CCEE", "#CC6677", "#DDCC77","#44AA99", "#117733", "#332288", "#AA4499", "#999933", "#661100"))
  

### Arrange Pannel and save figure!
g=arrangeGrob(g_sed, g_water, nrow=1, widths=c(2.5,3.6))

ggsave("Fig.1_MFclass_histogram_FINAL_08-22.png", dpi=300, height = 5, width=8, g)

###################### End #############################

#PS:
#Checking percentage of each class to describe in the results section

##Sediment
ggplot(cross_merge1, aes(x = cs.flag.emergent_sed, fill = Class)) +theme_classic()+
  geom_bar(position = "fill") +labs(title="Sediment", X=" ", y="Relative Contribution (%)")+
  theme(aspect.ratio=1, strip.text.x = element_text(size = 12), axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title.x = element_blank(), axis.text.x = element_text(size = 12))+
  scale_fill_manual(values = c("#88CCEE", "#CC6677", "#DDCC77","#44AA99", "#117733", "#332288", "#AA4499", "#999933", "#661100"))+
  theme(legend.position = "none")+
  geom_text(
    aes(label=signif(..count.. / tapply(..count.., ..x.., sum)[as.character(..x..)]*100, digits=3)),
    stat="count",
    position=position_fill(vjust=0.5))

##Water
ggplot(cross_merge1, aes(x = cs.flag.emergent_water, fill = Class)) +theme_classic()+
  geom_bar(position = "fill") +labs(title="Water", X=" ", y="Relative Contribution (%)")+
  theme(aspect.ratio=1, strip.text.x = element_text(size = 12), axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title.x = element_blank(), axis.text.x = element_text(size = 12))+
  scale_fill_manual(values = c("#88CCEE", "#CC6677", "#DDCC77","#44AA99", "#117733", "#332288", "#AA4499", "#999933", "#661100"))+
  theme(legend.position = "none")+
  geom_text(
    aes(label=signif(..count.. / tapply(..count.., ..x.., sum)[as.character(..x..)]*100, digits=3)),
    stat="count",
    position=position_fill(vjust=0.5))
