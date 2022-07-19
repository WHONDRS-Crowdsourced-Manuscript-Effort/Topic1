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

### Import data

#Path --> https://github.com/WHONDRS-Crowdsourced-Manuscript-Effort/Topic1/tree/main/4_gather.thresholds
setwd("")


# Tables used to generate figures

cross_merge1= read.csv('FTICR_crosstable_rep.merged1_all_em.thres_2022-05-05.csv')
names(cross_merge1)



### Remove rows with Flags=NA
cross_merge1=subset(cross_merge1, !is.na(cs.flag.emergent_sed))
cross_merge1=subset(cross_merge1, !is.na(cs.flag.emergent_water))



### PLot figures
g_sed=ggplot(cross_merge1, aes(x = cs.flag.emergent_sed, fill = Class)) +theme_classic()+
  geom_bar(position = "fill") +labs(title="Sediment", X=" ", y="Relative Contribution (%)")+
  theme(aspect.ratio=1, strip.text.x = element_text(size = 12), axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title.x = element_blank(), axis.text.x = element_text(size = 12))+
  scale_fill_viridis(discrete=TRUE)+theme(legend.position = "none")

g_water=ggplot(cross_merge1, aes(x = cs.flag.emergent_water, fill = Class)) +theme_classic()+
  geom_bar(position = "fill") +labs(title="Water", X=" ", y="Relative Contribution (%)")+
  theme(aspect.ratio=1, strip.text.x = element_text(size = 12),legend.title=element_text(12), legend.text=element_text(size=10), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank() , axis.title.x = element_blank(), axis.text.x = element_text(size = 12))+
  scale_fill_viridis(discrete=TRUE)


### Arrange Pannel and save figure!
g=arrangeGrob(g_sed, g_water, nrow=1, widths=c(2.5,3))

ggsave("Panel_Histogram_Class_merge1.png", dpi=300, height = 5, width=8, g)
