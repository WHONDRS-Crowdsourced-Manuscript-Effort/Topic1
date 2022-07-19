rm(list=ls());graphics.off()
library (vegan); library(ggrepel)
input.path = ("C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/Topic1/6_mulvar/PCA_loadings/")
setwd(input.path)
output.path = ("C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/Topic1/6_mulvar/PCA_loadings/")

# Read in the data

df = read.csv("FTICR_crosstable_rep.merged2_all_em.thres_2022-03-23.csv")

df = df[!is.na(df$cs.flag.emergent_sed),]

data = df[grep("Core|Satellite",df$cs.flag.emergent_sed),]

row.names(data) = data$Mass
data = data[,2:21] # Removing the column with core and Satellite specifications from the PCA and removing the first column since it is the mass

meta = cbind.data.frame(Mass = df$Mass,Category = df$cs.flag.emergent_sed)
pca = prcomp(data, center = T,scale = TRUE)


# Selecting plot objects for the PCA
scores.obj = data.frame(Mass = as.numeric(row.names(pca$x)), pca$x)%>% left_join(meta, by = "Mass")
arrow.obj = data.frame(Variable = row.names(pca$rotation), pca$rotation*15)
signif.obj = summary(pca)$importance

# Plotting the PCA
scores.obj %>%
  ggplot(aes(x = PC1, y = PC2))+
  geom_point(aes(fill = Category), shape = 21, size = 3, color = "black", alpha = 0.5)+
  geom_segment(data = arrow.obj, aes(x = 0, xend = PC1, y = 0, yend = PC2), 
               color = "black", lwd = 0.8, arrow = arrow(length = unit(0.3, "cm")))+
  geom_text_repel(data = arrow.obj, aes(x = PC1, y = PC2, 
                                        label = Variable), color = "black",size = 5)+
  xlab(paste0("PC1 (", round(signif.obj[2,1]*100, digits = 2), "%)"))+
  ylab(paste0("PC2 (", round(signif.obj[2,2]*100, digits = 2), "%)"))+
  scale_colour_viridis_c()+
  theme_bw() + theme(text = element_text(size = 12),
                     axis.title = element_text(color = "black", size = 14),
                     axis.text = element_text(color = "black"),
                     axis.ticks = element_line(color = "black"),
                     panel.border = element_rect(size = 1, color = "black"),
                     panel.background = element_blank(),
                     panel.grid = element_blank())

ggsave(paste0(output.path,"PCA_Molecular_characteristics_emergent_sediment_",Sys.Date(),".pdf"))