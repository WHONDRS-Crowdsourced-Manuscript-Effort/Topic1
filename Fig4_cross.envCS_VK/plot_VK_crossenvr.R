### Packages -------------------------------------------------------------------------------
pckgs <- list("plyr","tidyverse","data.table","doMC")

### Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
# install missing packages with following line. Please uncomment:
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)

### Load packages
invisible(lapply(pckgs, library, character.only = T))
rm(pckgs, miss.pckgs)

# Read in data-----------------------------------------------------------------------------
cross <- read.csv("./4_gather.thresholds/FTICR_crosstable_rep.merged1_all_em.thres_2022-05-05.csv",
                  sep = ",", dec = ".", stringsAsFactors = F) %>% setDT()

cross[, .(n = .N), by = .(cs.flag.emergent_sed)]
commat <- read.csv("./4_gather.thresholds/FTICR_commat_rep.merged1_2022-07-19.csv",
                   sep = ",", dec = ".", stringsAsFactors = F)

# only keep MF in cross table
commat <- commat[, which(colnames(commat) %in% cross$MolForm)]

cross[, .(n = .N), by = .(cs.flag.emergent_water)]
cross[, .(n = .N), by = .(cs.flag.emergent_sed)]

# add compound group lines
comp_groups <- data.frame(group = c("Saturated fatty acids",
                                    "Peptides and\nunsaturated\naliphatics",
                                    "Highly unsaturated\ncompounds",
                                    "Vascular plant-\nderived polyphenols\nand phenols"),
                          HtoC_ratio = c(2, 1.7, 1.25, 0.5),
                          cs.flag.emergent_overlap = "Global core") %>% setDT()

# re-label for plotting
comp_groups[, cs.flag.emergent_overlap := factor(cs.flag.emergent_overlap,
                                                 levels = c("Global core","Global in-between",
                                                            "Global satellite",
                                                            "Sed Core - Water Sat", "Sed Core - Water Inbetween", 
                                                            "Sed Sat - Water Core", "Sed Inbetween - Water Core",
                                                            "Sed Sat - Water Inbetween", "Sed Inbetween - Water Sat"),
                                                 labels = c("Global core","Global in-between",
                                                            "Global satellite",
                                                            "Sed Core - Water Sat/In-bet.", "Sed Core - Water Sat/In-bet.", 
                                                            "Water Core - Sed Sat/In-bet.", "Water Core - Sed Sat/In-bet.", 
                                                            "Sed Sat - Water Inbetween", "Water Sat - Sed Inbetween"))]

cross[, cs.flag.emergent_overlap := factor(cs.flag.emergent_overlap,
                                           levels = c("Global core","Global in-between",
                                                      "Global satellite",
                                                      "Sed Core - Water Sat", "Sed Core - Water Inbetween", 
                                                      "Sed Sat - Water Core", "Sed Inbetween - Water Core",
                                                      "Sed Sat - Water Inbetween", "Sed Inbetween - Water Sat"),
                                           labels = c("Global core","Global in-between",
                                                      "Global satellite",
                                                      "Sed Core - Water Sat/In-bet.", "Sed Core - Water Sat/In-bet.", 
                                                      "Water Core - Sed Sat/In-bet.", "Water Core - Sed Sat/In-bet.", 
                                                      "Sed Sat - Water Inbetween", "Water Sat - Sed Inbetween"))]

# Calculate proporiton of compound groups
cross[is.na(cs.flag.emergent_sed), cs.flag.emergent_overlap := "Water unique"]
cross[is.na(cs.flag.emergent_water), cs.flag.emergent_overlap := "Sediment unique"]

cross[, n.by.cross := .N, by = .(cs.flag.emergent_overlap)]
prop.df <- cross[, .(n.by.cross = unique(n.by.cross),
          n = .N), by = .(cs.flag.emergent_overlap, Class)] %>% arrange(cs.flag.emergent_overlap, Class) %>% setDT

prop.df[, cs.flag.emergent_overlap := factor(cs.flag.emergent_overlap,
                                           levels = c("Global core","Global in-between",
                                                      "Global satellite",
                                                      "Sed Core - Water Sat/In-bet.", 
                                                      "Water Core - Sed Sat/In-bet.",
                                                      "Sed Sat - Water Inbetween", "Water Sat - Sed Inbetween"),
                                           labels = c("Global core","Global in-between",
                                                      "Global satellite",
                                                      "Sed Core -\nWater Sat/In-bet.", 
                                                      "Water Core -\nSed Sat/In-bet.",
                                                      "Sed Sat -\nWater Inbetween", "Water Sat -\nSed Inbetween"))]

prop.df[, Class := factor(Class, levels = c("AminoSugar", "Carb", "ConHC", "Lignin", "Lipid", "Protein", "Tannin",
                                     "UnsatHC", "Other"),
                          labels = c("Amino Sugar", "Carbohydrates", "Condensed\naromatics", "Lignin",
                                     "Lipid", "Protein", "Tannin", "Unsaturated\nhydrocarbons", "Other"))]

prop.df[, prop := (n * 100) / n.by.cross]

# sanity check
prop.df[, .(sum = sum(prop)), by = .(cs.flag.emergent_overlap)]
prop.df[, .(sum = sum(n),
            n = unique(n.by.cross)), by = .(cs.flag.emergent_overlap)]
prop.df <- prop.df[cs.flag.emergent_overlap != "Water unique" & cs.flag.emergent_overlap != "Sediment unique", ]

library(ggpubr)

(prop.p <- ggplot(prop.df, aes(x = cs.flag.emergent_overlap, y = prop, fill = Class)) +
  theme_bw() +
  scale_fill_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#44AA99",
                               "#117733", "#332288", "#AA4499", "#999933", "#661100")) +
  geom_col() +
  labs(x = "", y = "Relative Contribution") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = 'black'),
        legend.position = "right", panel.grid = element_blank()))


# Island of stability
ios.df <- data.frame(HtoC_ratio = c(1.3,1.04), OtoC_ratio = c(0.62,0.42))

library(grid)
plot.df <- cross[cs.flag.emergent_overlap != "Water unique" & cs.flag.emergent_overlap != "Sediment unique",]


plot.ls <- list()
groups <- plot.df %>% arrange(cs.flag.emergent_overlap)
groups <- unique(groups$cs.flag.emergent_overlap)
cols <- c("#999999","#661100", "#0072B2","#FFB000", 
                   "#FE6100", "#648FFF","#DC267F")
for(i in 1:length(groups)){
  p <- ggplot(plot.df[cs.flag.emergent_overlap == groups[i],], 
              aes(x = OtoC_ratio, HtoC_ratio)) +
    theme_bw() +
    geom_point(aes(colour = cs.flag.emergent_overlap), alpha = 0.7) +
    scale_colour_manual(values = cols[i]) +
    facet_wrap(.~cs.flag.emergent_overlap, ncol =2) +
    geom_hline(yintercept = 1.5, linetype = "dashed") +
    geom_abline(intercept = 1.1, slope = -0.3, linetype = "dashed") +
    #stat_ellipse(data = ios.df, aes(x = OtoC_ratio, y = HtoC_ratio), type = "norm",
    #              colour = "tomato", linetype = "dashed") +
    labs(x = "O/C", y = "H/C") +
    lims(x = c(0,1.25), y = c(0,2)) +
    theme(panel.grid.minor = element_blank(), axis.title = element_blank()) +
    guides(colour = "none") +
    annotation_custom(grob=circleGrob(r=unit(1,"npc"),
                                      gp = gpar(col = 'black', lty = 3, fill = 'transparent')),
                      xmin=ios.df$OtoC_ratio[2], xmax=ios.df$OtoC_ratio[1],
                      ymin=ios.df$HtoC_ratio[2], ymax=ios.df$HtoC_ratio[1])
  
  if(i == 1){
    (p <- p + geom_text(data = comp_groups, aes(x = 0.9, y = HtoC_ratio, label = group), size = 3, 
                        hjust = 0, lineheight = 0.7) + annotate(geom = "text", x =0.67, y = 1.17, label = "IOS", size = 3))
  }
  
  plot.ls[[i]] <- p
}

library(gtable)


global <- ggarrange(plot.ls[[1]] + theme(axis.text.x = element_blank(),
                                         axis.ticks.x = element_blank()),
                    plot.ls[[2]]+ theme(axis.text.x = element_blank(),
                                        axis.ticks.x = element_blank()),
                    plot.ls[[3]],nrow = 3)

core <- ggarrange(plot.ls[[4]] +theme(axis.text.x = element_blank(),
                                      axis.ticks.x = element_blank()), plot.ls[[5]], nrow = 2)
sat <- ggarrange(plot.ls[[6]]+theme(axis.text = element_blank(),
                                    axis.ticks = element_blank()), 
                 plot.ls[[7]]+theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank()), nrow = 2)

# create a separate legend
leg.p <- get_legend(prop.p + theme(legend.position = "bottom",
                        legend.key.size = unit(0.25, 'cm'),
                        legend.text = element_text(size=8)))

temp <- ggarrange(core,sat, ncol = 2)
temp <- annotate_figure(temp, left = textGrob("H/C", rot = 90, vjust = 1, gp = gpar(cex = 0.8)),
                bottom = textGrob("O/C", gp = gpar(cex = 0.8)), fig.lab = "b.")
(out <- ggarrange(annotate_figure(global,left = textGrob("H/C", rot = 90, vjust = 1, gp = gpar(cex = 0.8)),
                          bottom = textGrob("O/C", gp = gpar(cex = 0.8)),
                          fig.lab = "a."), ggarrange(temp, annotate_figure(prop.p + guides(fill = "none"), 
                                                                                   fig.lab = "c."),
                            nrow = 2, legend.grob = leg.p, legend = "bottom"),
          ncol = 2))


ggsave("./Fig4_cross.envCS_VK/cross.envCS_VK_prop.bar.png", out, dpi = 300, width = 22, height = 20, units = "cm")

# p <- ggplot(plot.df, 
#             aes(x = OtoC_ratio, HtoC_ratio)) +
#   theme_bw() +
#   geom_point(aes(colour = cs.flag.emergent_overlap), alpha = 0.7) +
#   scale_colour_manual(values = c("#999999","#661100", "#0072B2","#FFB000", 
#                                  "#FE6100", "#648FFF","#DC267F")) +
#   facet_wrap(.~cs.flag.emergent_overlap, ncol =2) +
#   geom_hline(yintercept = 1.5, linetype = "dashed") +
#   geom_abline(intercept = 1.1, slope = -0.3, linetype = "dashed") +
#   #stat_ellipse(data = ios.df, aes(x = OtoC_ratio, y = HtoC_ratio), type = "norm",
#   #              colour = "tomato", linetype = "dashed") +
#   labs(x = "O/C", y = "H/C") +
#   lims(x = c(0,1.25), y = c(0,2)) +
#   theme(panel.grid.minor = element_blank()) +
#   annotation_custom(grob=circleGrob(r=unit(1,"npc"),
#                                     gp = gpar(col = 'black', lty = 3, fill = 'transparent')),
#                     xmin=ios.df$OtoC_ratio[2], xmax=ios.df$OtoC_ratio[1],
#                     ymin=ios.df$HtoC_ratio[2], ymax=ios.df$HtoC_ratio[1])
# 
# 
# (p <- p + geom_text(data = comp_groups, aes(x = 0.78, y = HtoC_ratio, label = group), size = 3, 
#                     hjust = 0, lineheight = 0.7))