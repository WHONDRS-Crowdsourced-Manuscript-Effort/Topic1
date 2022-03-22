# exploratory visualisations of differences between core and satellite molecules
# prepared by AJ Tanentzap (22 Mar 2022, ajt65 // @ // cam.ac.uk)

###############################################################################
# 1) load data
###############################################################################
#merge 1 = molecular formula present in at least 1 rep of the 3 reps per site
#merge 2 = at least 2 of the 3 reps
cross_merge1 = read.csv('FTICR_crosstable_rep.merged1_all_em.thres_2022-03-18.csv')
cross_merge2 = read.csv('FTICR_crosstable_rep.merged2_all_em.thres_2022-03-18.csv')


###############################################################################
# 2) Van Krevelen plots to visualise if core formulas are different than satellite ones
###############################################################################
#pdf("VK_emergent_merge2.pdf", width=6.3, height=3.9, useDingbats=F)
par(mfrow=c(1,2))
with(cross_merge2[which(is.na(cross_merge2$cs.flag.emergent_water)),], plot(OtoC_ratio,HtoC_ratio,pch=NA,main='a) water',xlab='H:C',ylab='O:C',las=1,col=alpha('gray', 0.1))) 
with(cross_merge2[which(cross_merge2$cs.flag.emergent_water == 'Satellite'),], points(OtoC_ratio,HtoC_ratio,pch=19,col=alpha('blue', 0.1))) 
with(cross_merge2[which(cross_merge2$cs.flag.emergent_water == 'Core'),], points(OtoC_ratio,HtoC_ratio,pch=19,col=alpha('red', 0.1))) 
points(min(cross_merge2$OtoC_ratio)*c(.9,.9),max(cross_merge2$HtoC_ratio)*c(.92,.99),pch=19,col=c('red','blue'))
text(c(.1,.1),max(cross_merge2$HtoC_ratio)*c(.92,.99),c('Core','Satellite'),adj=0)

with(cross_merge2[which(is.na(cross_merge2$cs.flag.emergent_sed)),], plot(OtoC_ratio,HtoC_ratio,pch=NA,main='b) sediment',xlab='H:C',ylab='O:C',las=1,col=alpha('gray', 0.1))) 
with(cross_merge2[which(cross_merge2$cs.flag.emergent_sed == 'Satellite'),], points(OtoC_ratio,HtoC_ratio,pch=19,col=alpha('blue', 0.1))) 
with(cross_merge2[which(cross_merge2$cs.flag.emergent_sed == 'Core'),], points(OtoC_ratio,HtoC_ratio,pch=19,col=alpha('red', 0.1))) 
#dev.off() 
 
 
###############################################################################
# 2) visualise if core formulas are comprised of different compound classes than satellite ones
###############################################################################
# compile aggregated compound class data
sed_emerge_merge2_comp_class <- with(cross_merge2, table(bs2_class,cs.flag.emergent_sed))
sed_emerge_merge2_comp_class <- t(apply(sed_emerge_merge2_comp_class,1,function(x){x/colSums(sed_emerge_merge2_comp_class)}))

wat_emerge_merge2_comp_class <- with(cross_merge2, table(bs2_class,cs.flag.emergent_water))
wat_emerge_merge2_comp_class <-  t(apply(wat_emerge_merge2_comp_class,1,function(x){x/colSums(wat_emerge_merge2_comp_class)}))

#pdf("compound_classes_emergent_merge2.pdf", width=9, height=3.5)
par(mfrow=c(1,3),oma = c(0,0,0,0),mar=c(4.5,4,4.5,0))
barplot(sed_emerge_merge2_comp_class*100, col = topo.colors(10),las=1,ylab='composition (%)')
mtext('a) sediment',3,adj=0)
barplot(wat_emerge_merge2_comp_class*100, col = topo.colors(10), yaxt='n')
mtext('b) water',3,adj=0)
par(oma = c(0,0,0,8), mar = c(0,0,0,0), new = TRUE)
legend("right", rownames(wat_emerge_merge2_comp_class), fill=topo.colors(10),bty="n",y.intersp=3)
#dev.off() 
 

###############################################################################
# 3) plot if occupancy of MFs in sediment and water is correlated and differently between C vs S
###############################################################################
#pdf("occupancy_by_type_source.pdf", width=6, height=6)
with(cross_merge2, plot(perc.occup_sed,perc.occup_water,pch=NA,las=1,xlab='sediment occupancy (%)',ylab='water occupancy (%)'))
abline(0,1,lwd=2)
with(cross_merge2[which(cross_merge2$cs.flag.emergent_sed == 'Core' & cross_merge2$cs.flag.emergent_water == 'Core'),], points(perc.occup_sed,perc.occup_water,pch=19,col=alpha('red', 0.1)))
with(cross_merge2[which(cross_merge2$cs.flag.emergent_sed == 'Satellite' & cross_merge2$cs.flag.emergent_water == 'Satellite'),], points(perc.occup_sed,perc.occup_water,pch=19,col=alpha('blue', 0.1)))
with(cross_merge2[which(cross_merge2$cs.flag.emergent_sed == 'In-between' & cross_merge2$cs.flag.emergent_water == 'In-between'),], points(perc.occup_sed,perc.occup_water,pch=19,col=alpha('green', 0.1)))
#dev.off() 
 
