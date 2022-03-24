PURPOSE:
visualize characteristics of core and satellite species. This is done using swarm plots for a number of properties derived from the molecular formulae, such as double bonds, aromaticity, nominal oxidation state of carbon etc. Figures distinguish water and porewater (sediment) samples (i.e. how do core and satellite properties compare in the water, in the sediment?) and (some) show these properties for individual compound classes (how does the distribution of properties compare between core and satellite speices in a given compound class)


SCRIPTS:
vis_swarm_withPCA_classes.m:
as previous, but now it also generates the swarm plots for the individual compound classes
So far, this is the latest version, so you can ignore the others

vis_swarm_withPCA_SatelliteFixed.m: 
as previous, but tapping into the data file that has the Satellite annotation fixed

vis_swarm_withPCA: 
original script, with PCA added. The PCA will probably get superseded by work in the multivariate analsyis folder

vis_swarm.m: original script


RESULTS:
I added a few screenshots for some of the results, in case you don't have matlab and the required toolboxes, or don't want to rerun but just glance at the results

swarm_screenshot...: 
swarm plots for properties in core, satellite (and in-between) separated by water/sediment and by method (emergent = second derivate; pca-based; random forest based)

merge1/2_compoundclass:
swarm plots as for swarm_screenshot, but now the plots are done separately for individual compound classes

PCA_Screenshot...: 
PCA biplots 


rep1/rep2 and merge1/merge2 refer to different data sets used: 
The former uses data sets merged1 (considering MFs present in any sample), the latter merged2 (MF present in at least 2 samples)
