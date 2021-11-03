# Merge replicates

### Packages -------------------------------------------------------------------------------
pckgs <- list("plyr","tidyverse","data.table")

### Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
# install missing packages with following line. Please uncomment:
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)

### Load packages
invisible(lapply(pckgs, library, character.only = T))
rm(pckgs, miss.pckgs)

# Read in data ----------------------------------------------------------------------------
# Pre-processed data:

# Community matrix
com.mat <- read.csv("./Dataset/FTICR_commat_2021-09-29.csv",
            sep = ",", stringsAsFactors = F) %>% setDT()
com.mat[1:5,1:5]

# cross table
cross <- read.csv("./Dataset/FTICR_crosstable_2021-09-29.csv",
            sep = ",", stringsAsFactors =  F) %>% setDT()
head(cross)

meta <- read.csv("./Dataset/FTICR_meta_all_2021-09-29.csv",
                 sep = ",", stringsAsFactors = F) %>% setDT()
head(meta)


# Remove rare peaks ---------------------------------------------------------------------
# Check if there are several peaks within the same molecular formula (= Mass)

mult.peaks <- cross[, .(n.peaks = .N), by = .(MolForm)]
mult.peaks[n.peaks > 1,]
# no, good we can merge peaks by molecular formula
rm(mult.peaks)

# Remove MF that only appear in <= 2 samples
# melt community matrix for easy classification
melt.commat <- melt(com.mat, id.vars = "ID", value.name = "PA", variable.name = "MolForm")
# enumerate frequency of each molecular formulae
melt.commat[, freq := nrow(.SD[PA > 0,]), by = .(MolForm)]
length(unique(melt.commat[freq <= 2,]$MolForm))
# 18225 molecular formulae of...
length(unique(melt.commat$MolForm))
# 37528 molecular formulae appear <= 2 times in the whole dataset, which means...
length(unique(melt.commat[freq <= 2,]$MolForm)) * 100/ length(unique(melt.commat$MolForm))
# We will loose 48.56 % of the data

# Ok, we will make one dataset with this 2 sample rarity cutoff
rar2 <- melt.commat[freq > 2,]

# Another rarity cutoff could be, 1 sample
length(unique(melt.commat[freq < 2,]$MolForm))
# 13720 molecular formulae of...
length(unique(melt.commat$MolForm))
# 37528 molecular formulae appear <= 2 times in the whole dataset, which means...
length(unique(melt.commat[freq < 2,]$MolForm)) * 100/ length(unique(melt.commat$MolForm))
# We will loose 36.56 % of the data

# We will make one dataset with a 1 sample rarity cutoff
rar1 <- melt.commat[freq >= 2,]

# So, we do the replicate merging for each site/sample type across three datasets:
# 1) Whole dataset (no rarity cutoff, "all")
# 2) Remove MF that are only present in 1 site ("rar1")
# 3) Remove MF that are only present in 2 sites ("rar2")

# Put the three dataset into list for easy processing
com.ls <- list(melt.commat, rar1, rar2)

# Write function to merge replicates

# Make function, could not fine one in base R
# Create function that can merge binary (presence absence) results
as.binary <- function(x){
  x <- ifelse(sum(x) >= 1, 1, 0)
  return(x)
}

merged.ls <- llply(com.ls, function(x){
   # first we split the ID column
  # tidyr::separate is too heavy, let's do a data.table solution
  #x <- x %>% separate(col = "ID", into = c("sample.type","sample.name", "location.id", "extra"), sep = "_", remove = F)
  x <- x[,c("sample.type","sample.name", "location.id", "extra") := 
      list(sapply(str_split(ID, pattern = "_"),"[[",1),
           sapply(str_split(ID, pattern = "_"),"[[",2),
           sapply(str_split(ID, pattern = "_"),"[[",3),
           sapply(str_split(ID, pattern = "_"),"[[",4))]
  # Sediments = replicates are in col 'location.id',
  # Surface waters = replicates are in col 'extra'
  merged <- x[, .(PA = as.binary(PA)), by = .(sample.type, sample.name, MolForm)]
  # sanity check
  #length(unique(merged$MolForm)) == length(unique(x$MolForm))
  
  # create new ID
  merged[, ID := paste(sample.type, sample.name, sep = "_")]
  # cast into wide format
  merged.commat <- dcast(merged, ID ~ MolForm, value.var = "PA")
  return(merged.commat)
})

# next, get cross table for each subset
cross.ls <- llply(merged.ls, function(x){
  subcross <- cross[MolForm %in% colnames(x),]
  return(subcross)
})

# Sanity check
# melt.commat <- melt(merged.ls[[1]], id.vars = "ID", value.name = "PA", variable.name = "MolForm")
# # any not binary?
# melt.commat[PA > 1, ] # no

# Revert to original meta data ---------------------------------------------------------------------------------------
# Get original meta data
header <- colnames(read.csv("./Dataset/Surface/WHONDRS_S19S_Metadata_v2.csv",
                   sep = ",", stringsAsFactors = F, header = T))
surface <- read.csv("./Dataset/Surface/WHONDRS_S19S_Metadata_v2.csv",
                    sep = ",", stringsAsFactors = F, header = F, skip = 2) %>% setDT()
colnames(surface) <- header

header <- colnames(read.csv("./Dataset/Sediment/WHONDRS_S19S_Metadata_v3.csv",
                            sep = ",", stringsAsFactors = F, header = T))
sed <- read.csv("./Dataset/Sediment/WHONDRS_S19S_Metadata_v3.csv",
                    sep = ",", stringsAsFactors = F, skip = 2, header = F) %>% setDT()
colnames(sed) <- header

# merge
colnames(surface) %in% colnames(sed)
colnames(sed)[!(colnames(sed) %in% colnames(surface))] # Strahler order missing in surface

# clean
# replace "Not_Provided" with NA to allow correct column types

# merge the stream order from sediment meta
surface <- surface[sed,c("Stream_Order") := list(i.Stream_Order), on = .(Sample_ID)]
meta <- surface
rm(surface, sed)

allcols <- colnames(meta)
meta[, (allcols) := lapply(.SD, function(x) ifelse(x == "Not_Provided", NA, x)), .SD = allcols]

# rename Sample_ID to river.id and replace _ with .
#meta[, river.id := str_replace(Sample_ID, "_", ".")]

# pH, temp, DO are character, clean
# take first value in pH, except when there is a ">"
meta[str_detect(SW_pH, ">"), new.pH := str_extract(SW_pH, pattern = "\\d+\\.*\\d*")]
meta[!str_detect(SW_pH, ">"), new.pH := str_extract(SW_pH, pattern = "^\\d+\\.*\\d*")]
meta[, SW_pH := as.numeric(new.pH)][, new.pH := NULL]

# For temperature, there are ranges, take first value
meta[, SW_Temp_degC := as.numeric(str_extract(SW_Temp_degC, pattern = "^\\d+\\.*\\d*"))]

# comments in DO, save in separate column
meta[str_detect(DO_perc.sat, pattern = "calibrated"), DO_comment := "may not be correctly calibrated"]
meta[, DO_perc.sat := as.numeric(sapply(str_split(DO_perc.sat, pattern = " "), "[[",1))]
meta[, DO_mg.per.L := as.numeric(sapply(str_split(DO_mg.per.L, pattern = " "), "[[",1))]

# duplicate for both SED and SW
meta <- bind_rows(meta[, sample.type := "SW"],meta[, sample.type := "SED"])

# rename Sample_ID to river.id and replace _ with .
meta[, river.id := str_replace(Sample_ID, "_", ".")]
meta[, ID := paste(sample.type, river.id, sep = "_")]

# re-order
meta <- meta %>% dplyr::select(Study_Code, ID, sample.type:river.id, Date:DO_comment)

# Save -------------------------------------------------------------------------------------------------------
# Save meta data
write.table(meta,
            paste0("./Dataset/FTICR_meta_eachriver_", Sys.Date(),".csv"), sep = ",", dec = ".", row.names = F)

# Save different rarity thresholds
# create list with names
save.vec <- c(paste0("./Dataset/FTICR_commat_rep.merged_all_",Sys.Date(),".csv"),
  paste0("./Dataset/FTICR_commat_rep.merged_rar1_",Sys.Date(),".csv"),
  paste0("./Dataset/FTICR_commat_rep.merged_rar2_",Sys.Date(),".csv"))

# Save community matrix
mapply(write.table, merged.ls, save.vec, MoreArgs = list(sep = ",", dec = ".", row.names = F))

# Save cross table
# create list with names
save.vec <- c(paste0("./Dataset/FTICR_cross.table_rep.merged_all_",Sys.Date(),".csv"),
              paste0("./Dataset/FTICR_cross.table_rep.merged_rar1_",Sys.Date(),".csv"),
              paste0("./Dataset/FTICR_cross.table_rep.merged_rar2_",Sys.Date(),".csv"))

mapply(write.table, cross.ls, save.vec, MoreArgs = list(sep = ",", dec = ".", row.names = F))

# Done.
