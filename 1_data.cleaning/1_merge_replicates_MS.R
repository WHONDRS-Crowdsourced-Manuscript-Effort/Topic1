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
com.mat <- read.csv("./1_data.cleaning/output/FTICR_commat_2021-09-29.csv",
            sep = ",", stringsAsFactors = F) %>% setDT()
com.mat[1:5,1:5]

# cross table
# cross <- read.csv("./1_data.cleaning/output/FTICR_crosstable_2021-09-29.csv",
#             sep = ",", stringsAsFactors =  F) %>% setDT()
# head(cross)

meta <- read.csv("./1_data.cleaning/output/FTICR_meta_all_2021-09-29.csv",
                 sep = ",", stringsAsFactors = F) %>% setDT()
head(meta)

# Add peak data
peaks <- read.csv("./1_data.cleaning/output/peaks/FTICR_raw.peaks_commat_2022-01-18.csv",
                  sep = ",", stringsAsFactors = F) %>% setDT()

# Merge replicates ----------------------------------------------------------------------
# melt community matrix
melt.commat <- melt(com.mat, id.vars = "ID", value.name = "PA", variable.name = "MolForm")

melt.commat <- melt.commat[,c("sample.type","sample.name", "location.id", "extra") := 
         list(sapply(str_split(ID, pattern = "_"),"[[",1),
              sapply(str_split(ID, pattern = "_"),"[[",2),
              sapply(str_split(ID, pattern = "_"),"[[",3),
              sapply(str_split(ID, pattern = "_"),"[[",4))]

# Sediments = replicates are in col 'location.id',
# Surface waters = replicates are in col 'extra'
colnames(melt.commat)[1] <- "old.ID"
# Add new ID
melt.commat[, ID := paste(sample.type, sample.name, sep = "_")]
# calculate the number of replicates for each site
replicates <- melt.commat %>% dplyr::select(-MolForm, -PA) %>% distinct()
replicates <- replicates[, .(n.replicate = .N), by = .(ID)]

# Create replicate column
# Change D, M, U to replicates
melt.commat[sample.type == "SED", replicate.id := factor(location.id, levels = c("D","M","U"),
                                                         labels = c("r1","r2","r3"))]
# Add surface water
melt.commat[sample.type == "SW", replicate.id := factor(extra, levels = c("1","2","3"),
                                                        labels = c("r1","r2","r3"))]

# sanity checks
melt.commat[is.na(replicate.id),] # any rows without a replicate ID? -> no
levels(factor(melt.commat$replicate.id)) # ok, three replicates

# Cast replicates into columns
rep.dt <- dcast(melt.commat, ID + MolForm ~ replicate.id, value.var = "PA")

# Calculate how many replicates observed a MF
rep.dt[, rep.sum := rowSums(.SD, na.rm = T), .SDcols = c("r1","r2","r3")]
rep.dt[replicates, n.rep := i.n.replicate, on = .(ID)]

# Keep two data frames
# 1. Keep observation if MF was found at least in one replicate
# 2. Keep observation if MF was found at least in two replicates
# BUT if there is only one replicate, we will keep the observation

one.rep <- rep.dt[n.rep > 1, PA := ifelse(rep.sum >= 1, 1, 0)] # at least in one replicate
one.rep <- one.rep[n.rep == 1, PA := ifelse(rep.sum == 1, 1, 0)]
# Cast into community matrix
one.commat <- dcast(one.rep, ID ~ MolForm, value.var = "PA")

two.rep <- rep.dt[, PA := NULL]
two.rep <- rep.dt[n.rep > 1, PA := ifelse(rep.sum >= 2, 1, 0)] # at least in two replicates
two.rep <- rep.dt[n.rep == 1, PA := ifelse(rep.sum == 1, 1, 0)]
# Cast into community matrix
two.commat <- dcast(two.rep, ID ~ MolForm, value.var = "PA")

# some sanity check
any(rowSums(one.commat[,-1]) != rowSums(two.commat[,-1])) # TRUE = ok

# Save
write.table(one.commat,
            paste0("./1_data.cleaning/output/FTICR_commat_rep.merged1_", Sys.Date(), ".csv"), sep = ",",
            row.names = F)
write.table(two.commat,
            paste0("./1_data.cleaning/output/FTICR_commat_rep.merged2_", Sys.Date(), ".csv"), sep = ",",
            row.names = F)

# Merge replicate data with meta data
meta[, merge.id := paste(sample.type, river.id, sep = "_")]
meta[replicates, no.replicates := i.n.replicate, on = c("merge.id==ID")]

# Overwrite
write.table(meta, paste0("./1_data.cleaning/output/FTICR_meta_all_", Sys.Date(),".csv"),
            row.names = F, sep = ",")

rm(two.rep, one.rep, melt.commat, rep.dt)

# Merge replicates of peaks -------------------------------------------------------------
peaks <- peaks[,c("sample.type","sample.name", "location.id", "extra") := 
                             list(sapply(str_split(ID, pattern = "_"),"[[",1),
                                  sapply(str_split(ID, pattern = "_"),"[[",2),
                                  sapply(str_split(ID, pattern = "_"),"[[",3),
                                  sapply(str_split(ID, pattern = "_"),"[[",4))]

# Sediments = replicates are in col 'location.id',
# Surface waters = replicates are in col 'extra'
colnames(peaks)[1] <- "old.ID"
# Add new ID
peaks[, ID := paste(sample.type, sample.name, sep = "_")]

# Create replicate column
# Change D, M, U to replicates
peaks[sample.type == "SED", replicate.id := factor(location.id, levels = c("D","M","U"),
                                                         labels = c("r1","r2","r3"))]
# Add surface water
peaks[sample.type == "SW", replicate.id := factor(extra, levels = c("1","2","3"),
                                                        labels = c("r1","r2","r3"))]

# sanity checks
peaks[is.na(replicate.id),] # any rows without a replicate ID? -> no
levels(factor(peaks$replicate.id)) # ok, three replicates

# remove unneeded columns
peaks[,c("old.ID","sample.type", "sample.name","location.id", "extra") := 
        list(NULL, NULL, NULL, NULL, NULL)]

# melt community matrix
melt.commat <- melt(peaks, id.vars = c("ID","replicate.id"), value.name = "rel.abun", variable.name = "mz")

# Make into presence absence
melt.commat[, PA := ifelse(rel.abun > 0, 1, 0)]

# Cast replicates into columns
rep.dt <- dcast(melt.commat, ID + mz ~ replicate.id, value.var = "PA")

# Calculate how many replicates observed a MF
rep.dt[, rep.sum := rowSums(.SD, na.rm = T), .SDcols = c("r1","r2","r3")]
rep.dt[replicates, n.rep := i.n.replicate, on = .(ID)]

# Keep two data frames
# 1. Keep observation if MF was found at least in one replicate
# 2. Keep observation if MF was found at least in two replicates
# BUT if there is only one replicate, we will keep the observation

one.rep <- rep.dt[n.rep > 1, PA := ifelse(rep.sum >= 1, 1, 0)] # at least in one replicate
one.rep <- one.rep[n.rep == 1, PA := ifelse(rep.sum == 1, 1, 0)]
# Cast into community matrix
one.commat <- dcast(one.rep, ID ~ mz, value.var = "PA")

two.rep <- rep.dt[, PA := NULL]
two.rep <- rep.dt[n.rep > 1, PA := ifelse(rep.sum >= 2, 1, 0)] # at least in two replicates
two.rep <- rep.dt[n.rep == 1, PA := ifelse(rep.sum == 1, 1, 0)]
# Cast into community matrix
two.commat <- dcast(two.rep, ID ~ mz, value.var = "PA")

# some sanity check
any(rowSums(one.commat[,-1]) != rowSums(two.commat[,-1])) # TRUE = ok

# Save
write.table(one.commat,
            paste0("./1_data.cleaning/output/peaks/FTICR_peaks_commat_rep.merged1_", Sys.Date(), ".csv"), sep = ",",
            row.names = F)
write.table(two.commat,
            paste0("./1_data.cleaning/output/peaks/FTICR_peaks_commat_rep.merged2_", Sys.Date(), ".csv"), sep = ",",
            row.names = F)
