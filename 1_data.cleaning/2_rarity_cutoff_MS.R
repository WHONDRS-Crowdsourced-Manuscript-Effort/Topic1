# Apply rarity cut-offs

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

meta <- read.csv("./1_data.cleaning/output/FTICR_meta_all_2022-01-18.csv",
                 sep = ",", stringsAsFactors = F) %>% setDT()
head(meta)

# cross table
cross <- read.csv("./1_data.cleaning/output/FTICR_crosstable_2021-09-29.csv",
            sep = ",", stringsAsFactors =  F) %>% setDT()
head(cross)

# Loading four data frames
# Molecular formulae vs peaks data
# Two replicate merging thresholds (present in at least one replicate, present in at least two replicates)
mf.rep1 <- read.csv("./1_data.cleaning/output/FTICR_commat_rep.merged1_2022-01-18.csv",
                    sep = ",", stringsAsFactors = F) %>% setDT()
mf.rep2 <- read.csv("./1_data.cleaning/output/FTICR_commat_rep.merged2_2022-01-18.csv",
                    sep = ",", stringsAsFactors = F) %>% setDT()

mz.rep1 <- read.csv("./1_data.cleaning/output/peaks/FTICR_peaks_commat_rep.merged1_2022-01-18.csv",
                    sep = ",", stringsAsFactors = F) %>% setDT()
mz.rep2 <- read.csv("./1_data.cleaning/output/peaks/FTICR_peaks_commat_rep.merged2_2022-01-18.csv",
                    sep = ",", stringsAsFactors = F) %>% setDT()


# Remove rare peaks ---------------------------------------------------------------------
# Check if there are several peaks within the same molecular formula

mult.peaks <- cross[, .(n.peaks = .N), by = .(MolForm)]
mult.peaks[n.peaks > 1,]
# no, good we can merge peaks by molecular formula
rm(mult.peaks)

df.list <- list(mf.rep1, mf.rep2, mz.rep1, mz.rep2)

rar.ls <- llply(df.list, function(x){
  # melt community matrix for easy classification
  melt.commat <- melt(x, id.vars = "ID", value.name = "PA", variable.name = "MolForm")
  # enumerate frequency of each molecular formulae
  melt.commat[, freq := nrow(.SD[PA > 0,]), by = .(MolForm)]
  
  # How many MF are kept if we only include MF that appear in more than 2 samples
  stats.df.rar2 <- data.frame(no.mf.kept = length(unique(melt.commat[freq > 2,]$MolForm)),
                         no.mf = length(unique(melt.commat$MolForm))) %>%
    mutate(perc.kept = no.mf.kept * 100 / no.mf,
           rarity.cutoff = "More than 2 samples (rar2)")
  
  stats.df.rar1 <- data.frame(no.mf.kept = length(unique(melt.commat[freq >= 2,]$MolForm)),
             no.mf = length(unique(melt.commat$MolForm))) %>%
    mutate(perc.kept = no.mf.kept * 100 / no.mf,
           rarity.cutoff = "More than 1 sample (rar1)")
  
  # One dataset with this 2 sample rarity cutoff
  rar2 <- melt.commat[freq > 2,]
  rar2 <- dcast(rar2, ID ~ MolForm, value.var = "PA") %>% setDF()
  
  # One dataset with a 1 sample rarity cutoff
  rar1 <- melt.commat[freq >= 2,]
  rar1 <- dcast(rar1, ID ~ MolForm, value.var = "PA") %>% setDF()
  
  rar.ls <- list(rar1 = list(rar1, stats.df.rar1),
                 rar2 = list(rar2, stats.df.rar2))
  
  return(rar.ls)
})

# Add list names to merge in to data frame
names(rar.ls) <- c("Molecular formulae (mf)-MF in at least 1 replicate",
                   "Molecular formulae (mf)-MF in at least 2 replicates",
                   "Peaks (mz)-MZ in at least 1 replicate",
                   "Peaks (mz)-MZ in at least 2 replicates")
# Merge stats into table
stats.tb <- ldply(rar.ls, function(x){
  out <- bind_rows(lapply(x, "[[", 2))
  return(out)
})

# Clean data frame to export
stats.tb <- stats.tb %>% select(.id, rarity.cutoff, no.mf, everything()) %>%
  separate(.id, into = c("Dataset","Replicate merging"), sep = "-")
colnames(stats.tb)[3:6] <- c("Used rarity cutoff",
                             "Original number of mf/mz",
                             "Number of mf/mz kept",
                             "Percentage retained")

# Upload to Github
knitr::kable(stats.tb)

# Save
write.table(stats.tb, "./1_data.cleaning/output/cleaning_summary.csv", sep = ",", row.names = F)

# Extract community matrix in each list
commat.ls <- llply(rar.ls, function(x){
  out <- lapply(x, "[[",1)
  return(out)
})

# Flatten list into a one-level list
flattened.ls <- purrr::flatten(commat.ls)

# Some sanity check
sapply(flattened.ls, dim)

# Saving paths
# create list with names
save.vec <- c(paste0("./1_data.cleaning/output/FTICR_commat_rep.merged1_rar1_",Sys.Date(),".csv"),
              paste0("./1_data.cleaning/output/FTICR_commat_rep.merged1_rar2_",Sys.Date(),".csv"),
              paste0("./1_data.cleaning/output/FTICR_commat_rep.merged2_rar1_",Sys.Date(),".csv"),
              paste0("./1_data.cleaning/output/FTICR_commat_rep.merged2_rar2_",Sys.Date(),".csv"),
              paste0("./1_data.cleaning/output/peaks/FTICR_peaks_commat_rep.merged1_rar1_",Sys.Date(),".csv"),
              paste0("./1_data.cleaning/output/peaks/FTICR_peaks_commat_rep.merged1_rar2_",Sys.Date(),".csv"),
              paste0("./1_data.cleaning/output/peaks/FTICR_peaks_commat_rep.merged2_rar1_",Sys.Date(),".csv"),
              paste0("./1_data.cleaning/output/peaks/FTICR_peaks_commat_rep.merged2_rar2_",Sys.Date(),".csv"))

# Save community matrix
mapply(write.table, flattened.ls, save.vec, MoreArgs = list(sep = ",", dec = ".", row.names = F))

# So we end up with the following data sets:
# A)) Molecular formulae
  # 1) Replicate merged in at least 1 sample - all MF
  # -- FTICR_commat_rep.merged1_DATE.csv
  # 2) Replicate merged in at least 1 sample - MF that are present in more than 1 site (rar1)
  # -- FTICR_commat_rep.merged1_rar1_DATE.csv
  # 3) Replicate merged in at least 1 sample - MF that are present in more than 2 sites (rar2)
  # -- FTICR_commat_rep.merged1_rar2_DATE.csv
  # 4) Replicate merged in at least 2 samples - all MF
  # -- FTICR_commat_rep.merged2_DATE.csv
  # 5) Replicate merged in at least 2 samples - MF that are present in more than 1 site (rar1)
  # -- FTICR_commat_rep.merged2_rar1_DATE.csv
  # 6) Replicate merged in at least 2 samples - MF that are present in more than 2 sites (rar2)
  # -- FTICR_commat_rep.merged2_rar2_DATE.csv
# B)) Peaks
  # 7) Replicate merged in at least 1 sample - all peaks
  # -- ./peaks/FTICR_peaks_commat_rep.merged1_DATE.csv
  # 8) Replicate merged in at least 1 sample - peaks that are present in more than 1 site (rar1)
  # -- ./peaks/FTICR_peaks_commat_rep.merged1_rar1_DATE.csv
  # 9) Replicate merged in at least 1 sample - peaks that are present in more than 2 sites (rar2)
  # -- ./peaks/FTICR_peaks_commat_rep.merged1_rar2_DATE.csv
  # 10) Replicate merged in at least 2 samples - all peaks
  # -- ./peaks/FTICR_peaks_commat_rep.merged2_DATE.csv
  # 11) Replicate merged in at least 2 samples - peaks that are present in more than 1 site (rar1)
  # -- ./peaks/FTICR_peaks_commat_rep.merged2_rar1_DATE.csv
  # 12) Replicate merged in at least 2 samples - peaks that are present in more than 2 sites (rar2)
  # -- ./peaks/FTICR_peaks_commat_rep.merged2_rar2_DATE.csv


# Put the three dataset into list for easy processing
com.ls <- list(mf.rep1,
               flattened.ls[[1]],  flattened.ls[[2]],
               mf.rep2,
               flattened.ls[[3]], flattened.ls[[4]])

# next, get cross table for each subset
cross.ls <- llply(com.ls, function(x){
  subcross <- cross[MolForm %in% colnames(x),]
  return(subcross)
})

# Save different rarity thresholds' cross tables
# create list with names
save.vec <- c(paste0("./1_data.cleaning/output/FTICR_crosstable_rep.merged1_",Sys.Date(),".csv"),
              paste0("./1_data.cleaning/output/FTICR_crosstable_rep.merged1_rar1_",Sys.Date(),".csv"),
              paste0("./1_data.cleaning/output/FTICR_crosstable_rep.merged1_rar2_",Sys.Date(),".csv"),
              paste0("./1_data.cleaning/output/FTICR_crosstable_rep.merged2_",Sys.Date(),".csv"),
              paste0("./1_data.cleaning/output/FTICR_crosstable_rep.merged2_rar1_",Sys.Date(),".csv"),
              paste0("./1_data.cleaning/output/FTICR_crosstable_rep.merged2_rar2_",Sys.Date(),".csv"))

# Save cross table
mapply(write.table, cross.ls, save.vec, MoreArgs = list(sep = ",", dec = ".", row.names = F))


# Revert to original by river meta data --------------------------------------------------------------------------------
# Keep meta data of replicates and measurements
rep.meta <- meta

# Get original meta data
header <- colnames(read.csv("./1_data.cleaning/raw.data/Surface/WHONDRS_S19S_Metadata_v2.csv",
                            sep = ",", stringsAsFactors = F, header = T))
surface <- read.csv("./1_data.cleaning/raw.data/Surface/WHONDRS_S19S_Metadata_v2.csv",
                    sep = ",", stringsAsFactors = F, header = F, skip = 2) %>% setDT()
colnames(surface) <- header

header <- colnames(read.csv("./1_data.cleaning/raw.data/Sediment/WHONDRS_S19S_Metadata_v3.csv",
                            sep = ",", stringsAsFactors = F, header = T))
sed <- read.csv("./1_data.cleaning/raw.data/Sediment/WHONDRS_S19S_Metadata_v3.csv",
                sep = ",", stringsAsFactors = F, skip = 2, header = F) %>% setDT()
colnames(sed) <- header

# merge
colnames(surface) %in% colnames(sed)
colnames(sed)[!(colnames(sed) %in% colnames(surface))] # Strahler order missing in surface

# clean
# replace "Not_Provided" with NA to allow correct column types

# merge the stream order from sediment meta
surface <- surface[sed,c("Stream_Order") := list(i.Stream_Order), on = .(Sample_ID)]

# rename Sample_ID to river.id and replace _ with .
surface[, river.id := str_replace(Sample_ID, "_", ".")]
surface[, sample.type := "SW"]
surface[, ID := paste(sample.type, river.id, sep = "_")]

sed[, river.id := str_replace(Sample_ID, "_", ".")]
sed[, sample.type := "SED"]
sed[, ID := paste(sample.type, river.id, sep = "_")]

meta <- bind_rows(surface, sed)
rm(surface, sed, header)

allcols <- colnames(meta)
meta[, (allcols) := lapply(.SD, function(x) ifelse(x == "Not_Provided", NA, x)), .SD = allcols]

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

# re-order
meta <- meta %>% dplyr::select(Study_Code, ID, sample.type, river.id, Date:Stream_Order, DO_comment)

# Add other meta data
colnames(rep.meta)[!(colnames(rep.meta) %in% colnames(meta))]
# We added the last 13 columns in script 0
colnames(rep.meta)[(ncol(rep.meta)-13):ncol(rep.meta)]

# select only columns needed and summarise
rep.means <- left_join(rep.meta %>% dplyr::select(ID = merge.id, resprate_mg.L.h:F_mgL) %>%
  group_by(ID) %>%
  summarise_all(.funs = mean, na.rm = T),
  rep.meta %>% dplyr::select(ID = merge.id, no.replicates) %>%
  group_by(ID) %>%
  summarise_all(.funs = unique),
  by = "ID")
# Replace NaN with NA
rep.means <- rep.means %>%
  mutate_at(.vars = vars(-ID), .funs = ~ifelse(is.nan(.), NA, .)) %>% setDT()

# Merge
meta$ID[!(meta$ID %in% rar1$ID)]
# some rivers don't have Sediments -> NA, one doesn't have surface water -> NA
# Keep only rows that have FT data
meta <- meta[rep.means, , on = .(ID)]

# Save meta data ---------------------------------------------------------------------------------------------------------------
write.table(meta,
            paste0("./1_data.cleaning/output/FTICR_meta_eachriver_", Sys.Date(),".csv"), sep = ",", dec = ".", row.names = F)


# Done.