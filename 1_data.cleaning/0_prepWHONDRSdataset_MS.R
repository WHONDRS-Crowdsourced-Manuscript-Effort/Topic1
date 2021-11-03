# Read in processed data of FT-ICR-MS
# Data were processed by Kayla with the R package (https://github.com/EMSL-Computing/ftmsRanalysis)
# and are the processed versions of the raw FT-ICR-MS data available in the EMS drive.
# Non-isotopic metabolites with a mass between 200-900 m/z were included in this processed dataset.

# Set-up -----------------------------------------------------------------------------------

### Packages -------------------------------------------------------------------------------
pckgs <- list("tidyverse","data.table")

### Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
# install missing packages with following line. Please uncomment:
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)

### Load packages
invisible(lapply(pckgs, library, character.only = T))
rm(pckgs, miss.pckgs)

# Read in data ----------------------------------------------------------------------------
proc.ft <- read.csv("./1_data.cleaning/raw.data/Processed_S19S_Sediments_Water_2-2_newcode.csv",
                    sep = ",", stringsAsFactors = F) %>% setDT()
head(proc.ft)[1:5,1:20]
dim(proc.ft)
# 41822 formulae as rows
# samples start at column 39, until that we have variables related to molecular formulae categorization
ncol(proc.ft)- 39
# 504 samples

# Separate into community matrix of abundances
# check if molecular formulae are unique
nrow(proc.ft) == length(unique(proc.ft$MolForm)) # no
# There are duplicated molecular formulae
proc.ft$MolForm[duplicated(proc.ft$MolForm)] # quite a lot
# seems like duplicates are usually almost identical, minimal mz difference
# how many duplicates do we have per molecular formula?
proc.ft[, n.mf := .N, by = .(MolForm)]
hist(proc.ft$n.mf) # ~ 10000 have 2, very few have 3

# 1. merge duplicated molecular formulae for presence-absence matrix
# 2. take always the first entry for the remaining variables
var.dt <- proc.ft[!duplicated(proc.ft$MolForm),] %>%
  dplyr::select(!S19S_0004_ICR.2_p05:S19S_0004_ICR.1_p05)

# Data is in presence absence?
t<-proc.ft[, .SD , .SDcols = S19S_0004_ICR.2_p05:S19S_0004_ICR.1_p05]
any(t[,-1] > 1) # FALSE, data is in presence absence
rm(t)
# Merge duplicates, take the unique peaks

# Make function, could not fine one in base R
# Create function that can merge duplicate entries and keep binary (presence absence) result
as.binary <- function(x){
  x <- ifelse(sum(x) >= 1, 1, 0)
  return(x)
}

abund.dt <- proc.ft[, lapply(.SD, as.binary), by = .(MolForm), .SDcols = S19S_0004_ICR.2_p05:S19S_0004_ICR.1_p05]

# merge back
proc.ft <- merge(var.dt, abund.dt, by = "MolForm")

# do a long fromat of var.dt to decipher sample names
melt.abund <- melt(abund.dt, id.vars = c("MolForm"), variable.name = "Sample", value.name = "peak.int")
# separate "Sample" column
# naming strategy different between surface water and sediment -> separate
sed <- melt.abund[str_detect(Sample,"Field"),]
sur <- melt.abund[!(Sample %in% sed$Sample),]

# First separate surface water
sur <- sur[, c("project.name", "sample.name", "location.id", "unknown") := 
      tstrsplit(Sample, "_", fixed = T)]

levels(factor(sur$location.id)) # Replicate ID
levels(factor(sur$unknown)) # p05 all the same

# Extract only replicate from "location.id"
# and overwrite original column with "U" for upstream so that it matches the sediment ID
# according to the methods, surface water was only sampled upstream
sur[, replicate.id := str_extract(location.id, "[:digit:]$")]
sur[, location.id := "U"] # what does ICR stand for? Necessary? -> remove
sur[, sample.type := "SW"] # add another sample type identifier

# create new IDs
# For specific rivers
sur[, river.id := paste(project.name, sample.name, sep = ".")]
# true unique (cleaned) identifier
sur[, ID := paste(sample.type, river.id, location.id, replicate.id, sep = "_")]
# check if all samples are uniquely identofied
length(unique(sur$Sample)) == length(unique(sur$ID)) # yes

# Next, sediment:
sed <- sed[, c("project.name", "sample.name", "sample.type", "sampling.type", "location.id", "unknown") := 
             tstrsplit(Sample, "_", fixed = T)]

levels(factor(sed$location.id)) # location in river, U, M, D = Upstream, Midstream, Downstream
levels(factor(sed$unknown)) # p1, p15, p2?
levels(factor(sed$sample.type)) # Only "Sed"
levels(factor(sed$sampling.type)) # Only "field", probably other samples are "incubation" not included here

# Extract only letter for location ID
sed[, location.id := str_extract(location.id, "[:alpha:]$")]
sed[, sample.type := "SED"] # overwrite to all capitals

# create new IDs
# For specific rivers
sed[, river.id := paste(project.name, sample.name, sep = ".")]
# true unique (cleaned) identifier
sed[, ID := paste(sample.type, river.id, location.id, sampling.type, sep = "_")]
# check if all samples are uniquely identified
length(unique(sed$Sample)) == length(unique(sed$ID)) # yes

# clean columns for merging back
sed <- sed %>%
  mutate(replicate.id = NA) %>%
  select(ID, river.id, sample.type, location.id, sampling.type, replicate.id,
         original.id = Sample, MolForm, peak.int)
sur <- sur %>%
  mutate(sampling.type = NA) %>%
  select(ID, river.id, sample.type, location.id, sampling.type, replicate.id,
         original.id = Sample, MolForm, peak.int)

# Merge
cleaned.dt <- bind_rows(sed,sur)
# last check
nrow(melt.abund) == nrow(cleaned.dt) # yes
# remove obsolete objects
rm(melt.abund, sed, sur)

# save as long format comes in handy sometimes
#write.table(cleaned.dt, paste0("./1_data.cleaning/output/SED.SW_FTICR_MFpeakint_long_", Sys.Date(),".csv"),
#            sep = ",", row.names = F)


# Meta data -----------------------------------------------------------------------------------------------------
#surface <- read.csv("./1_data.cleaning/raw.data/Surface/WHONDRS_S19S_Metadata_v2.csv",
#                    sep = ",", stringsAsFactors = F) %>% setDT()

#sed <- read.csv("./1_data.cleaning/raw.data/Sediment/WHONDRS_S19S_Metadata_v3.csv",
#                    sep = ",", stringsAsFactors = F) %>% setDT()

# are they the same?
#colnames(surface)[!(colnames(surface) %in% colnames(sed))]
#colnames(sed)[!(colnames(sed) %in% colnames(surface))] # Stream_Order
#surface$Sample_ID %in% sed$Sample_ID 
#sed$Sample_ID %in% surface$Sample_ID # all true

# take only v3

meta <- read.csv("./1_data.cleaning/raw.data/Sediment/WHONDRS_S19S_Metadata_v3.csv",
                 sep = ",", stringsAsFactors = F) %>% setDT()
# 66 cols

# clean
# remove first row, replace "Not_Provided" with NA to allow correct column types
col.description <- meta[1,]
meta <- meta[-1,]
allcols <- colnames(meta)
meta[, (allcols) := lapply(.SD, function(x) ifelse(x == "Not_Provided", NA, x)), .SD = allcols]

# rename Sample_ID to river.id and replace _ with .
meta[, river.id := str_replace(Sample_ID, "_", ".")]

# create same structure as FT-ICR-MS data
# US, MS, DS are corresponding with location.id
sed.meta <- melt(meta, id.vars = "river.id",
           measure.vars = patterns("^US_|^MS_|^DS_"))
sed.meta[, location.id := str_extract(variable, "^[:alpha:]")]
sed.meta[, variable := sub(".*(^US_|^MS_|^DS_)", "", variable)]
sed.meta <- dcast(sed.meta, formula = river.id + location.id ~ variable, value.var = "value")
sed.meta[,c("sampling.type", "sample.type") := list("Field","SED")]
# re-order
sed.meta <- sed.meta %>%
  select(river.id, sample.type, location.id, sampling.type, Algal.Mat.Coverage:Water.Column.Height_cm)

# merge back
meta <- meta[, Sample_ID := NULL][, .SD, .SDcols = !(colnames(meta)[str_detect(colnames(meta), "^US_|^MS_|^DS_")])]

# separate "universal" meta, applicable for both Sediment and surface water
uni.meta <- meta %>% 
  select(Study_Code, river.id, Date:Distance_MS.and.US_meters,
         Approx.Distance.From.Gauge_meters:Notes)

# separate meta data for surface water
sw.meta <- meta %>%
  select(river.id, pH = SW_pH, DO_perc.sat, DO_mg.per.L, Temp_degC = SW_Temp_degC) %>%
  mutate(sample.type = "SW", location.id = "U") %>%
  select(river.id, sample.type, location.id, pH, Temp_degC, DO_perc.sat:DO_mg.per.L)

# create a main meta file matching to the FT data
ft.meta <- cleaned.dt %>%
  select(ID:original.id) %>%
  unique()

# add surface water data
ft.meta <- merge(ft.meta, sw.meta,
                 by = c("river.id","sample.type", "location.id"), all.x = T)
# add sediment data
ft.meta <- merge(ft.meta, sed.meta,
                 by = c("river.id","sample.type","location.id","sampling.type"), all.x = T)

# add a selection of universal meta data
ft.meta <- merge(ft.meta, uni.meta,
                 by = "river.id", all.x = T)

# re-order
ft.meta<-ft.meta %>%
  select(river.id:original.id,
         Date:Time_zone,
         Stream_Name:Distance_MS.and.US_meters,
         Primary.Sources.Flow.Variation, pH:Water.Column.Height_cm,
         Sampler_Organization:Sampler_Name,
         Approx.Distance.From.Gauge_meters:Gauge_Longitude_dec.deg, Met.Station.Nearby:Notes, Study_Code)

# add respiration ---------------------------------------------------------------------------------------
resp <- read.csv("./1_data.cleaning/raw.data/Sediment/WHONDRS_S19S_Sediment_Incubations_Respiration_Rates.csv",
                 sep = ",", stringsAsFactors = F) %>% setDT()
resp[,c("project.name", "sample.name", "sample.type", "location.id") := 
       tstrsplit(Sample_ID, "_", fixed = T)]
resp[, c("river.id", "location.id") := list(paste(project.name, sample.name, sep = "."),
                                            str_extract(location.id, "[:alpha:]$"))]
# re-order and take a selection of columns
resp <- resp %>% select(river.id, location.id, sample.type, resprate_mg.L.h = rate_mg_per_L_per_h,
                resp.inc.time.min = total_incubation_time_min, resp.r2 = R_squared, resp.pval = p_value)

# merge with rest
ft.meta <- merge(ft.meta, resp, by = c("river.id", "location.id", "sample.type"), all.x = T)


# add NPOC Sediment ---------------------------------------------------------------------------------------
# Includes incubation data
sed.npoc <- read.csv("./1_data.cleaning/raw.data/Sediment/WHONDRS_S19S_Sediment_NPOC.csv",
                 sep = ",", stringsAsFactors = F) %>% setDT()

sed.npoc[, c("project.name", "sample.name", "sample.type", "sampling.type", "location.id") := 
             tstrsplit(Sample_ID, "_", fixed = T)]
sed.npoc[, c("river.id", "location.id", "sample.type") := list(paste(project.name, sample.name, sep = "."),
                                            str_extract(location.id, "[:alpha:]$"),
                                            "SED")]
levels(factor(sed.npoc$sampling.type))
sed.npoc[sampling.type == "Inc" | sampling.type == "INC", sampling.type := "Inc"]

# re-order and take a selection of columns
sed.npoc <- sed.npoc %>%
  select(river.id, location.id, sample.type, sampling.type, NPOC_mg.L.asC = X00681_NPOC_mg_per_L_as_C) %>%
  filter(sampling.type == "Field") %>% setDT()
sed.npoc[str_detect(NPOC_mg.L.asC, "Above_Range"), NPOC_mg.L.asC := 22]
sed.npoc[, NPOC_mg.L.asC := as.numeric(NPOC_mg.L.asC)]

# There are some duplicates most are easily resolvable, as they are exactly the same
# However S19S.0067_U_SED_Field has two values 29.74 and 46.59 ?
# Calculate mean...
sed.npoc <- sed.npoc[, .(NPOC_mg.L.asC = mean(NPOC_mg.L.asC)),
                     by = .(river.id, location.id, sample.type, sampling.type)]

# merge with rest
ft.meta <- merge(ft.meta, sed.npoc, by = c("river.id", "location.id", "sample.type","sampling.type"), all.x = T)

# add NPOC surface ---------------------------------------------------------------------------------------
sur.npoc <- read.csv("./1_data.cleaning/raw.data/Surface/WHONDRS_S19S_SW_NPOC.csv",
                     sep = ",", stringsAsFactors = F) %>% setDT()

sur.npoc <- sur.npoc[, c("project.name", "sample.name", "replicate.id") := 
             tstrsplit(Sample_ID, "_", fixed = T)]
sur.npoc[, c("river.id", "replicate.id", "sample.type") := list(paste(project.name, sample.name, sep = "."),
                                                               str_extract(replicate.id, "[:digit:]$"),
                                                               "SW")]

# re-order and take a selection of columns
sur.npoc <- sur.npoc %>%
  select(river.id, replicate.id, sample.type, NPOC_mg.L.asC = X000681_NPOC_mg_per_L_as_C)
sur.npoc[str_detect(NPOC_mg.L.asC, "Below_Range"), NPOC_mg.L.asC := 0.45]
sur.npoc[, NPOC_mg.L.asC := as.numeric(NPOC_mg.L.asC)]

# merge with rest
ft.meta <- ft.meta[sur.npoc, c("NPOC_mg.L.asC") := list(i.NPOC_mg.L.asC), on = .(river.id, replicate.id, sample.type)]

# add isotopes surface ---------------------------------------------------------------------------------------
iso <- read.csv("./1_data.cleaning/raw.data/Surface/WHONDRS_S19S_SW_Isotopes.csv",
                sep = ",", stringsAsFactors = F) %>% setDT()

iso <- iso[, c("project.name", "sample.name", "replicate.id") := 
                       tstrsplit(Sample_ID, "_", fixed = T)]
iso[, c("river.id", "replicate.id", "sample.type") := list(paste(project.name, sample.name, sep = "."),
                                                                str_extract(replicate.id, "[:digit:]$"),
                                                                "SW")]

# re-order and take a selection of columns
iso <- iso %>%
  select(river.id, replicate.id, sample.type,
         del2H_permil = del_2H_permil,
         del18O_permil = del_18O_permil)

# merge with rest
ft.meta <- merge(ft.meta, iso, by = c("river.id", "replicate.id","sample.type"), all.x = T)

# add anions surface --------------------------------------------------------------------------------------
ani <- read.csv("./1_data.cleaning/raw.data/Surface/S19S_Merged_Anions.csv",
                sep = ",", stringsAsFactors = F) %>% setDT()

ani <- ani[, c("project.name", "sample.name", "replicate.id") := 
             tstrsplit(Sample_ID, "_", fixed = T)]
ani[, c("river.id", "replicate.id", "sample.type") := list(paste(project.name, sample.name, sep = "."),
                                                           str_extract(replicate.id, "[:digit:]$"),
                                                           "SW")]

# re-order and take a selection of columns
ani <- ani %>%
  select(river.id, replicate.id, sample.type,
         Cl_mgL = X00940_Cl_mg_per_L,
         SO4_mgL = X00945_SO4_mg_per_L_as_SO4,
         NO3_mgL = X71851_NO3_mg_per_L_as_NO3,
         NO2_mgL = X71856_NO2_mg_per_L_as_NO2,
         F_mgL = X00950_F_mg_per_L)

# correct below detection with minimum value
ani[str_detect(Cl_mgL,"Below_Range"), Cl_mgL := 0.16] #Below_Range_Less_Than_0.16
ani[str_detect(Cl_mgL,"Above_Range"), Cl_mgL := 132] #Above_Range_Greater_Than_132
ani[str_detect(NO3_mgL,"Below_Range"), NO3_mgL := 0.12] #Below_Range_Less_Than_0.12
ani[str_detect(NO2_mgL,"Below_Range"), NO2_mgL := 0.04] #Below_Range_Less_Than_0.04
ani[str_detect(F_mgL,"Below_Range"), F_mgL := 0.02] #Below_Range_Less_Than_0.02
#update col type
ani <- ani[, c("Cl_mgL","SO4_mgL",
             "NO3_mgL","NO2_mgL","F_mgL") := lapply(.SD, as.numeric), .SDcols = c("Cl_mgL","SO4_mgL",
                                           "NO3_mgL","NO2_mgL","F_mgL")]

# merge with rest
ft.meta <- merge(ft.meta, ani, by = c("river.id", "replicate.id","sample.type"), all.x = T)

# reorder
ft.meta <- ft.meta %>% select(ID, river.id, sample.type, location.id, sampling.type, replicate.id,
                   everything())

# save
# write.table(ft.meta, paste0("./1_data.cleaning/output/FTICR_meta_all_",Sys.Date(),".csv"),
#             sep = ",", row.names = F)


# Export "community" matrix
# MolForm cols and ID rows
com.mat <- dcast(cleaned.dt, ID ~ MolForm, value.var = "peak.int")
write.table(com.mat, paste0("./1_data.cleaning/output/FTICR_commat_",Sys.Date(),".csv"),
            sep = ",", row.names = F)

# cross table
write.table(var.dt, paste0("./1_data.cleaning/output/FTICR_crosstable_",Sys.Date(),".csv"),
            sep = ",", row.names = F)

