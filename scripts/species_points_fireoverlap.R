## Fire overlap analysis: foreach{}

## Points overlap for all species


## Set working environment ####
rm(list = ls())
gc()
# system("ps")
# system("pkill -f R")
x <- c("data.table", "sp", "raster", "rgdal", 
       "gdalUtils", "rgeos", "doMC", "foreach")
lapply(x, require, character.only = TRUE)
rm(x)

## File paths and folders
fungi_dir = "/tempdata/research-cifs/uom_data/fungi_fire_data"
spdata_dir = file.path(fungi_dir ,"spdata_rds")
output_dir = file.path(fungi_dir ,"outputs")
overlap_dir = file.path(output_dir, "points_overlap")
  ## Remove existing overlap folder
  # file.remove(file.path(overlap_dir, dir(path = overlap_dir)))
  # unlink(overlap_dir, recursive = TRUE)
if(!dir.exists(overlap_dir)){dir.create(overlap_dir)}

## Points overlap function
source("~/gsdms_r_vol/tempdata/workdir/nesp_bugs/scripts/points_overlap.R")

## Get species data
spfiles <- list.files(spdata_dir, pattern= ".rds$", full.names = TRUE)
length(spfiles)

## Load in fire severity raster (re-classed) and get unique classes
fire_severity <- raster(file.path("~/gsdms_r_vol/tempdata/research-cifs/uom_data/nesp_bugs_data/outputs/fire/severity5_eqar250_native_paa.tif"))
fire_classes <- sort(unique(na.omit(fire_severity[])))

## Function parameters
wgs_crs  <-  "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
eqarea_crs <- "+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=134 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"



## Run overlap analysis in parallel: doMC ####
registerDoMC(future::availableCores())
system.time(log <- foreach(species_dat = spfiles, 
                           .combine = rbind,
                           .errorhandling = "pass",
                           .packages = c('sp', 'raster', 'rgdal', 'data.table')) %dopar%{
                             
                             points_overlap(data_rds = species_dat, 
                                            crs_org = wgs_crs, 
                                            crs_new = eqarea_crs, 
                                            fire_severity = fire_severity,
                                            fire_classes = fire_classes,
                                            outdir = overlap_dir)
                           })




## Error checking ####
log
csvfiles <- list.files(overlap_dir, pattern = ".csv$",
                       full.names = TRUE, all.files = TRUE)
message(cat("Number of input species: "),
        length(spfiles))
message(cat("Number of output files: "),
        length(csvfiles))


## Output table ####
## >> Merge csv files ####
out <- do.call("rbind", lapply(csvfiles, fread)); dim(out)
names(out)[1] <- "spfile"
setDT(out, key = "spfile")
write.csv(out, file = file.path(output_dir, "species_points_fireoverlap.csv"), row.names = FALSE)

## >> Checks ####
## >> Total overlapped points should be <= Total Occurrence points
message(cat("Check if #points overlapping with fire is always <= Total # points for species: "),
        all(rowSums(out[,.(Fire_Class_1, Fire_Class_2, Fire_Class_3, Fire_Class_4, Fire_Class_5)], na.rm = TRUE) <= out$Occurrence_Points))

## >> Look for NAs in table - shouldn't be any
message(cat("Check for NAs: "),
        sum(is.na(out)))

## >> Look for Occurrence_Points = 0 in table - shouldn't be any
message(cat("Number of species with Occurrence_Points = 0: "),
        sum(out$Occurrence_Points == 0))


## >> Add percentage overlap columns ####
out$Overlap_Points_Fire345_GEEBAM2_as_unburnt <- 
  ((out$Fire_Class_3 + out$Fire_Class_4 + out$Fire_Class_5)/
     (out$Occurrence_Points - out$Fire_Class_1)) * 100

out$Overlap_Points_Fire2345_GEEBAM2_as_burnt <- 
  ((out$Fire_Class_2 + out$Fire_Class_3 + out$Fire_Class_4 + out$Fire_Class_5)/
     (out$Occurrence_Points - out$Fire_Class_1)) * 100

out$Overlap_Points_Severe_Fire45 <- 
  ((out$Fire_Class_4 + out$Fire_Class_5)/
     (out$Occurrence_Points - out$Fire_Class_1)) * 100

sum(is.na(out$Overlap_Points_Fire345_GEEBAM2_as_unburnt))
sum(is.na(out$Overlap_Points_Fire2345_GEEBAM2_as_burnt))
sum(is.na(out$Overlap_Points_Severe_Fire45))

## Check if NAs are same in all three columns
all(which(is.na(out$Overlap_Points_Fire345_GEEBAM2_as_unburnt)) == 
      which(is.na(out$Overlap_Points_Fire2345_GEEBAM2_as_burnt)))
all(which(is.na(out$Overlap_Points_Fire2345_GEEBAM2_as_burnt)) == 
      which(is.na(out$Overlap_Points_Severe_Fire45)))

## Look at rows with NAs: happens because we get a 0 in the denominator
naidx <- which(is.na(out$Overlap_Points_Severe_Fire45))
out[naidx]
## 0/0 divisions
out[naidx]$spfile


## >> Add class/family information to output table ####
tax <- fread(file = file.path(output_dir, "fungi_data.csv"))
tax <- setDT(tax, key = "spfile")[, .SD[1L] ,.(scientificName, class, order, family, spfile)]
tax <- tax[,.(scientificName, class, order, family, spfile)]
dim(tax)


## Checks
all(out$spfile %in% tax$spfile)
dim(out); length(unique(out$spfile))
dim(tax); length(unique(tax$spfile)); length(unique(tax$scientificName))

out <- merge(out, tax, by = "spfile")
setDT(out, key = "spfile")
write.csv(out, file = file.path(output_dir, "species_points_fireoverlap.csv"), 
          row.names = FALSE)


## >> Add EOO and AOO information from species_EOO_AOO_ahullareas.csv ####
polyareas <- fread(file.path(output_dir, "species_EOO_AOO_ahullareas.csv"))
polyareas <- setDT(polyareas, key = "spfile")[spfile %in% out$spfile][, .(spfile, EOO, AOO, Nbe_unique_occ., scientificName)]

dim(out); length(unique(out$spfile)); length(unique(out$scientificName))
dim(polyareas); length(unique(polyareas$spfile)); length(unique(polyareas$scientificName))

## Check
out2 <- merge(out, polyareas, by = "spfile")
all(out2$scientificName.x == out2$scientificName.y)
which(!(out2$scientificName.x == out2$scientificName.y))
sum(!(out2$scientificName.x == out2$scientificName.y))
out2[which(!(out2$scientificName.x == out2$scientificName.y))][,.(scientificName.x, scientificName.y)]
## out2$scientificName.x & out2$scientificName.y are the same
rm(out2)

polyareas[, scientificName := NULL]
out <- merge(out, polyareas, by = "spfile")
write.csv(out, file = file.path(output_dir, "species_points_fireoverlap.csv"),
          row.names = FALSE)

## Checks
all(out$Nbe_unique_occ <= out$Occurrence_Points)
  

## Final formatting ####
## Column names
names(out)[which(names(out) == "Nbe_unique_occ.")] <- "Nbe_unique_occ"
names(out)[grep("Fire_Class_", names(out))] <- paste0(names(out)[grep("Fire_Class_", names(out))], "_Points")


## Save table
setDT(out, key = "spfile")
write.csv(out, file = file.path(output_dir, "species_points_fireoverlap.csv"),
          row.names = FALSE)


## Summarize outputs ####
message(cat("NA in scientificName: "),
        length(which(is.na(out$scientificName))))

message(cat("Total number of species: "),
        nrow(out))

message(cat("Proportion of species showing any overlap: "),
        nrow(out[Overlap_Points_Fire2345_GEEBAM2_as_burnt != 0])/nrow(out))
nrow(out[Overlap_Points_Fire2345_GEEBAM2_as_burnt != 0])

## Unique class-family summaries in 100% overlap
## Any fire severity
message(cat("Unique species showing 100% overlap with any fire:"))
unique(out[Overlap_Points_Fire2345_GEEBAM2_as_burnt == 100]$scientificName)

message(cat("Unique classes showing 100% overlap with any fire:"))
out[Overlap_Points_Fire2345_GEEBAM2_as_burnt == 100][, .SD[1L] ,.(class)][,.(class)][, .N, class]$class

message(cat("Unique families showing 100% overlap with any fire:"))
out[Overlap_Points_Fire2345_GEEBAM2_as_burnt == 100][, .SD[1L] ,.(family)][,.(family)][, .N, family]

message(cat("Unique classes & families showing 100% overlap with any fire:"))
out[Overlap_Points_Fire2345_GEEBAM2_as_burnt == 100][, .SD[1L] ,.(class, family)][,.(class,family)]


## High fire severity
message(cat("Unique species showing 100% overlap with severe fire:"))
unique(out[Overlap_Points_Severe_Fire45 == 100]$scientificName)

message(cat("Unique classes showing 100% overlap with severe fire:"))
out[Overlap_Points_Severe_Fire45 == 100][, .SD[1L] ,.(class)][,.(class)][, .N, class]$class

message(cat("Unique families showing 100% overlap with any fire:"))
out[Overlap_Points_Severe_Fire45 == 100][, .SD[1L] ,.(family)][,.(family)][, .N, family]

message(cat("Unique classes & families showing 100% overlap with any fire:"))
out[Overlap_Points_Severe_Fire45 == 100][, .SD[1L] ,.(class, family)][,.(class,family)]



