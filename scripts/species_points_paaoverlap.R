## PAA region overlap analysis: foreach{}


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
output_dir = file.path(fungi_dir ,"outputs")
spdata_dir = file.path(fungi_dir ,"spdata_rds")

source("/tempdata/workdir/nesp_bugs/scripts/region_overlap.R")

## Load species data rds files
spfiles <- list.files(spdata_dir, pattern= ".rds$", full.names = TRUE)
length(spfiles)



## Load PAA layer ####
## NOTE: Using states clipped by PAA layer; PAA.tif file to be created
region <- raster(file.path("~/gsdms_r_vol/tempdata/research-cifs/uom_data/nesp_bugs_data/outputs/regions/auslands_wgs84_p_paa.tif"))
region[which(!is.na(region[]))] <- 1
region_vals <- region[]
region_classes <- 1


## Specify overlap folder ####
overlap_dir = file.path(output_dir, "paa_overlap")
# ## Remove existing overlap folder
# file.remove(file.path(overlap_dir, dir(path = overlap_dir)))
# unlink(overlap_dir, recursive = TRUE)
if(!dir.exists(overlap_dir)){dir.create(overlap_dir)}

## Function parameters
wgs_crs  <-  "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
eqarea_crs <- "+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=134 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

## Run overlap analysis in parallel: doMC ####
## 'log' only useful when running small number of species
registerDoMC(future::availableCores()-2)
system.time(log <- foreach(species_dat = spfiles, 
                           .combine = rbind,
                           .errorhandling = "pass",
                           .packages = c('sp', 'raster', 'rgdal', 'data.table')) %dopar%{
                             
                             region_overlap(data_rds = species_dat, 
                                            crs_org = wgs_crs, 
                                            crs_new = eqarea_crs, 
                                            region_raster = region,
                                            region_classes = region_classes,
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
## Merge csv files
out <- do.call("rbind", lapply(csvfiles, fread)); dim(out)
message(cat("Check for NAs: "),
        sum(is.na(out)))
names(out) <- c("spfile", "paa_points")

## Save output table
setDT(out, key = "spfile")
write.csv(out, file = file.path(output_dir, "species_points_in_paa.csv"), row.names = FALSE)

## List of species in/out of PAA
message(cat("Number of species with 0 records inside PAA: "),
        nrow(out[paa_points == 0]))
message(cat("Number of species with at least 1 record inside PAA: "),
        nrow(out[paa_points != 0]))

write.csv(out[paa_points == 0]$spfile, 
          file = file.path(output_dir, "PAA_out_species.csv"),
          row.names = FALSE)

write.csv(out[paa_points > 0]$spfile, 
          file = file.path(output_dir, "PAA_in_species.csv"),
          row.names = FALSE)

