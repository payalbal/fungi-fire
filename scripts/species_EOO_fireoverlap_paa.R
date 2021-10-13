## Fire overlap analysis: foreach{}

## Polygon overlap for species with 3+ records


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

shapefile_dir = file.path(output_dir, "species_shapefiles")
overlap_dir = file.path(output_dir, "polygon_overlap")
# ## Remove exisitng overlap and shapefiles folder
# unlink(shapefile_dir, recursive = TRUE, force = TRUE)
# file.remove(file.path(overlap_dir, dir(path = overlap_dir)))
# unlink(overlap_dir, recursive = TRUE)
if(!dir.exists(shapefile_dir)){dir.create(shapefile_dir)}
if(!dir.exists(overlap_dir)){dir.create(overlap_dir)}

## Load spdf data for species with EOO ####
species_maps <- readRDS(file.path(output_dir, "species_ahullEOOspdf.rds"))
polygon_list <- names(species_maps)
message(cat("Number of species with polygons: "),
        length(polygon_list))


# ## To only run overlaps for species with > 0 points in PAA
# paa_species <- fread(file.path(output_dir, "PAA_in_species.csv"))$x
# polygon_list <- polygon_list[polygon_list %in% paa_species]
# message(cat("Number of species with polygons and > 0 points in PAA: "),
#         length(polygon_list))

# region_overlap()

## Run overlap analysis in parallel: doMC ####
job_script <- file.path("/tempdata/workdir/fungi-fire", 
                        "scripts", "species_EOO_fireoverlap_paa_job.R")
rstudioapi::jobRunScript(job_script, encoding = "unknown", 
                         workingDir = "/tempdata/workdir/fungi-fire",
                         importEnv = FALSE, exportEnv = "")


## Error reruns for individual species runs if required ####
## Check files
csvfiles <- list.files(overlap_dir, pattern = ".csv$",
                       full.names = TRUE, all.files = TRUE)
message(cat("Number of input species: "),
        length(polygon_list))
message(cat("Number of output files: "),
        length(csvfiles))
length(polygon_list) - length(csvfiles)

## Find missing species from outputs
csvnames <- basename(tools::file_path_sans_ext(csvfiles))
error_list <- polygon_list[!polygon_list %in% csvnames]
message(cat("Number of species in error list: "),
        length(error_list))

## Run polygon_paa_overlap.R stepwise if needed
polys = error_list[1]
species_name = polys
species_poly = species_maps[[polys]]
outdir = overlap_dir
...

## Check output is created
fread(grep(polys, csvfiles, value = TRUE))


## Output table ####
## Merge csv files
out <- do.call("rbind", lapply(csvfiles, fread)); dim(out)
names(out)[1] <- "spfile"
setDT(out, key = "spfile")
write.csv(out, file = file.path(output_dir, "species_polygon_fireoverlap.csv"), row.names = FALSE)

## Checks
## >> Total overlapped area should be <= Species Polygon
message(cat("Check if area overlapping with fire is always <= Total polygon area for species: "),
        all(rowSums(out[,.(Fire_Class_1, Fire_Class_2, Fire_Class_3, Fire_Class_4, Fire_Class_5)]) <= out$Species_Polygon))

## >> Look for NAs in table - shoudln't be any
message(cat("Check for NAs: "),
        sum(is.na(out)))

## >> Look for Species_Polygon = 0 in table
message(cat("Number of species with Species_Polygon area = 0: "),
        sum(out$Species_Polygon == 0))

## Add percentage overlap columns
out$Overlap_Polygons_Fire345_GEEBAM2_as_unburnt <- 
  ((out$Fire_Class_3 + out$Fire_Class_4 + out$Fire_Class_5)/
     (out$Species_Polygon - out$Fire_Class_1)) * 100

out$Overlap_Polygons_Fire2345_GEEBAM2_as_burnt <- 
  ((out$Fire_Class_2 + out$Fire_Class_3 + out$Fire_Class_4 + out$Fire_Class_5)/
     (out$Species_Polygon - out$Fire_Class_1)) * 100

out$Overlap_Polygons_Severe_Fire45 <- 
  ((out$Fire_Class_4 + out$Fire_Class_5)/
     (out$Species_Polygon - out$Fire_Class_1)) * 100


sum(is.na(out$Overlap_Polygons_Fire345_GEEBAM2_as_unburnt))
sum(is.na(out$Overlap_Polygons_Fire2345_GEEBAM2_as_burnt))
sum(is.na(out$Overlap_Polygons_Severe_Fire45))


## Final formatting ####
names(out)[grep("Fire_Class_", names(out))] <- paste0(names(out)[grep("Fire_Class_", names(out))], "_Area")
setDT(out, key = "spfile")
write.csv(out, file = file.path(output_dir, "species_polygon_fireoverlap.csv"), 
          row.names = FALSE)



## Summarize outputs ####
message(cat("NA in scientificName: "),
        length(which(is.na(out$scientificName))))

message(cat("Total number of species: "),
        nrow(out))




