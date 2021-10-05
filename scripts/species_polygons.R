## EOO and AOO polygons

## Set working environment ####
rm(list = ls())
gc()
# system("ps")
# system("pkill -f R")
x <- c("data.table", "sp", "raster", "rgdal", "rgeos",
       "quickPlot", "fastshp",
       "alphahull", "ConR", "rnaturalearthdata", 
       "future", "future.apply", "parallel")
lapply(x, require, character.only = TRUE)
# options(rgl.useNULL=TRUE) ## to suppress warnings when using library(red)
# x <- c("red", "rCAT")
# lapply(x, require, character.only = TRUE)
rm(x)

## File paths and folders
fungi_dir = "/tempdata/research-cifs/uom_data/fungi_fire_data"
spdata_dir = file.path(fungi_dir ,"spdata_rds")
output_dir = file.path(fungi_dir ,"outputs")
polygons_dir = file.path(output_dir,"species_polygons")
if (!dir.exists(polygons_dir)) {dir.create(polygons_dir)}

# ## Remove existing polygons folder
# file.remove(file.path(polygons_dir, dir(path = polygons_dir)))
# unlink(polygons_dir, recursive = TRUE)

working_dir <- paste0("~/gsdms_r_vol", polygons_dir)


## Using ConR package - minimum convex polygon or alpha hulls
## https://cran.r-project.org/web/packages/ConR/index.html
## Notes: 
##  WGS84 required for data
##  data format: latitude, longitude (in decimal degrees), and taxon name
##  outputs for EOO and AOO same as from {red} when method.range = "convex.hull"
##  ** country_map + exclude.area can be used for cropping to prelimiunary analysis area
spfiles <- list.files(spdata_dir, pattern= ".rds$", full.names = TRUE)
message(cat("Total number of species in cleaned ALA data: "),
        length(spfiles))

## Run IUCN.eval in parallel ####
source("~/gsdms_r_vol/tempdata/workdir/nesp_bugs/scripts/conr_iucn_eval.R")
basemap_file <- "~/gsdms_r_vol/tempdata/research-cifs/uom_data/nesp_bugs_data/outputs/masks/auslands_1poly_wgs84.shp"
file.exists(basemap_file)

hull.method <- "alpha.hull" # "convex.hull" 

## Package: future - for catching errors
plan(multiprocess, workers = future::availableCores()-2)
options(future.globals.maxSize = +Inf) ## CAUTION: Set this to a value, e.g. availablecores-1?/RAM-10?
errorlog <- paste0(output_dir, "/errorlog_species_polygons_", gsub("-", "", Sys.Date()), ".txt")
# if(file.exists(errorlog)){unlink(errorlog)}
writeLines(c(""), errorlog)

system.time(
  suppressWarnings(
    future.apply::future_lapply(
      spfiles,
      function(x){
        tmp <- tryCatch(expr = conr_iucn_eval(species_filename = x,
                                              hull.method = hull.method,
                                              exclude.by.map = TRUE,
                                              basemap_path = basemap_file,
                                              working_dir = working_dir,
                                              iucn_outpath = polygons_dir),
                        error = function(e) {
                          cat(
                            paste(as.character(x), "\n"),
                            file = errorlog,
                            append = TRUE)
                        }
        )
      }, future.seed = TRUE)))

## Check files
csvfiles <- list.files(polygons_dir, pattern = ".csv$",
                       full.names = TRUE, all.files = TRUE)
rdsfiles <- list.files(polygons_dir, pattern = ".rds$",
                       full.names = TRUE, all.files = TRUE)
message(cat("Total number of input species files from cleaned ALA data: "),
        length(spfiles))
message(cat("Number of .csv output files created from IUCN.eval(): "),
        length(csvfiles))
message(cat("Number of .rds output files created from IUCN.eval(): "),
        length(rdsfiles))
# ## Check when files were created
# unique(lubridate::date(file.info(csvfiles, extra_cols = TRUE)$mtime))



## Error runs - I ####
## List species not in output files
errorlog <- file.path(output_dir, "errorlog_species_polygons_20211005.txt")
errorfiles <- trimws(readLines(errorlog)[-1])



## Error runs - II ####
output_sp <- basename(tools::file_path_sans_ext(rdsfiles))
input_sp <- basename(tools::file_path_sans_ext(spfiles))
input_sp[!input_sp %in% output_sp]
errorfiles <- spfiles[!input_sp %in% output_sp]

## Checks & reurns ####
message(cat("Number of species showing errors: "),
        length(errorfiles))

length(csvfiles)+length(errorfiles) == length(spfiles)
all(errorfiles %in% spfiles)

message(cat("error species == species not found in rds output files: "),
        all(basename(tools::file_path_sans_ext(errorfiles))
            %in% input_sp[!input_sp %in% output_sp]))

## Run IUCN.eval in parallel ####
if (length(errorfiles) > future::availableCores()) {
  mc.cores = future::availableCores()-2
} else {
  mc.cores = length(errorfiles)
}

set.seed(1, kind = "L'Ecuyer-CMRG")
system.time(invisible(mclapply(errorfiles,
                               conr_iucn_eval,
                               hull.method = hull.method,
                               exclude.by.map = TRUE,
                               basemap_path = basemap_file,
                               working_dir = working_dir,
                               iucn_outpath = polygons_dir,
                               mc.cores = mc.cores)))

csvfiles <- list.files(polygons_dir, pattern = ".csv$",
                       full.names = TRUE, all.files = TRUE)
rdsfiles <- list.files(polygons_dir, pattern = ".rds$",
                       full.names = TRUE, all.files = TRUE)
message(cat("Total number of input species files from cleaned ALA data: "),
        length(spfiles))
message(cat("Number of .csv output files created from IUCN.eval(): "),
        length(csvfiles))
message(cat("Number of .rds output files created from IUCN.eval(): "),
        length(rdsfiles))



  # Error runs - III ####
  # output_sp <- basename(tools::file_path_sans_ext(rdsfiles))
  # input_sp[!input_sp %in% output_sp]
  # errorfiles <- spfiles[!input_sp %in% output_sp]
  # 
  # setwd(working_dir)
  # basemap <- readOGR(basemap_file)
  # 
  # dat <- as.data.table(readRDS(grep("desognaphosa_yabbra_1043", spfiles,
  #                                   value = TRUE)))
  # dim(dat)
  # spname <- unique(dat$spfile)
  # message(cat("Processing species... ",
  #             spname))
  # dat <- dat[ , .(latitude, longitude, scientificName, family, year)]
  # names(dat) <- c("latitude", "longitude", "tax", "family", "coly")
  # 
  # ## Run ConR function
  # out <- IUCN.eval(dat,
  #                  method.range = hull.method,
  #                  alpha = 2,
  #                  Cell_size_AOO = 2,
  #                  Cell_size_locations = 2,
  #                  country_map = basemap,
  #                  exclude.area = TRUE,
  #                  write_file_option = "csv",
  #                  file_name = spname,
  #                  export_shp = TRUE, ## to get SpatialPolygonsDataFrame in output
  #                  write_results = TRUE,
  #                  write_shp = FALSE, ## to write shapefile files to folder
  #                  SubPop = FALSE,
  #                  DrawMap = FALSE)
  # 
  # saveRDS(out, file = paste0(polygons_dir, "/", spname, ".rds"))
  # rm(dat, spname, out)
  # 
  # ## >> Visualise species with unresolved errors ####
  # rdsfiles <- list.files(polygons_dir, pattern = ".rds$",
  #                        full.names = TRUE, all.files = TRUE)
  # output_sp <- basename(tools::file_path_sans_ext(rdsfiles))
  # input_sp[!input_sp %in% output_sp]
  # errorfiles <- spfiles[!input_sp %in% output_sp]
  # 
  # dat <- do.call("rbind", lapply(errorfiles, readRDS))
  # wgs_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  # dat.spdf <- sp::SpatialPointsDataFrame(coords = dat[, .(longitude, latitude)],
  #                                        data = dat, proj4string = CRS(wgs_crs))
  # 
  # plot(basemap, col = "wheat")
  # plot(dat.spdf, add = TRUE, pch = 18, cex = 1, col = "tomato3")
  # 
  # par(mfrow= c(2,3), mar=c(0,0,1,0))
  # for (e in errorfiles){
  #   dat <- as.data.table(readRDS(e))
  #   sp <- SpatialPoints(dat[, .(longitude, latitude)],
  #                       proj4string = CRS(wgs_crs))
  #   plot(basemap, col = "wheat", main = basename(tools::file_path_sans_ext(e)))
  #   plot(sp, add = TRUE, pch = 8, cex = 3, col = "tomato3")
  #   
  # }
  # rm(dat, dat.spdf, sp, wgs_crs)


## Save IUCN.eval() outputs ####
## >> AOO & EOO areas ####
## >> Save combined .csv output files from IUCN.eval() as data.table
csvfiles <- list.files(polygons_dir, pattern = ".csv$",
                       full.names = TRUE, all.files = TRUE)
out <- do.call("rbind", lapply(csvfiles , read.csv)); dim(out)

out$spfile <- basename(tools::file_path_sans_ext(csvfiles))
names(out)[1] <- "scientificName"

setDT(out, key = "spfile")
out <- out[,Nbe_subPop := NULL]

  # ## Add rows for species with errors
  # rdsfiles <- list.files(polygons_dir, pattern = ".rds$",
  #                        full.names = TRUE, all.files = TRUE)
  # output_sp <- basename(tools::file_path_sans_ext(rdsfiles))
  # input_sp <- basename(tools::file_path_sans_ext(spfiles))
  # errorsp <- input_sp[!input_sp %in% output_sp]
  # errorfiles <- spfiles[!input_sp %in% output_sp]
  # 
  # dat <- do.call("rbind", lapply(errorfiles, readRDS))
  # 
  # temp <- data.frame(matrix(NA, length(errorfiles), dim(out)[2]))
  # names(temp) <- names(out)
  # temp$scientificName <- dat[, .N, scientificName]$scientificName
  # temp$Nbe_unique_occ. <- dat[, .N, scientificName]$N
  # temp$spfile <- basename(tools::file_path_sans_ext(errorfiles))
  # setDT(temp)
  # 
  # out <- rbind(out, temp)

## Save outputs
setDT(out, key = "spfile"); dim(out)
write.csv(out, file = file.path(output_dir, "species_EOO_AOO_ahullareas.csv"),
          row.names = FALSE)

message(cat("Number of species with EOOs: "),
        nrow(out[!is.na(EOO)]))
message(cat("Number of species without EOOs (no error species): "),
        nrow(out[is.na(EOO)]))
message(cat("max number of records for species without EOOs: "),
        max(out[is.na(EOO)]$Nbe_unique_occ.))
  # message(cat("max number of records for species without EOOs (excluding species with errors): "),
  #         max(out[!(scientificName %in% temp$scientificName) & is.na(EOO)]$Nbe_unique_occ.))



## >> EOO polygons ####
## >> Save combined .rds output files from IUCN.eval() as a list
rdsfiles <- list.files(polygons_dir, pattern = ".rds$",
                       full.names = TRUE, all.files = TRUE)
polynames <- basename(tools::file_path_sans_ext(rdsfiles))
temp <- lapply(rdsfiles, readRDS)
names(temp) <- polynames
length(temp)
saveRDS(temp, file = file.path(output_dir, "species_ahull_outputs.rds"))


## >> Save list of SPDF from .rds output files
temp2 <- lapply(temp, "[[", 1)
temp2 <- lapply(temp2, "[[", 2)
length(temp2)
saveRDS(temp2, file = file.path(output_dir, "species_ahullspdf.rds"))

## Create list of non-NULL SPDF for species with EOOs
na.eooIDX <- sapply(temp2, length)
na.eooIDX <- which(na.eooIDX == 0)
temp3 <- temp2[-na.eooIDX]
length(temp3)
saveRDS(temp3, file = file.path(output_dir, "species_ahullEOOspdf.rds"))

## Create list of species names without EOOs (no error species)
length(sort(names(na.eooIDX)))
write.csv(sort(names(na.eooIDX)),
          file = file.path(output_dir, "species_ahullnoEOO.csv"),
          row.names = FALSE)
  # ## Create list of species names without EOOs (including error species)
  # length(sort(c(names(na.eooIDX), errorsp)))
  # write.csv(sort(c(names(na.eooIDX), errorsp)),
  #           file = file.path(output_dir, "species_ahullnoEOO.csv"),
  #           row.names = FALSE)

## Checks
length(names(na.eooIDX)) + length(temp3) == length(spfiles)
message(cat("Total number of species: "),
        length(names(na.eooIDX)) + length(temp3))


