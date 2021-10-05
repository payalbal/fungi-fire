## Data preparation

## Set working environment ####

rm(list = ls())
gc()
# system("ps")
# system("pkill -f R")

x <- c("data.table", 
       "sp", "raster", "sf", "rgdal",
       "rje", "stringr", 
       "lubridate", 
       "future", "future.apply")
lapply(x, require, character.only = TRUE)
rm(x)


## File paths and folders ####
fungi_dir = "/tempdata/research-cifs/uom_data/fungi_fire_data"
spdata_dir = file.path(fungi_dir ,"spdata_rds")
if(!dir.exists(spdata_dir)) {
  dir.create(spdata_dir)
}
output_dir = file.path(fungi_dir ,"outputs")
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}

## Load data ####
dat_all <- dat <- fread(file.path(fungi_dir, "AAFcleaned 2021 sep.csv"))

## Plot
basemap_file <- file.path("/tempdata/research-cifs/uom_data/nesp_bugs_data/outputs", "masks", "auslands_1poly_wgs84.shp")
basemap <- readOGR(basemap_file)
wgs_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

sp <- SpatialPoints(dat[, .(LONGITUDE, LATITUDE)],
                    proj4string = CRS(wgs_crs))
plot(basemap, col = "wheat")
plot(sp, add = TRUE, pch = 8, cex = 0.5, col = "tomato3")




## Filter data ####

## >> Basis of record - no filter applied ####
dat[, .N, BASIS.OF.RECORD]


## >> Year - no filter applied ####
str(dat$DATE)
sum(is.na(dat$DATE))

dat$year <- year(dat$DATE) ## NOTE: 'DATE' is in IDate format
dat[, .N, year]
print(setorder(dat[, .N, by = year], year))

# ## >> Find records without year information
# range(dat$DATE, na.rm = TRUE)
# plot(dat[!is.na(year)][, .N, year], xaxp = c(1761, 2021, 10), pch = 20)
# 
# message(cat("Number of records with year = NA: "),
#         sum(is.na(dat$year))); sum(is.na(dat$year))/nrow(dat)
# 
# message(cat("Number of records with year < 1990: "),
#         nrow(dat[year < 1990])); nrow(dat[year < 1990])/nrow(dat)
# 
# message(cat("Number of records with year as NA or < 1990: "),
#         nrow(dat[is.na(year) | year < 1990])); nrow(dat[is.na(year) | year < 1990])/nrow(dat)
# 
# message(cat("Number of records with year >= 1990: "),
#         nrow(dat[year >= 1990])); nrow(dat[year >= 1990])/nrow(dat)
# 
# length(unique(dat$SPECIES))
# 
# ## >> Number of records lost with year-filter
# t1 <- dat[, .N, SPECIES]
# t2 <- dat[!is.na(year) & year >= 1990][, .N, SPECIES]
# t <- merge(t1, t2, by = "SPECIES", all.x = TRUE)
# dim(t)
# names(t)[2:3] <- c("n.all", "n.sub")
# t[which(is.na(n.sub))]$n.sub = 0
# sum(is.na(t$n.sub))
# setorder(t, n.sub)
# # t$p.lost <- (t$n.all - t$n.sub)/t$n.all
# 
# message(cat("Total number of species in data:"),
#         nrow(t))
# message(cat("Number of species with < 3 records prior to applying year filter: "),
#         nrow(t[n.all < 3])); nrow(t[n.all < 3])/nrow(t)
# message(cat("Number of species with >= 3 records prior to applying year filter: "),
#         nrow(t[n.all >= 3])); nrow(t[n.all >= 3])/nrow(t)
# 
# message(cat("Number of species with < 3 records after filter for 1990 and NAs: "),
#         nrow(t[n.sub < 3])); nrow(t[n.sub < 3])/nrow(t)
# message(cat("Number of species with >= 3 records after filter for 1990 and NAs: "),
#         nrow(t[n.sub >= 3])); nrow(t[n.sub >= 3])/nrow(t)
# message(cat("Number of species with < 3 records after filter that had >=3 records before filter: "),
#         nrow(t[n.sub < 3 & n.all >=3]))
# 
# ## >> Remove records based on year-filter rule: 
# ##  if >= 3 records after filter, remove NA and <1990
# ##  if < 3 records after filter, keep NA and <1990
# sp_applyfilter <- t[n.sub >= 3]$SPECIES
# length(sp_applyfilter)
# 
# dim(dat)
# dat0 <- dat[SPECIES %in% sp_applyfilter][!is.na(year) & year >= 1990]
# dat1 <- dat[!(SPECIES %in% sp_applyfilter)]
# 
# ## Checks
# sum(dat0[, .N, SPECIES]$N < 3)
# dim(dat0)[1] + dim(dat1)[1]
# length(unique(dat0$SPECIES)) + 
#   length(unique(dat1$SPECIES)) == length(unique(dat$SPECIES))
# 
# n <- dim(dat)[1]
# dat <- rbind(dat0, dat1)
# points(dat[!is.na(year)][, .N, year], xaxp = c(1761, 2021, 10), pch = 20, col = "tomato3")
# 
# message(cat("Number of records lost with year-filer: "),
#         n-dim(dat)[1])
# message(cat("Proportion of data lost with year-filer: "),
#         (n-dim(dat)[1])/n)
# message(cat("Number of unique species in data: "),
#         length(unique(dat$SPECIES)))



## >> Geographic accuracy ####
sum(!is.na(dat$POSITIONAL.ACCURACY))
range(dat$POSITIONAL.ACCURACY, na.rm = TRUE)
tail(sort(unique(dat$POSITIONAL.ACCURACY)))

hist(dat$POSITIONAL.ACCURACY, prob = TRUE)
plot(density(dat$POSITIONAL.ACCURACY, na.rm = TRUE, bw = 10))

## Filter if recorded accuracy is > 5k
message(cat("Proportion of records with POSITIONAL.ACCURACY recorded as > 5k: "),
        nrow(dat[POSITIONAL.ACCURACY > 5000])/nrow(dat)); nrow(dat[POSITIONAL.ACCURACY > 5000])

dat <- dat[POSITIONAL.ACCURACY <= 5000 | is.na(POSITIONAL.ACCURACY)]; dim(dat)

## Check
nrow(dat_all) - nrow(dat_all[POSITIONAL.ACCURACY > 5000]) == nrow(dat)

message(cat("Number of species lost: "),
        length(unique(dat_all$SPECIES)) - length(unique(dat$SPECIES)))

## >> Outlier status - no filter applied ####
unique(dat$OUTLIER.STATUS)
dat[, .N, OUTLIER.STATUS]
message(cat("Number of records labelled as outliers: "),
        sum(dat[, .N, OUTLIER.STATUS][OUTLIER.STATUS != "NONE"]$N))





## Subset data ####
dim(dat_all)
length(unique(dat_all$ID.IN.THIS.DATABASE))
length(unique(dat_all$ID.IN.ORIGINAL.DATABASE))
## Note: Not the same. Use ID.IN.THIS.DATABASE.

names(dat)
dat <- dat[,.(ID.IN.THIS.DATABASE, LATITUDE, LONGITUDE, SPECIES, order, class, family, year, ORIGINAL.DATABASE, `ALBERS X`, `ALBERS Y`)]
names(dat) <- c("id", "latitude", "longitude", "scientificName", 
                "order", "class", "family", "year", "data_source",
                "X", "Y")



## Create new column to give unique ID by scientificName ####
## Create spfile column
y <- dat$scientificName
grep('\\d', y, value = TRUE) ## looks for digits

y <- tolower(str_replace_all(y, " ", "_"))
grep('\\s', y, value = TRUE) ## looks for spaces

dat[, spfile := y]
names(dat)

## Create unique ID by "scientificName", "class", "family"
setDT(dat)[, new_id := .GRP, by = c("scientificName", "class", "family")]; names(dat)
length(unique(dat$new_id))
range(unique(dat$new_id))
length(unique(dat$scientificName))

## Merge spfile and new_ID (overwrite existing spfile column)
dat$spfile <- paste0(dat$spfile, "_", dat$new_id)
message(cat("Number of unique scientificName in new dataset: "),
        length(unique(dat$scientificName)))
message(cat("Number of unique spfile in new dataset: "),
        length(unique(dat$spfile)))
dat[,new_id := NULL]



## Save as rds file
## Save rds files for species ####
save_spdata2 <- function(species_uid, data, data_dir){
  
  temp <- dat[spfile == species_uid]
  spfile <- unique(temp$spfile)
  
  if (length(spfile) > 1){
    stop("Error: More than 1 unique spfile for naming species file...")
  }
  
  saveRDS(as.data.table(temp),
          file = file.path(data_dir, paste0(spfile, ".rds")))
}


## Run function in parallel
errorlog <- paste0(output_dir, "/errorlog_rdsdata_", gsub("-", "", Sys.Date()), ".txt")
writeLines(c(""), errorlog)

# data <- fread(file.path(output_dir, "data_ALAnonALA_wgs84_corrected.csv"))
dat <- dat
all_species <- unique(dat$spfile)

plan(multiprocess, workers = future::availableCores()-2)
options(future.globals.maxSize = +Inf)

system.time(invisible(future.apply::future_lapply(all_species,
                            function(x){
                              tmp <- tryCatch(expr = save_spdata2(species_uid = x, 
                                                                  data = dat, 
                                                                  data_dir = spdata_dir),
                                              error = function(e){ 
                                                print(paste("\nError: More than 1 unique spfile for naming species file for...", x))
                                                cat(paste(x, "\n"),
                                                    file = errorlog, 
                                                    append = TRUE)
                                              })
                            })))

## Check files
length(list.files(spdata_dir, pattern = ".rds$"))
length(all_species)