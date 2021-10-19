## Combining fire overlap results


## Set working environment ####
rm(list = ls())
gc()
# system("ps")
# system("pkill -f R")
x <- c("data.table", "readxl") 
lapply(x, require, character.only = TRUE)
rm(x)

## File paths and folders
fungi_dir = "/tempdata/research-cifs/uom_data/fungi_fire_data"
output_dir = file.path(fungi_dir ,"outputs")


## Combine output tables ####
poly <- fread(file.path(output_dir, "species_polygon_fireoverlap.csv")); dim(poly)
point <- fread(file.path(output_dir, "species_points_fireoverlap.csv")); dim(point)

## Merge rows for species with polygon overlap innfo (n = 29029)
out1 <- merge(point, poly, by = "spfile")
dim(out1)

## Get rows for species without polygon informatiom (n = 29948)
sum(!(point$spfile %in% out1$spfile))
out2 <- point[!(point$spfile %in% out1$spfile)]
dim(out2)

## Add empty polygon information columns to out2 
dt <- setNames(data.table(matrix(nrow = nrow(out2), ncol = length(names(poly)[-1]))), names(poly)[-1])
out2 <- cbind(out2, dt)
dim(out2)
rm(dt)

## Check
dim(out1)[1] + dim(out2)[1] == dim(point)[1]
ncol(out1) == ncol(out2)

## Combine both tables
out <- rbind(out1, out2); dim(out)
rm(out1, out2, point, poly)
setDT(out, key = "spfile")


## Additional info in output table ####
## Total records within PAA
paa <- fread(file = file.path(output_dir, "species_points_in_paa.csv"))

dim(out)[1] == dim(paa)[1]
all(out$spfile == paa$spfile)

setDT(paa, key = "spfile")
out <- merge(out, paa, by = "spfile")
rm(paa)

message(cat("Number of species with 0 records inside PAA: "),
        nrow(out[paa_points == 0]))
message(cat("Number of species with at least 1 record inside PAA: "),
        nrow(out[paa_points != 0]))


## Total number of points overlaping GEEBAM
grep("_Points$", names(out), value = TRUE)[1:5]
num <- grep("_Points$", names(out), value = TRUE)[1:5]
out$Total_fire_points <- rowSums(out[, ..num])
nrow(out[Total_fire_points == 0])

## Total polygon area overlaping GEEBAM
grep("_Area$", names(out), value = TRUE)[1:5]
num <- grep("_Area$", names(out), value = TRUE)[1:5]
out$Total_fire_polygon <- rowSums(out[, ..num])
nrow(out[Total_fire_polygon == 0])
rm(num)

## Save output table ####
## Reorder columns
names(out)[c(1, 11:14, 2:6, 8:10, 28, 27, 7, 17, 18:22, 24:26, 29, 23, 15, 16)]
out <- out[, c(1, 11:14, 2:6, 8:10, 28, 27, 7, 17, 18:22, 24:26, 29, 23, 15, 16)]

setDT(out, key = "spfile")
write.csv(out, file = file.path(output_dir, paste0("fungi_overlap_", Sys.Date(), ".csv")), 
          row.names = FALSE)


## Checks ####
out <- fread(file.path(output_dir, "fungi_overlap_2021-10-12.csv"))

## >> Check that species with and without polygons add up to total species  ####
sum(is.na(out$Species_Polygon)) + sum(!is.na(out$Species_Polygon)) == nrow(out)

## Check that unique occurrence points < total occurrence points
all(out$Nbe_unique_occ <= out$Occurrence_Points)

## Check that points with paa < total occurrence points
all(out$paa_points <= out$Occurrence_Points)


## >> Check for Species_Polygon = 0 ####
message(cat("Number of species with Species_Polygon area = 0: "),
        sum(out$Species_Polygon == 0))


## >> Check if total points in PAA = 0, the fire overlap points = 0 ####
## All TRUE
all(out[paa_points == 0]$Total_fire_points == 0)


## >> Check if total points in PAA = 0, then polygon fire overlap = 0 ####
## All TRUE
all(out[paa_points == 0]$Total_fire_polygon == 0, na.rm = TRUE)


## >> Check if total points in PAA = 0, species polygon (raster) is also = 0 ####
## Species polygon (raster, cipped to PAA) 
## This can be FALSE because species polygon clipped to PAA may or may not overlap with fire within PAA
## All TRUE
all(out[paa_points == 0 & !is.na(Total_fire_polygon)]$Species_Polygon != 0)


## >> Check if polygon fire overlap < species polygon clipped to PAA ####
## Area of polygon fire overlap will always be LESS THAN OR EQUAL TO area of species polygon clipped to PAA
## All TRUE
all(out[!is.na(Species_Polygon)]$Total_fire_polygon <= out[!is.na(Species_Polygon)]$Species_Polygon)
out[!is.na(Species_Polygon) & Total_fire_polygon == Species_Polygon]
summary(out$Overlap_Polygons_Fire2345_GEEBAM2_as_burnt)
summary(out$Total_fire_polygon)

## >> Check if Species_Polygon (raster) > EOO (shapefile) ####
## Area of species polygon raster clipped to PAA compared to area of species polygon shapefile from ConR
all(out[!is.na(Species_Polygon)]$Species_Polygon > out[!is.na(Species_Polygon)]$EOO)

plot(density(out$Species_Polygon/out$EOO, na.rm = TRUE))
range(out$Species_Polygon/out$EOO, na.rm = TRUE)
## some values are marginally > 1
## This might be a problem because of species raster (Species_Polygon) ...

## Check for Species_Polygon/EOO > 1 
nrow(out[Species_Polygon/EOO > 1])
out[Species_Polygon/EOO > 1][,.(spfile, Total_fire_points, Occurrence_Points, Total_fire_polygon, Species_Polygon, EOO, paa_points)]



## Save table truncated to species with > 0 records in PAA ####
## List of species in/out of PAA
dim(out); out1 <- out[paa_points > 0]; dim(out1)
setDT(out1, key = "spfile")
write.csv(out1, file = file.path(output_dir, 
                                 paste0("fungi_overlap_PAAonly_", 
                                        Sys.Date(), ".csv")), 
          row.names = FALSE)


## Table summary ####
message(cat("Total number of species with outputs: "),
        nrow(out))


message(cat("Number of species with polygon overlaps: "),
        sum(!is.na(out$Species_Polygon)))
message(cat("Number of species without polygon overlaps: "),
        sum(is.na(out$Species_Polygon)))


message(cat("Number of species with 0 records inside PAA: "),
        nrow(out[paa_points == 0]))
message(cat("Number of species with at least 1 record inside PAA: "),
        nrow(out[paa_points != 0]))
