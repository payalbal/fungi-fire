# Fire overlap analysis for fungi data (based on nesp_bugs)

## Workflow:
1. Data preparation and save individual species data as .rds files: https://github.com/payalbal/fungi-fire/blob/main/scripts/data_prep.R  


2. Create polygons using IUCN.eval() function from ConR package for all data: https://github.com/payalbal/fungi-fire/blob/main/scripts/species_polygons.R  

  - Function based on IUCN.eval() function from ConR package: https://github.com/payalbal/nesp_bugs/blob/master/scripts/conr_iucn_eval.R 


3. Convert .rds to .shp files: https://github.com/payalbal/nesp_bugs/blob/master/scripts/rds_to_shp_conversion.R   

4. PAA overlap for all species with point data: https://github.com/payalbal/fungi-fire/blob/main/scripts/species_points_paaoverlap.R  


5. Fire overlap for all species with point data: https://github.com/payalbal/fungi-fire/blob/main/scripts/species_points_fireoverlap.R   

  - Points overlap function: https://github.com/payalbal/nesp_bugs/blob/master/scripts/points_overlap.R 


6. Fire overlap for species with polygons: https://github.com/payalbal/fungi-fire/blob/main/scripts/species_EOO_fireoverlap_paa.R  

  - Polygon overlap function: https://github.com/payalbal/nesp_bugs/blob/master/scripts/polygon_paa_overlap.R  


7. Combine polygon and point output table: https://github.com/payalbal/fungi-fire/blob/main/scripts/fireoverlap_table_PAA.R  


