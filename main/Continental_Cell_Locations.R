# For input in MCME Model
# Get cell values of continents
library(terra)
library(raster)

# Convert grid boxes to the lat lon values
lon_centers <- c(0.0000, 3.7500,7.5000,11.2500,15.0000,18.7500,22.5000,26.2500,
                 30.0000,33.7500,37.5000,41.2500,45.0000,48.7500,52.5000,56.2500,
                 60.0000,63.7500,67.5000,71.2500,75.0000,78.7500,82.5000,86.2500,
                 90.0000,93.7500,97.5000,101.2500,105.0000,108.7500,112.5000,116.2500,
                 120.0000,123.7500,127.5000,131.2500,135.0000,138.7500,142.5000,146.2500,
                 150.0000,153.7500,157.5000,161.2500,165.0000,168.7500,172.5000,176.2500,
                 180.0000,183.7500,187.5000,191.2500,195.0000,198.7500,202.5000,206.2500,
                 210.0000,213.7500,217.5000,221.2500,225.0000,228.7500,232.5000,236.2500,
                 240.0000,243.7500,247.5000,251.2500,255.0000,258.7500,262.5000,266.2500,
                 270.0000,273.7500,277.5000,281.2500,285.0000,288.7500,292.5000,296.2500,
                 300.0000,303.7500,307.5000,311.2500,315.0000,318.7500,322.5000,326.2500,
                 330.0000,333.7500,337.5000,341.2500,345.0000,348.7500,352.5000,356.2500)
lat_centers <- c(90.0000,87.5000,85.0000,82.5000,80.0000,77.5000,75.0000,72.5000,
                 70.0000,67.5000,65.0000,62.5000,60.0000,57.5000,55.0000,52.5000,
                 50.0000,47.5000,45.0000,42.5000,40.0000,37.5000,35.0000,32.5000,
                 30.0000,27.5000,25.0000,22.5000,20.0000,17.5000,15.0000,12.5000,
                 10.0000,7.5000,5.0000,2.5000,0.0000,-2.5000,-5.0000,-7.5000,
                 -10.0000,-12.5000,-15.0000,-17.5000,-20.0000,-22.5000,-25.0000,-27.5000,
                 -30.0000,-32.5000,-35.0000,-37.5000,-40.0000,-42.5000,-45.0000,-47.5000,
                 -50.0000,-52.5000,-55.0000,-57.5000,-60.0000,-62.5000,-65.0000,-67.5000,
                 -70.0000,-72.5000,-75.0000,-77.5000,-80.0000,-82.5000,-85.0000,-87.5000,
                 -90.0000)

# Import Rasters
Tmax_Import_4kyr = raster("/Users/wyattpetryshen/Documents/Nature Paper 2024/Climate/Proc_Bricks/Tmax.mon.4kyr.tif")

test_map <-  terra::rotate(rast(Tmax_Import_4kyr), left=F)
ext(test_map) <- c(-180,180,-90,90) 
world_cont <- sf::read_sf("/Users/wyattpetryshen/Documents/Nature Paper 2024/Climate/Maps/World_Continents_-8398826466908339531/World_Continents.shp")
crs_in = "+proj=longlat +lon_0=0 +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
Africa <- sf::st_transform(world_cont["CONTINENT"][1,], crs(crs_in))
Asia <- sf::st_transform(world_cont["CONTINENT"][2,], crs(crs_in))
Australia <- sf::st_transform(world_cont["CONTINENT"][3,], crs(crs_in))
Oceania <- sf::st_transform(world_cont["CONTINENT"][4,], crs(crs_in))
SouthAmerica <- sf::st_transform(world_cont["CONTINENT"][5,], crs(crs_in))
Antartica <- sf::st_transform(world_cont["CONTINENT"][6,], crs(crs_in))
Europe <- sf::st_transform(world_cont["CONTINENT"][7,], crs(crs_in))
NorthAmerica <- sf::st_transform(world_cont["CONTINENT"][8,], crs(crs_in))

# Template Raster
template_raster = test_map[[1]]

# Create Raster to Store Results
NA_r = template_raster
AF_r = template_raster
AS_r = template_raster
AU_r = template_raster
OC_r = template_raster
SA_r = template_raster
AN_r = template_raster
EU_r = template_raster

# Get coordinate cell index of continents
NA_r[][-extract(template_raster,NorthAmerica,cell=TRUE)[,3]] = NA
AF_r[][-extract(template_raster,Africa,cell=TRUE)[,3]] = NA
AS_r[][-extract(template_raster,Asia,cell=TRUE)[,3]] = NA
AU_r[][-extract(template_raster,Australia,cell=TRUE)[,3]] = NA
OC_r[][-extract(template_raster,Oceania,cell=TRUE)[,3]] = NA
SA_r[][-extract(template_raster,SouthAmerica,cell=TRUE)[,3]] = NA
AN_r[][-extract(template_raster,Antartica,cell=TRUE)[,3]] = NA
EU_r[][-extract(template_raster,Europe,cell=TRUE)[,3]] = NA

# Convert to array
NA_arr = as.array(NA_r)
AF_arr = as.array(AF_r)
AS_arr = as.array(AS_r)
AU_arr = as.array(AU_r)
OC_arr = as.array(OC_r)
SA_arr = as.array(SA_r)
AN_arr = as.array(AN_r)
EU_arr = as.array(EU_r)


# Convert cell index to coordinates
NA_arr_vec = which(!is.na(NA_arr), TRUE)[,1:2]
AF_arr_vec = which(!is.na(AF_arr), TRUE)[,1:2]
AS_arr_vec = which(!is.na(AS_arr), TRUE)[,1:2]
AU_arr_vec = which(!is.na(AU_arr), TRUE)[,1:2]
OC_arr_vec = which(!is.na(OC_arr), TRUE)[,1:2]
SA_arr_vec = which(!is.na(SA_arr), TRUE)[,1:2]
AN_arr_vec = which(!is.na(AN_arr), TRUE)[,1:2]
EU_arr_vec = which(!is.na(EU_arr), TRUE)[,1:2]

# Rename Columns
colnames(NA_arr_vec) = c("Lat","Lon")
colnames(AF_arr_vec) = c("Lat","Lon")
colnames(AS_arr_vec) = c("Lat","Lon")
colnames(AU_arr_vec) = c("Lat","Lon")
colnames(OC_arr_vec) = c("Lat","Lon")
colnames(SA_arr_vec) = c("Lat","Lon")
colnames(AN_arr_vec) = c("Lat","Lon")
colnames(EU_arr_vec) = c("Lat","Lon")

# Assign coordinates
NA_latlon = data.frame(matrix(ncol=5,nrow=dim(NA_arr_vec)[1], dimnames=list(NULL, c("Lat", "Lon", "Lat_d", "Lon_d","Cont"))))
NA_latlon[,"Lat"] = NA_arr_vec[,1]
NA_latlon[,"Lon"] = NA_arr_vec[,2]
NA_latlon[,"Cont"] = "North America"

for(i in 1:dim(NA_latlon)[1]){
  NA_latlon[i, "Lat_d"] = lat_centers[as.numeric(NA_arr_vec[i,1])]
  NA_latlon[i, "Lon_d"] = lon_centers[as.numeric(NA_arr_vec[i,2])]
}

AS_latlon = data.frame(matrix(ncol=5,nrow=dim(AS_arr_vec)[1], dimnames=list(NULL, c("Lat", "Lon", "Lat_d", "Lon_d","Cont"))))
AS_latlon[,"Lat"] = AS_arr_vec[,1]
AS_latlon[,"Lon"] = AS_arr_vec[,2]
AS_latlon[,"Cont"] = "Asia"

for(i in 1:dim(AS_latlon)[1]){
  AS_latlon[i, "Lat_d"] = lat_centers[as.numeric(AS_arr_vec[i,1])]
  AS_latlon[i, "Lon_d"] = lon_centers[as.numeric(AS_arr_vec[i,2])]
}

AU_latlon = data.frame(matrix(ncol=5,nrow=dim(AU_arr_vec)[1], dimnames=list(NULL, c("Lat", "Lon", "Lat_d", "Lon_d","Cont"))))
AU_latlon[,"Lat"] = AU_arr_vec[,1]
AU_latlon[,"Lon"] = AU_arr_vec[,2]
AU_latlon[,"Cont"] = "Australia"

for(i in 1:dim(AU_latlon)[1]){
  AU_latlon[i, "Lat_d"] = lat_centers[as.numeric(AU_arr_vec[i,1])]
  AU_latlon[i, "Lon_d"] = lon_centers[as.numeric(AU_arr_vec[i,2])]
}

OC_latlon = data.frame(matrix(ncol=5,nrow=dim(OC_arr_vec)[1], dimnames=list(NULL, c("Lat", "Lon", "Lat_d", "Lon_d","Cont"))))
OC_latlon[,"Lat"] = OC_arr_vec[,1]
OC_latlon[,"Lon"] = OC_arr_vec[,2]
OC_latlon[,"Cont"] = "Oceania"

for(i in 1:dim(OC_latlon)[1]){
  OC_latlon[i, "Lat_d"] = lat_centers[as.numeric(OC_arr_vec[i,1])]
  OC_latlon[i, "Lon_d"] = lon_centers[as.numeric(OC_arr_vec[i,2])]
}

SA_latlon = data.frame(matrix(ncol=5,nrow=dim(SA_arr_vec)[1], dimnames=list(NULL, c("Lat", "Lon", "Lat_d", "Lon_d","Cont"))))
SA_latlon[,"Lat"] = SA_arr_vec[,1]
SA_latlon[,"Lon"] = SA_arr_vec[,2]
SA_latlon[,"Cont"] = "South America"

for(i in 1:dim(SA_latlon)[1]){
  SA_latlon[i, "Lat_d"] = lat_centers[as.numeric(SA_arr_vec[i,1])]
  SA_latlon[i, "Lon_d"] = lon_centers[as.numeric(SA_arr_vec[i,2])]
}

AN_latlon = data.frame(matrix(ncol=5,nrow=dim(AN_arr_vec)[1], dimnames=list(NULL, c("Lat", "Lon", "Lat_d", "Lon_d","Cont"))))
AN_latlon[,"Lat"] = AN_arr_vec[,1]
AN_latlon[,"Lon"] = AN_arr_vec[,2]
AN_latlon[,"Cont"] = "Antartica"

for(i in 1:dim(AN_latlon)[1]){
  AN_latlon[i, "Lat_d"] = lat_centers[as.numeric(AN_arr_vec[i,1])]
  AN_latlon[i, "Lon_d"] = lon_centers[as.numeric(AN_arr_vec[i,2])]
}

EU_latlon = data.frame(matrix(ncol=5,nrow=dim(EU_arr_vec)[1], dimnames=list(NULL, c("Lat", "Lon", "Lat_d", "Lon_d","Cont"))))
EU_latlon[,"Lat"] = EU_arr_vec[,1]
EU_latlon[,"Lon"] = EU_arr_vec[,2]
EU_latlon[,"Cont"] = "Europe"

for(i in 1:dim(EU_latlon)[1]){
  EU_latlon[i, "Lat_d"] = lat_centers[as.numeric(EU_arr_vec[i,1])]
  EU_latlon[i, "Lon_d"] = lon_centers[as.numeric(EU_arr_vec[i,2])]
}

AF_latlon = data.frame(matrix(ncol=5,nrow=dim(AF_arr_vec)[1], dimnames=list(NULL, c("Lat", "Lon", "Lat_d", "Lon_d","Cont"))))
AF_latlon[,"Lat"] = AF_arr_vec[,1]
AF_latlon[,"Lon"] = AF_arr_vec[,2]
AF_latlon[,"Cont"] = "Africa"

for(i in 1:dim(AF_latlon)[1]){
  AF_latlon[i, "Lat_d"] = lat_centers[as.numeric(AF_arr_vec[i,1])]
  AF_latlon[i, "Lon_d"] = lon_centers[as.numeric(AF_arr_vec[i,2])]
}

# Put into Single Data Frame
Cont_Coords = rbind(NA_latlon, AS_latlon)
Cont_Coords = rbind(Cont_Coords, AU_latlon)
Cont_Coords = rbind(Cont_Coords, OC_latlon)
Cont_Coords = rbind(Cont_Coords, SA_latlon)
Cont_Coords = rbind(Cont_Coords, AN_latlon)
Cont_Coords = rbind(Cont_Coords, EU_latlon)
Cont_Coords = rbind(Cont_Coords, AF_latlon)

write.csv(Cont_Coords, "/Users/wyattpetryshen/Documents/GitHub/Metacommunity-Models-In-Python/main/World_Continent_Coords.csv", row.names=FALSE)

