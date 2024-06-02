# For input in MCME Model
# Get cell values of continents
library(terra)
library(raster)

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

# Get Values
Af_vals <- extract(test_map,Africa,cell=TRUE)
As_vals <- extract(test_map,Asia,cell=TRUE)
Au_vals <- extract(test_map,Australia,cell=TRUE)
OC_vals <- extract(test_map,Oceania,cell=TRUE)
SA_vals <- extract(test_map,SouthAmerica,cell=TRUE)
An_vals <- extract(test_map,Antartica,cell=TRUE)
EU_vals <- extract(test_map,Europe,cell=TRUE)
NA_vals <- extract(test_map,NorthAmerica,cell=TRUE)

# Get index values of cells
cell_index = expand.grid(1:73,1:96)
colnames(cell_index) = c("Lat", "Lon")

# NA
NA_latlon = data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("Lat", "Lon", "Cont"))))
for(i in 1:dim(NA_vals)[1]){
  NA_latlon[i,] <- c(unlist(cell_index[NA_vals$cell[i],]),"North America")
}
# AS
AS_latlon = data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("Lat", "Lon", "Cont"))))
for(i in 1:dim(As_vals)[1]){
  AS_latlon[i,] <- c(unlist(cell_index[As_vals$cell[i],]),"Asia")
}
# AU
AU_latlon = data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("Lat", "Lon", "Cont"))))
for(i in 1:dim(Au_vals)[1]){
  AU_latlon[i,] <- c(unlist(cell_index[Au_vals$cell[i],]),"Australia")
}
# OC
OC_latlon = data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("Lat", "Lon", "Cont"))))
for(i in 1:dim(OC_vals)[1]){
  OC_latlon[i,] <- c(unlist(cell_index[OC_vals$cell[i],]),"Oceania")
}
# SA
SA_latlon = data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("Lat", "Lon", "Cont"))))
for(i in 1:dim(SA_vals)[1]){
  SA_latlon[i,] <- c(unlist(cell_index[SA_vals$cell[i],]),"South America")
}
# An
AN_latlon = data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("Lat", "Lon", "Cont"))))
for(i in 1:dim(An_vals)[1]){
  AN_latlon[i,] <- c(unlist(cell_index[An_vals$cell[i],]),"Antartica")
}
# EU
EU_latlon = data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("Lat", "Lon", "Cont"))))
for(i in 1:dim(EU_vals)[1]){
  EU_latlon[i,] <- c(unlist(cell_index[EU_vals$cell[i],]),"Europe")
}
# AF
AF_latlon = data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("Lat", "Lon", "Cont"))))
for(i in 1:dim(Af_vals)[1]){
  AF_latlon[i,] <- c(unlist(cell_index[Af_vals$cell[i],]),"Africa")
}

Cont_Coords = rbind(NA_latlon, AS_latlon)
Cont_Coords = rbind(Cont_Coords, AU_latlon)
Cont_Coords = rbind(Cont_Coords, OC_latlon)
Cont_Coords = rbind(Cont_Coords, SA_latlon)
Cont_Coords = rbind(Cont_Coords, AN_latlon)
Cont_Coords = rbind(Cont_Coords, EU_latlon)
Cont_Coords = rbind(Cont_Coords, AF_latlon)

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

Copy_Cont_Coords = Cont_Coords

for(i in 1:dim(Copy_Cont_Coords)[1]){
  Copy_Cont_Coords[i, "Lat_d"] = lat_centers[as.numeric(Cont_Coords[i,"Lat"])]
  Copy_Cont_Coords[i, "Lon_d"] = lon_centers[as.numeric(Cont_Coords[i,"Lon"])]
}


# Save as text file for import into python
write.csv(Copy_Cont_Coords, "/Users/wyattpetryshen/Documents/GitHub/Metacommunity-Models-In-Python/main/Continental_Coords.csv", row.names=FALSE)


