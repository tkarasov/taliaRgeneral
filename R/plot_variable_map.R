
# tutorial on sf package: https://geocompr.robinlovelace.net/spatial-operations.html#spatial-ras
library(sf)
library(raster)
library(spData)
library(spDataLarge)
library(tmap)

# convert data to sf
path = "/ebio"
OTU5_tab = read.table(paste(path, "/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/all_runs/demult_python/OTU5_RA_table.txt", sep =""), sep="\t", header = T )
OTU5_tab$OTU5 = as.numeric(as.character(OTU5_tab$OTU5))
coordinates(OTU5_tab) = ~Latitude + Longitude
my.sf = st_as_sf(x = OTU5_tab, coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84")
my.sp = as_Spatial(my.sf)
#my.sp = spTransform(my.sp, CRS("+proj=laea"))
my.df = data.frame(my.sp)[,c("coords.x2", "coords.x1", "OTU5")]
my.raster = rasterFromXYZ(my.df)


gridded(my.sp) = TRUE
my.raster = raster(my.sp) #raster(extent(my.sp))
projection(my.raster) = "+proj=longlat +datum=WGS84"
values(my.raster) = OTU5_tab$OTU5
crs(my.raster) = coordinates(my.sp)


# subset world map. Still need to remove canary islands and northern part of Norway
europe = world[world$continent == "Europe" | world$continent == "Asia" ,]
#europe = europe[europe$iso_a2 != "RU",]
#europe = europe[europe$name_long != "Greenland",]
#europe = europe[europe$name_long != "Iceland",]

#Plot world map with data points
plot(st_geometry(europe), xlim = c(-10, 40.08079), ylim = c(28.07888, 70.65714))

plot(my.sf, add = TRUE)

#Interpolate over the region












# Calculate density over region
density <-getdensityraster(coords$longitude, coords$latitude, refraster=baseraster,dilute = 50,method='bilinear')
densityall <- focal(my.raster, w=matrix(1, 31, 31), mean)

