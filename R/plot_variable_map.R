
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
proj4string(OTU5_tab) = CRS("+proj=longlat +datum=WGS84")
my.sf = st_as_sf(x = OTU5_tab, coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84")
my.sp = as_Spatial(my.sf)
#my.sp = spTransform(my.sp, CRS("+proj=laea"))
#my.df = data.frame(my.sp)[,c("coords.x1", "coords.x2", "OTU5")]
my.df = data.frame(my.sf)

# create an empty raster object to the extent of the points
e = extent(my.sp)
r <- raster(e, ncol=100, nrow=200)

# rasterize your irregular points 
my.raster <- rasterize(my.df[, 1:2], r, my.df[,3], fun=mean)



# subset world map. Still need to remove canary islands and northern part of Norway
europe = world[world$continent == "Europe" | world$continent == "Asia" ,]
#europe = europe[europe$iso_a2 != "RU",]
#europe = europe[europe$name_long != "Greenland",]
#europe = europe[europe$name_long != "Iceland",]

#Plot world map with data points
plot(st_geometry(europe), xlim = c(-10, 40.08079), ylim = c(28.07888, 70.65714))

plot(my.sf, add = TRUE)

plot(my.raster, add = TRUE)

#Interpolate over the region












# Calculate density over region
density <-getdensityraster(coords$longitude, coords$latitude, refraster=baseraster,dilute = 50,method='bilinear')
densityall <- focal(my.raster, w=matrix(1, 31, 31), mean)

