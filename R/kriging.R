# @input_data is a spatial dataframe from sp package
# @
# Need automap

kriging_world<-function(input_data, large.map, input_variable){

####################################################################################
# Smooth over the region with Kriging. 
####################################################################################

#I'm going to use the function krige from gstat
  my.transform <- "+proj=longlat +datum=WGS84 +no_defs"
  #"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  laea <- "+init=epsg:3857"  #laea for europe
  
  my.sp <- spTransform(input_data, CRS(my.transform))
  
  #I need to create a grid of points I would like to predict the values for.
  lon <- seq(extent(my.sp)[1] - 5, extent(my.sp)[2] + 5, length.out = 500)
  lat <- seq(extent(my.sp)[3] -5, extent(my.sp)[4] + 5, length.out = 500)
  grd <- expand.grid(lon = lon, lat = lat)
  grd_sf  <-  st_as_sf(grd, coords = c("lon", "lat"), crs = my.transform, agr = "constant")
  grd_sp <- as_Spatial(grd_sf)
  crs(grd_sp) = my.transform

  #automap seems a little faster
  sp.laea <- spTransform(my.sp, CRS(laea))
  map.laea <-spTransform(large.map, CRS(laea))
  grd.laea <- spTransform(grd_sp, CRS(laea))
  kriging_result = autoKrige(formula = OTU5~1, input_data = sp.laea, new_data = grd.laea)
  
  #kriging_result1 = autoKrige(formula = OTU5~1, input_data = sp.laea, new_data = europe.laea)
  Krig = kriging_result$krige_output
  # 
  # # take only the points falling in polygons
  Krig = kriging_result$krige_output[!is.na(over(kriging_result$krige_output,as(map.laea,"SpatialPolygons"))),]  
  # Krig.fin = spTransform(Krig, my.transform)
  # Krig_df = as.data.frame(Krig.fin)
  # names(Krig_df) = c(  "APPT_pred","APPT_var","APPT_stdev", "longitude","latitude")
  # my.df = as.data.frame(my.sp)

  return(c(kriging_result, Krig))
}




# ####################################################################################
# # Read in data and project data and europe to laea centered on europe
# ####################################################################################
# path = "/ebio"
# OTU5_tab = read.table(paste(path, "/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/all_runs/demult_python/OTU5_RA_table.txt", sep =""), 
#                       sep="\t", 
#                       header = T )
# OTU5_tab$OTU5 = as.numeric(as.character(OTU5_tab$OTU5))
# 
# #First add some wiggle to the OTU5 points then reproject the points
# OTU5_tab$Longitude = OTU5_tab$Longitude + runif(n=length(OTU5_tab$Longitude), min=0, max=0.000001)
# OTU5_tab$Latitude = OTU5_tab$Latitude + runif(n=length(OTU5_tab$Latitude), min=0, max=0.000001)
# 
# #Now make sure everything has the same projection
# coordinates(OTU5_tab) = c("Longitude", "Latitude")
# proj4string(OTU5_tab) = CRS(my.transform)
# my.sf = st_as_sf(x = OTU5_tab, coords = c("Longitude", "Latitude"), crs = my.transform)
# my.sp = as_Spatial(my.sf)
# my.sp = spTransform(my.sp, CRS(my.transform))
# my.df = data.frame(my.sp)
# 
# # subset world map. Still need to remove canary islands and northern part of Norway
# europe = world[world$continent == "Europe" | world$continent == "Asia" ,]
# europe = as_Spatial(europe, )
# crs(europe) = my.transform
# europe.sf = st_as_sf(europe)
# 
