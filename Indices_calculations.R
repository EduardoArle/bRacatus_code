library(rgdal);library(raster);library(rgeos)

#set paths to files that will be used (species list available in Supplementary Information)

wd_ranges <- "C:/Users/ca13kute/Documents/bRacatus/bRacatus_data/Range_maps"
wd_sim_pts <- "C:/Users/ca13kute/Documents/bRacatus/bRacatus_data/GBIF_occurrence"


#load sps point records

setwd(wd_sim_pts)

sps_names <- list.files()

gbif_rpts_status <- lapply(sps_names,readRDS)

#load the iucn shp

setwd(wd_ranges)

a <- list.files(pattern=".shp$")
b <- gsub(".shp","",a)

ranges <- readOGR(dsn=wd_ranges,b)  #shp terrestrial mammals


#create one raster per sps with cell ID values

setwd("C:\\Users\\ca13kute\\Documents\\Workshop_21_06_18") #load raster with ID_values for land
ID_raster <- raster("ID_raster.img")

range_simp <- lapply(range_sps,function(x){
              gSimplify(x,.05)})

range_buffer <- lapply(range_simp,function(x){
              gBuffer(x,width=.25)})

range_raster <- lapply(range_buffer,function(x){
              mask(ID_raster,mask=x)})

#obtain point records' cell ID values

gbif_sp <- lapply(gbif_rpts_status,as.data.frame)

i=1
for(i in 1:length(gbif_sp))
{
  coordinates(gbif_sp[[i]]) <- ~decimalLongitude+decimalLatitude
  proj4string(gbif_sp[[i]]) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  print(i)
}

ID_points <- lapply(gbif_sp,function(x){
            extract(ID_raster,x)})

gbif_rpts_status2 <- list()
i=1
for(i in 1:length(gbif_rpts_status))
{
  gbif_rpts_status2[[i]] <- cbind(gbif_rpts_status[[i]],ID_points[[i]])
}

## import file for distances to each ID wich points


ID_raster_comp <- which(!is.na(ID_raster[])) #ID of cells on land

gbif_rpts_status3 <- gbif_rpts_status2
i=1
for(i in 1:length(gbif_rpts_status2))
{
  dist <- numeric()
  j=1
  for(j in 1:nrow(gbif_rpts_status2[[i]]))
  {
    sps_range_ID <- which(!is.na(range_raster[[i]][]))
    if(is.na(ID_points[[i]][j])){
      dist[j] <- NA
    }else{
      setwd("C:\\Users\\ca13kute\\Documents\\Workshop_21_06_18\\Distance_between_cells_compact")
      a <- readRDS(paste(ID_points[[i]][j]))  #load file with distances from point
      b <- as.data.frame(cbind(ID_raster_comp,a))  #bind ID_cells and distance
      c <- b[which(b$ID_raster_comp %in% sps_range_ID),] #distance of every cell in the range to each point
      dist[j] <- min(c$a)  #minimum distace to the polygon
      print(j)
    }
  }
  gbif_rpts_status3[[i]] <- cbind(gbif_rpts_status2[[i]],dist*100)
  setwd("C:\\Users\\ca13kute\\Documents\\GIT_projects\\GBIF_record_status_distance\\TERRESTRIAL_MAMMALS9T")
  saveRDS(gbif_rpts_status3[[i]],sps_names[i])
  print(i)
}
