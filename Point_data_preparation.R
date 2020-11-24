#load necessary packages

library(raster);library(rgdal);library(rgeos)

#set paths to files that will be used (species list available in Supplementary Information)

wd_points <- "C:/Users/ca13kute/Documents/bRacatus/bRacatus_data/GBIF_occurrence_not_process"

#import GBIF point records

setwd(wd_points)

sps_gbif <- list.files()

gbif <- lapply(sps_gbif,readRDS) #load gbif data

gbif2 <- lapply(gbif,function(x){ #select only records with lat and lon
  x[complete.cases(x[,c("decimalLatitude","decimalLongitude")]),]})

gbif3 <- lapply(gbif2,function(x){    #select only the relevant col
  x[,c("species","decimalLongitude","decimalLatitude")]})

gbif4 <- lapply(gbif3,function(x){   #add a col with record type
  cbind(x,type="GBIF")})

gbif_sp <- lapply(gbif4,as.data.frame)  #create spatial points from gbif records

i=1
for(i in 1:length(gbif_sp))
{
  coordinates(gbif_sp[[i]]) <- ~decimalLongitude+decimalLatitude
  proj4string(gbif_sp[[i]]) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
}

#rarefy points, keeping one per grid cell (0.5 degree)

setwd("C:/Users/ca13kute/Documents/bRacatus/bRacatus_data")

ID_raster <- raster("ID_raster.img")  #load raster with cell IDs

gbif_ID <- lapply(gbif_sp,function(x){extract(ID_raster,x)}) #extract cell ID from each gbif record

gbif5 <- list()
i=1
for(i in 1:length(gbif4))
{
  gbif5[[i]] <- cbind(gbif4[[i]],gbif_ID[[i]])
  names(gbif5[[i]])[5] <- "cell_ID"
}

gbif6 <- lapply(gbif5,function(x){x[!duplicated(x$cell_ID),]})

### save rarefied version of the datasets

setwd("C:/Users/ca13kute/Documents/bRacatus/bRacatus_data/GBIF_occurrence_unique_cell")

for(i in 1:length(gbif6))
{
  saveRDS(gbif6[[i]],sps_gbif[i])
}

#load suitability maps for each species in its native range

setwd("/gpfs1/data/idiv_meyer/suitable_habitat/Esa_10km_NoElev_suit_unWeight/mamm/2015")

suit_lst <- list.files()

# load maps that are available for the selected species

suit_maps <- list()
for(i in 1:length(gbif6))
{
  a <- grep(names(gbif6)[i],suit_lst)
  if(length(a)!=0){
    suit_maps[[length(suit_maps)+1]] <- raster(suit_lst[a])
    names(suit_maps)[length(suit_maps)] <- names(gbif6)[i]
  }
  print(i)
}

suit_hab2 <- suit_maps

i=1
for(i in 1:length(suit_hab2))   #transform 0 values in NA
{
  suit_hab2[[i]][] <- ifelse(suit_hab2[[i]][]==0,NA,suit_hab2[[i]][])
  suit_hab2[[i]] <- aggregate(suit_hab2[[i]],fact=6)
  print(i)
}

### save simplified mammal suitability maps

setwd("/gpfs1/data/idiv_meyer/01_projects/eduardo/Suitability_selected_mammals_half_degree")

mapply(function(x,y) writeRaster(x,filename=y,datatype="img",overwrite=F),
       suit_hab2,names(suit_hab2))

# suit_hab3 <- lapply(suit_hab2,function(x){aggregate(x,fact=6)})   #aggregate cells to 0.5 resolution


#load, simplify and save range maps per group (in parallel sessions)

wd_anura <- "/gpfs1/data/idiv_meyer/01_projects/eduardo/Range_maps/AMPHIBIANS/ANURA/Anura_ranges_sps"
wd_anura2 <- "/gpfs1/data/idiv_meyer/01_projects/eduardo/Range_maps/AMPHIBIANS/ANURA/Anura_ranges_sps_simp"

range_anura <- list()
for(i in 1:length(gbif6))
{
  a <- try(readOGR(dsn=wd_anura,layer=gsub("_"," ",names(gbif6)[i])),silent=T)
  if(class(a)!="try-error"){
    b <- gSimplify(a,.1,topologyPreserve=T)
    c <- SpatialPolygonsDataFrame(b, data=a@data)
    range_anura[[length(range_anura)+1]] <- c
    names(range_anura)[length(range_anura)] <- names(gbif6)[i]
    writeOGR(c,dsn=wd_anura2,layer=names(gbif6)[i],driver="ESRI Shapefile")
  }
  print(i)
}

wd_mammal <- "/gpfs1/data/idiv_meyer/01_projects/eduardo/Range_maps/TERRESTRIAL_MAMMALS/Terrestrial_mammals_ranges_sps"
wd_mammal2 <- "/gpfs1/data/idiv_meyer/01_projects/eduardo/Range_maps/TERRESTRIAL_MAMMALS/Terrestrial_mammals_ranges_sps_simp"

range_mammal <- list()
for(i in 1:length(gbif6))
{
  a <- try(readOGR(dsn=wd_mammal,layer=gsub("_"," ",names(gbif6)[i])),silent=T)
  if(class(a)!="try-error"){
    b <- gSimplify(a,.1,topologyPreserve=T)
    c <- SpatialPolygonsDataFrame(b, data=a@data)
    range_mammal[[length(range_mammal)+1]] <- c
    names(range_mammal)[length(range_mammal)] <- names(gbif6)[i]
    writeOGR(c,dsn=wd_mammal2,layer=names(gbif6)[i],driver="ESRI Shapefile")
  }
  print(i)
}

wd_aves <- "/gpfs1/data/idiv_meyer/01_projects/eduardo/Range_maps/AVES/Aves_ranges_sps"
wd_aves2 <- "/gpfs1/data/idiv_meyer/01_projects/eduardo/Range_maps/AVES/Aves_ranges_sps_simp"

range_aves <- list()
for(i in 1:length(gbif6))
{
  a <- try(readOGR(dsn=wd_aves,layer=gsub("_"," ",names(gbif6)[i])),silent=T)
  if(class(a)!="try-error"){
    b <- gSimplify(a,.1,topologyPreserve=T)
    c <- SpatialPolygonsDataFrame(b, data=a@data)
    range_aves[[length(range_aves)+1]] <- c
    names(range_aves)[length(range_aves)] <- names(gbif6)[i]
    writeOGR(c,dsn=wd_aves2,layer=names(gbif6)[i],driver="ESRI Shapefile")
  }
  print(i)
}



range_mam <- "C:\\Users\\ca13kute\\Documents\\GIT_projects\\Range_maps\\Mammals_corrected"

setwd(range_mam)
range_lst <- gsub(".shp","",list.files(pattern="shp$"))

ranges <- lapply(range_lst,FUN=function(x){
  readOGR(dsn=range_mam,layer=x)})

sps <- unique(gsub("(^.*)_.*$","\\1",range_lst)) #list unique sps names


#create one map per sps with only native range

ranges_sps_native <- lapply(range_lst[grep(pattern="_native$",range_lst)],function(x){
  readOGR(dsn=range_mam,layer=x)})

ranges_sps_native2 <- lapply(ranges_sps_native,function(x){
  gSimplify(x,.1)})    #simplify the polygon lines

#create one map per sps (native,alien and unknown)

ranges_sps <- list()
i=1
for(i in 1:length(sps))
{
  a <- ranges[grep(sps[i],range_lst)]
  if(length(a)>1){
    ranges_sps[[i]] <- a[[1]]
    for(j in 2:length(a))
    {
      ranges_sps[[i]] <- gUnion(ranges_sps[[i]],a[[j]])
      ranges_sps[[i]]$data <- 0
    }
  }else{
    ranges_sps[[i]] <- a[[1]]
  }
  print(i)
}

#fix bug due to holes in shp

ranges_sps[[1]] <- gBuffer(ranges_sps[[1]],width=.1)
ranges_sps[[1]]$data <- 0

ranges_sps[[10]] <- gBuffer(ranges_sps[[10]],width=.1)
ranges_sps[[10]]$data <- 0

#############################################

#####suit hab comparison###

suit_hab2

setwd(paste0(wd_local,"\\01_projects\\eduardo\\Suit_hab_2"))

a <- list.files()

suit_hab_ruben <- lapply(a,raster)

par(mfrow=c(1,2))

plot(suit_hab2[[19]])
plot(ranges_sps[[19]],add=T)

plot(suit_hab_ruben[[5]])
plot(ranges_sps[[19]],add=T)


##############################

wd_suit <- paste0(wd,"/03_tests/ESHM_edu/habitat_maps_area")

setwd(wd_suit)

files <- list.files()


maps <- list()
for(i in 1:length(files))
{
  maps[[i]] <- raster(files[i])
  print(i)
}

range <- list()
for(i in 1:length(maps))
{
  range[[i]] <- range(maps[[i]][],na.rm=T)
  print(i)
}

range[[100]]
range(maps2[[10]][],na.rm=T)
range(maps[[1]][])

which(!is.na(maps2[[10]][]))

a <- as.numeric(which(sapply(range,class)=="integer"))
maps2 <- maps[a]

plot(maps2[[8]])


b <- c(1,4,5)
class(a)

plot(maps[[1]])

class(test3)
class(test)

plot(test)
plot(test2)
plot(test3)

big <- character()
for(i in 1:ncell(test3))
{
  if(!is.na(test3[i])){
    big[length(big)+1] <- test3[1]}
  print(i)
}


#extract suitability index for GBIF points

gbif_sp2 <- lapply(gbif6,as.data.frame)  #create spatial points from gbif records

i=1
for(i in 1:length(gbif_sp2))
{
  coordinates(gbif_sp2[[i]]) <- ~decimalLongitude+decimalLatitude
  proj4string(gbif_sp2[[i]]) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
}

gbif_suit <- mapply(function(x,y) extract(x,y), suit_hab3, gbif_sp2)

gbif7 <- list()
i=1
for(i in 1:length(gbif6))
{
  gbif7[[i]] <- cbind(gbif6[[i]],gbif_suit[[i]])
  names(gbif7[[i]])[6] <- "suitability"
}


