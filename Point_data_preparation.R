#load necessary packages
#
library(raster);library(rgdal);library(rgeos);library(rworldmap)

#set paths to files that will be used (species list available in Supplementary Information)

wd_points <- "C:/Users/ca13kute/Documents/bRacatus/bRacatus_data/GBIF_occurrence_not_process"
wd_suit_layers <- "C:/Users/ca13kute/Documents/bRacatus/bRacatus_data/Suitability_layers"
wd_ranges <- "C:/Users/ca13kute/Documents/bRacatus/bRacatus_data/Range_maps"
wd_ecorregions <- "C:/Users/ca13kute/Documents/bRacatus/bRacatus_data/Ecoregions"
wd_sim_pts <- "C:/Users/ca13kute/Documents/bRacatus/bRacatus_data/GBIF_occurrence"

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

eco_reg <- readOGR(dsn=wd_ecorregions,layer="wwf_terr_ecos") #load terrestrial ecoregions shp

for(i in 1:length(gbif6))
{
  pts <- gbif6[[i]]

  setwd(wd_suit_layers) #load sps suitability map
  suit <- raster(paste0(sps_gbif[i],".tif"))

  range <- readOGR(sps_names[i],dsn=wd_ranges) #load range map for the species

  #select ecoregions out of ranges
  eco_reg_abs <- eco_reg[which(eco_reg$OBJECTID %ni% eco_reg[subset(range),]$OBJECTID),]
  #seed half the number of presence gbif points out of the ecoregions
  random_pts <- spsample(eco_reg_abs,n=ceiling(nrow(pts_pr)/2),type="random")
  #prepare and add the random points to gbif data
  coords_rpts <- as.data.table(coordinates(random_pts))
  names(coords_rpts) <- c("decimalLongitude","decimalLatitude")
  random_pts_ID <- extract(ID_raster,random_pts) #extract cell ID from each gbif record
  final_rpts <- cbind(species=gsub("_"," ",sps_names[i]),coords_rpts,occurrenceStatus="simulated",
                      type="Random_false",cell_ID=random_pts_ID)
  pts2 <- rbind(pts[,-1],final_rpts)

  #crop the suitability maps by the range plus one degree buffer

  range_buffer <- gBuffer(range,width=1) #create a buffer of 1 degree on the range map

  a <- extend(suit2,extent(range_buffer),snap="out")
  b <- rasterize(range_buffer,a) #create a mask to crop the raster
  suit3 <- mask(a,mask=b)

  #get coordinates and values of cells with suitability for the species

  suit4 <- as.data.table(rasterToPoints(suit3))

  #create max of 3/12 GBIF number as pseudo false records within cells with 0% suitability
  if(length(which(suit4[,3]==0))==0){
    rows <- 0
  }else{
    rows <-as.numeric(unique(sample(as.character(which(suit4[,3]==0)),
                                    ceiling(nrow(pts_pr)*3/12),replace=T)))
  }
  pseudo_false1 <- suit4[rows,]
  names(pseudo_false1) <- c("decimalLongitude","decimalLatitude","suitability")

  #create max of 2/12 GBIF number as pseudo false records within cells with less than 25% suitability
  if(length(which(suit4[,3]<=25))==0){
    rows <- 0
  }else{
    rows <-as.numeric(unique(sample(as.character(which(suit4[,3]<=25)),
                                    ceiling(nrow(pts_pr)*2/12),replace=T)))
  }
  pseudo_false2 <- suit4[rows,]
  names(pseudo_false2) <- c("decimalLongitude","decimalLatitude","suitability")

  #create max of 1/12 GBIF number as pseudo false records within cells with between 25 and 50% suitability
  #create max of 1/12 GBIF number as pseudo true records within cells with between 25% and 50% suitability
  if(length(which(suit4[,3]>25 & suit4[,3]<=50))==0){
    rows <- 0
  }else{
    rows <- as.numeric(unique(sample(as.character(which(suit4[,3]>25 & suit4[,3]<=50)),
                                     ceiling(nrow(pts_pr)*2/12),replace=T)))
  }
  a <- sample(seq_len(length(rows)),size=ceiling(1/2*length(rows)))
  pseudo_false3 <- suit4[as.numeric(rows[-a]),]
  names(pseudo_false3) <- c("decimalLongitude","decimalLatitude","suitability")
  pseudo_true1 <- suit4[as.numeric(rows[a]),]
  names(pseudo_true1) <- c("decimalLongitude","decimalLatitude","suitability")

  #create max of 2/12 GBIF number as pseudo true records within cells with between 50% and 75% suitability
  if(length(which(suit4[,3]>50 & suit4[,3]<=75))==0){
    rows <- 0
  }else{
    rows <- as.numeric(unique(sample(as.character(which(suit4[,3]>50 & suit4[,3]<=75)),
                                     ceiling(nrow(pts_pr)*2/12),replace=T)))
  }
  pseudo_true2 <- suit4[rows,]
  names(pseudo_true2) <- c("decimalLongitude","decimalLatitude","suitability")

  #create max of 3/12 GBIF number as pseudo true records within cells with more than 75% suitability
  if(length(which(suit4[,3]>75))==0){
    rows <- 0
  }else{
    rows <-as.numeric(unique(sample(as.character(which(suit4[,3]>75)),
                                    ceiling(nrow(pts_pr)*3/12),replace=T)))
  }
  pseudo_true3 <- suit4[rows,]
  names(pseudo_true3) <- c("decimalLongitude","decimalLatitude","suitability")

  #rbind all pseudo_true and all pseudo_false points per species
  pseudo_true <- rbind(pseudo_true1,pseudo_true2,pseudo_true3)
  pseudo_false <- rbind(pseudo_false1,pseudo_false2,pseudo_false3)
  #create spatial points from pseudo_true and pseudo_false
  pseudo_true_sp <- as.data.frame(pseudo_true)
  if(nrow(pseudo_true_sp)>0){
    coordinates(pseudo_true_sp) <- ~decimalLongitude+decimalLatitude
    proj4string(pseudo_true_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    #get cell_ID of pseudo_true points
    pseudo_true_ID <- extract(ID_raster,pseudo_true_sp)
    #standardise pseudo_true points to gbif
    final_pseudo_true <- cbind(species=sps_names[i],pseudo_true[,c(1,2),],occurrenceStatus="simulated",type="Pseudo_true",cell_ID=pseudo_true_ID)
    #join gbif and pseudo_true points
    pts2 <- rbind(pts2,final_pseudo_true)
  }
  pseudo_false_sp <- as.data.frame(pseudo_false)
  if(nrow(pseudo_false_sp)>0){
    coordinates(pseudo_false_sp) <- ~decimalLongitude+decimalLatitude
    proj4string(pseudo_false_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    #get cell_ID of pseudo_false points
    pseudo_false_ID <- extract(ID_raster,pseudo_false_sp)
    #standardise pseudo_false points to gbif
    final_pseudo_false <- cbind(species=sps_names[i],pseudo_false[,c(1,2)],occurrenceStatus="simulated",type="Pseudo_false",cell_ID=pseudo_false_ID)
    #join gbif and pseudo_false points
    pts2 <- rbind(pts2,final_pseudo_false)
  }

  setwd(wd_sim_pts)
  saveRDS(pts2,sps_names[i])
  print(i)
}
