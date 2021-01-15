
fitModels <- function(x){
  models <- list()
  models[[1]] <- glm(obs ~ 1, data=x, family=binomial(link="logit")) #null model
  models[[2]] <- glm(obs ~ dist, data=x, family=binomial(link="logit"))
  models[[3]] <- glm(obs ~ dist, data=x, family=binomial(link="probit"))
  models[[4]] <- glm(obs ~ dist, data=x, family=binomial(link="cauchit"))
  models[[5]] <- glm(obs ~ dist, data=x, family=binomial(link="cloglog"))
  models[[6]] <- glm(obs ~ dist+range_size, data=x, family=binomial(link="logit"))
  models[[7]] <- glm(obs ~ dist_log, data=x, family=binomial(link="logit"))
  models[[8]] <- glm(obs ~ dist_log+range_size, data=x, family=binomial(link="logit"))

  return(models)
}

dataSet <- function(x){
  dist <- x$dist
  dist2 <- dist
  dist2[which(dist2==0)] <- 1   #distance with 0 values transformed to 1 to work with log
  dist_log <- log(dist2)
  obs <- x$binary
  range_size <- x$range_size
  data <- as.data.frame(cbind(obs,dist,dist_log,range_size))
  return(data)
}

fitTest <- function(wd){
  library(rworldmap);library(rgdal);library(data.table);library(rlist)
  setwd(wd)
  sps <- list.files()  #list files in the folder
  points <- lapply(sps,readRDS)  #load point records

  i=1
  for(i in 1:length(points))    #change the name of the cols cell_ID and dist
  {
    names(points[[i]])[c(8,9)] <- c("cell_ID","dist")
  }

  points2 <- list()   #add a binary col for likely_true and likely_false
  i=1
  for(i in 1:length(points))
  {
    binary <- ifelse(points[[i]]$pseudo_status=="likely_true",1,0)
    points2[[i]] <- cbind(points[[i]],binary)
  }

  points3 <- lapply(points2,function(x){   #eliminate duplicated lines
    x[!duplicated(x),]})

  points4 <- lapply(points3,function(x){  #eliminate lines with NA in distance (points in water)
    x[!is.na(x$dist),]})

  train_sps <- list()     #subset each species' points in train and test
  test_sps <- list()
  i=1
  for(i in 1:length(points4))
  {
    a <- sample(seq_len(nrow(points4[[i]])),size=floor(.5*nrow(points4[[i]])))
    train_sps[[i]] <- points4[[i]][a,]
    test_sps[[i]] <- points4[[i]][-a,]
  }

  train_sps2 <- rbindlist(train_sps)

  train_data <- dataSet(train_sps2)

  #fit the models with all sps

  models <- fitModels(train_data)

  #apply all models to train data

  m.preds_train <- lapply(models,function(x){predict(x,train_data, type="response")})

  #prepare test dataset

  test_sps2 <- rbindlist(test_sps)

  test_data <- dataSet(test_sps2)

  #apply all models to test data

  m.preds_test <- lapply(models,function(x){predict(x,test_data, type="response")})

  dir.create(paste(getwd(),"\\Model_tests",sep="")) #create folder for results
  wd2 <- paste(getwd(),"\\Model_tests\\Trained_All_sps",sep="")
  dir.create(wd2) #create subfolder for results for all species
  setwd(wd2)

  write.csv(data.frame(AIC=sapply(models,AIC),
                       BIC=sapply(models,BIC),
                       mse_train=sapply(m.preds_train,function(x){mean((train_data$obs - x)^2)}),
                       mse_test=sapply(m.preds_test,function(x){mean((test_data$obs - x)^2)})),
            "Stats_all_sps.csv")   #save csv with AIC values

  #boxplot with models' performance

  jpeg('Boxplot_All_sps.jpg')

  par(mfrow=c(round(sqrt(length(m.preds_test))),ceiling(sqrt(length(m.preds_test)))),mar=c(2,2,2,2),bg="white")
  lapply(m.preds_test,function(x){boxplot(x~obs,test_data)})
  dev.off()

  #train the models with all sps and test to individual sps test datasets

  all_points <- as.data.table(rbind(train_sps2,test_sps2))

  train_data <- dataSet(all_points)

  #fit the models with all sps

  models <- fitModels(train_data)

  #apply all models to train data

  m.preds_train <- lapply(models,function(x){predict(x,train_data, type="response")})

  #prepare all points per species as test_datasets

  all_points_sps <- list()
  i=1
  for(i in 1:length(train_sps))
  {
    all_points_sps[[i]] <- rbind(train_sps[[i]], test_sps[[i]])
  }

  test_data <- lapply(all_points_sps,dataSet)

  #apply the model to each test dataset

  i=1
  for(i in 1:length(sps))
  {
    m.preds_test <- lapply(models,function(x){predict(x,test_data[[i]], type="response")})

    write.csv(data.frame(AIC=sapply(models,AIC),
                         BIC=sapply(models,BIC),
                         mse_train=sapply(m.preds_train,function(x){mean((train_data$obs - x)^2)}),
                         mse_test=sapply(m.preds_test,function(x){mean((test_data[[i]]$obs - x)^2)})),
              paste("Stats_",sps[i],".csv",sep=""))   #save csv with AIC values

    #boxplot with models' performance

    jpeg(paste("Boxplot_",sps[i],".jpg",sep=""))

    par(mfrow=c(round(sqrt(length(m.preds_test))),ceiling(sqrt(length(m.preds_test)))),mar=c(2,2,2,2),bg="white")
    lapply(m.preds_test,function(x){boxplot(x~obs,test_data[[i]])})
    dev.off()
  }

  #fit the and test the models for each species individually

  wd2 <- paste(wd,"\\Model_tests\\Trained_tested_each_sps",sep="")
  dir.create(wd2) #create subfolder for results for all species
  setwd(wd2)

  i=1
  for(i in 1:length(sps))
  {
    train_data <- dataSet(train_sps[[i]])
    models <- fitModels(train_data)

    #apply all models to train data

    m.preds_train <- lapply(models,function(x){predict(x,train_data, type="response")})

    #apply all models to train data

    test_data <- dataSet(test_sps[[i]])

    m.preds_test <- lapply(models,function(x){predict(x,test_data, type="response")})

    write.csv(data.frame(AIC=sapply(models,AIC),
                         BIC=sapply(models,BIC),
                         mse_train=sapply(m.preds_train,function(x){mean((train_data$obs - x)^2)}),
                         mse_test=sapply(m.preds_test,function(x){mean((test_data$obs - x)^2)})),
              paste("Stats_",sps[i],".csv",sep=""))   #save csv with AIC values

    #boxplot with models' performance

    jpeg(paste("Boxplot_",sps[i],".jpg",sep=""))

    par(mfrow=c(round(sqrt(length(m.preds_test))),ceiling(sqrt(length(m.preds_test)))),mar=c(2,2,2,2),bg="white")
    lapply(m.preds_test,function(x){boxplot(x~obs,test_data)})
    dev.off()
  }

  #fit the models with each sps data and test for all and each sps

  ind_sps_data <- lapply(all_points_sps,dataSet)

  i=1
  for(i in 1:length(sps))
  {

    models <- fitModels(ind_sps_data[[i]])

    #apply each sps'model to train data

    m.preds_train <- lapply(models,function(x){predict(x,ind_sps_data[[i]], type="response")})

    #apply each sps' model to full data

    test_data <- dataSet(all_points)

    m.preds_test <- lapply(models,function(x){predict(x,test_data, type="response")})

    wd3 <- paste(wd,"\\Model_tests","\\Trained_",sps[[i]],sep="")
    dir.create(wd3)
    setwd(wd3)

    write.csv(data.frame(AIC=sapply(models,AIC),
                         BIC=sapply(models,BIC),
                         mse_train=sapply(m.preds_train,function(x){mean((ind_sps_data[[i]]$obs - x)^2)}),
                         mse_test=sapply(m.preds_test,function(x){mean((test_data$obs - x)^2)})),
              paste("Stats_all_sps.csv",sep=""))   #save csv with AIC values

    #boxplot with models' performance

    jpeg(paste("Boxplot_all_sps.jpg",sep=""))

    par(mfrow=c(round(sqrt(length(m.preds_test))),ceiling(sqrt(length(m.preds_test)))),mar=c(2,2,2,2),bg="white")
    lapply(m.preds_test,function(x){boxplot(x~obs,test_data)})
    dev.off()

    #apply models trained with each individual sps to each other as test dataset

    j=1
    for(j in 1:length(sps))
    {
      test_data <- ind_sps_data[[j]]

      m.preds_test <- lapply(models,function(x){predict(x,test_data, type="response")})

      write.csv(data.frame(AIC=sapply(models,AIC),
                           BIC=sapply(models,BIC),
                           mse_train=sapply(m.preds_train,function(x){mean((ind_sps_data[[i]]$obs - x)^2)}),
                           mse_test=sapply(m.preds_test,function(x){mean((test_data$obs - x)^2)})),
                paste("Stats_",sps[j],".csv",sep=""))   #save csv with AIC values

      #boxplot with models' performance

      jpeg(paste("Boxplot_",sps[j],".jpg",sep=""))

      par(mfrow=c(round(sqrt(length(m.preds_test))),ceiling(sqrt(length(m.preds_test)))),mar=c(2,2,2,2),bg="white")
      lapply(m.preds_test,function(x){boxplot(x~obs,test_data)})
      dev.off()
    }
  }

  #fit and test the models to the whole world and to each continent

  train_data <- dataSet(train_sps2)

  #fit the models with all sps

  models <- fitModels(train_data)

  #apply all models to train data

  m.preds_train <- lapply(models,function(x){predict(x,train_data, type="response")})

  #prepare test dataset

  test_sps2 <- rbindlist(test_sps)

  test_data <- dataSet(test_sps2)

  #apply all models to test data

  m.preds_test <- lapply(models,function(x){predict(x,test_data, type="response")})

  dir.create( paste(wd,"\\Model_tests\\Region",sep=""))
  wd4 <- paste(wd,"\\Model_tests\\Region\\Trained_Whole_world",sep="")
  dir.create(wd4) #create subfolder for results for the world
  setwd(wd4)

  write.csv(data.frame(AIC=sapply(models,AIC),
                       BIC=sapply(models,BIC),
                       mse_train=sapply(m.preds_train,function(x){mean((train_data$obs - x)^2)}),
                       mse_test=sapply(m.preds_test,function(x){mean((test_data$obs - x)^2)})),
            paste("Stats_whole_world.csv",sep=""))   #save csv with AIC values

  #boxplot with models' performance

  jpeg('Boxplot_Whole_world.jpg')

  par(mfrow=c(round(sqrt(length(m.preds_test))),ceiling(sqrt(length(m.preds_test)))),mar=c(2,2,2,2),bg="white")
  lapply(m.preds_test,function(x){boxplot(x~obs,test_data)})
  dev.off()

  #create spatial objects of all points

  all_points_sp <- all_points
  coordinates(all_points_sp) <- ~decimalLongitude+decimalLatitude
  proj4string(all_points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  map <- getMap() #get world map

  a <- as.character(unique(map$continent))
  b <- a[-c(which(is.na(a)),4)]    #create list of continents

  map_cont <- lapply(b,function(x){map[which(map$continent==x),]})  #create polygons for each continent
  pts_cont <- lapply(map_cont,function(x){over(all_points_sp,x)})  #extract continent info for each point
  pts_cont2 <- lapply(pts_cont,function(x){all_points[which(!is.na(x$continent)),]})  #separate datasets per continents

  #apply the model for the whole world to each continent

  i=1
  for(i in 1:length(b))
  {
    test_data <- dataSet(pts_cont2[[i]])

    m.preds_test <- lapply(models,function(x){predict(x,test_data, type="response")})

    write.csv(data.frame(AIC=sapply(models,AIC),
                         BIC=sapply(models,BIC),
                         mse_train=sapply(m.preds_train,function(x){mean((train_data$obs - x)^2)}),
                         mse_test=sapply(m.preds_test,function(x){mean((test_data$obs - x)^2)})),
              paste("Stats_",b[i],".csv",sep=""))   #save csv with AIC values

    jpeg(paste("Boxplot_",b[i],".jpg",sep=""))

    par(mfrow=c(round(sqrt(length(m.preds_test))),ceiling(sqrt(length(m.preds_test)))),mar=c(2,2,2,2),bg="white")
    lapply(m.preds_test,function(x){boxplot(x~obs,test_data)})
    dev.off()
  }

  #fit the models with each continent train data and test for all and each continent

  i=1
  for(i in 1:length(pts_cont2))
  {
    #fit the models with each continent

    train_data <- dataSet(pts_cont2[[i]])

    models <- fitModels(train_data)

    #apply all models to train data

    m.preds_train <- lapply(models,function(x){predict(x,train_data, type="response")})

    #prepare test dataset (whole world)

    test_data <- dataSet(all_points)

    #apply models for each continent to complete test data

    m.preds_test <- lapply(models,function(x){predict(x,test_data, type="response")})

    wd4 <- paste(wd,"\\Model_tests\\Region\\Trained_",b[i],sep="")
    dir.create(wd4) #create subfolder for results for the world
    setwd(wd4)

    write.csv(data.frame(AIC=sapply(models,AIC),
                         BIC=sapply(models,BIC),
                         mse_train=sapply(m.preds_train,function(x){mean((train_data$obs - x)^2)}),
                         mse_test=sapply(m.preds_test,function(x){mean((test_data$obs - x)^2)})),
              paste("Stats_whole_world.csv",sep=""))   #save csv with AIC values

    jpeg(paste("Boxplot_whole_world.jpg",sep=""))

    par(mfrow=c(round(sqrt(length(m.preds_test))),ceiling(sqrt(length(m.preds_test)))),mar=c(2,2,2,2),bg="white")
    lapply(m.preds_test,function(x){boxplot(x~obs,test_data)})
    dev.off()

    #apply models trained with each continent to each other test dataset

    j=1
    for(j in 1:length(b))
    {
      test_data <- dataSet(pts_cont2[[j]])

      m.preds_test <- lapply(models,function(x){predict(x,test_data, type="response")})

      write.csv(data.frame(AIC=sapply(models,AIC),
                           BIC=sapply(models,BIC),
                           mse_train=sapply(m.preds_train,function(x){mean((train_data$obs - x)^2)}),
                           mse_test=sapply(m.preds_test,function(x){mean((test_data$obs - x)^2)})),
                paste("Stats_",b[j],".csv",sep=""))   #save csv with AIC values

      #boxplot with models' performance

      jpeg(paste("Boxplot_",b[j],".jpg",sep=""))

      par(mfrow=c(round(sqrt(length(m.preds_test))),ceiling(sqrt(length(m.preds_test)))),mar=c(2,2,2,2),bg="white")
      lapply(m.preds_test,function(x){boxplot(x~obs,test_data)})
      dev.off()
    }
  }
}

wd <- "C:\\Users\\ca13kute\\Documents\\GIT_projects\\GBIF_record_status_distance\\TERRESTRIAL_MAMMALS11b"
setwd(wd)
getwd()
fitTest(wd)

#############################################

# modified for checklists

fitModels <- function(x){
  models <- list()
  models[[1]] <- glm(obs ~ 1, data=x, family=binomial(link="logit")) #null model
  models[[2]] <- glm(obs ~ ADI, data=x, family=binomial(link="logit"))
  models[[3]] <- glm(obs ~ ADI, data=x, family=binomial(link="probit"))
  models[[4]] <- glm(obs ~ ADI, data=x, family=binomial(link="cauchit"))
  models[[5]] <- glm(obs ~ ADI, data=x, family=binomial(link="cloglog"))
  models[[6]] <- glm(obs ~ ADI_logneg, data=x, family=binomial(link="logit"))
  models[[7]] <- glm(obs ~ ADI_1000, data=x, family=binomial(link="logit"))
  models[[8]] <- glm(obs ~ ADI_1000, data=x, family=binomial(link="probit"))
  models[[9]] <- glm(obs ~ ADI_1000, data=x, family=binomial(link="cauchit"))
  models[[10]] <- glm(obs ~ ADI_1000, data=x, family=binomial(link="cloglog"))
  models[[11]] <- glm(obs ~ ADI_1000_logneg, data=x, family=binomial(link="logit"))

  return(models)
}

dataSet <- function(x){
  dist <- x$dist
  dist2 <- dist
  dist2[which(dist2==0)] <- 1   #distance with 0 values transformed to 1 to work with log
  dist_log <- log(dist2)
  obs <- x$binary
  range_size <- x$range_size
  data <- as.data.frame(cbind(obs,dist,dist_log,range_size,
                              ADI=x$ADI,ADI_logneg=x$ADI_logneg,
                              ADI_1000=x$ADI_1000,ADI_1000_logneg=x$ADI_1000_logneg))
  return(data)
}

fitTest <- function(wd){
  library(rworldmap);library(rgdal);library(data.table);library(rlist)
  setwd(wd)
  sps <- list.files()  #list files in the folder
  points <- lapply(sps,readRDS)  #load point records

  i=1
  for(i in 1:length(points))    #change the name of the cols cell_ID and dist
  {
    names(points[[i]])[c(8,9)] <- c("cell_ID","dist")
  }

  points2 <- list()   #add a binary col for likely_true and likely_false
  i=1
  for(i in 1:length(points))
  {
    binary <- ifelse(points[[i]]$pseudo_status=="likely_true",1,0)
    points2[[i]] <- cbind(points[[i]],binary)
  }

  points3 <- lapply(points2,function(x){   #eliminate duplicated lines
    x[!duplicated(x),]})

  points4 <- lapply(points3,function(x){  #eliminate lines with NA in distance (points in water)
    x[!is.na(x$dist),]})

  train_sps <- list()     #subset each species' points in train and test
  test_sps <- list()
  i=1
  for(i in 1:length(points4))
  {
    a <- sample(seq_len(nrow(points4[[i]])),size=floor(.5*nrow(points4[[i]])))
    train_sps[[i]] <- points4[[i]][a,]
    test_sps[[i]] <- points4[[i]][-a,]
  }

  train_sps2 <- rbindlist(train_sps)

  train_data <- dataSet(train_sps2)

  #fit the models with all sps

  models <- fitModels(train_data)

  #apply all models to train data

  m.preds_train <- lapply(models,function(x){predict(x,train_data, type="response")})

  #prepare test dataset

  test_sps2 <- rbindlist(test_sps)

  test_data <- dataSet(test_sps2)

  #apply all models to test data

  m.preds_test <- lapply(models,function(x){predict(x,test_data, type="response")})

  dir.create(paste(getwd(),"\\Model_tests",sep="")) #create folder for results
  wd2 <- paste(getwd(),"\\Model_tests\\Trained_All_sps",sep="")
  dir.create(wd2) #create subfolder for results for all species
  setwd(wd2)

  write.csv(data.frame(AIC=sapply(models,AIC),
                       BIC=sapply(models,BIC),
                       mse_train=sapply(m.preds_train,function(x){mean((train_data$obs - x)^2,na.rm=T)}),
                       mse_test=sapply(m.preds_test,function(x){mean((test_data$obs - x)^2,na.rm=T)})),
            "Stats_all_sps.csv")   #save csv with AIC values

  #boxplot with models' performance

  jpeg('Boxplot_All_sps.jpg')

  par(mfrow=c(round(sqrt(length(m.preds_test))),ceiling(sqrt(length(m.preds_test)))),mar=c(2,2,2,2),bg="white")
  lapply(m.preds_test,function(x){boxplot(x~obs,test_data)})
  dev.off()

  #train the models with all sps and test to individual sps test datasets

  all_points <- as.data.table(rbind(train_sps2,test_sps2))

  train_data <- dataSet(all_points)

  #fit the models with all sps

  models <- fitModels(train_data)

  #apply all models to train data

  m.preds_train <- lapply(models,function(x){predict(x,train_data, type="response")})

  #prepare all points per species as test_datasets

  all_points_sps <- list()
  i=1
  for(i in 1:length(train_sps))
  {
    all_points_sps[[i]] <- rbind(train_sps[[i]], test_sps[[i]])
  }

  test_data <- lapply(all_points_sps,dataSet)

  #apply the model to each test dataset

  i=1
  for(i in 1:length(sps))
  {
    m.preds_test <- lapply(models,function(x){predict(x,test_data[[i]], type="response")})

    write.csv(data.frame(AIC=sapply(models,AIC),
                         BIC=sapply(models,BIC),
                         mse_train=sapply(m.preds_train,function(x){mean((train_data$obs - x)^2,na.rm=T)}),
                         mse_test=sapply(m.preds_test,function(x){mean((test_data[[i]]$obs - x)^2,na.rm=T)})),
              paste("Stats_",sps[i],".csv",sep=""))   #save csv with AIC values

    #boxplot with models' performance

    jpeg(paste("Boxplot_",sps[i],".jpg",sep=""))

    par(mfrow=c(round(sqrt(length(m.preds_test))),ceiling(sqrt(length(m.preds_test)))),mar=c(2,2,2,2),bg="white")
    lapply(m.preds_test,function(x){
      if(length(grep(F,is.na(x)))==0){
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, paste("NA"),cex = 6, col = "black")
      }else{
        boxplot(x~obs,test_data[[i]])
      }
    })
    dev.off()
  }

  #fit the and test the models for each species individually

  wd2 <- paste(wd,"\\Model_tests\\Trained_tested_each_sps",sep="")
  dir.create(wd2) #create subfolder for results for all species
  setwd(wd2)

  i=1
  for(i in 1:length(sps))
  {
    train_data <- dataSet(train_sps[[i]])

    if(length(grep(F,is.na(train_data$ADI)))!=0)
    {
      models <- fitModels(train_data)

      #apply all models to train data

      m.preds_train <- lapply(models,function(x){predict(x,train_data, type="response")})

      #apply all models to train data

      test_data <- dataSet(test_sps[[i]])

      m.preds_test <- lapply(models,function(x){predict(x,test_data, type="response")})

      write.csv(data.frame(AIC=sapply(models,AIC),
                           BIC=sapply(models,BIC),
                           mse_train=sapply(m.preds_train,function(x){mean((train_data$obs - x)^2,na.rm=T)}),
                           mse_test=sapply(m.preds_test,function(x){mean((test_data$obs - x)^2,na.rm=T)})),
                paste("Stats_",sps[i],".csv",sep=""))   #save csv with AIC values

      #boxplot with models' performance

      jpeg(paste("Boxplot_",sps[i],".jpg",sep=""))

      par(mfrow=c(round(sqrt(length(m.preds_test))),ceiling(sqrt(length(m.preds_test)))),mar=c(2,2,2,2),bg="white")
      lapply(m.preds_test,function(x){boxplot(x~obs,test_data)})
      dev.off()
    }else{
      write.csv(data.frame("No checklists availabe"),
                paste("Stats_",sps[i],".csv",sep=""))   #save csv with AIC values

      jpeg(paste("Boxplot_",sps[i],".jpg",sep=""))

      par(mfrow=c(1,1),mar=c(2,2,2,2),bg="white")
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.5, paste("No checklists availabe"),cex = 3, col = "black")
      dev.off()
    }
  }

  #fit the models with each sps data and test for all and each sps

  ind_sps_data <- lapply(all_points_sps,dataSet)

  i=1
  for(i in 1:length(sps))
  {
    if(length(grep(F,is.na(ind_sps_data[[i]]$ADI)))!=0){
      models <- fitModels(ind_sps_data[[i]])

      #apply each sps'model to train data

      m.preds_train <- lapply(models,function(x){predict(x,ind_sps_data[[i]], type="response")})

      #apply each sps' model to full data

      test_data <- dataSet(all_points)

      m.preds_test <- lapply(models,function(x){predict(x,test_data, type="response")})

      wd3 <- paste(wd,"\\Model_tests","\\Trained_",sps[[i]],sep="")
      dir.create(wd3)
      setwd(wd3)

      write.csv(data.frame(AIC=sapply(models,AIC),
                           BIC=sapply(models,BIC),
                           mse_train=sapply(m.preds_train,function(x){mean((ind_sps_data[[i]]$obs - x)^2,na.rm=T)}),
                           mse_test=sapply(m.preds_test,function(x){mean((test_data$obs - x)^2,na.rm=T)})),
                paste("Stats_all_sps.csv",sep=""))   #save csv with AIC values

      #boxplot with models' performance

      jpeg(paste("Boxplot_all_sps.jpg",sep=""))

      par(mfrow=c(round(sqrt(length(m.preds_test))),ceiling(sqrt(length(m.preds_test)))),mar=c(2,2,2,2),bg="white")
      lapply(m.preds_test,function(x){boxplot(x~obs,test_data)})
      dev.off()

      #apply models trained with each individual sps to each other as test dataset

      j=1
      for(j in 1:length(sps))
      {
        test_data <- ind_sps_data[[j]]

        m.preds_test <- lapply(models,function(x){predict(x,test_data, type="response")})

        write.csv(data.frame(AIC=sapply(models,AIC),
                             BIC=sapply(models,BIC),
                             mse_train=sapply(m.preds_train,function(x){mean((ind_sps_data[[i]]$obs - x)^2,na.rm=T)}),
                             mse_test=sapply(m.preds_test,function(x){mean((test_data$obs - x)^2,na.rm=T)})),
                  paste("Stats_",sps[j],".csv",sep=""))   #save csv with AIC values

        #boxplot with models' performance

        jpeg(paste("Boxplot_",sps[j],".jpg",sep=""))

        par(mfrow=c(round(sqrt(length(m.preds_test))),ceiling(sqrt(length(m.preds_test)))),mar=c(2,2,2,2),bg="white")
        lapply(m.preds_test,function(x){
          if(length(grep(F,is.na(x)))==0){
            plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
            text(x = 0.5, y = 0.5, paste("NA"),cex = 6, col = "black")
          }else{
            boxplot(x~obs,test_data)
          }
        })
        dev.off()
      }
    }else{
      wd3 <- paste(wd,"\\Model_tests","\\Trained_",sps[[i]],sep="")
      dir.create(wd3)
      setwd(wd3)

      write.csv(data.frame("No checklists availabe"),
                "No checklists availabe.csv")   #save csv with AIC values
    }
  }

  #fit and test the models to the whole world and to each continent

  train_data <- dataSet(train_sps2)

  #fit the models with all sps

  models <- fitModels(train_data)

  #apply all models to train data

  m.preds_train <- lapply(models,function(x){predict(x,train_data, type="response")})

  #prepare test dataset

  test_sps2 <- rbindlist(test_sps)

  test_data <- dataSet(test_sps2)

  #apply all models to test data

  m.preds_test <- lapply(models,function(x){predict(x,test_data, type="response")})

  dir.create( paste(wd,"\\Model_tests\\Region",sep=""))
  wd4 <- paste(wd,"\\Model_tests\\Region\\Trained_Whole_world",sep="")
  dir.create(wd4) #create subfolder for results for the world
  setwd(wd4)

  write.csv(data.frame(AIC=sapply(models,AIC),
                       BIC=sapply(models,BIC),
                       mse_train=sapply(m.preds_train,function(x){mean((train_data$obs - x)^2,na.rm=T)}),
                       mse_test=sapply(m.preds_test,function(x){mean((test_data$obs - x)^2,na.rm=T)})),
            paste("Stats_whole_world.csv",sep=""))   #save csv with AIC values

  #boxplot with models' performance

  jpeg('Boxplot_Whole_world.jpg')

  par(mfrow=c(round(sqrt(length(m.preds_test))),ceiling(sqrt(length(m.preds_test)))),mar=c(2,2,2,2),bg="white")
  lapply(m.preds_test,function(x){boxplot(x~obs,test_data)})
  dev.off()

  #create spatial objects of all points

  all_points_sp <- all_points
  coordinates(all_points_sp) <- ~decimalLongitude+decimalLatitude
  proj4string(all_points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  map <- getMap() #get world map

  a <- as.character(unique(map$continent))
  b <- a[-c(which(is.na(a)),4)]    #create list of continents

  map_cont <- lapply(b,function(x){map[which(map$continent==x),]})  #create polygons for each continent
  pts_cont <- lapply(map_cont,function(x){over(all_points_sp,x)})  #extract continent info for each point
  pts_cont2 <- lapply(pts_cont,function(x){all_points[which(!is.na(x$continent)),]})  #separate datasets per continents

  #apply the model for the whole world to each continent

  i=1
  for(i in 1:length(b))
  {
    test_data <- dataSet(pts_cont2[[i]])

    m.preds_test <- lapply(models,function(x){predict(x,test_data, type="response")})

    write.csv(data.frame(AIC=sapply(models,AIC),
                         BIC=sapply(models,BIC),
                         mse_train=sapply(m.preds_train,function(x){mean((train_data$obs - x)^2,na.rm=T)}),
                         mse_test=sapply(m.preds_test,function(x){mean((test_data$obs - x)^2,na.rm=T)})),
              paste("Stats_",b[i],".csv",sep=""))   #save csv with AIC values

    jpeg(paste("Boxplot_",b[i],".jpg",sep=""))

    par(mfrow=c(round(sqrt(length(m.preds_test))),ceiling(sqrt(length(m.preds_test)))),mar=c(2,2,2,2),bg="white")
    lapply(m.preds_test,function(x){boxplot(x~obs,test_data)})
    dev.off()
  }

  #fit the models with each continent train data and test for all and each continent

  i=1
  for(i in 1:length(pts_cont2))
  {
    #fit the models with each continent

    train_data <- dataSet(pts_cont2[[i]])

    models <- fitModels(train_data)

    #apply all models to train data

    m.preds_train <- lapply(models,function(x){predict(x,train_data, type="response")})

    #prepare test dataset (whole world)

    test_data <- dataSet(all_points)

    #apply models for each continent to complete test data

    m.preds_test <- lapply(models,function(x){predict(x,test_data, type="response")})

    wd4 <- paste(wd,"\\Model_tests\\Region\\Trained_",b[i],sep="")
    dir.create(wd4) #create subfolder for results for the world
    setwd(wd4)

    write.csv(data.frame(AIC=sapply(models,AIC),
                         BIC=sapply(models,BIC),
                         mse_train=sapply(m.preds_train,function(x){mean((train_data$obs - x)^2)}),
                         mse_test=sapply(m.preds_test,function(x){mean((test_data$obs - x)^2)})),
              paste("Stats_whole_world.csv",sep=""))   #save csv with AIC values

    jpeg(paste("Boxplot_whole_world.jpg",sep=""))

    par(mfrow=c(round(sqrt(length(m.preds_test))),ceiling(sqrt(length(m.preds_test)))),mar=c(2,2,2,2),bg="white")
    lapply(m.preds_test,function(x){boxplot(x~obs,test_data)})
    dev.off()

    #apply models trained with each continent to each other test dataset

    j=1
    for(j in 1:length(b))
    {
      test_data <- dataSet(pts_cont2[[j]])

      m.preds_test <- lapply(models,function(x){predict(x,test_data, type="response")})

      write.csv(data.frame(AIC=sapply(models,AIC),
                           BIC=sapply(models,BIC),
                           mse_train=sapply(m.preds_train,function(x){mean((train_data$obs - x)^2)}),
                           mse_test=sapply(m.preds_test,function(x){mean((test_data$obs - x)^2)})),
                paste("Stats_",b[j],".csv",sep=""))   #save csv with AIC values

      #boxplot with models' performance

      jpeg(paste("Boxplot_",b[j],".jpg",sep=""))

      par(mfrow=c(round(sqrt(length(m.preds_test))),ceiling(sqrt(length(m.preds_test)))),mar=c(2,2,2,2),bg="white")
      lapply(m.preds_test,function(x){boxplot(x~obs,test_data)})
      dev.off()
    }
  }
}

wd <- "C:\\Users\\ca13kute\\Documents\\GIT_projects\\GBIF_record_status_distance\\TERRESTRIAL_MAMMALS12"

#test the models for all the sizes

wd <- list()
setwd("C:\\Users\\ca13kute\\Documents\\GIT_projects\\Check_lists\\Raster_size")
a <- list.files()
for(i in 1:length(a))
{
  wd[[i]] <- paste(getwd(),a[[i]],"Tables",sep="\\")
}

lapply(wd,fitTest)

fitTest(wd)

setwd(wd[[1]])
setwd("C:\\Users\\ca13kute\\Documents\\GIT_projects\\Package\\Models_checklists")

i=1
for(i in 1:length(models))
{
  saveRDS(models[[i]],paste("Model",i,sep="_"))
}

sapply(models,AIC)
