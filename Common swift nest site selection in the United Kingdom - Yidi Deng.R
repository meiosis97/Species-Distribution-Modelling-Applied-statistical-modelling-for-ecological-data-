nests <- read.csv('Symbols_Map_data.csv')
setwd("C:/Users/yidid/Downloads")

#load data
load('data.RData')
plot(nests$Lon, nests$Lat)

blank <- theme(axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="none")
require(ggplot2)
library(sp)
library(raster)
library(rgdal)
require(elevatr)
require(tiff)
require(ggmap)
require(ncdf4)
require(ggfortify)
require(ggpubr)



##### scatter plot
ggplot(nests, aes(x = Lon, y = Lat)) +
  stat_density2d(aes(alpha =..level.., fill =..level..), geom = 'polygon', bins = 50)+
  scale_fill_gradientn(colours = rainbow(100)[1:20])+
  scale_alpha_continuous(range = c(0,0.3), guide = F)+
  geom_density2d(aes(alpha = 0.1),color = 'grey', size = 1, bins = 10)+
  geom_point(aes(alpha = 0.05), color = 'black', size = 0.5)
lat.range <- range(nests$Lat)
lon.range <- range(nests$Lon)



((lat.range[2] - lat.range[1] )/0.1) * ((lon.range[2] - lon.range[1] )/0.1)



##### sampling range of swift survey
a <- seq(lat.range[1],lat.range[2] , by = 0.1)
b <- seq(lon.range[1],lon.range[2] , by = 0.1)




##### generate bivariate-uniform RV
runif2d <- function(a,b,c,d){
  x <- runif(1,a,b)
  y <- runif(1,c,d)
  return(c(x,y))
}




##### sample pseudo absent points
df <- data.frame(x = NA, y = NA)
for(i in 1:(length(a) - 1)){
  for(j in 1:(length(b) - 1)){
    df <- rbind(df,runif2d(a[i],a[i+1], b[j], b[j+1]))
  }
}
df <- df[-1,]




##### get the map of the UK
uk <- getData('GADM', country='GBR', level = 2)  

coordinates(df) <- ~ y + x
proj4string(df) <- CRS(proj4string(uk))
is.uk <- !is.na(over(df, uk)[,2])




##### remove pseudo absent point that is not inside the UK, set projection
noise <- data.frame(df[is.uk,])
colnames(noise) <- c('x', 'y')
remove <- c()
for(i in 1:nrow(nests)){
  remove <- c(remove,which(sqrt((noise[,'x'] - nests[i,c('Lon')])^2 + (noise[,'y'] - nests[i,c('Lat')])^2) < 0.1))
}

df <- noise[-unique(remove),]
coordinates(df) <- ~x + y
proj4string(df) <- CRS(proj4string(uk))
nest <- nests[,c('Lon', 'Lat')]
coordinates(nest) <- ~Lon + Lat
proj4string(nest) <- CRS(proj4string(uk))




##### visulization figure 1
ggplot(nests, aes(x = Lon, y = Lat)) +
  stat_density2d(aes(alpha =..level.., fill =..level..), geom = 'polygon', bins = 50)+
  scale_fill_gradientn(colours = rainbow(100)[1:20])+
  scale_alpha_continuous(range = c(0,0.3), guide = F)+
  geom_density2d(aes(alpha = 0.1),color = 'grey', size = 1, bins = 10)+
  geom_point(aes(alpha = 0.05), color = 'black', size = 0.5)+
  geom_point(data =data.frame(df) ,aes(x,y, alpha = 0.05), size = 0.5, col='steelblue')



map <- get_stamenmap(c(lon.range[1],lat.range[1], lon.range[2], lat.range[2]), zoom = 7)
ggmap(map) + geom_point(data =data.frame(df) ,aes(x,y, alpha = 0.1), size = 0.5, col='steelblue')+
  stat_density2d(data = nests, aes(x = Lon, y = Lat, alpha =..level.., fill =..level..), geom = 'polygon', bins = 20)+
  scale_fill_gradientn(colours = rainbow(100)[1:20])+
  scale_alpha_continuous(range = c(0,0.25), guide = F) +
  geom_point(data = nests, aes(x = Lon, y = Lat, alpha = 0.05), color = 'black', size = 0.5)+
  geom_density2d(data = nests, aes(x = Lon, y = Lat, alpha = 0.1),color = 'grey', size = 1, bins = 10)+
  theme(legend.position = "none") 




##### extract raster data
bio10=raster('wc2.0_bio_30s_10.tif')
bio18=raster('wc2.0_bio_30s_18.tif')
f.bio10 = raster('cc85bi5010.tif')
f.bio18 = raster('cc85bi5018.tif')
hfp = raster('wildareas-v3-2009-human-footprint.tif')
land = raster('gpw_v4_land_water_area_rev11_landareakm_2pt5_min.tif')
water = raster('gpw_v4_land_water_area_rev11_waterareakm_2pt5_min.tif')

bio10 <- crop(bio10,extent(uk))
bio10 <- mask(bio10, uk)
bio18 <- crop(bio18,extent(uk))
bio18 <- mask(bio18, uk)
f.bio10 <- crop(f.bio10,extent(uk))
f.bio10 <- mask(f.bio10, uk)
f.bio18 <- crop(f.bio18,extent(uk))
f.bio18 <- mask(f.bio18, uk)
elevation <- get_elev_raster(uk, z = 6)
elevation <- crop(elevation,extent(uk))
elevation <- mask(elevation, uk)
try <- spTransform(uk, CRS(projection(hfp)))
hfp <- crop(hfp,extent(try))
hfp <- projectRaster(hfp, crs = CRS(projection(uk)))
hfp <- mask(hfp, uk)
values(hfp)[values(hfp) > 60] <- NA
elev.sd <- terrain(elevation)
land <- crop(land,extent(uk))
land <- mask(land, uk)
water <- crop(water,extent(uk))
water <- mask(water, uk)

ggmap(map) + geom_point(data =data.frame(df) ,aes(x,y, alpha = 0.1, col=extract(bio10, df)), size = 0.5)

ggplot() + geom_raster(data = data.frame(rasterToPoints(bio10)), aes(x,y, fill = wc2.0_bio_30s_10)) + theme_bw()+
  scale_fill_viridis_c()
ggmap(map) + geom_point(data = data.frame(rasterToPoints(bio10)), aes(x,y, col = wc2.0_bio_30s_10)) + theme_bw() +
  scale_colour_viridis_c('')

ggplot() + geom_raster(data = data.frame(rasterToPoints(bio18)), aes(x,y, fill = wc2.0_bio_30s_18)) + theme_bw() +
  scale_fill_viridis_c()
ggmap(map) + geom_point(data = data.frame(rasterToPoints(bio18)), aes(x,y, col = wc2.0_bio_30s_18)) + theme_bw() +
  scale_colour_viridis_c('')

ggplot() + geom_raster(data = data.frame(rasterToPoints(elevation)), aes(x,y, fill = layer)) + 
  scale_fill_viridis_c('Elevation')+ xlim(-10,3) + xlab('Longitude') + ylab('Latitude') + theme_classic()
ggmap(map) + geom_point(data = data.frame(rasterToPoints(elevation)), aes(x,y, col = layer)) + theme_bw() +
  scale_colour_viridis_c('')

ggplot() + geom_raster(data = data.frame(rasterToPoints(f.bio10)), aes(x,y, fill = cc85bi5010/10)) + theme_bw() +
  scale_fill_viridis_c()
ggmap(map) + geom_point(data = data.frame(rasterToPoints(f.bio10)), aes(x,y, col = cc85bi5010/10)) + theme_bw() +
  scale_colour_viridis_c('')

ggplot() + geom_raster(data = data.frame(rasterToPoints(f.bio18)), aes(x,y, fill = cc85bi5018)) + theme_bw() +
  scale_fill_viridis_c()
ggmap(map) + geom_point(data = data.frame(rasterToPoints(f.bio18)), aes(x,y, col = cc85bi5018)) + theme_bw() +
  scale_colour_viridis_c('')

ggplot() + geom_raster(data = data.frame(rasterToPoints(hfp)), aes(x,y, fill = wildareas.v3.2009.human.footprint)) +
  scale_fill_gradientn(colors = heat.col) + blank
ggmap(map) + geom_point(data = data.frame(rasterToPoints(hfp)), aes(x,y, col = wildareas.v3.2009.human.footprint), size = 0.5) + theme_bw() +
  scale_colour_viridis_c('')+
  geom_point(data = nests, aes(x = Lon, y = Lat, alpha = 0.05), color = 'tomato1', size = 0.5, shape = 88)
  

##### contrast plot now and future: figure 2
heat.col <- rev(rainbow(100))[c(30:50,80:100)]

contrast.bio10 <- extract(f.bio10, rasterToPoints(bio10)[,1:2])
contrast.bio10 <- contrast.bio10/10 - rasterToPoints(bio10)[,3]
ggplot() + geom_raster(data = data.frame(rasterToPoints(bio10)), aes(x,y, fill = contrast.bio10))  +
          scale_fill_gradientn('BIO10 Difference',colors = heat.col) + xlim(-10,3) + xlab('Longitude') + ylab('Latitude') +
  theme_classic()
          

heat.col <- rev(rainbow(100))[c(18:50,80:100)]

contrast.bio18 <- extract(f.bio18, rasterToPoints(bio18)[,1:2])
contrast.bio18 <- contrast.bio18 - rasterToPoints(bio18)[,3]
ggplot() + geom_raster(data = data.frame(rasterToPoints(bio18)), aes(x,y, fill = contrast.bio18))  +
  scale_fill_gradientn('Precipitation Difference',colors = rev(heat.col))+ xlim(-10,3) + xlab('Longitude') + ylab('Latitude') +
  theme_classic()



##### average temperature in midland england
mean(rasterToPoints(bio10)[rasterToPoints(bio10)[,2] < 53,3])
mean(rasterToPoints(f.bio10)[rasterToPoints(f.bio10)[,2] < 53,3])

mean(rasterToPoints(bio18)[rasterToPoints(bio18)[,2] < 53,3])
mean(rasterToPoints(f.bio18)[rasterToPoints(f.bio18)[,2] < 53,3])




###### dataset aggregation
data.noise <- data.frame(is.nest = FALSE, bio10 = extract(bio10, df), bio18 = extract(bio18, df), 
                   elev = extract(elevation, df), elev.sd = extract(elev.sd, df), land = extract(land, df),
                   water = extract(water, df), hfp = extract(hfp, df))

data.real <- data.frame(is.nest = TRUE, bio10 = extract(bio10, nest), bio18 = extract(bio18, nest), 
                   elev = extract(elevation, nest), elev.sd = extract(elev.sd, nest), land = extract(land, nest),
                   water = extract(water, nest), hfp = extract(hfp, nest))

data.full <- rbind(data.real, data.noise)
data.full <- data.full[-which(rowSums(is.na(data.full))>0),]

pca <- prcomp(scale(data.full[,-1]))

autoplot(pca, data = data.full, colour = 'is.nest', loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 5) + theme_classic()

data.full$is.nest <- as.factor(data.full$is.nest)


######
require(caret)
require(pROC)
require(psych)
require(mgcv)
require(randomForest)
require(e1071)
require(mda)
require(earth)
require(xgboost)
require(nnet)
require(tree)


###### data to be predicted
cv <- createFolds(data.full$is.nest, k = 10)

newpoints <- rasterToPoints(bio10)[,1:2]
newpoints <- data.frame(newpoints)
coordinates(newpoints) <- ~x + y
proj4string(newpoints) <- CRS(proj4string(uk))

new.data <- data.frame(data.frame(newpoints), bio10 = extract(bio10, newpoints), bio18 = extract(bio18, newpoints), 
                       elev = extract(elevation, newpoints), elev.sd = extract(elev.sd, newpoints), land = extract(land, newpoints),
                       water = extract(water, newpoints), hfp = extract(hfp, newpoints))
new.data <- new.data[-which(rowSums(is.na(new.data))>0),]

future.data <- data.frame(data.frame(newpoints), bio10 = raster::extract(f.bio10, newpoints), bio18 = raster::extract(f.bio18, newpoints), 
                          elev = raster::extract(elevation, newpoints), elev.sd = raster::extract(elev.sd, newpoints), land = raster::extract(land, newpoints),
                          water = raster::extract(water, newpoints), hfp = raster::extract(hfp, newpoints))
future.data <- future.data[-which(rowSums(is.na(future.data))>0),]

future.data[,'bio10'] <- future.data[,'bio10']/10


###### model evaluation
modEva <- function(data){
  use = data[,"use"]; prb = data[,"p"]
  roc1 = roc(use, prb , percent = T, auc = T, plot = T) # for AUC
  SST = coords(roc1, 'best', ret = c('spec', 'sens', 'threshold')) # for specificity, sensitivity, threshold
  ACC = (SST[1]*sum(use) + SST[2]*(nrow(data) - sum(use))) / nrow(data) # for accuracy
  names(ACC) = "accuracy(%)"
  AUC = roc1$auc # AUC
  names(AUC) = "AUC(%)"
  Pred <- ifelse(prb >= SST[3], 1, 0) # classification based on best threshold
  data = cbind(data, Pred)
  kappa = cohen.kappa(data[, c("use", "Pred")])$kappa
  names(kappa)="kappa"
  L = c(AUC, ACC, kappa , SST) ; L = round(L, 3)
  return(L)
}

###### glm #####
mod.glm <- glm(is.nest~., data.full, family = binomial())
predict.glm <- predict(mod.glm, new.data[,-(1:2)], type = 'response',se.fit = TRUE)


p.glm <- ggplot(data = new.data) + geom_raster(aes(x,y,fill = predict.glm, alpha = predict.glm))+
  scale_fill_gradientn('Probability, GLM',colors = heat.col) + xlim(-10,3) + xlab('Longitude') + ylab('Latitude')+
  blank + ggtitle('GLM')

ggplot(data = new.data) + geom_raster(aes(x,y,fill = log(predict.glm$se.fit)))+
  scale_fill_gradientn('Probability, GLM',colors = rev(heat.col)) + xlim(-10,3) + xlab('Longitude') + ylab('Latitude')+
  blank + ggtitle('GLM')

ass.result <- data.frame(matrix(rep(NA,60), ncol = 6))
for(i in 1:10){
  mod <- glm(is.nest~., data.full[-cv[[i]],], family = binomial())
  ass <- data.frame(use = as.numeric(data.full[cv[[i]],]$is.nest)-1, p = predict(mod, data.full[cv[[i]],], type = 'response'))
  ass.result[i,] <- modEva(ass)
}
colnames(ass.result) <- names(modEva(ass))
colMeans(ass.result)

imp <- varImp(mod.glm)/max(varImp(mod.glm))
round(imp$Overall *100,2)



##### gam #####
mod.gam<-  gam(is.nest~ s(bio10) + s(bio18) + s(elev) + s(elev.sd) + s(hfp)
               +s(land) + s(water), data.full, family = binomial())
predict.gam <- predict(mod.gam, new.data[,-(1:2)], type = 'response')
p.gam <- ggplot(data = new.data) + geom_raster(aes(x,y,fill = predict.gam, alpha = predict.gam))+
  scale_fill_gradientn('Probability, GAM',colors = heat.col) + xlim(-10,3) + xlab('Longitude') + ylab('Latitude') +
  blank + ggtitle('GAM')


ass.result <- data.frame(matrix(rep(NA,60), ncol = 6))
for(i in 1:10){
  mod <- gam(is.nest~ s(bio10) + s(bio18) + s(elev) + s(elev.sd) + s(hfp)
             +s(land) + s(water), data.full[-cv[[i]],], family = binomial())
  ass <- data.frame(use = as.numeric(data.full[cv[[i]],]$is.nest)-1, p = predict(mod, data.full[cv[[i]],], type = 'response'))
  ass.result[i,] <- modEva(ass)
}
colnames(ass.result) <- names(modEva(ass))
colMeans(ass.result)


imp <- varImp(mod.gam)/max(varImp(mod.gam))
round(imp$Overall *100,2)



##### RF #####
mod.rf <- randomForest(is.nest~., data.full, ntree = 1000)
predict.rf <- predict(mod.rf, new.data[,-(1:2)],type = 'prob')[,2]
p.rf <- ggplot(data = new.data) + geom_raster(aes(x,y,fill = predict.rf, alpha = predict.rf))+
  scale_fill_gradientn('Probability, RF',colors = heat.col) + xlim(-10,3) + xlab('Longitude') + ylab('Latitude') +
  blank + ggtitle('RF')

predict.rf.future <- predict(mod.rf, future.data[,-(1:2)],type = 'prob')[,2]

partials <- list()
for(i in 1:7){
  var <- names(data.full)[1+i]
  partial <- eval(parse(text = paste('partialPlot(mod.rf, data.full, x.var = ', var ,', which.class = TRUE)')))
  partials[[i]] <- ggplot() + geom_point(aes(x = partial$x, y = partial$y))
}
ggarrange(plotlist=partials)


##### Partial effect plot figure 7
require(plotly)
means <- colMeans(data.full[data.full$is.nest,-1])
range(data.real$bio18, na.rm = T)

bio10.int <- seq(12,20,length.out = 200)
bio18.int <- seq(130,356,length.out = 200)
hfp.int <- seq(1,50, length.out = 100)

int.df.1 <- data.frame(expand.grid(bio10.int,hfp.int), matrix(rep(means[c(-1,-7)],20000),nrow = 20000, byrow = T))
int.df.2 <- data.frame(expand.grid(bio18.int,hfp.int), matrix(rep(means[c(-2,-7)],20000),nrow = 20000, byrow = T))

colnames(int.df.1) <- c('bio10','hfp',names(means[c(-1,-7)]))
colnames(int.df.2) <- c('bio18','hfp',names(means[c(-2,-7)]))


pred.int.1 <- predict(mod.rf, int.df.1,type = 'prob')[,2]
plot_ly(type  = 'mesh3d',x = int.df.1$bio10, y = int.df.1$hfp, z= pred.int.1, color = pred.int.1,
        colors = colorRamp(heat.col), intensity = pred.int.1)%>%
layout(
  title = "Marginal effect: BIO10 and HFI",
  scene = list(
    xaxis = list(title = "temperature"),
    yaxis = list(title = "human footprint index"),
    zaxis = list(title = "probability")
  ))


pred.int.2 <- predict(mod.rf, int.df.2,type = 'prob')[,2]
plot_ly(type  = 'mesh3d',x = int.df.2$bio18, y = int.df.2$hfp, z= pred.int.2, color = pred.int.2,
        colors = colorRamp(heat.col), intensity = pred.int.2) %>%
        layout(
          title = "Marginal effect: BIO18 and HFI",
          scene = list(
            xaxis = list(title = "precipitation"),
            yaxis = list(title = "human footprint index"),
            zaxis = list(title = "probability")
          ))




ass.result <- data.frame(matrix(rep(NA,60), ncol = 6))
for(i in 1:10){
  mod <-  randomForest(is.nest~., data.full[-cv[[i]],], ntree = 1000)
  ass <- data.frame(use = as.numeric(data.full[cv[[i]],]$is.nest)-1, p = predict(mod, data.full[cv[[i]],], type = 'prob')[,2])
  ass.result[i,] <- modEva(ass)
}
colnames(ass.result) <- names(modEva(ass))
colMeans(ass.result)


imp <- varImp(mod.rf)/max(varImp(mod.rf))
round(imp$Overall *100,2)



###### figure 5
sum(include)
include <- runif(length(predict.rf)) > (1-predict.rf)^(1/15)
p1 <- ggmap(map) + geom_point(data = new.data[include,],aes(x,y, col = predict.rf[include,]), size = 0.5, alpha = 0.1, col = 'tomato2')+
   xlab('Longitude') + ylab('Latitude') + 
   blank + ggtitle('2000')
include.future <- runif(length(predict.rf.future)) > (1-predict.rf.future)^(1/15)
p2 <- ggmap(map) + geom_point(data = future.data[include.future,],aes(x,y, col = predict.rf.future[include,]), size = 0.5, alpha = 0.1, col = 'tomato2')+
  xlab('Longitude') + ylab('Latitude') +
  blank + ggtitle('2050')
p4 <- ggplot(data = future.data) + geom_raster(aes(x,y,fill = predict.rf.future, alpha = predict.rf.future))+
  scale_fill_gradientn('Probability, RF',colors = heat.col) + xlim(-10,3) + xlab('Longitude') + ylab('Latitude') +
  blank + ggtitle('2050')
p3 <- ggplot(data = new.data) + geom_raster(aes(x,y,fill = predict.rf, alpha = predict.rf))+
  scale_fill_gradientn('Probability, RF',colors = heat.col) + xlim(-10,3) + xlab('Longitude') + ylab('Latitude') +
  blank + ggtitle('2000')

ggarrange(p1,p2)
ggarrange(p3,p4)


##### nb #####
mod.nb <- naiveBayes(is.nest~., data.full)
predict.nb <- predict(mod.nb, new.data[,-(1:2)],type = 'raw')[,2]
p.nb <- ggplot(data = new.data) + geom_raster(aes(x,y,fill = predict.nb, alpha = predict.nb))+
  scale_fill_gradientn('Probability, RF',colors = heat.col) + xlim(-10,3) + xlab('Longitude') + ylab('Latitude') +
  blank + ggtitle('NBC')



ass.result <- data.frame(matrix(rep(NA,60), ncol = 6))
for(i in 1:10){
  mod <-  naiveBayes(is.nest~., data.full[-cv[[i]],])
  ass <- data.frame(use = as.numeric(data.full[cv[[i]],]$is.nest)-1, p = predict(mod, data.full[cv[[i]],], type = 'raw')[,2])
  ass.result[i,] <- modEva(ass)
}
colnames(ass.result) <- names(modEva(ass))
colMeans(ass.result)

imp <- varImp(train(is.nest~., data = data.full, method = 'nb', trControl = trainControl(method = "cv")))$importance[,1]
imp <- imp/sum(imp)
round(imp *100,2)



##### FDA ######
mod.fda <- fda(is.nest~., data.full)
predict.fda <- predict(mod.fda, new.data[,-(1:2)],type = 'posterior')[,2]
p.fda <- ggplot(data = new.data) + geom_raster(aes(x,y,fill = predict.fda, alpha = predict.fda))+
  scale_fill_gradientn('Probability, RF',colors = heat.col) + xlim(-10,3) + xlab('Longitude') + ylab('Latitude') +
  blank + ggtitle('FDA')


ass.result <- data.frame(matrix(rep(NA,60), ncol = 6))
for(i in 1:10){
  mod <-  fda(is.nest~., data.full[-cv[[i]],])
  ass <- data.frame(use = as.numeric(data.full[cv[[i]],]$is.nest)-1, p = predict(mod, data.full[cv[[i]],],type = 'posterior')[,2])
  ass.result[i,] <- modEva(ass)
}
colnames(ass.result) <- names(modEva(ass))
colMeans(ass.result)

varImp(train(is.nest~., data = data.full, method = 'bagFDA', trControl = trainControl(method = "cv")))
imp <- c(100.0000, 38.1226,11.6524,5.7870,2.4853,0.5043,0.0000)
imp <- imp/sum(imp)
round(imp *100,2)
imp



##### MARS ######
mod.mars <- earth(is.nest~., data = data.full, glm=list(family=binomial))
predict.mars <- abs(predict(mod.mars, new.data[,-(1:2)], type = 'response'))
p.mars <- ggplot(data = new.data) + geom_raster(aes(x,y,fill = predict.mars, alpha = predict.mars))+
  scale_fill_gradientn('Probability, RF',colors = heat.col) + xlim(-10,3) + xlab('Longitude') + ylab('Latitude') +
  blank + ggtitle('MARS')

predict.mars.future <- predict(mod.mars, future.data[,-(1:2)], type = 'response')


ass.result <- data.frame(matrix(rep(NA,60), ncol = 6))
for(i in 1:10){
  mod <-  earth(is.nest~., data = data.full[-cv[[i]],], glm=list(family=binomial))
  ass <- data.frame(use = as.numeric(data.full[cv[[i]],]$is.nest)-1, p = as.numeric(predict(mod, data.full[cv[[i]],],type = 'response')))
  ass.result[i,] <- modEva(ass)
}
colnames(ass.result) <- names(modEva(ass))
colMeans(ass.result)

imp <- evimp(mod.mars)[,'rss']
imp <- imp/sum(imp)
round(imp *100,2)


###### xgboost ######
mod.gbm <- xgb.DMatrix(data = as.matrix(data.full[,-1]), label = as.numeric(data.full[,1]) - 1, 
                   nthread = 2, nrounds = 1000, objective = "binary:logistic",  verbose = 0, early_stopping_rounds = 10)


predict.gbm <- predict(mod.gbm, as.matrix(new.data[,-(1:2)]))
p.gbm <- ggplot(data = new.data) + geom_raster(aes(x,y,fill = predict.gbm, alpha = predict.gbm))+
  scale_fill_gradientn('Probability, RF',colors = heat.col) + xlim(-10,3) + xlab('Longitude') + ylab('Latitude') +
  blank + ggtitle('GBM')


ass.result <- data.frame(matrix(rep(NA,60), ncol = 6))
for(i in 1:10){
  mod <-  xgboost(data = as.matrix(data.full[-cv[[i]],-1]), label = as.numeric(data.full[-cv[[i]],1]) - 1, nthread = 2, nrounds = 1000, objective = "binary:logistic", verbose = 0)
  ass <- data.frame(use = as.numeric(data.full[cv[[i]],]$is.nest)-1, p = predict(mod, as.matrix(data.full[cv[[i]],-1])))
  ass.result[i,] <- modEva(ass)
}


colnames(ass.result) <- names(modEva(ass))
colMeans(ass.result)

xgb.importance(model = mod.gbm)
imp <- xgb.importance(model = mod.gbm)[,'Gain']
imp <- imp/sum(imp)
round(imp *100,2)


##### ann ######
s <- scale(data.full[,-1])
center <- attributes(s)$`scaled:center`
scale <- attributes(s)$`scaled:scale`
s.new <- new.data
for(i in 1:7){
  s.new[,2+i]  <-  (s.new[,2+i] - center[i])/scale[i]
}

mod.ann <- nnet(data.full$is.nest~., data = s, size = 4)
predict.ann <- predict(mod.ann, s.new[,-(1:2)])
p.ann <- ggplot(data = new.data) + geom_raster(aes(x,y,fill = predict.ann, alpha = predict.ann))+
  scale_fill_gradientn('Probability, RF',colors = heat.col) + xlim(-10,3) + xlab('Longitude') + ylab('Latitude') +
  blank + ggtitle('ANN')

ass.result <- data.frame(matrix(rep(NA,60), ncol = 6))
for(i in 1:10){
  mod <-  nnet(data.full[-cv[[i]],]$is.nest~., data = s[-cv[[i]],],  size = 4)
  ass <- data.frame(use = as.numeric(data.full[cv[[i]],]$is.nest)-1, p = predict(mod, s[cv[[i]],]))
  ass.result[i,] <- modEva(ass)
}
colnames(ass.result) <- names(modEva(ass))
colMeans(ass.result)


imp <- varImp(mod.ann)
imp <- imp/sum(imp)
round(imp *100,2)


###### tree ######
mod.tree <- rpart(is.nest~., data = data.full)
predict.tree <- abs(predict(mod.tree, new.data[,-(1:2)])[,2])
p.tree <- ggplot(data = new.data) + geom_raster(aes(x,y,fill = predict.tree, alpha = predict.tree))+
  scale_fill_gradientn('Probability, RF',colors = heat.col) + xlim(-10,3) + xlab('Longitude') + ylab('Latitude') +
  blank + ggtitle('CDF')


ass.result <- data.frame(matrix(rep(NA,60), ncol = 6))
for(i in 1:10){
  mod <-  rpart(is.nest~., data = data.full[-cv[[i]],])
  ass <- data.frame(use = as.numeric(data.full[cv[[i]],]$is.nest)-1, p = as.numeric(predict(mod, data.full[cv[[i]],])[,2]))
  ass.result[i,] <- modEva(ass)
}
colnames(ass.result) <- names(modEva(ass))
colMeans(ass.result)

imp <- varImp(mod.tree)
imp <- imp/sum(imp)
round(imp *100,2)

require(rpart)
predict(rpart(is.nest~., data = data.full), new.data[,-(1:2)])
ggarrange(p.glm, p.gam, p.mars,p.nb,p.rf,p.fda,p.gbm, p.tree, p.ann)


