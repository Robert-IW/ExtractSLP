## Requires the EOT stacks created by '2.0_getEOT.R' (stored on external)
## Plots the locations of the first and second EOF for each month
## Creates a trasect from the base location to the coast

library(raster)
library(rasterVis)
library(maptools)
library(ncdf4)

plotBasePoints <-  TRUE

setwd("/media/robert/Seagate Backup Plus Drive/SLP/NCEP/L4/EOT/")

# setup mapy
bathy <- raster("/home/robert/R/x86_64-pc-linux-gnu-library/3.3/bathy-GEBCO/gebco-southernAfrica_1min.nc")

locbb <- matrix(c(-30,-50,80,-10),nrow=2,ncol=2)
x <- extent(locbb)
y <- crop(bathy,x)      # crop bathy map
z <- y >= 0
depz <- c(-100, -500, -2000)

file.list <- list.files(path = "data/", pattern = "eot.*.Rdata")

eot.names <- month.name[1:12]
eot.list <- sapply(eot.names,function(x) NULL)

# get a list of the monthly base point coords for modes 1 & 2
for (i in 1:12){
  
  temp <- readRDS(paste0("data/",file.list[grep(month.name[i],file.list)]))
  eot.list[[i]] <- rbind(temp@modes$mode_01@coords_bp,temp@modes$mode_02@coords_bp,temp@modes$mode_03@coords_bp)
  rm(temp)
}

# plot the base points
mode01_bp <- do.call(rbind, lapply(eot.list,'[',1,))
mode02_bp <- do.call(rbind, lapply(eot.list,'[',2,))
mode03_bp <- do.call(rbind, lapply(eot.list,'[',3,))

random <- seq(-1,1,length.out = 12)meanPoint01 <- list(eot.list$June[3,],eot.list$August[3,])

mode01_txt <- mode01_bp
mode01_txt <- mode01_txt+random
mode02_txt <- mode02_bp
mode02_txt <- mode02_txt+random
mode03_txt <- mode03_bp
mode03_txt <- mode03_txt+random

mode01.pts <- SpatialPoints(mode01_bp)
mode02.pts <- SpatialPoints(mode02_bp)
mode03.pts <- SpatialPoints(mode03_bp)

if (plotBasePoints){

df.spatial <- SpatialPointsDataFrame(rbind(mode01.pts,mode02.pts,mode03.pts),
                                     data.frame("month"=rep(rownames(mode01_bp),3),
                                                "mode"=matrix(cbind(rep("mode_01",12),
                                                                    rep("mode_02",12),
                                                                    rep("mode_03",12)),36,1),
                                                "modeID"=matrix(cbind(rep("cornflowerblue",12),
                                                                     rep("darkorange1",12),
                                                                     rep("chartreuse4",12)),36,1),
                                                stringsAsFactors = FALSE))


myColors <- c('cornflowerblue', 'darkorange1', 'chartreuse4')
myKey <- list(text=list(lab=c("1st mode","2nd mode","3rd mode")),
              rectangles=list(col = myColors),
              space='inside', padding.text = 5,
              size= 3, height = 0.75,
              rows=4)

# plot EOT base points
png("images/EOTbasepoints.png",width = 8.5,height = 5.8,units = 'in',res = 300)

levelplot(z, margin = FALSE, colorkey = FALSE, key = myKey,col.regions = c("white","gray"),
          main=list("EOT Modes 1-3 Base Points",cex=2.4),
          scales=list(x=list(cex=1.2),y=list(cex=1.2)),
          xlab=list("Latitude",cex=1.8),ylab=list("Longitude",cex=1.8)) +
  contourplot(y,at=c(0,-200,-500,-2000), margin = F, labels = F) +
  layer(sp.points(df.spatial,col=df.spatial$modeID,fill=TRUE,pch=19,cex=1.5)) +
  layer(sp.pointLabel(df.spatial, df.spatial$month, cex=1.2))
dev.off()
}

###############################################################################################################3
# Get the groups of Base Points and their mean locations
groupPoint <- list()
groupPoint[[1]] <- rbind(eot.list$June[3,],eot.list$August[3,])
groupPoint[[2]] <- rbind(eot.list$October[1,],eot.list$February[1,],eot.list$November[1,],eot.list$May[2,],eot.list$January[1,])
groupPoint[[3]] <- rbind(eot.list$March[1,],eot.list$September[2,],eot.list$December[2,],eot.list$April[2,])
groupPoint[[4]] <- rbind(eot.list$April[2,],eot.list$June[2,],eot.list$July[2,],eot.list$August[2,])
groupPoint[[5]] <- rbind(eot.list$January[3,],eot.list$May[3,])
groupPoint[[6]] <- rbind(eot.list$November[2,],eot.list$February[2,])
groupPoint[[7]] <- rbind(eot.list$September[1,],eot.list$December[1,],eot.list$October[2,])
groupPoint[[8]] <- rbind(eot.list$June[1,],eot.list$July[1,],eot.list$January[2,])
groupPoint[[9]] <- rbind(eot.list$April[1,],eot.list$May[1,],eot.list$March[2,])

meanPoint <- lapply(groupPoint, function(x) colMeans(x))

w.base <- (do.call("rbind",meanPoint))
w.base.pts <- SpatialPoints(w.base)

# get the extension points for the transect  
u <- rasterToPoints(z, function(x) x == TRUE)
v <- do.call("rbind",meanPoint)

# function to return the nearest land coords to a point
getLand <- function(pnt, landpnts){
  temp <- spDists(landpnts[,1:2], matrix(pnt,nrow=1,ncol=2))
  coords <- landpnts[which(temp == min(temp)),1:2]
  # make sure the points are within SA
  i <- 0
  while (coords[1]<10 | coords[1] > 35 | coords[2] > -18 | coords[2] < -36){
    i <- i+1
    n <- length(unique(temp))
    coords <- landpnts[which(temp == sort(unique(temp),decreasing=TRUE)[n-i]),1:2]
  }
  return(coords)
}

w.ext <-  apply(v,1,getLand,u)
w.ext <- t(w.ext)
# point #9 is on land so needs to be extended seaward
w.ext[,9] <- c(41,-27)
w.ext.pts <- SpatialPoints(w.ext)

w.lines <- vector("list", length(w.ext.pts))
for (i in seq_along(w.lines)) {
  w.lines[[i]] <- Lines(list(Line(rbind(w.base.pts[i, ], w.ext.pts[i,]))), as.character(i))
}
w.lines <- SpatialLines(w.lines)

# plot the base and extension points
locbb <- matrix(c(-15,-50,60,-15),nrow=2,ncol=2)
x <- extent(locbb)
y <- crop(bathy,x)      # crop bathy map
z <- y >= 0

png("images/EOTtransects.png",width = 8.5,height = 5.8,units = 'in',res = 300)

levelplot(z, margin=FALSE, main=list("Locations of SLP Trasects",cex=2.2),
     colorkey = FALSE, col.regions = c("white","gray"),
     xlab=list("Longitude",cex=1.8),ylab=list("Latitude",cex=1.8),
     scales=list(x=list(cex=1.2),y=list(cex=1.2)),
     xlim = c(-10,60),ylim = c(-50,-15)) +
  contourplot(y,at=c(0,-200,-500,-2000), margin = F, labels = F) +
  layer(sp.points(w.base.pts, pch=19, col="red",cex=1.4)) +
  layer(sp.points(w.ext.pts, pch=19, col="blue",cex=1.4)) +
  layer(sp.lines(w.lines,lty=2))
dev.off()
print("Image saved")

save(w.lines,file="~/R_projects-CS/ame-temporalbabe/CS-ExtractFeatures/Output/SLPtransects.Rdata")




  
