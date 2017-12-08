# Uses the spatial lines to extract a SLP gradiant along each transe


library(reshape2)
library(plyr)
library(dplyr)
library(raster)
library(rasterVis)

doGraphic <- F

load("SLPtransects.Rdata")    # a list of spatial lines it the start and end coords
load("SLPstack.Rdata")

# extract each raster layer from the stack and extract SLP along each track
# record the time and date
grad.lines <- extract(r.stack,w.lines,df=TRUE,cellnumbers=T)
temp <- split(grad.lines,f=grad.lines$ID)

# get the lat and lon from the second column of temp (cellnumbers)
#station.loc <- xyFromCell(r.stack[[1]][[1]], wave.trans[[1]][[1]][,'cell'])   # a matrix
station.loc <- lapply(temp, function(x) xyFromCell(r.stack[[1]], x[,'SLP.1']))

# for each element in the list reorder "date","lon","lat","St1"
reshapeDF <- function(a,b){
  #a <- cbind(a[,1],data.frame("station"=seq(1:nrow(a))),a[,2:ncol(a)]) # add a stations column
  u <- cbind(a[,1],b,a[,3:ncol(a)])
  names(u)[1:3] <- c("ID","lon","lat")      # reset the column name
  v  <- melt(u, id = c("ID","lon","lat"))
  colnames(v)[4] <- "date"
  #d  <- reshape(c, idvar = c("ID","date"), timevar = "station", direction = "wide")
  n.row <- nrow(v)
  date.fill <- rep(idx,n.row)
  c$date <- idx
  return(c)
}

# extract the data from the brick along the transect at each station
grad.station <- lapply(temp, reshapeDF)

# plot results
st <- 1:9         # number of transects
days <- 1:10      # number of consecutive days to plot

for (j in 1:length(st)){
  
  # save each item in the list an daily transect for the time series
  df <- data.frame(grad.station[[j]])
  
  if (doGraphic){
    n.st = length(grad.station[[j]][1,])-2  # number of stations along the transect
    
    for (i in 1:length(days)){
      
      filename1 <- paste0("Station_0",j,"_",grad.station[[j]][i,2],".png")
      
      filename2 <- paste0("/media/robert/KELP-HDD-Portable/SLP/NOAA/L4/NCEP/EOT/images/",filename1)
      
      png(filename = filename2, width = 5,height = 5, units = "in", res = 300)
      
      plot(1:n.st,grad.station[[j]][i,3:(n.st+2)],
           type='o',pch=19,col="blue",xlab="Stations",ylab="SLP (hPa)",
           ylim = c(992,1032),
           main=paste0("Station ",grad.station[[j]][1,1],"\n",as.character(grad.station[[j]][i,2])))
      dev.off
    }
  }
}



