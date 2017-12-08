## Requires the object r.stack from '1.0_getDataStack.R'
## Uses Empirical Orthogobnal Teleconnections

## r,stack is a raster stack of all the avialable data for a product

library(remote)

setwd("~/R_projects-CS/ame-temporalbabe/CS-ExtractFeatures/")
load("Output/stack_SLP.Rdata")

r.stack <- setZ(r.stack, idx)

## replace NA values
indNA <- which(is.na(r.stack@data@values))
# get number of pixels in each layer
no.pixel <- nrow(r.stack)*ncol(r.stack)
# replace NA with mean of day before and after
r.stack@data@values[indNA] <- rowMeans(cbind(r.stack@data@values[indNA-no.pixel],r.stack@data@values[indNA+no.pixel]))
rm(indNA,no.pixel)

############################################### Monthly EOT using 'remote'
# idxM is the months of the year
for (i in 1:12){

  sub.stack <- subset(r.stack, which(idxM==i))
  eot.stack <- eot(sub.stack,y=NULL,n=5,verbose=TRUE)
  
  # for (j in 1:5){
  #   filename.png <- paste0("Output/SLP_EOT",j,"_",month.name[i],".png")
  #   png(filename.png,width = 8,height = 8,units = 'in',res = 300)
  #   plot(eot.stack, y=j, show.bp=TRUE, arrange="wide",anomalies=TRUE)
  #   dev.off()
  # }

  filename1 <- paste0("Output/substack",month.name[i],".Rdata")
  filename2 <- paste0("Output/eotstack",month.name[i],".Rdata")
  saveRDS(sub.stack,filename1)
  saveRDS(eot.stack,filename2)
}
