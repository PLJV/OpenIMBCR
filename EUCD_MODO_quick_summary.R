s <- imbcrTableToShapefile(recursiveFindFile(root="/home/ktaylora/Incoming/",name="RawData_PLJV_IMBCR_20161024.csv"))

transects <- unique(s[grepl(s$birdcode,pattern="MODO|EUCD"),]$transectnum)
  transects <- as.vector(transects[grepl(transects,pattern="^TX-")])

counts <- matrix(0,nrow=length(transects),ncol=2)

for(i in 1:length(transects)){
  s_focal <- s[s$transectnum == transects[i],]
  counts[i,1] <- sum(s_focal[s_focal$birdcode == "EUCD",]$cl_count)
  counts[i,2] <- sum(s_focal[s_focal$birdcode == "MODO",]$cl_count)
}

counts <- data.frame(counts)
  counts <- cbind(trans=as.vector(transects),counts)
    names(counts) <- c("TRANSECT","EUCD","MODO")

outliers <- (counts$EUCD > quantile(counts$EUCD,p=0.99))+(counts$MODO > quantile(counts$MODO,p=0.99)) > 0 # check for outliers

if(sum(outliers)>0){
  cat(paste(" -- removing :",sum(outliers)," outliers"))
  counts <- counts[!outliers,]
}

plot_max <- max(counts[,c(2,3)])

plot(counts[,c(2,3)],xlim=c(0,plot_max),ylim=c(0,plot_max),pch=22,cex=1.1,col="white",main=paste(sep="","Texas Transects (N=",length(transects),")"))
  grid();grid();
    #points(runif(nrow(counts),min=0,max=plot_max),runif(nrow(counts),min=0,max=plot_max),col="DarkGrey")
      points(counts[,c(2,3)],xlim=c(0,plot_max),ylim=c(0,plot_max))
        

x <- counts$EUCD
y <- counts$MODO
m <- lm(y ~ poly(x,2))

p <- predict(m,newdata=data.frame(x=0:plot_max),interval='confidence',
                               level=0.99)

lines(y=p[,1],x=0:plot_max,col="red",lwd=2)
lines(y=p[,2],x=0:plot_max,col="pink",lwd=2)
lines(y=p[,3],x=0:plot_max,col="pink",lwd=2)

lines(x=0:plot_max,y=plot_max:0,lwd=1.1,lty=4)

s_out <- s[s$transectnum %in% counts$TRANSECT,]

out <- NA
transects <- as.vector(counts$TRANSECT) # different if we removed outliers

for(i in 1:length(transects)){
  focal <- s_out[s_out$transectnum == as.vector(counts$TRANSECT)[i],]
    focal <- focal[!duplicated(focal$point),]
      focal <- rgeos::gCentroid(focal)
      
  focal <- SpatialPointsDataFrame(focal,data=data.frame(EUCD=counts[i,"EUCD"],MODO=counts[i,"MODO"]))

  if(is.na(out)){
    out <- focal
  } else {
    out <- rbind(out,focal)
  }
}

writeOGR(out,"/tmp","tx_eucd_modo",driver="ESRI Shapefile",overwrite=T)
