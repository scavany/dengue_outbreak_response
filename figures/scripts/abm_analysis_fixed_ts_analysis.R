#======================================================================================
# Author: Sean Cavany
# project: Outbreak response
# Year: 2019
# 
# Code to analyze the output of the agent based model following space spraying
#
# Files needed:
#======================================================================================
library(stringr)
library(vioplot)
directory <- ".."
moz.outname <- "immature"
spr.outname <- "spray"
foi.outname <- "foi"
none.folder <- 'spray_none'
folder <- "spray_threshold_fixed"
suffixes <- c("/", "_2/")

n.sims <- 1000
n.batches <- length(suffixes)
n.cats <- 2
n.outputs <- 3
n.years <- 11
n.days <- 4018
n.burn.years <- 0
if (n.burn.years > 0) {
    burnin.period <- seq(1,n.burn.years*365)
} else {
    burnin.period <- NULL
}
drop.days <- seq(365*n.years+1, n.days)
drops <- c(burnin.period, drop.days)

## parms <- read.csv(paste0(folder, "/parameters.csv"))
## thresholds <- parms[,"outbreak_response_threshold"]
## foi.names <- colnames(read.csv(paste0(folder, "/foi/sim_00000_foi.csv")))
## inf.parms <- grep("Infectious", foi.names)
## moz.parms <- grep("Moz", foi.names)
## inf.hum.parms <- setdiff(inf.parms, moz.parms)

## ##spraying strategies
## print(folder)
## setwd(directory)
## spr.files <- list.files(path=paste0(folder,"/",spr.outname),
##                         pattern=paste0(spr.outname, ".csv"),
##                         full.names=T)
## moz.files <- list.files(path=paste0(folder,"/",moz.outname),
##                         pattern=paste0(moz.outname, ".csv"),
##                         full.names=T)
## foi.files <- list.files(path=paste0(folder,"/",foi.outname),
##                         pattern=paste0(foi.outname, ".csv"),
##                         full.names=T)
## all.files <- data.frame(spr=spr.files, moz=moz.files,
##                         foi=foi.files)
## tot.output <- array(dim=c(length(all.files),
##                           length(foi.files),
##                           n.days-length(drops)))
## for (i in seq(nrow(all.files))) {
##     print(i)
##     ##spr output
##     spr.out <- read.csv(as.character(all.files$spr[i]))
##     tot.output[1,i,] <- (spr.out[-drops,2] > 0)
##     ##moz output
##     moz.out <- read.csv(as.character(all.files$moz[i]))
##     tot.output[2,i,] <- moz.out[-drops,"Adults"]
##     ##foi output
##     foi.out <- read.csv(as.character(all.files$foi[i]))
##     tot.output[3, i,]  <- apply(foi.out[-drops,inf.hum.parms], MARGIN=1, FUN=sum)
## }
    

##save output
setwd(directory)
#save(thresholds, tot.output, file="fixed_ts_analysis_output.RData")

load("scripts/fixed_ts_analysis_output.RData")

divisions <- 9

monthly.max <- 1000
weekly.max <- 1000*12/52

monthly.increment <- 1000/divisions
weekly.increment <- monthly.increment*12/52

monthly.lims <- seq(0,monthly.max, monthly.increment)
weekly.lims <- seq(0,weekly.max, weekly.increment)

monthly.lims <- c(0,400, 400, 1000)
weekly.lims <- c(0, 130, 130, 250) 

monthly.stratified <- list()
weekly.stratified <- list()

for (i in seq(length(monthly.lims)-1)) {
    monthly.stratified[[i]] <- tot.output[,which(thresholds[1:1000]>monthly.lims[i]
                                                 & thresholds[1:1000]<monthly.lims[i+1]),]
    weekly.stratified[[i]] <- tot.output[,1000+which(thresholds[1001:2000]>weekly.lims[i]
                                                     & thresholds[1001:2000]<weekly.lims[i+1]),]
}

monthly.sprays <- sapply(monthly.stratified, function(x) apply(x[1,,], MARGIN=2, FUN=sum))
monthly.med.cases <- sapply(monthly.stratified, function(x) apply(x[3,,], MARGIN=2, FUN=median))
monthly.upp.cases <- sapply(monthly.stratified, function(x) apply(x[3,,], MARGIN=2, FUN=function(y) quantile(y,0.75)))
monthly.low.cases <- sapply(monthly.stratified, function(x) apply(x[3,,], MARGIN=2, FUN=function(y) quantile(y,0.25)))
monthly.med.mozzes <- sapply(monthly.stratified, function(x) apply(x[2,,], MARGIN=2, FUN=median))
monthly.upp.mozzes <- sapply(monthly.stratified, function(x) apply(x[2,,], MARGIN=2, FUN=function(y) quantile(y,0.75)))
monthly.low.mozzes <- sapply(weekly.stratified, function(x) apply(x[2,,], MARGIN=2, FUN=function(y) quantile(y,0.25)))

weekly.sprays <- sapply(weekly.stratified, function(x) apply(x[1,,], MARGIN=2, FUN=sum))
weekly.med.cases <- sapply(weekly.stratified, function(x) apply(x[3,,], MARGIN=2, FUN=median))
weekly.upp.cases <- sapply(weekly.stratified, function(x) apply(x[3,,], MARGIN=2, FUN=function(y) quantile(y,0.75)))
weekly.low.cases <- sapply(weekly.stratified, function(x) apply(x[3,,], MARGIN=2, FUN=function(y) quantile(y,0.25)))
weekly.med.mozzes <- sapply(weekly.stratified, function(x) apply(x[2,,], MARGIN=2, FUN=median))
weekly.upp.mozzes <- sapply(weekly.stratified, function(x) apply(x[2,,], MARGIN=2, FUN=function(y) quantile(y,0.75)))
weekly.low.mozzes <- sapply(weekly.stratified, function(x) apply(x[2,,], MARGIN=2, FUN=function(y) quantile(y,0.25)))

library(viridis)
cols <- viridis(5)
cols.alpha <- viridis(5, alpha=0.4)
y.up <- max(max(weekly.upp.cases[,c(1,3)]), max(monthly.upp.cases[,c(1,3)]))

tiff("figures/Fig7.tif", res=600,
     width=4152, height=4152, compression="lzw")
par(mfrow=c(2,2), mar=c(1.1,4.1,1.1,1.1), ps=10,
    oma=c(4,0,3,3), xpd=NA)
for (i in c(1,3)) {
    plot(2000+seq(4015)/365, monthly.sprays[,i]/dim(monthly.stratified[[i]])[2],
         type='l', col=cols.alpha[1], xlab=NA, ylab=NA, axes=F)
    axis(side=4, col.axis=cols[1])
    par(new=T)
    plot(2000+seq(4015)/365, monthly.med.cases[,i], col=cols[3],
         ylim=c(0,y.up),xlab=NA, ylab="Infections", type='l',
         lwd=1, col.lab=cols[3], axes=F)
    axis(side=1)
    axis(side=2, col.lab=cols[3], col.axis=cols[3])
    polygon(c(2000+seq(4015)/365, rev(2000+seq(4015)/365)),
            c(monthly.low.cases[,i], rev(monthly.upp.cases[,i])),
            col=cols.alpha[3], border=NA)
    if (i==1) {
        title(main="Monthly")
    } else if (i==length(monthly.stratified)) {
        title(xlab="Year")
    }
    plot(2000+seq(4015)/365, weekly.sprays[,i]/dim(weekly.stratified[[i]])[2],
         type='l', col=cols.alpha[1], xlab=NA, ylab=NA, axes=F)
    axis(side=4, col.axis=cols[1])
    mtext(side=4, line=2, "Proportion simulations sprayed",
          col=cols[1])
    par(new=T)
    plot(2000+seq(4015)/365, weekly.med.cases[,i], col=cols[3],
         ylim=c(0,y.up),xlab=NA, ylab=NA, type='l', lwd=1, axes=F)
    axis(side=1)
    axis(side=2, col.lab=cols[3], col.axis=cols[3])
    polygon(c(2000+seq(4015)/365, rev(2000+seq(4015)/365)),
            c(weekly.low.cases[,i], rev(weekly.upp.cases[,i])),
            col=cols.alpha[3], border=NA)
    if (i==1) {
        title(main="Weekly")
    }
    else if (i==length(monthly.stratified)) {
        title(xlab="Year")
    }
}
dev.off()
