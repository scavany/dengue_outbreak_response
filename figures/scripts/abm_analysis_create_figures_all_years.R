#======================================================================================
# Author: Sean Cavany
# project: Outbreak response
# Year: 2019
# 
# Code to analyze the output of the agent based model following space spraying
#
# Files needed:
#======================================================================================
library(lattice)
library(stringr)
library(vioplot)
library(viridis)
library(gridExtra)
library(fields)
directory <- ".."
moz.outname <- "immature"
spr.outname <- "spray"
foi.outname <- "foi"
none.folder <- 'spray_none_mean'
analysis.folders <- c("spray_threshold_mean",
                      "spray_yearly_all",
                      "spray_twice_yearly_all")
suffixes <- c("/", "_2/")


n.sims <- 200
n.batches <- length(suffixes)
threshold.cats <- 4
yearly.cats <- 12
twice.yearly.cats <- 12*11/2
n.cats <- c(threshold.cats, yearly.cats, twice.yearly.cats)
names(n.cats) <- analysis.folders
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


##analysis - load output
setwd(directory)
load("scripts/spray_analysis_paired_all_years.RData")

#obtain quantities of interest
med.none <- apply(none.output, MARGIN=c(2), FUN=median, na.rm=T)
low.none <- apply(none.output, MARGIN=c(2),
                  FUN=function(x){quantile(x, 0.25, na.rm=T)})
upp.none <- apply(none.output, MARGIN=c(2),
                  FUN=function(x){quantile(x, 0.75, na.rm=T)})
mean.none <- apply(none.output, MARGIN=c(2), FUN=mean, na.rm=T)
sd.none <- apply(none.output, MARGIN=c(2), FUN=sd, na.rm=T)
cv.none <- sd.none/mean.none

med.out <- lapply(tot.output, function(x){apply(x, MARGIN=c(2,3),
                                                FUN=median,
                                                na.rm=T)})
low.out <- lapply(tot.output,
                     function(x){apply(x,
                                       MARGIN=c(2,3),
                                       FUN=function(x){
                                           quantile(x,
                                                    probs=0.25,
                                                    na.rm=T)})}) 
upp.out <- lapply(tot.output,
                     function(x){apply(x,
                                       MARGIN=c(2,3),
                                       FUN=function(x){
                                           quantile(x,
                                                    0.75,
                                                    na.rm=T)})})
mean.out <- lapply(tot.output, function(x){apply(x,
                                                 MARGIN=c(2,3),
                                                 FUN=mean,
                                                 na.rm=T)})
sd.out <- lapply(tot.output, function(x){apply(x,
                                               MARGIN=c(2,3),
                                               FUN=sd, na.rm=T)})
cv.out <- Map('/',sd.out,mean.out)
## mean.out <- lapply(tot.output, function(x){apply(x,
##                                                  MARGIN=c(2,3),
##                                                  FUN=mean,
##                                                  na.rm=T)})
## sum.out <- lapply(tot.output, function(x){apply(x,
##                                                 MARGIN=c(2,3),
##                                                 FUN=sum,
##                                                 na.rm=T)})
## n.out <- Map('/',sum.output,mean.output)
## max.out <- lapply(tot.output, function(x){apply(x,
##                                                 MARGIN=c(2,3),
##                                                 FUN=max,
##                                                 na.rm=T)})
## min.out <- lapply(tot.output, function(x){apply(x,
##                                                 MARGIN=c(2,3),
##                                                 FUN=min,
##                                                 na.rm=T)})
## sd.out <- lapply(tot.output, function(x){apply(x,
##                                                MARGIN=c(2,3),
##                                                FUN=sd, na.rm=T)})
## se.out <- Map('/', sd.out, lapply(n.out, sqrt))

best.locs <- lapply(med.out, function(x) max.col(-x))
best.values <- lapply(med.out, function(x) apply(x,1,min))

case.cat <- 3
best.total.cases <- list(none.folder=none.output[,case.cat])
for (folder in analysis.folders) {
    temp.loc <- best.locs[[folder]][case.cat]
    best.total.cases[[folder]] <- tot.output[[folder]][,
                                                       case.cat,
                                                       temp.loc]
}

##setup plots ========================================
if(!file.exists("figures")){dir.create("figures")}

threshold.mean.list <- list()
for (i in seq(n.cats["spray_threshold_mean"])) {
    threshold.mean.list[[i]] <- tot.output[["spray_threshold_mean"]][,3,i]
}

##yearly violin plot ===================================
months <- c("January", "February", "March", "April",
            "May", "June", "July", "August",
            "September", "October", "November", "December")
months.short <- c("J", "F", "M", "A", "M", "J",
                  "J", "A", "S", "O", "N", "D")
spray.yearly.list <- list()
for (i in seq(n.cats["spray_yearly_all"])) {
    spray.yearly.list[[months[i]]] <- tot.output[["spray_yearly_all"]][,3,i]
}

##twice yearly violin plot ===================================
pairs.by.month <- matrix(data=NA, nrow=12, ncol=12)
for (i in seq(11)) {
    for (j in seq(i+1, 12)) {
        pairs.by.month[i,j] <- -0.5*(i^2-25*i+24)+(j-i)
        pairs.by.month[j,i] <- -0.5*(i^2-25*i+24)+(j-i)
    }
}
load("scripts/best_timeseries_all_years.RData", verbose = T)

## plot timeseries =======================
med.none.ts <- apply(none.ts, c(1,3), median)
low.none.ts <- apply(none.ts, c(1,3),
                     function(x) quantile(x, 0.05))
upp.none.ts <- apply(none.ts, c(1,3),
                     function(x) quantile(x, 0.95))

med.tot.ts <- lapply(tot.ts,
                     function(x) apply(x, c(1,3), median))
low.tot.ts <- lapply(tot.ts,
                     function(x) apply(x, c(1,3),
                                       function(y) quantile(y, 0.05)))
upp.tot.ts <- lapply(tot.ts,
                     function(x) apply(x, c(1,3),
                                       function(y) quantile(y, 0.95)))

none.cum.ts <- apply(none.ts, MARGIN=c(2,3),FUN=cumsum)
tot.cum.ts <- lapply(tot.ts, function(x) apply(x,c(2,3),cumsum))
med.none.cum.ts <- apply(none.cum.ts, c(1,3), median)
low.none.cum.ts <- apply(none.cum.ts, c(1,3),
                         function(x) quantile(x, 0.05))
upp.none.cum.ts <- apply(none.cum.ts, c(1,3),
                         function(x) quantile(x, 0.95))

med.tot.cum.ts <- lapply(tot.cum.ts,
                         function(x) apply(x, c(1,3), median))
low.tot.cum.ts <- lapply(tot.cum.ts,
                         function(x) apply(x, c(1,3),
                                           function(y) quantile(y, 0.05)))
upp.tot.cum.ts <- lapply(tot.cum.ts,
                         function(x) apply(x, c(1,3),
                                           function(y) quantile(y, 0.95)))

october.days <- 274+seq(0,n.years-1)*365
march.days <- 60+seq(0,n.years-1)*365
september.days <- 244+seq(0,n.years-1)*365
november.days <- 305+seq(0,n.years-1)*365
tvec <- 2000+seq(n.days)[-drops]/365

##moz timeseries ==========================
cols <- c(red="#AE1A19", black="#252525", grey="#B1B1A6",
          brown="#604B4A", lgrey="#D0D0D0", white="#FFFFFF")
dayvec <- seq(as.Date("2000-01-01"), as.Date("2011-01-01"), by="1 day")


##case timeseries ================================
tiff("figures/Fig5.tif", res=600, units="in",
     width=4152/600, height=4982/600, compression="lzw",
     pointsize=10)
par(mfrow=c(3,1), #ps=10,
    mar=c(4.1,4.1,2.1,1.1),
    mgp=c(2,1,0))
cols <- viridis(5)
cols.alpha <- viridis(5, alpha=c(0.3, 0.3))
for (folder in analysis.folders) {
    plot(tvec,med.none.ts[,2], type='l', col=cols[1],
         ylim=c(0,1.1*max(med.none.ts[,2])), ylab="Infections",
         xlab="Year")
    CI.col <- cols.alpha[1] 
    polygon(c(tvec, rev(tvec)),
            c(upp.none.ts[,2], rev(low.none.ts[,2])),
            col=CI.col, border=NA)
    lines(tvec,med.tot.ts[[folder]][,2], type='l', col=cols[3], lwd=1.5)
    CI.col <- cols.alpha[3] 
    polygon(c(tvec, rev(tvec)),
            c(upp.tot.ts[[folder]][,2],
              rev(low.tot.ts[[folder]][,2])),
            col=CI.col, border=NA)
    if (folder == "spray_yearly_all" ) {
        abline(v=tvec[september.days], lty=2)
        mtext("(b) Once yearly strategy", side = 3, adj = -0.02, line = 1)
    } else if (folder == "spray_twice_yearly_all") {
        abline(v=tvec[september.days], lty=2)
        abline(v=tvec[november.days], lty=2)
        mtext("(c) Twice yearly strategy", side = 3, adj = -0.02, line = 1)
    } else {
        mtext("(a) Adaptive threshold strategy", side = 3, adj = -0.02, line = 1)
        legend("top", legend=c("Not sprayed", "Sprayed"),
               lty=rep("solid",2), col=cols[c(1,3)],
               cex=1)
    }
}
dev.off()

m.ulv <- med.tot.ts
u.ulv <- upp.tot.ts
l.ulv <- low.tot.ts
save(m.ulv, u.ulv, l.ulv, file="ulv_timeseries.RData")

##surface of median cases =================================================
median.matrix <- matrix(NA, nrow=12, ncol=12)
for (i in seq(12)) {
    median.matrix[i,i] <- median(tot.output[["spray_yearly_all"]][,3,i])
    for (j in seq(12)[-i]) {
        median.matrix[i,j] <- median(tot.output[["spray_twice_yearly_all"]][,3,pairs.by.month[i,j]])
    }
}


##Summary figure ==========================================================
tiff("figures/Fig4.tif", res=600,
     width=4152, height=4152, compression="lzw")
layout(matrix(c(1,1,1,2,2,2,2,3,3,3,4,5,5,5),2,7, byrow=T))
par(mgp=c(2,1,0),
    #mfrow=c(2,2), ps=10,
    ps=10, mar=c(3.0, 2.7, 2.6, 2.3),
    oma=c(1,1,1,0.11))
vioplot(threshold.mean.list,
        names=c(expression(paste("1",sigma," monthly")),
                expression(paste("2",sigma," monthly")),
                expression(paste("1",sigma," weekly")),
                expression(paste("2",sigma," weekly"))),
        col=viridis(4)[3],
        ylab="",
        xlab="",
        asp=1)
title(ylab = "Total infections", xlab = "Spraying strategy",
      cex.lab = 1)
mtext("(a) Adaptive threshold strategies", side = 3, adj = -0.1, line = 1)
vioplot(spray.yearly.list,
        names=months.short,
        col=viridis(4)[2],
        ylab="",
        xlab="",
        cex.axis = 1,
        cex.lab = 1)
title(ylab = "Total infections", xlab = "Spraying strategy",
      cex.lab = 1)
mtext("(b) Once yearly strategies", side = 3, adj = -0.1, line = 1)
med.matrix.new <- median.matrix
for (i in seq(ncol(median.matrix))) {
    med.matrix.new[i,i] <- NA
}
image(z=med.matrix.new, x=seq(12), y=seq(12), xlab="Month", ylab="Month",
      col=colorRampPalette(c(viridis(4)[1], "white"))(100)[1:81],
      asp=1, axes=F)
axis(2, at=seq(12),labels=months.short)
axis(1, at=seq(12),labels=months.short)
mtext("(c) Twice yearly strategies", side = 3, adj = -0.1, line = 0.6)
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)
##    dev.new(width=1.75, height=5)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
}
color.bar(colorRampPalette(c(viridis(4)[1], "white"))(100)[1:81],
          min=min(med.matrix.new, na.rm=T), max=max(med.matrix.new, na.rm=T),
          nticks=length(seq(180000,440000, 20000)), ticks=seq(180000,440000, 20000))
vioplot(best.total.cases,
        names=c("None", "Adaptive\n threshold",
                "Once\n yearly", "Twice\n yearly"),
        col=rev(viridis(4)),
        ylab="",
        xlab="",
        cex.axis = 1,
        cex.lab = 1)
title(ylab = "Total infections", xlab = "Spraying strategy",
      cex.lab = 1)
mtext("(d) Optimum strategies by type", side = 3, adj = -0.1, line = 0.6)
dev.off()

btc.ulv <- best.total.cases
seasonal.ulv <- spray.yearly.list
save(btc.ulv, seasonal.ulv, file="best_total_cases_ulv.RData")
