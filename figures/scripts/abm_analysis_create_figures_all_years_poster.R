#======================================================================================
# Author: Sean Cavany
# project: Outbreak response
# Year: 2019
# 
# Code to analyze the output of the agent based model following space spraying
#
# Files needed:
#======================================================================================
library(scales)
library(stringr)
library(viridis)
library(vioplot)
directory <- ".."
moz.outname <- "immature"
spr.outname <- "spray"
foi.outname <- "foi"
none.folder <- 'spray_none_mean'
analysis.folders <- c("irs_threshold_mean",
                      "irs_yearly_all")
suffixes <- c("/", "_2/")


n.sims <- 200
n.batches <- length(suffixes)
threshold.cats <- 4
yearly.cats <- 12
twice.yearly.cats <- 12*11/2
n.cats <- c(threshold.cats, yearly.cats)
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
load("scripts/irs_analysis_paired_all_years.RData", verbose=T)

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

if(!file.exists("figures")){dir.create("figures")}

threshold.mean.list <- list()
for (i in seq(n.cats["irs_threshold_mean"])) {
    threshold.mean.list[[i]] <- tot.output[["irs_threshold_mean"]][,3,i]
}
months <- c("January", "February", "March", "April",
            "May", "June", "July", "August",
            "September", "October", "November", "December")
months.short <- c("J", "F", "M", "A", "M", "J",
                  "J", "A", "S", "O", "N", "D")
spray.yearly.list <- list()
for (i in seq(n.cats["irs_yearly_all"])) {
    spray.yearly.list[[months[i]]] <- tot.output[["irs_yearly_all"]][,3,i]
}
load("scripts/best_irs_timeseries_all_years.RData", verbose = T)

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

september.days <- 244+seq(0,n.years-1)*365
march.days <- 60+seq(0,n.years-1)*365
tvec <- 2000+seq(n.days)[-drops]/365

##moz timeseries ==========================
m.irs <- med.tot.ts
u.irs <- upp.tot.ts
l.irs <- low.tot.ts


##case timeseries ================================
tiff("figures/Fig9.tif", res=600,
     width=4152, height=3320, compression="lzw")
par(mfrow=c(2,1), ps=10,
    mar=c(4.1,4.1,2.1,1.1),
    mgp=c(2,1,0))
cols <- viridis(5)
cols.alpha <- viridis(5, alpha=0.3)
for (folder in analysis.folders) {
    plot(tvec,med.none.ts[,2], type='l', col=cols[1],
         ylim=c(0,1.1*max(med.none.ts[,2])), ylab="Infections",
         xlab="Year", cex.lab=1, cex.axis=1)
    CI.col <- cols.alpha[1]
    polygon(c(tvec, rev(tvec)),
            c(upp.none.ts[,2], rev(low.none.ts[,2])),
            col=CI.col, border=NA)
    lines(tvec,med.tot.ts[[folder]][,2], type='l', col=cols[3])
    CI.col <- cols.alpha[3]
    polygon(c(tvec, rev(tvec)),
            c(upp.tot.ts[[folder]][,2],
              rev(low.tot.ts[[folder]][,2])),
            col=CI.col, border=NA)
    if (folder == "irs_yearly_all") {
        abline(v=tvec[september.days], lty=2)
        mtext("(b) Once yearly strategy", side = 3, adj = -0.02, line = 1)
    } else {
        mtext("(a) Adaptive threshold strategy", side = 3, adj = -0.02, line = 1)
        legend("top", legend=c("Not sprayed", "Sprayed"),
               lty=rep("solid",2), col=cols[c(1,3)],
               cex=1)
    }
##    title(folder, cex.lab=1.5)
}
dev.off()


##summary figure
load("scripts/best_total_cases_ulv.RData")
btc.irs <- best.total.cases
tiff("figures/Fig8.tif", res=600,
     width=4152, height=4152, compression="lzw")
layout(matrix(c(1,2,3,3),2,2, byrow=T))
par(ps=10,
    mar=c(3.0, 2.7, 2.6, 2.3),
    mgp=c(2,1,0),
    oma=c(1,1,1,0.11))
vioplot(threshold.mean.list,
        names=c(expression(paste("1",sigma," monthly")),
                expression(paste("2",sigma," monthly")),
                expression(paste("1",sigma," weekly")),
                expression(paste("2",sigma," weekly"))),
        col=viridis(4)[3],
        ylab="",
        xlab="",
        cex.axis = 1,
        cex.lab = 1)
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

vioplot(c(btc.ulv, list("1"=btc.irs[[2]], "2"=btc.irs[[3]])),
        col=c(rev(viridis(4)), rev(viridis(4))[2:3]),
        names=c("\n None", "Adaptive\n threshold",
                "Once\n yearly", "Twice\n yearly",
                "Adaptive\n threshold", "Once\n yearly"),
        ylog=TRUE)
title(ylab = "Total infections", xlab = "Spraying strategy",
      cex.lab = 1)
mtext("(c) Comparison of best strategies", side = 3, adj = -0.1, line = 0.6)
text(3, log(55000), labels="ULV", cex=1.2)
text(5.5, log(55000), labels="TIRS", cex=1.2)
abline(v=4.5, lty="dashed")
abline(v=1.5, lty="dashed")
abline(h=0)
dev.off()
