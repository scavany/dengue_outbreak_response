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
library(viridis)
directory <- ".."
moz.outname <- "immature"
spr.outname <- "spray"
foi.outname <- "foi"
folder <- "irs_threshold_fixed"
suffixes <- c("/", "_2/")

n.sims <- 200
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


##analysis - load output
setwd(directory)
load("scripts/irs_analysis_fixed.RData", verbose=T)
parms <- read.csv("scripts/irs_fixed_parms.csv")

##obtain quantities of interest
med.none <- apply(none.output, MARGIN=c(2), FUN=median, na.rm=T)
low.none <- apply(none.output, MARGIN=c(2),
                  FUN=function(x){quantile(x, 0.25, na.rm=T)})
upp.none <- apply(none.output, MARGIN=c(2),
                  FUN=function(x){quantile(x, 0.75, na.rm=T)})
mean.none <- apply(none.output, MARGIN=c(2), FUN=mean, na.rm=T)
sd.none <- apply(none.output, MARGIN=c(2), FUN=sd, na.rm=T)
cv.none <- sd.none/mean.none

med.out <- apply(tot.output, MARGIN=c(2,3),
                 FUN=median,
                 na.rm=T)
low.out <- apply(tot.output, MARGIN=c(2,3),
                 FUN=function(x){quantile(x,
                                          probs=0.25,
                                          na.rm=T)})
upp.out <- apply(tot.output, MARGIN=c(2,3),
                 FUN=function(x){quantile(x,
                                          probs=0.75,
                                          na.rm=T)}) 
mean.out <- apply(tot.output, MARGIN=c(2,3),
                  FUN=mean, na.rm=T)
sd.out <- apply(tot.output, MARGIN=c(2,3), FUN=sd, na.rm=T)
cv.out <- sd.out/mean.out


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

library(mgcv)
output.monthly <- as.data.frame(cbind(threshold=parms$outbreak_response_threshold[1:1000],
                                      sprays=tot.output[,1,1],
                                      mozzes=tot.output[,2,1],
                                      cases=tot.output[,3,1]))
output.weekly <- as.data.frame(cbind(threshold=parms$outbreak_response_threshold[1001:2000],
                                      sprays=tot.output[,1,2],
                                      mozzes=tot.output[,2,2],
                                      cases=tot.output[,3,2]))

## plot(sprays~threshold, data=output.monthly,
##      type='p', xlab=NA, ylab="number of days sprayed")
## plot(sprays~threshold, data=output.weekly,
##      type='p', xlab=NA, ylab=NA)
## plot(mozzes~threshold, data=output.monthly,
##      type='p', xlab=NA, ylab="average number mosquitoes")
## plot(mozzes~threshold, data=output.weekly,
##      type='p', xlab=NA, ylab=NA)
## plot(cases~threshold, data=output.monthly,
##      type='p', xlab="threshold (cases/month)", ylab="total cases")
## plot(cases~threshold, data=output.weekly,
##      type='p', xlab="threshold (cases/week)", ylab=NA)

##gam
axis.size <- 1
label.size <- 1
spray.ylim <- c(min(c(output.monthly$sprays, output.weekly$sprays)),
                max(c(output.monthly$sprays, output.weekly$sprays)))
moz.ylim <- c(min(c(output.monthly$mozzes, output.weekly$mozzes)),
              max(c(output.monthly$mozzes, output.weekly$mozzes)))
case.ylim <- c(min(c(output.monthly$cases, output.weekly$cases)),
                max(c(output.monthly$cases, output.weekly$cases)))

load("scripts/ulv_gam_outputs.RData", verbose=T)

sp.ylim <- c(min(c(spray.ylim, s.ylim)),
             max(c(spray.ylim, s.ylim)))
mo.ylim <- c(min(c(moz.ylim, m.ylim)),
             max(c(moz.ylim, m.ylim)))
ca.ylim <- c(min(c(case.ylim, c.ylim)),
             max(c(case.ylim, c.ylim)))

tiff("figures/Fig6.tif", res=600,
     width=4152, height=4152, compression="lzw")
par(mfrow=c(3,2), ps=10,
    mar=c(2.1, 2.1, 1.1, 1.1),
    oma=c(3,2,3,1), xpd=NA)
plot(gms.ulv, residuals=T, rug=F, shift=coef(gms.ulv)[1],
     se=F, xlab=NA, ylab="# days sprayed", cex.axis=axis.size, cex.lab=label.size,
     ylim=sp.ylim)
title(main="ULV")

gam.monthly.sprays <- gam(sprays~s(threshold), data=output.monthly)
gam.weekly.sprays <- gam(sprays~s(threshold), data=output.weekly)
plot(gam.monthly.sprays, residuals=T, rug=F, shift=coef(gam.monthly.sprays)[1],
     se=F, xlab=NA, ylab=NA, cex.axis=axis.size, cex.lab=label.size,
     ylim=sp.ylim)
title(main="TIRS")

plot(gmm.ulv, residuals=T, rug=F, shift=coef(gmm.ulv)[1],
     se=F, xlab=NA, ylab="average mosquitoes", cex.axis=axis.size, cex.lab=label.size,
     ylim=mo.ylim)

gam.monthly.mozzes <- gam(mozzes~s(threshold), data=output.monthly)
gam.weekly.mozzes <- gam(mozzes~s(threshold), data=output.weekly)
plot(gam.monthly.mozzes, residuals=T, rug=F, shift=coef(gam.monthly.mozzes)[1],
     se=F, xlab=NA, ylab=NA, cex.axis=axis.size, cex.lab=label.size,
     ylim=mo.ylim)

plot(gmc.ulv, residuals=T, rug=F, shift=coef(gmc.ulv)[1],
     se=F, xlab="threshold (cases/month)", ylab="total infections", cex.axis=axis.size, cex.lab=label.size,
     ylim=ca.ylim)

gam.monthly.cases <- gam(cases~s(threshold), data=output.monthly)
gam.weekly.cases <- gam(cases~s(threshold), data=output.weekly)
plot(gam.monthly.cases, residuals=T, rug=F, shift=coef(gam.monthly.cases)[1],
     se=F, xlab="threshold (cases/month)", ylab=NA, cex.axis=axis.size, cex.lab=label.size,
     ylim=ca.ylim)

dev.off()
