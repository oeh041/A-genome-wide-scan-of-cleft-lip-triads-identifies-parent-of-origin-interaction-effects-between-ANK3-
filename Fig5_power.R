## Power calculations in Haplin
if(T){
  # setwd("W:/Jobb/Skriblerier/Postdoc/Haplin/GxE/")
  setwd("D:/Haplin/Artikkel_POOxE_CL")
  rm(list=ls(all=T))
  gc()
  cat("\014")
  ls(all=T)
}

## Installing Haplin
install.packages("Haplin", dependencies = T)
## loading new version
library("Haplin")
## Needs to be verson 7 or later
packageVersion("Haplin")

## Testing that all is well
hapPowerAsymp(n.strata = 2,                             ## smoke, alcohol, vitamin
              cases = list(c(mfc=1500), c(mfc=300)),    ## CL/P
              haplo.freq = list(c(.8,.2)),              ## MAF=0.3
              RRcm = list(c(1,2),c(1,1)),               ## RR=2 for mom in strata 1
              RRcf = list(c(1,1),c(1,2)),               ## RR=0.5 for dad in strata 1
              alpha = 0.05)

## In general
## MAF=0.2
## RRR from 1 to 3
## 1400 cases and 300 controls

## Panel A
## Vary number of trios from 500 unexposed and 300 exposed to 1400 unexposed and 600 controls
nexp <- seq(from = 500, to = 1400, by = 300);nexp
nexp2 <- seq(from = 300, to = 600, by = 100);nexp
RR <- seq(from = 1, to = 3, length.out = 9);RR
maf <- .2
RRRn <- matrix(0, nrow = length(nexp), ncol = length(RR))
rownames(RRRn) <- nexp
colnames(RRRn) <- RR
RRRn
for(i in seq_along(nexp)){
  for(j in seq_along(RR)){
    RRRn[i,j] <- as.numeric(hapPowerAsymp(n.strata = 2,                      ## smoke, alcohol, vitamin
                               cases = list(c(mfc=nexp[i]), c(mfc=nexp2[i])),     ## CL/P
                               haplo.freq = list(c(1-maf,maf)),              ## MAF=0.3
                               RRcm = list(c(1,1),c(1,RR[j])),               ## RR=2 for mom in strata 1
                               RRcf = list(c(1,1),c(1,1)),                   ## RR=1 for dad in strata 1
                               alpha = 0.05)$haplo.power[2,"RRcm_cf.power"])
  }
}
## Panel B
## Vary MAF from 0.1 to 0.4 (steps of 0.1)
maf <- seq(from = .1, to = .4, by = .1);maf
RRRmaf <- matrix(0, nrow = length(maf), ncol = length(RR))
rownames(RRRmaf) <- maf
colnames(RRRmaf) <- RR
RRRmaf
for(i in seq_along(maf)){
  for(j in seq_along(RR)){
    RRRmaf[i,j] <- as.numeric(hapPowerAsymp(n.strata = 2,                      ## smoke, alcohol, vitamin
                                 cases = list(c(mfc=1100), c(mfc=500)),        ## CL/P
                                 haplo.freq = list(c(1-maf[i],maf[i])),        ## MAF=0.3
                                 RRcm = list(c(1,1),c(1,RR[j])),               ## RR=2 for mom in strata 1
                                 RRcf = list(c(1,1),c(1,1)),                   ## RR=1 for dad in strata 1
                                 alpha = 0.05)$haplo.power[2,"RRcm_cf.power"])
  }
}

## Figure
# windows(NA,width = 20, height = 10)
xlim <- range(RR)
ylim <- 0:1
ltyn <- c(4,3,1,2)
coln <- c(4,3,1,2)
pchn <- c(4,3,1,2)
ltymaf <- c(5,1,6,7)
colmaf <- c(5,1,6,7)
pchmaf <- c(5,1,6,7)
lwd <- 3
cex <- 2.5
tiff(filename = "FigPower.tiff",
     width = 20, height = 10, units = "in", res = 450)
par(mfrow = c(1,2), mar = c(6, 6, 4.1, 2.1))
## Ns
plot(NA, main = "Varying N", xlab = "RRR", ylab = "Power", xlim = xlim, ylim = ylim, cex.lab = cex, cex.axis = cex, cex.main = cex)
abline(h = seq(from = 0, to = 1, by = .1), lty = 3, col = "grey")
for(i in 1:nrow(RRRn)){
  lines(RR, RRRn[i,], col = coln[i], lty = ltyn[i], lwd = lwd)
  points(RR, RRRn[i,], col = coln[i], pch = pchn[i], cex = cex, lwd = lwd)
}
legend("bottomright", legend = paste(nexp,"-",nexp2,sep=""), lty = ltyn, col = coln, pch = pchn, bg = "white", cex = cex, lwd = lwd)
## MAFs
plot(NA, main = "Varying MAF", xlab = "RRR", ylab = "Power", xlim = xlim, ylim = ylim, cex.lab = cex, cex.axis = cex, cex.main = cex)
abline(h = seq(from = 0, to = 1, by = .1), lty = 3, col = "grey")
for(i in 1:nrow(RRRmaf)){
  lines(RR, RRRmaf[i,], col = colmaf[i], lty = ltymaf[i], lwd = lwd)
  points(RR, RRRmaf[i,], col = colmaf[i], pch = pchmaf[i], cex = cex, lwd = lwd)
}
legend("bottomright", legend = paste("MAF=",maf,sep=""), lty = ltymaf, col = colmaf, pch = pchmaf, bg = "white", cex = cex, lwd = lwd)
dev.off()
shell("start FigPower.tiff")