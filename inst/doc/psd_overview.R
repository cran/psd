
## ----eval=TRUE, echo=FALSE-----------------------------------------------
par(las=1)


## ----eval=TRUE, echo=TRUE, label="Load library."-------------------------
library(psd)


## ----eval=TRUE, eval=TRUE, label="Load Project MAGNET data."-------------
data(magnet)


## ----eval=TRUE, eval=TRUE, label="Show contents of Project MAGNET."------
names(magnet)


## ----eval=TRUE, echo=TRUE, label="Outliers."-----------------------------
subset(magnet, abs(mdiff)>0)


## ----eval=TRUE, echo=TRUE, label=MAGPSDS---------------------------------
psdr <- pspectrum(magnet$raw)
psdc <- pspectrum(magnet$clean)


## ----echo=TRUE, eval=TRUE------------------------------------------------
psdc_recovered <- psd_envGet("final_psd")
all.equal(psdc, psdc_recovered)


## ----eval=TRUE, echo=TRUE, label=RAWvCLEAN-------------------------------
plot(psdc, log="dB", main="Raw and Clean Project MAGNET power spectral density", 
       lwd=3, ci.col=NA, ylim=c(0,32), yaxs="i")
#plot(psdc_ar, log="dB", add=TRUE, lwd=3, col="red")
plot(psdr, log="dB", add=TRUE, lwd=3, lty=5)
text(c(0.25,0.34), c(11,24), c("Clean","Raw"), cex=1)


## ----eval=FALSE, echo=TRUE, label="Naive spectrum estimation."-----------
## spec.pgram(X, pad=1, taper=0.2, detrend=FALSE, demean=FALSE, plot=F)


## ----eval=TRUE, echo=TRUE,  label=MAGNETNAIVE----------------------------
ntap <- psdc$taper
psdcore(magnet$clean, ntaper=ntap, refresh=TRUE, plotpsd=TRUE)


## ----eval=TRUE, echo=TRUE, label="Load RSEIS package."-------------------
require(RSEIS)
dt=1 # km
# prewhiten the data after adding a linear trend + offset
summary(prewhiten(mc <- (ts(magnet$clean+1e3,deltat=dt)+seq_along(magnet$clean)), plot=FALSE))


## ----eval=TRUE, echo=TRUE, label="AR prewhiten"--------------------------
summary(atsar <- prewhiten(mc, AR.max=100, plot=FALSE))
print(atsar$ardfit)
ats_lm <- atsar$prew_lm
ats_ar <- atsar$prew_ar


## ----eval=TRUE, echo=TRUE, fig.height=5, fig.width=5.5, label=ARFITPLT----
plot(ts.union(orig.plus.trend=mc, linear=ats_lm, ar=ats_ar), yax.flip=TRUE, 
     main=sprintf("Prewhitened Project MAGNET series"))
mtext(sprintf("linear and linear+AR(%s)", atsar$ardfit$order), line=1.1)


## ----eval=TRUE, echo=TRUE, label="Sampling rate versus interval."--------
a <- rnorm(32)
all.equal(psdcore(a,1)$spec, psdcore(a,-1)$spec)


## ----eval=TRUE, echo=TRUE, label="Compute PSD with mtapspec."------------
tapinit <- 10
Mspec <- mtapspec(ats_lm, deltat(ats_lm), MTP=list(kind=2, inorm=3, nwin=tapinit, npi=0))


## ----eval=TRUE, echo=TRUE------------------------------------------------
str(Mspec)


## ----eval=TRUE, echo=TRUE, label="Comparative spectra."------------------
Xspec <- spec.pgram(ats_lm, pad=1, taper=0.2, detr=TRUE, dem=TRUE, plot=FALSE)
Pspec <- psdcore(ats_lm, dt, tapinit)
Aspec <- pspectrum(ats_lm, dt, tapinit, plot=FALSE)
# Correct for double-sidedness of spectrum and mtapspec results
class(Mspec)
Mspec <- normalize(Mspec, dt, "spectrum")
nt <- 1:Mspec$numfreqs
mspec <- Mspec$spec[nt]
class(Xspec)
Xspec <- normalize(Xspec, dt, "spectrum")


## ----eval=TRUE, echo=TRUE, fig.width=6.0, fig.height=5.4, label=RSEIS----
require(RColorBrewer)
cols <- c("dark grey", brewer.pal(8, "Set1")[c(5:4,2)])
lwds <- c(1,2,2,5)
par(las=1)
plot(Xspec, log="dB", ylim=40*c(-0.4,1), ci.col=NA, 
       col=cols[1], lwd=lwds[1], main="PSD Comparisons") 
pltf <- Mspec$freq
lines(pltf, pltp <- dB(mspec), col=cols[2], lwd=lwds[2]) 
plot(Pspec, log="dB",  add=TRUE, col=cols[3], lwd=lwds[3]) 
plot(Aspec, log="dB", add=TRUE, col=cols[4], lwd=lwds[4]) 
legend("topright", 
  c("spec.pgram","RSEIS::mtapspec","psdcore","pspectrum"), 
  title="Estimator", lwd=3, cex=1.1, col=cols)


## ----eval=TRUE, echo=TRUE, label="Interpolate results."------------------
require(signal)
pltpi <- interp1(pltf, pltp, Pspec$freq)


## ----eval=TRUE, echo=TRUE, label="Summarize regression statistics."------
df <- data.frame(x=dB(Pspec$spec), y=pltpi, tap=unclass(Aspec$taper))
summary(dflm <- lm(y ~ x + 0, df))
df$res <- residuals(dflm)


## ----eval=TRUE, echo=TRUE, fig.width=6, fig.height=2.5, label=RSEISvsRLP2----
require(ggplot2)
gr <- ggplot(df, aes(x=x, y=res)) + geom_abline(intercept=0, slope=0, size=2, color="salmon") + 
geom_point(aes(color=tap))
print(gr + theme_bw() + 
ggtitle("Regression residuals, colored by optimized tapers")+
xlab("Power levels, dB") + ylab("")
)


## ----eval=TRUE, echo=TRUE, label=BSPEC-----------------------------------
require(bspec)
print(Bspec <- bspec(ts(magnet$clean)))


## ----eval=TRUE, echo=TRUE, fig.width=6, fig.height=5, label=BSPECFIG-----
Bspec_plt <- plot(Bspec)
lines(Pspec$freq, Pspec$spec, col="red", lwd=2)


## ----eval=TRUE, echo=TRUE, label="AR spectrum"---------------------------
ntap <- 7
psd_ar <- psdcore(ats_ar, ntaper=ntap, refresh=TRUE)
dB(mean(psd_ar$spec))


## ----eval=TRUE, echo=TRUE, fig.width=6.0, fig.height=5.4, label=MAGPSDAR----
pilot_spec(ats_lm, ntap=ntap, remove.AR=100, plot=TRUE)
plot(Aspec, log="dB", add=TRUE, col="grey", lwd=4) 
plot(Aspec, log="dB", add=TRUE, lwd=3, lty=3)
spec.ar(ats_lm, log="dB", add=TRUE, lwd=2, col="grey40")


## ----eval=TRUE, echo=TRUE, fig.width=5, fig.height=4.5, label=SPECERR----
sp <- spectral_properties(as.tapers(1:50), p=0.95, db.ci=TRUE)
par(las=1)
plot(stderr.chi.upper ~ taper, sp, type="s", 
       ylim=c(-10,20), yaxs="i", xaxs="i",
       xlab=expression("number of tapers ("* nu/2 *")"), ylab="dB",
       main="Spectral uncertainties")
mtext("(additive factor)", line=.3)
lines(stderr.chi.lower ~ taper, sp, type="s")
lines(stderr.chi.median ~ taper, sp, type="s", lwd=2)
lines(stderr.chi.approx ~ taper, sp, type="s", col="red",lwd=2)
# to reach 3 db width confidence interval at p=.95
abline(v=33, lty=3)
legend("topright",
	c(expression("Based on "* chi^2 *"(p,"*nu*") and (1-p,"*nu*")"),
	  expression(""* chi^2 *"(p=0.5,"*nu*")"), 
	  "approximation"),
lwd=c(1,3,3), col=c("black","black","red"), bg="white")


## ----eval=TRUE, echo=TRUE, label="Compute spectral properties."----------
spp <- spectral_properties(Pspec$taper, db.ci=TRUE)
spa <- spectral_properties(Aspec$taper, db.ci=TRUE)
str(spa)
create_poly <- function(x, y, dy, from.lower=FALSE){
  xx <- c(x, rev(x))
  if (from.lower){
    yy <- c(y, rev(y+dy))
  } else {
    yy <- c(y+dy, rev(y-dy))
  }
  return(data.frame(xx=xx, yy=yy))
}
psppu <- create_poly(Pspec$freq, dB(Pspec$spec), spp$stderr.chi.upper)
pspau <- create_poly(Aspec$freq, dB(Aspec$spec), spa$stderr.chi.upper)
# and the Bayesian spectrum 95% limits
pspb <- create_poly(Bspec_plt$freq, Bspec_plt$spectrum[,1], Bspec_plt$spectrum[,3], from.lower=TRUE)


## ----eval=TRUE, echo=TRUE, fig.width=7, fig.height=4.5, label=MAGERR-----
plot(c(0,0.5),c(-5,40),col="white", 
       main="Project MAGNET Spectral Uncertainty (p > 0.95)",
       ylab="", xlab="spatial frequency, 1/km", yaxt="n", frame.plot=FALSE)
lines(c(2,1,1,2)*0.01,c(0,0,7,7))
text(.04, 3.5, "7 dB")
polygon(pspb$xx, dB(pspb$yy), col="light blue", border=NA)
text(0.26, 37, "Bayesian (bspec)", col="#0099FF", cex=cx<-0.9)
polygon(psppu$xx, psppu$yy, col="dark grey", border="black", lwd=0.2)
text(0.15, 6, "Light: adaptive\ntaper refinement\n(pspectrum)", cex=cx)
polygon(pspau$xx, pspau$yy, col="light grey", border="black", lwd=0.2)
text(0.40, 22, "Dark: Uniform\ntapering (psdcore)", cex=cx)


## ----eval=TRUE, echo=TRUE, fig.width=6, fig.height=3.5, label=MAGRES-----
frq <- Aspec$freq
relp <- (spa$resolution - spp$resolution) / spp$resolution
par(las=1, oma=rep(0,4), omi=rep(0,4), mar=c(4,3,2,0))
layout(matrix(c(1,2),1,2,byrow=TRUE), heights=c(2,2), widths=c(3,0.5), respect=TRUE)
plot(frq, relp,
     main="Percent change in spectral resolution",
     col="light grey", 
     ylim=yl<-c(0,35),
     type="h", xaxs="i", yaxs="i", 
     ylab="dB", xlab="frequency, 1/km")
lines(frq, relp)
text(0.25, 45, "Adaptive relative to fixed", cex=0.9)
par(mar=c(4,0,2,2))
# empirical distribution of values
boxplot(relp, range=0, main=sprintf("%.01f",median(relp)), axes=FALSE, ylim=yl, yaxs="i", notch=TRUE)
axis(4)


## ----eval=TRUE, echo=TRUE, label="Get adaptive history."-----------------
pspectrum(ats_lm, niter=4, plot=FALSE)
str(AH <- get_adapt_history())


## ----eval=TRUE, echo=TRUE, label="Some manipulation."--------------------
Freqs <- (AH$freq)
Dat <- AH$stg_psd
numd <- length(Freqs)
numit <- length(Dat)
StgPsd <- dB(matrix(unlist(Dat), ncol=numit))
Dat <- AH$stg_kopt
StgTap <- matrix(unlist(Dat), ncol=numit)
rm(Dat, AH)


## ----eval=TRUE, echo=FALSE, fig.width=4.8, fig.height=2.3, label=HIST1----
seqcols <- seq_len(numit)
itseq <- seqcols - 1
toadd <- matrix(rep(itseq, numd), ncol=numit, byrow=T)
par(xpd=TRUE, oma=rep(0,4), mar=c(1,4,3,2))
matplot(Freqs, StgPsd + (sc<-6)*toadd, type="l", lty=1, lwd=2, col="black",
             main="Adaptive estimation history", ylab="", xlab="",
             yaxt="n", frame.plot=FALSE)
text(0.52, 1.06*sc*itseq, itseq, cex=0.9)
lines(-1*c(1.5,1,1,1.5)*0.02,c(0,0,7,7))
text(-.06, 3.5, "7 dB", cex=0.8)
mtext("(a)", font=2, adj=0, line=0.6)
mtext("PSDs by stage", line=-0.4)


## ----eval=TRUE, echo=FALSE, fig.width=4.8, fig.height=2.1, label=HIST2----
par(xpd=TRUE, las=1, oma=rep(0,4), mar=c(1,4,2,2))
Cols <- rev(rev(brewer.pal(9, "PuBuGn"))[seqcols])
invisible(lapply(rev(seqcols), FUN=function(mcol, niter=numit, Frq=Freqs, Dat=StgTap, cols=Cols){
  iter <- (niter+1)-mcol
  y <- Dat[,mcol]
  icol <- Cols[mcol]
  if (iter==1){
    plot(Frq, y, type="h", col=icol, 
           main="", ylab="", 
           xlab="", #Spatial frequency",
           ylim=c(0,650), yaxs="i", frame.plot=FALSE)
  } else {
    lines(Frq, y, type="h", col=icol)
  }
  if (iter >= mcol){
    yf <- Dat[,(mcol+2)]
    lcol <- Cols[(mcol+2)]
    lines(Frq, yf, lty=3)
  }
  lines(Frq, y)
  #print(c(iter,mcol)) #1 5, 2 4, 3 3, 4 2, 5 1.
  x <- (c(0,1)+mcol)*.05+0.075
  y <- c(600,600,655,655,600)
  text(mean(x),max(y)-1.65*diff(range(y)), mcol-1, cex=0.9)
  polygon(c(x,rev(x),x[1]),y,border="black",col=icol)
}))
mtext("(b)", font=2, adj=0, line=0.5)
mtext("Tapers by stage", line=0.5)


## ----eval=TRUE, echo=FALSE, fig.width=4.8, fig.height=2.7, label=HIST3----
par(xpd=TRUE, las=1, oma=rep(0,4), mar=c(3.5,4,2,2))
#Cols <- rev(rev(brewer.pal(9, "PuBuGn"))[seqcols])
invisible(lapply(rev(seqcols), FUN=function(mcol, niter=numit, Frq=Freqs, Tap=StgTap, cols=Cols){
  iter <- (niter+1)-mcol
  tap <- Tap[,mcol]
  icol <- Cols[mcol]
  spp <- spectral_properties(as.tapers(tap), db.ci=TRUE)
  psppu <- create_poly(Frq, tap*0-(iter*0.76)**2, spp$stderr.chi.upper)
  if (iter==1){
	 plot(psppu$xx, psppu$yy, type="l", col=NA,
	 	main="", ylab="", xlab="", yaxt="n",
	    ylim=18*c(-1,0), 
	    yaxs="i", frame.plot=FALSE)
  }
  polygon(psppu$xx, psppu$yy, col=icol, border = "black") #, lwd = 0.2)
}))
mtext("(c)", font=2, adj=0, line=0.6)
lines(-1*c(1.5,1,1,1.5)*0.02, -1*c(0,0,7,7)-10)
text(-0.06, -3.5-10, "7 dB", cex=0.8)
mtext("Uncertainties by stage", line=0.6)
mtext("Spatial frequency, 1/km", side=1, line=2.3)
text(0.25, -14.5, "(uniform tapers)", font=3, cex=0.7)


## ----eval=TRUE, echo=TRUE, label=SYMCORT---------------------------------
suppressWarnings(symnum( cT <- cor(StgTap) ))


## ----eval=TRUE, echo=TRUE, label=SYMCORP---------------------------------
suppressWarnings(symnum( cP <- cor(StgPsd) ))


## ----eval=TRUE, echo=TRUE, label=SI--------------------------------------
sessionInfo()


