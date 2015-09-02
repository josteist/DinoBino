# Making the rate figure. We want two panels
# (A) estimated interval sampling rates for all 4 groups
# (B) estimated binomial sampling probability for all 4 groups.

# Two figures, with and w/o confidence intervals.
docis = TRUE;

# Loading color library.
library(RColorBrewer)
color.list <- brewer.pal(4,"Set1")
color.list2 <- color.list
for (ii in 1:4){
  color.list2[ii] = paste(color.list[ii],'44',sep=""); 
  # Making one slightly transparent.
}
# Also printed as Pngs
# jpeg("SamplingRates_170815_ci.jpg",width=18,height=12,units='cm',res=2400)
par(mfrow=c(2,1),cex=0.6)
#(A) midpoints p_int_median
plbins = 1:27;
par(mar=c(4,4,4,4))
#   plot(midpoints[plbins],N_est_glob[plbins,1,ii], type='l',col=color.list[1],ylim=c(0,max(c(N_est_glob[plbins,,ii],N_est_per[plbins,,ii],N_est_epo[plbins,,ii]),na.rm=TRUE)*1.1),xlim = rev(range(midpoints[plbins])),xlab="Ma",ylab="Eestimated true richness")
# docis = TRUE
ylimits=c(0,max(p_int_median,na.rm=TRUE)*1.1)
ylimits = c(-.1,1.3)
xlimits = rev(range(midpoints[plbins]))
# Need to loop over polygons here since some are NA's.
plot(midpoints[plbins],p_int_median[plbins,ii,1], type='o',col=color.list[1],ylim=ylimits,xlim = rev(range(midpoints[plbins])),xlab="Ma",ylab=expression(paste(lambda, " - sampling rate")))
# for (ii in 1:4){
for (ii in 1:4){
  lines(midpoints[plbins],p_int_median[plbins,ii,1], type='o',col=color.list[ii],ylim=ylimits,xlim = rev(range(midpoints[plbins])),xlab="Ma")
  if (docis==TRUE){
    for (ss in 2:27){
      polygon(c(rev(midpoints[plbins[(ss-1):ss]]),midpoints[plbins[(ss-1):ss]]),c(rev(p_int_median[plbins[(ss-1):ss],ii,2]),p_int_median[plbins[(ss-1):ss],ii,3]),col=color.list2[ii],border=NA)
    }
  }
}
for (jj in plbins){
  abline(v=Bins[jj,1],col="lightgray",lwd=.5)
  #ep(Bins[ii,1],2),c(yplmin,0))
}
abline(v=Bins[1,2],col="lightgray",lwd=.5)

for (jj in 1:6){
  abline(v=max(Bins[Bins[,6]==jj,1]),col="grey10",lwd=0.5)
}
abline(v=min(Bins[Bins[,6]==6,2]),col="grey10",lwd=0.5)
# legend("topright",Dinogroups,pch='l',col=color.list[c(1,2,3,4)],bg="white")
legend(125,par("usr")[4],Dinogroups,pch='l',col=color.list[c(1,2,3,4)],bg="white")

for (ii in 1:27){
  # rect(Bins[ii,1],1,Bins[ii,2],2)
  rect(Bins[ii,1],-.1,Bins[ii,2],0,lwd=0.3)
}
# rect(xleft, ybottom, xright, ytop,
for (ii in 1:27){
  if (Bins[ii,3]>4){
    text(midpoints[ii],-.05,substr(interval.names[ii],1,2),cex=.65)
  }
}

# Probability plot
ylimits = c(-.1,1)
plot(midpoints[plbins],p_binos[plbins,1,1], type='o',col=color.list[1],ylim=ylimits,xlim = rev(range(midpoints[plbins])),xlab="Ma",ylab="Binomial sampling probability")
# for (ii in 1:4){
for (ii in 1:4){
  lines(midpoints[plbins],p_binos[plbins,ii,1], type='o',col=color.list[ii],ylim=ylimits,xlim = rev(range(midpoints[plbins])),xlab="Ma")
  if (docis==TRUE){
    for (ss in 2:27){
      polygon(c(rev(midpoints[plbins[(ss-1):ss]]),midpoints[plbins[(ss-1):ss]]),c(rev(p_binos[plbins[(ss-1):ss],ii,2]),p_binos[plbins[(ss-1):ss],ii,3]),col=color.list2[ii],border=NA)
    }
  }
}
for (jj in plbins){
  abline(v=Bins[jj,1],col="lightgray",lwd=.5)
  #ep(Bins[ii,1],2),c(yplmin,0))
}
abline(v=Bins[1,2],col="lightgray",lwd=.5)

for (jj in 1:6){
  abline(v=max(Bins[Bins[,6]==jj,1]),col="grey10",lwd=0.5)
}
abline(v=min(Bins[Bins[,6]==6,2]),col="grey10",lwd=0.5)
# legend("topright",Dinogroups,pch='l',col=color.list[c(1,2,3,4)],bg="white")
# trying to draw abbreviations for all intervals in the richness plots.

figlims = par("usr"); # getting figure limits. Try to box the abbreviations from miny to 0

for (ii in 1:27){
  # rect(Bins[ii,1],1,Bins[ii,2],2)
  rect(Bins[ii,1],-.1,Bins[ii,2],0,lwd=0.3)
}
# rect(xleft, ybottom, xright, ytop,
for (ii in 1:27){
  if (Bins[ii,3]>4){
    text(midpoints[ii],-.05,substr(interval.names[ii],1,2),cex=.65)
  }
}

# dev.off()

