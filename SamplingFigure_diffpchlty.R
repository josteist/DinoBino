# Samp figs. Make them possible to interpret in black/white, i.e. change lines and patches
# Making the rate figure. We want two panels
# (A) estimated interval sampling rates for all 4 groups
# (B) estimated binomial sampling probability for all 4 groups.

#Trying to standardize figs: 17 cm wide (ca a4 page) and 4 cm height per panel

# Two figures, with and w/o confidence intervals. Just change this to false/true
# docis = FALSE;

# Loading color library.
library(RColorBrewer)
color.list <- brewer.pal(4,"Set1")
color.list2 <- color.list
for (ii in 1:4){
  color.list2[ii] = paste(color.list[ii],'44',sep=""); 
  # Making one slightly transparent.
}
patch.list = c(15,16,17,18)
lines.list = c(1,2,3,4)
# 
# jpeg("SamplingRates_250815_genera_noci.jpg",width=18,height=8,units='cm',res=2400)

# jpeg(printfilename,width=18,height=8,units='cm',res=1200)
par(mfrow=c(2,1),cex=0.5)
#(A) midpoints p_int_median
plbins = 1:27;
if (docis==FALSE){
ylimits = c(0,max(p_int_median[,,1],na.rm=TRUE))
} else {
  ylimits = c(0,max(p_int_median,na.rm=TRUE))
}
xlimits = rev(range(midpoints[plbins]))
# Need to loop over polygons here since some are NA's.
par(mar=(c(0,4,4,2)))
ii=1
# plot.new()
plot(midpoints[plbins],p_int_median[plbins,ii,1], type='o',col=color.list[1],ylim=ylimits,xlim = rev(range(midpoints[plbins])),xlab="Ma",ylab=expression(paste(lambda, " - sampling rate")),xaxt="n",pch=patch.list[ii],lty=lines.list[ii])
# for (ii in 1:4){
abline(h = 0,col="black",lwd=1)
for (ii in 1:4){
  lines(midpoints[plbins],p_int_median[plbins,ii,1], type='o',col=color.list[ii],ylim=ylimits,xlim = rev(range(midpoints[plbins])),xlab="Ma",pch=patch.list[ii],lty=lines.list[ii])
  if (docis==TRUE){
    for (ss in 2:27){
      polygon(c(rev(midpoints[plbins[(ss-1):ss]]),midpoints[plbins[(ss-1):ss]]),c(rev(p_int_median[plbins[(ss-1):ss],ii,2]),p_int_median[plbins[(ss-1):ss],ii,3]),col=color.list2[ii],border=NA)
    }
  }
}
for (jj in plbins){
  lines(rep(Bins[jj,1],2),c(-.2,ylimits[2]*1.2),col="lightgray",lwd=1)
}
lines(rep(Bins[1,2],2),c(-.2,ylimits[2]*1.2),col="lightgray",lwd=1)
for (jj in 1:6){
  lines(rep(max(Bins[Bins[,6]==jj,1]),2),c(-.2,ylimits[2]*1.2),col="grey10",lwd=1)
}
lines(rep(min(Bins[Bins[,6]==6,2]),2),c(-.2,ylimits[2]*1.2),col="grey10",lwd=1)

# legend(125,par("usr")[4],Dinogroups,pch='l',col=color.list[c(1,2,3,4)],bg="white")

# Probability plot
par(mar=(c(4,4,0,2)))
ylimits = c(-.2,1)
ii = 1;
plot(midpoints[plbins],p_binos[plbins,1,1], type='o',col=color.list[1],ylim=ylimits,xlim = rev(range(midpoints[plbins])),yaxt="n",ylab = "Binomial sampling probability",xlab="Ma",pch=patch.list[ii],lty=lines.list[ii],xaxt="n")
# supressing axis at first with yaxt="n", then defining with axis(side)
# ,xlab="Ma",ylab="Binomial sampling probability"
axis(side=2,at=c(0.0,0.2,0.4,0.6,.8,1.0))
# title()
# for (ii in 1:4){
for (ii in 1:4){
  lines(midpoints[plbins],p_binos[plbins,ii,1], type='o',col=color.list[ii],ylim=ylimits,xlim = rev(range(midpoints[plbins])),xlab="Ma",pch=patch.list[ii],lty=lines.list[ii])
  if (docis==TRUE){
    for (ss in 2:27){
      polygon(c(rev(midpoints[plbins[(ss-1):ss]]),midpoints[plbins[(ss-1):ss]]),c(rev(p_binos[plbins[(ss-1):ss],ii,2]),p_binos[plbins[(ss-1):ss],ii,3]),col=color.list2[ii],border=NA)
    }
  }
}
for (jj in plbins){
  lines(rep(Bins[jj,1],2),c(0,1.2),col="lightgray",lwd=1)
}
lines(rep(Bins[1,2],2),c(0,1.2),col="lightgray",lwd=1)
for (jj in 1:6){
  lines(rep(max(Bins[Bins[,6]==jj,1]),2),c(0,1.2),col="grey10",lwd=1)
}
lines(rep(min(Bins[Bins[,6]==6,2]),2),c(0,1.2),col="grey10",lwd=1)

figlims = par("usr"); # getting figure limits. Try to box the abbreviations from miny to 0

# xtx = array(NA,c(7,1))
# for (ii in 1:6){
  # xtx[ii] = max(Bins[Bins[,6]==ii,1])
  # 
# }
# xtx[7] = min(Bins[Bins[,6]==6],1)
axis(1,at=c(240,200,150,100,66)) #xtx)

for (ii in 1:27){
  # rect(Bins[ii,1],1,Bins[ii,2],2)
  rect(Bins[ii,1],-.1,Bins[ii,2],0,lwd=1)
}
color.list.eps = c("grey70","grey90")
for (ii in 1:27){
  rect(Bins[ii,1],-.1,Bins[ii,2],0,lwd=1,col = color.list.eps[Bins[ii,6] %% 2])
}

for (ii in 1:6){
  rect(max(Bins[Bins[,6]==ii,1]),-.2,min(Bins[Bins[,6]==ii,2]),-.1,lwd=1,col = color.list.eps[ii %% 2])
}
for (ii in 1:27){
  text(midpoints[ii],-.05,substr(interval.names[ii],1,3),cex=.5,srt=90)
}
for (ii in 1:6){
  text((min(Bins[Bins[,6]==ii,2])+max(Bins[Bins[,6]==ii,1]))/2,-.15,epoch.names[ii],cex=0.6)
}

# dev.off()

