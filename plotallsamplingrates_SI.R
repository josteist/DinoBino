# Samp figs. Make them possible to interpret in black/white, i.e. change lines and patches
# Making the rate figure. We want two panels
# (A) estimated interval sampling rates for all 4 groups
# (B) estimated binomial sampling probability for all 4 groups.

# FOr CI, plotting all the 100 replicated sampling rates.

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
docis = FALSE
# jpeg(printfilename,width=18,height=8,units='cm',res=1200)
par(mfrow=c(4,1),cex=0.5)
#(A) midpoints p_int_median
for (dd in 1:4){
  if (dd==1){
    p_int_plot = p_interval_dinos;
  } else if (dd==2){
    p_int_plot = p_interval_ornits;
  } else if (dd==3){
    p_int_plot = p_interval_sauros;
  } else if (dd==4){
    p_int_plot = p_interval_theros;
  }
  
  plbins = 1:27;
  
  ylimits = c(0,max(p_int_plot[,,1],na.rm=TRUE))
  xlimits = rev(range(midpoints[plbins]))
  if (dd==4){
    ylimits = c(-.2,max(p_int_plot[,,1],na.rm=TRUE))
  }
  # Need to loop over polygons here since some are NA's.
  par(mar=(c(0,4,0,2)))
  
  # plot.new()
  plot(midpoints[plbins],p_int_plot[1,plbins,1], type='l',col=color.list[dd],ylim=ylimits,xlim = rev(range(midpoints[plbins])),xlab="Ma",ylab=expression(paste(lambda, " - sampling rate")),xaxt="n",pch=patch.list[dd],lty=lines.list[dd])
  # for (ii in 1:4){
  # abline(h = 0,col="black",lwd=1)
  for (rr in 2:100){
    lines(midpoints[plbins],p_int_plot[rr,plbins,1], type='l',col=color.list[dd],ylim=ylimits,xlim = rev(range(midpoints[plbins])),xlab="Ma",pch=patch.list[dd],lty=lines.list[dd])
  }
  lines(midpoints[plbins],p_int_median[plbins,dd,1],type='l',col='black')
  mtext(Dinogroups[dd],side=2,adj=.6,padj=1.6,cex=0.5)
  for (jj in plbins){
    lines(rep(Bins[jj,1],2),c(-.2,ylimits[2]*1.2),col="lightgray",lwd=1)
  }
  lines(rep(Bins[1,2],2),c(-.2,ylimits[2]*1.2),col="lightgray",lwd=1)
  for (jj in 1:6){
    lines(rep(max(Bins[Bins[,6]==jj,1]),2),c(-.2,ylimits[2]*1.2),col="grey10",lwd=1)
  }
  lines(rep(min(Bins[Bins[,6]==6,2]),2),c(-.2,ylimits[2]*1.2),col="grey10",lwd=1)
  
  # legend(125,par("usr")[4],Dinogroups,pch='l',col=color.list[c(1,2,3,4)],bg="white")
  
}
figlims = par("usr"); # getting figure limits. Try to box the abbreviations from miny to 0
axis(1,at=c(240,200,150,100,66)) #xtx)

for (ii in 1:27){
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

