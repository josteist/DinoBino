
# testing our figureprinting stuff.

# pdf("LogRichness_4gr_wholmeso.pdf")
png("Richness_WholeMesozoicII_withints.png",width=12,height=12,units='cm',res=2400)
par(mfrow=c(2,2),cex=0.4) #making six plots in one
# cex scales all text and symbols [NOT mtext]
# Trying to plot the CI's
par(mar=c(3.5,3.5,3.5,3.5))
labs = c("(a)","ign","(b)","(c)","(d)")
doints = TRUE
plbins = 1:27;
for (ii in c(1,3,4,5)){
  #   plot(midpoints[plbins],N_est_glob[plbins,1,ii], type='l',col=color.list[1],ylim=c(0,max(c(N_est_glob[plbins,,ii],N_est_per[plbins,,ii],N_est_epo[plbins,,ii]),na.rm=TRUE)*1.1),xlim = rev(range(midpoints[plbins])),xlab="Ma",ylab="Eestimated true richness")
  if (doints==TRUE){  
    ylimits=c(.1,max(c(N_est_int[plbins,,ii],N_est_per[plbins,,ii],N_est_epo[plbins,,ii]),na.rm=TRUE)*1.1)
    xlimits = rev(range(midpoints[plbins]))
    plot(midpoints[plbins],N_est_glob[plbins,1,ii], type='o',col=color.list[1],ylim=c(.1,max(c(N_est_int[plbins,,ii],N_est_per[plbins,,ii],N_est_epo[plbins,,ii]),na.rm=TRUE)*1.1),xlim = rev(range(midpoints[plbins])),log="y",xlab="",ylab="")
  } else {
    ylimits=c(0,max(c(N_est_glob[plbins,,ii],N_est_per[plbins,,ii],N_est_epo[plbins,,ii]),na.rm=TRUE)*1.1)
    xlimits=rev(range(midpoints[plbins]))
    plot(midpoints[plbins],N_est_glob[plbins,1,ii], type='o',col=color.list[1],ylim=c(0,max(c(N_est_glob[plbins,,ii],N_est_per[plbins,,ii],N_est_epo[plbins,,ii]),na.rm=TRUE)*1.1),xlim = rev(range(midpoints[plbins])),xlab="",ylab="")
    
  }
  
  # Not interval lines, but epochs/periods
  for (jj in 1:3){
    abline(v=max(Bins[Bins[,5]==jj,1]),lwd=0.75)
  }
  abline(v=min(Bins[Bins[,5]==3,2]),lwd=0.75)
  for (jj in 1:6){
    abline(v=max(Bins[Bins[,6]==jj,1]),lwd=0.5)
  }
  abline(v=min(Bins[Bins[,6]==6,2]),lwd=0.5)
  
  polygon(c(rev(midpoints[plbins]),midpoints[plbins]),c(rev(N_est_glob[plbins,2,ii]),N_est_glob[plbins,3,ii]),col=color.list2[1],border=NA)
  lines(midpoints[plbins],N_est_glob[plbins,1,ii],type='o',col=color.list[1],lwd=1)
  
  polygon(c(rev(midpoints[plbins]),midpoints[plbins]),c(rev(N_est_per[plbins,2,ii]),N_est_per[plbins,3,ii]),col=color.list2[2],border=NA)
  lines(midpoints[plbins],N_est_per[plbins,1,ii],type='o',col=color.list[2],lwd=1)
  
  polygon(c(rev(midpoints[plbins]),midpoints[plbins]),c(rev(N_est_epo[plbins,2,ii]),N_est_epo[plbins,3,ii]),col=color.list2[3],border=NA)
  lines(midpoints[plbins],N_est_epo[plbins,1,ii],type='o',col=color.list[3],lwd=1)
  
  lines(midpoints[plbins],MeanSpec[plbins,ii,2],type='o')

    # For intervals, not laaaarge CIs
  if (doints==TRUE){
    for (jj in 1:(length(plbins)-1)){
      polygon(c(rev(midpoints[plbins[jj:(jj+1)]]),midpoints[plbins[jj:(jj+1)]]),c(rev(N_est_int[plbins[jj:(jj+1)],2,ii]),N_est_int[plbins[jj:(jj+1)],3,ii]),col=color.list2[4],border=NA)
      lines(midpoints[plbins[jj:(jj+1)]],N_est_int[plbins[jj:(jj+1)],1,ii],type='o',col=color.list[4],lwd=1)
      
    }
  }
  title(Dinogroups[ii])
  # 
  if (ii==1){
    mtext("Species richness",side=2,padj=-3,cex=0.5)
  }
  if (ii==4){
    mtext("Species richness",side=2,padj=-3,cex=0.5);
    mtext("Myr",side=1,padj=2.8,cex=.5);
  }
  if (ii==5){
    mtext("Myr",side=1,padj=2.8,cex=.5)
  }
  mtext(labs[ii],side=3,adj=0,padj=-0.3,cex=0.5)
}
# legend("topleft",samp.names[c(1,2,3)],pch='l',col=color.list[c(1,2,3)],bg="white")
legend("topleft",c(samp.names[c(1,2,3)],"Observed species count"),pch='l',col=c(color.list[c(1,2,3)],"black"),bg="white")


dev.off()
