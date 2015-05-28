
# 4 groups (see abvoe for 6 groups) Log richness
samp.names[1] = "Mesozoic sampling rate"
Dinogroups[3] = "Ornithischia"
Dinogroups[4] = "Sauropodomorpha"
pdf("LogRichness_4gr_wholmeso.pdf")

par(mfrow=c(2,2)) #making six plots in one
# Trying to plot the CI's
par(mar=c(3,3,3,3))

doints = FALSE
plbins = 1:27;
for (ii in c(1,3,4,5)){
  #   plot(midpoints[plbins],N_est_glob[plbins,1,ii], type='l',col=color.list[1],ylim=c(0,max(c(N_est_glob[plbins,,ii],N_est_per[plbins,,ii],N_est_epo[plbins,,ii]),na.rm=TRUE)*1.1),xlim = rev(range(midpoints[plbins])),xlab="Ma",ylab="Eestimated true richness")
  if (doints==TRUE){  
  plot(midpoints[plbins],N_est_glob[plbins,1,ii], type='l',col=color.list[1],ylim=c(.1,max(c(N_est_int[plbins,,ii],N_est_per[plbins,,ii],N_est_epo[plbins,,ii]),na.rm=TRUE)*1.1),xlim = rev(range(midpoints[plbins])),log="y",xlab="Ma",ylab="Eestimated true richness")
  } else {
    plot(midpoints[plbins],N_est_glob[plbins,1,ii], type='l',col=color.list[1],ylim=c(0,max(c(N_est_glob[plbins,,ii],N_est_per[plbins,,ii],N_est_epo[plbins,,ii]),na.rm=TRUE)*1.1),xlim = rev(range(midpoints[plbins])),xlab="Ma",ylab="Eestimated true richness")
    
  }
#   for (jj in plbins){
#     abline(v=Bins[jj,1],col="lightgray",lwd=.5)
#     #ep(Bins[ii,1],2),c(yplmin,0))
#   }
#   abline(v=Bins[1,2],col="lightgray",lwd=.5)
#   
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
  lines(midpoints[plbins],N_est_glob[plbins,1,ii],type='l',col=color.list[1],lwd=1)
  
    polygon(c(rev(midpoints[plbins]),midpoints[plbins]),c(rev(N_est_per[plbins,2,ii]),N_est_per[plbins,3,ii]),col=color.list2[2],border=NA)
    lines(midpoints[plbins],N_est_per[plbins,1,ii],type='l',col=color.list[2],lwd=1)
  
  polygon(c(rev(midpoints[plbins]),midpoints[plbins]),c(rev(N_est_epo[plbins,2,ii]),N_est_epo[plbins,3,ii]),col=color.list2[3],border=NA)
  lines(midpoints[plbins],N_est_epo[plbins,1,ii],type='l',col=color.list[3],lwd=1)
  
  lines(midpoints[plbins],MeanSpec[plbins,ii,2])
  
  #   
  #   
  # For intervals, not laaaarge CIs
  if (doints==TRUE){
  for (jj in 1:(length(plbins)-1)){
    polygon(c(rev(midpoints[plbins[jj:(jj+1)]]),midpoints[plbins[jj:(jj+1)]]),c(rev(N_est_int[plbins[jj:(jj+1)],2,ii]),N_est_int[plbins[jj:(jj+1)],3,ii]),col=color.list2[4],border=NA)
    lines(midpoints[plbins[jj:(jj+1)]],N_est_int[plbins[jj:(jj+1)],1,ii],type='l',col=color.list[4],lwd=1)
    
  }
  }
  title(Dinogroups[ii])
# 
# if (ii==1){ylabel("Species richness")}
# if (ii==4){ylabel("Species richness"),xlabel("Myr")}
# if (ii==5){xlabel("Myr")}

}
legend("topleft",samp.names[c(1,2,3)],pch='l',col=color.list[c(1,2,3)],bg="white")
dev.off()








# 4 groups (see abvoe for 6 groups) Log richness

pdf("Richness_4gr_int_cretaceous.pdf")

par(mfrow=c(2,2)) #making six plots in one
# Trying to plot the CI's
par(mar=c(3,3,3,3))


plbins = 1:6;
for (ii in c(1,3,4,5)){
  #   plot(midpoints[plbins],N_est_glob[plbins,1,ii], type='l',col=color.list[1],ylim=c(0,max(c(N_est_glob[plbins,,ii],N_est_per[plbins,,ii],N_est_epo[plbins,,ii]),na.rm=TRUE)*1.1),xlim = rev(range(midpoints[plbins])),xlab="Ma",ylab="Eestimated true richness")
  
  plot(midpoints[plbins],N_est_glob[plbins,1,ii], type='l',col=color.list[1],ylim=c(5,max(c(N_est_int[plbins,,ii],N_est_per[plbins,,ii],N_est_epo[plbins,,ii]),na.rm=TRUE)*1.1),xlim = rev(range(midpoints[plbins])),log="y",xlab="Ma",ylab="Eestimated true richness")
  for (jj in plbins){
    abline(v=Bins[jj,1],col="lightgray",lwd=.5)
    #ep(Bins[ii,1],2),c(yplmin,0))
  }
  abline(v=Bins[1,2],col="lightgray",lwd=.5)
  polygon(c(rev(midpoints[plbins]),midpoints[plbins]),c(rev(N_est_glob[plbins,2,ii]),N_est_glob[plbins,3,ii]),col=color.list2[1],border=NA)
  lines(midpoints[plbins],N_est_glob[plbins,1,ii],type='l',col=color.list[1],lwd=1)
  
  #   polygon(c(rev(midpoints[plbins]),midpoints[plbins]),c(rev(N_est_per[plbins,2,ii]),N_est_per[plbins,3,ii]),col=color.list2[2],border=NA)
  #   lines(midpoints[plbins],N_est_per[plbins,1,ii],type='l',col=color.list[2],lwd=1)
  
  polygon(c(rev(midpoints[plbins]),midpoints[plbins]),c(rev(N_est_epo[plbins,2,ii]),N_est_epo[plbins,3,ii]),col=color.list2[3],border=NA)
  lines(midpoints[plbins],N_est_epo[plbins,1,ii],type='l',col=color.list[3],lwd=1)
  
  lines(midpoints[plbins],MeanSpec[plbins,ii,2])
  
  #   
  #   
  # For intervals, not laaaarge CIs
  for (jj in 1:(length(plbins)-1)){
    polygon(c(rev(midpoints[plbins[jj:(jj+1)]]),midpoints[plbins[jj:(jj+1)]]),c(rev(N_est_int[plbins[jj:(jj+1)],2,ii]),N_est_int[plbins[jj:(jj+1)],3,ii]),col=color.list2[4],border=NA)
    lines(midpoints[plbins[jj:(jj+1)]],N_est_int[plbins[jj:(jj+1)],1,ii],type='l',col=color.list[4],lwd=1)
    
  }
  title(Dinogroups[ii])
}
legend("bottomright",samp.names[c(1,3,4)],pch='l',col=color.list[c(1,3,4)],bg="white")
dev.off()






# Main figure plots for paper
# Plotting sampling rates (not for interval) for 4 groups


# Plotting the sampling rates per interval for early cretaceous. 
# rate plotting ---- 4 groups
pdf("SamplingRates_4gr_wholemeso.pdf")
par(mfrow=c(2,2)) #making six plots in one
# Trying to plot the CI's
par(mar=c(3,3,3,3))
plbins <- seq(1,27)
doints=FALSE

for (ii in c(1,3,4,5)){ # for each dinosaur group
  plot(midpoints[plbins],pois_glob[plbins,1,ii], type='l',col=color.list[1],ylim=c(0,.8),xlim = rev(range(midpoints[plbins])),xlab="Ma",ylab="Binomial sampling probability")
  
  for (jj in 1:3){
    abline(v=max(Bins[Bins[,5]==jj,1]),lwd=0.75)
  }
  abline(v=min(Bins[Bins[,5]==3,2]),lwd=0.75)
  
  for (jj in 1:6){
    abline(v=max(Bins[Bins[,6]==jj,1]),lwd=0.5)
  }
  abline(v=min(Bins[Bins[,6]==6,2]),lwd=0.5)
  
  # Adding interval lines

#   
#   for (jj in plbins){
#     abline(v=Bins[jj,1],col="lightgray",lwd=.5)
#     #ep(Bins[ii,1],2),c(yplmin,0))
#   }
#   abline(v=Bins[1,2],col="lightgray",lwd=.5)
#   #   
  #   # Adding period/epoch lines
  #   perlims = c(Bins[c(1,13,24),2],Bins[27,1])
  #   epolims = c(Bins[c(1,7,13,16,20,24),2],Bins[27,1])
  #   for (jj in 1:4){
  #     abline(v=perlims[jj],col="gray25")
  #   }
  #   for (jj in 1:7){
  #     abline(v=epolims[jj],col="gray14")
  #   }
  
  
  polygon(c(rev(midpoints[plbins]),midpoints[plbins]),c(rev(pois_glob[plbins,2,ii]),pois_glob[plbins,3,ii]),col=color.list2[1],border=NA)
  lines(midpoints[plbins],pois_glob[plbins,1,ii],type='l',col=color.list[1],lwd=1)
  
  
  # plot(midpoints,p_interval[,1,1], type='l',col=color.list[5],ylim=c(0,1))
  polygon(c(rev(midpoints[plbins]),midpoints[plbins]),c(rev(pois_period[plbins,2,ii]),pois_period[plbins,3,ii]),col=color.list2[2],border=NA)
  lines(midpoints[plbins],pois_period[plbins,1,ii],type='l',col=color.list[2],lwd=1)
  
  # plot(midpoints,p_epoch[,1,1], type='l',col=color.list[3],ylim=c(0,1))
  polygon(c(rev(midpoints[plbins]),midpoints[plbins]),c(rev(pois_epoch[plbins,2,ii]),pois_epoch[plbins,3,ii]),col=color.list2[3],border=NA)
  lines(midpoints[plbins],pois_epoch[plbins,1,ii],type='l',col=color.list[3],lwd=1)
  
  
  # We need something here for the missing values, the polygon doesnn't work for the whole period.
  # perhaps loop over them and make a polygon for each bin, with check for Na?
  # below commented out for the intervals
  if (doints==TRUE){
    for (jj in 1:(length(plbins)-1)){
      polygon(c(rev(midpoints[plbins[jj:(jj+1)]]),midpoints[plbins[jj:(jj+1)]]),c(rev(pois_interval[plbins[jj:(jj+1)],2,ii]),pois_interval[plbins[jj:(jj+1)],3,ii]),col=color.list2[4],border=NA)
      lines(midpoints[plbins[jj:(jj+1)]],pois_interval[plbins[jj:(jj+1)],1,ii],type='l',col=color.list[4],lwd=1)
      
    }
  }
  
  title(Dinogroups[ii])
  
  # Adding legend
  if (ii==5){
    if (doints==TRUE){
      legend("topleft",samp.names[1:4],pch='l',col=color.list[1:4],bg="white")
    } else {
      legend("topleft",samp.names[1:3],pch='l',col=color.list[1:3],bg="white")
    }
  }
  #,bty="n")
  # bty="n" to make legend box not transparent?
  
  
}
dev.off()# closing pdf print



# Plotting the sampling rates per interval for early cretaceous. 
# rate plotting ---- 4 groups
pdf("SamplingRates_4gr_int_lateCret.pdf")
par(mfrow=c(2,2)) #making six plots in one
# Trying to plot the CI's
par(mar=c(3,3,3,3))
plbins <- seq(1,6)
doints=TRUE

for (ii in c(1,3,4,5)){ # for each dinosaur group
  plot(midpoints[plbins],pois_glob[plbins,1,ii], type='l',col=color.list[1],ylim=c(0,.8),xlim = rev(range(midpoints[plbins])),xlab="Ma",ylab="Binomial sampling probability")
  # Adding interval lines
  for (jj in plbins){
    abline(v=Bins[jj,1],col="lightgray",lwd=.5)
    #ep(Bins[ii,1],2),c(yplmin,0))
  }
  abline(v=Bins[1,2],col="lightgray",lwd=.5)
  #   
  #   # Adding period/epoch lines
  #   perlims = c(Bins[c(1,13,24),2],Bins[27,1])
  #   epolims = c(Bins[c(1,7,13,16,20,24),2],Bins[27,1])
  #   for (jj in 1:4){
  #     abline(v=perlims[jj],col="gray25")
  #   }
  #   for (jj in 1:7){
  #     abline(v=epolims[jj],col="gray14")
  #   }
  
  
  polygon(c(rev(midpoints[plbins]),midpoints[plbins]),c(rev(pois_glob[plbins,2,ii]),pois_glob[plbins,3,ii]),col=color.list2[1],border=NA)
  lines(midpoints[plbins],pois_glob[plbins,1,ii],type='l',col=color.list[1],lwd=1)
  
  
  # plot(midpoints,p_interval[,1,1], type='l',col=color.list[5],ylim=c(0,1))
  polygon(c(rev(midpoints[plbins]),midpoints[plbins]),c(rev(pois_period[plbins,2,ii]),pois_period[plbins,3,ii]),col=color.list2[2],border=NA)
  lines(midpoints[plbins],pois_period[plbins,1,ii],type='l',col=color.list[2],lwd=1)
  
  # plot(midpoints,p_epoch[,1,1], type='l',col=color.list[3],ylim=c(0,1))
  polygon(c(rev(midpoints[plbins]),midpoints[plbins]),c(rev(pois_epoch[plbins,2,ii]),pois_epoch[plbins,3,ii]),col=color.list2[3],border=NA)
  lines(midpoints[plbins],pois_epoch[plbins,1,ii],type='l',col=color.list[3],lwd=1)
  
  
  # We need something here for the missing values, the polygon doesnn't work for the whole period.
  # perhaps loop over them and make a polygon for each bin, with check for Na?
  # below commented out for the intervals
  if (doints==TRUE){
    for (jj in 1:(length(plbins)-1)){
      polygon(c(rev(midpoints[plbins[jj:(jj+1)]]),midpoints[plbins[jj:(jj+1)]]),c(rev(pois_interval[plbins[jj:(jj+1)],2,ii]),pois_interval[plbins[jj:(jj+1)],3,ii]),col=color.list2[4],border=NA)
      lines(midpoints[plbins[jj:(jj+1)]],pois_interval[plbins[jj:(jj+1)],1,ii],type='l',col=color.list[4],lwd=1)
      
    }
  }
  
  title(Dinogroups[ii])
  
  # Adding legend
  if (ii==5){
    if (doints==TRUE){
      legend("topleft",samp.names[1:4],pch='l',col=color.list[1:4],bg="white")
    } else {
      legend("topleft",samp.names[1:3],pch='l',col=color.list[1:3],bg="white")
    }
  }
  #,bty="n")
  # bty="n" to make legend box not transparent?
  
  
}
dev.off()# closing pdf print

