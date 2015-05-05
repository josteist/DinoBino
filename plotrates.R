# Generating large arrays for poisson rates to be plotted
pois_glob = array(NA,c(27,3,6))
pois_period = array(NA,c(27,3,6))
pois_epoch  = array(NA,c(27,3,6))
pois_interval = poisrates_interval;
for (ii in 1:27){
  for (jj in 1:6){
  pois_glob[ii,,jj]<-poisrates[jj,]
  }
  pois_period[ii,,]<-poisrates_period[Bins[ii,5],,]
  pois_epoch[ii,,] <- poisrates_epoch[Bins[ii,6],,]
  
}


# rate plotting ----
pdf("SamplingRates_cret_int.pdf")
par(mfrow=c(3,2)) #making six plots in one
# Trying to plot the CI's
par(mar=c(2,2,2,2))
plbins <- seq(1,6)
doints=TRUE

for (ii in 1:6){ # for each dinosaur group
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
  if (ii==6){
    if (doints==TRUE){
      legend("bottomright",samp.names[1:4],pch='l',col=color.list[1:4],bg="white")
    } else {
    legend("bottomright",samp.names[1:3],pch='l',col=color.list[1:3],bg="white")
    }
  }
  #,bty="n")
  # bty="n" to make legend box not transparent?
  
  
}
dev.off()# closing pdf print




## Below are binomial sampling probabilities ----
# Cretaceous only.
plbins = which(Bins[,5]==3)
# if all;
plbins = 1:27;
# pdf("SamplingProbs_Cret_all.pdf")
par(mfrow=c(3,2)) #making six plots in one
# Trying to plot the CI's
par(mar=c(2,2,2,2))



for (ii in 1:6){ # for each dinosaur group
  plot(midpoints[plbins],p_glob[plbins,1,ii], type='l',col=color.list[1],ylim=c(0,1),xlim = rev(range(midpoints[plbins])),xlab="Ma",ylab="Binomial sampling probability")
  
  
  
  polygon(c(rev(midpoints[plbins]),midpoints[plbins]),c(rev(p_glob[plbins,2,ii]),p_glob[plbins,3,ii]),col=color.list2[1],border=NA)
  lines(midpoints[plbins],p_glob[plbins,1,ii],type='l',col=color.list[1],lwd=1)
  
  
  # plot(midpoints,p_interval[,1,1], type='l',col=color.list[5],ylim=c(0,1))
  polygon(c(rev(midpoints[plbins]),midpoints[plbins]),c(rev(p_period[plbins,2,ii]),p_period[plbins,3,ii]),col=color.list2[2],border=NA)
  lines(midpoints[plbins],p_period[plbins,1,ii],type='l',col=color.list[2],lwd=1)
  
  # plot(midpoints,p_epoch[,1,1], type='l',col=color.list[3],ylim=c(0,1))
  polygon(c(rev(midpoints[plbins]),midpoints[plbins]),c(rev(p_epoch[plbins,2,ii]),p_epoch[plbins,3,ii]),col=color.list2[3],border=NA)
  lines(midpoints[plbins],p_epoch[plbins,1,ii],type='l',col=color.list[3],lwd=1)
  
  
  # We need something here for the missing values, the polygon doesnn't work for the whole period.
  # perhaps loop over them and make a polygon for each bin, with check for Na?
  for (jj in 1:(length(plbins)-1)){
  polygon(c(rev(midpoints[plbins[jj:(jj+1)]]),midpoints[plbins[jj:(jj+1)]]),c(rev(p_interval[plbins[jj:(jj+1)],2,ii]),p_interval[plbins[jj:(jj+1)],3,ii]),col=color.list2[4],border=NA)
  lines(midpoints[plbins[jj:(jj+1)]],p_interval[plbins[jj:(jj+1)],1,ii],type='l',col=color.list[4],lwd=1)
  
  }
  
  title(Dinogroups[ii])
  # Adding interval lines
  for (jj in plbins){
    abline(v=Bins[jj,1],col="lightgray")
    #ep(Bins[ii,1],2),c(yplmin,0))
  }
  abline(v=Bins[1,2],col="lightgray")
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
  # Adding legend
  if (ii==6){
    legend("bottomright",samp.names[1:4],pch='l',col=color.list[1:4],bg="white")}
  #,bty="n")
  # bty="n" to make legend box not transparent?
  
  
}
# dev.off()# closing pdf print
