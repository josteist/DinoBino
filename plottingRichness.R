
## Here we'll try to plot the MLE richnesses ----
# Though the assumptions of the methods are not really met, they do seem to work for 
# intervals of moderate length, with not so large rates of extinction/speciation.

# We have binomial sampling probs in arrays;
# interval by [mle,ci1,ci2] by dinogroup
# p_glob
# p_period
# p_epoch
# p_interval (should be called stage/age)

# We will use the mean/median/mode of species richnesses stored in
# Meanspec [age by group by mean/median/mode]

# We will for the mesozoiic plot not use interval/stage specific sampling rates. So for the big plot we will have
# six groups, three types of sampling rates, mle/ci/ci2/ci3/ci4 (possibly we'll skip ci1 and ci2)
# 
floor(MeanSpec[,1,2]/p_glob[,1,1])
floor(MeanSpec[,1,2]/p_glob[,2,1])
floor(MeanSpec[,1,2]/p_glob[,3,1])

N_est_glob = array(NA,c(27,3,6))
N_est_per  = array(NA,c(27,3,6))
N_est_epo  = array(NA,c(27,3,6))
N_est_int  = array(NA,c(27,3,6))

MeanSpec = round(MeanSpec)
dg = 1; #dinogroup
int = 1; #interval
for (dg in 1:6){
  for (int in 1:27){
    if (MeanSpec[int,dg,2]>0){
      N_est_glob[int,1,dg] =     estimatetrue(MeanSpec[int,dg,2],p_glob[int,1,dg])[1] #double MLE use first as this is MLE
      N_est_glob[int,2,dg] = max(estimatetrue(MeanSpec[int,dg,2],p_glob[int,2,dg])) #use the max here to get upper CI
      N_est_glob[int,3,dg] = min(estimatetrue(MeanSpec[int,dg,2],p_glob[int,3,dg])) #use the min here to get lower CI
      
      N_est_per[int,1,dg] =     estimatetrue(MeanSpec[int,dg,2],p_period[int,1,dg])[1] #double MLE use first as this is MLE
      N_est_per[int,2,dg] = max(estimatetrue(MeanSpec[int,dg,2],p_period[int,2,dg])) #use the max here to get upper CI
      N_est_per[int,3,dg] = min(estimatetrue(MeanSpec[int,dg,2],p_period[int,3,dg])) #use the min here to get lower CI
      
      N_est_epo[int,1,dg] =     estimatetrue(MeanSpec[int,dg,2],p_epoch[int,1,dg])[1] #double MLE use first as this is MLE
      N_est_epo[int,2,dg] = max(estimatetrue(MeanSpec[int,dg,2],p_epoch[int,2,dg])) #use the max here to get upper CI
      N_est_epo[int,3,dg] = min(estimatetrue(MeanSpec[int,dg,2],p_epoch[int,3,dg])) #use the min here to get lower CI
      if (is.na(p_interval[int,1,dg])){
      N_est_int[int,1,dg] =     estimatetrue(MeanSpec[int,dg,2],p_interval[int,1,1])[1] #double MLE use first as this is MLE
      N_est_int[int,2,dg] = max(estimatetrue(MeanSpec[int,dg,2],p_interval[int,2,1])) #use the max here to get upper CI
      N_est_int[int,3,dg] = min(estimatetrue(MeanSpec[int,dg,2],p_interval[int,3,1])) #use the min here to get lower CI
      } else {
      # If interval specific rate not available for the group, use the estimate for all dinosaurs
        N_est_int[int,1,dg] =     estimatetrue(MeanSpec[int,dg,2],p_interval[int,1,dg])[1] #double MLE use first as this is MLE
        N_est_int[int,2,dg] = max(estimatetrue(MeanSpec[int,dg,2],p_interval[int,2,dg])) #use the max here to get upper CI
        N_est_int[int,3,dg] = min(estimatetrue(MeanSpec[int,dg,2],p_interval[int,3,dg])) #use the min here to get lower CI
        
      }
    }
  }
}

# All 6 groups (see below for 4 groups)
par(mfrow=c(3,2)) #making six plots in one
# Trying to plot the CI's
par(mar=c(2,2,2,2))



plbins = 1:6;
for (ii in 1:6){
#   plot(midpoints[plbins],N_est_glob[plbins,1,ii], type='l',col=color.list[1],ylim=c(0,max(c(N_est_glob[plbins,,ii],N_est_per[plbins,,ii],N_est_epo[plbins,,ii]),na.rm=TRUE)*1.1),xlim = rev(range(midpoints[plbins])),xlab="Ma",ylab="Eestimated true richness")
  
  plot(midpoints[plbins],N_est_glob[plbins,1,ii], type='l',col=color.list[1],ylim=c(5,max(c(N_est_int[plbins,,ii],N_est_per[plbins,,ii],N_est_epo[plbins,,ii]),na.rm=TRUE)*1.1),xlim = rev(range(midpoints[plbins])),log="y",xlab="Ma",ylab="Eestimated true richness")
  
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




# 4 groups (see abvoe for 6 groups)
par(mfrow=c(2,2)) #making six plots in one
# Trying to plot the CI's
par(mar=c(2,2,2,2))



plbins = 1:6;
for (ii in c(1,3,4,5)){
  #   plot(midpoints[plbins],N_est_glob[plbins,1,ii], type='l',col=color.list[1],ylim=c(0,max(c(N_est_glob[plbins,,ii],N_est_per[plbins,,ii],N_est_epo[plbins,,ii]),na.rm=TRUE)*1.1),xlim = rev(range(midpoints[plbins])),xlab="Ma",ylab="Eestimated true richness")
  
  plot(midpoints[plbins],N_est_glob[plbins,1,ii], type='l',col=color.list[1],ylim=c(5,max(c(N_est_int[plbins,,ii],N_est_per[plbins,,ii],N_est_epo[plbins,,ii]),na.rm=TRUE)*1.1),xlim = rev(range(midpoints[plbins])),log="y",xlab="Ma",ylab="Eestimated true richness")
  
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



