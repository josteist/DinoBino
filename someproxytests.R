# Making some comparisons between our sampling rates and earlier used proxies.

# THe dinosaurbearing formation count of Barrett et al 
proxy <- c(19,30,33,38,27,28,44,41,33,33,23,23,21,21,14,14,17,20,32,32,24,24,65,70,62,61,55,60,97,102,107,106,113,118,109,105,112,120,120,123,107,100,80,81,85,86,88,93,135,153,156,172)
# Th
# These intervals are the other way around and do not have values for Ladinian, they are also late and early
proxy1 = array(NA,dim=c(27,1))
for (ss in 1:26){
  proxy1[ss] = mean(proxy[(52-(ss-1)*2):(53-(ss*2))])
}
sealevel <- c(202.5281489,226.3890431,216.7429194,225.7313194,230.9619407,227.4257238,181.2598273,154.1008864,166.3174695,121.83234,77.37530179,135.7678534,122.3555231,125.825266,89.41541861,51.91750712,69.34271044,94.25307352,80.69268013,59.00848256,49.10071958,29.28765987,2.330933794,-12.52399474,52.72477256,17.42634031,17.93230681)
# These data are taken from Butler, Benson, Carrano Mannion and Upchurch 2011 Proc.Roy.
# THey are Haq numbers and does not include Ladinian
sealevelM <- c(30.91627451,32.17527132,24.69173913,17.45714286,22.61785714,45.586137,7.097748974,-7.313862805,-1.043939233,-24.1136305,-18.30485144,-18.82420665,-8.762986,-1.180405815,-13.77864975,-21.63454761,-30.79817945,-31.53682994)
# These are Miller numbers also from Butler. They include bins 112:129
plot(p_int_median[,1,1],proxy1)
# Continental area of Butler et al 2011
contarea = c(109,110,118,118,121,118,129,131,125,125,127,127,130,130,132,128,133,138,138.5,139,142,146,150,134,134,134)

# I guess we could have a plot with the main message that we do capture something that is not in 
# - the richness counts themselves
# - three commonly used proxies
# DBF from Barrett
# DBC from our data
# Sea-level [which is also what others have found]

plot(contarea,p_binos[1:26,1,1])




# Thus it seems like our sampling rates capture an aspect of bias that has not been captured before.
par(mfrow=c(2,2),cex=0.8)
plot(log10(SpecRich_test[,1]),p_int_median[,1,1],ylab="Binomial sampling probability",xlab="Observed species richness")
plot(log10(proxy1),p_int_median[,1,1],ylab="Binomial sampling probability",xlab="Dinosaur bearing formations")
plot(log10(collections),p_int_median[,1,1],ylab="Binomial sampling probability",xlab="Dinosaur bearing collections")
plot((sealevel),p_int_median[,1,1],ylab="Binomial sampling probability",xlab="Sea level")

# We could compare this to plots of these 'proxies' agains observed species counts.
# If we also compare with our estimated diversities we can comment more on either the 'common-cause' (Peters)
# or the 'redundancy' (Benton) of these proxies.
# 

# I guess that if there is a relationship between, e.g. sea-level, and estimated diversities this indicates that
# sea-level has been a driver of species dynamics, since our estimates have taken sampling bias into account. 
# If the 'common-cause' (i.e. sealevel driving both sampling and richness), and we have reconstructed richnesses
# devoid of sampling bias then a remaining relationship indicates that sea-level has driven dynamics.

fit_ss <- lm(log10(N_est_int[,1,1])~sealevel)
plot(sealevel,log10(N_est_int[,1,1]))
plot(diff(sealevel),diff(log10(N_est_int[,1,1])))
fit_seatime <- lm(sealevel~midpoints)
plot(residuals(fit_seatime),log10(N_est_int[,1,1]))
plot(residuals(fit_seatime),log10(SpecRich_test[,1]))
plot(residuals(fit_seatime))

fit1 <- lm(log10(SpecRich_test[,1])~p_int_median[,1,1])
fit2 <- lm(p_int_median[,1,1]~log10(proxy1))
plot(log10(proxy1),p_int_median[,1,1])
abline(fit2)
fit3 <- lm(log10(SpecRich_test[,1])~p_int_median[,1,1])
fit4 <- lm(log10(SpecRich_test[,1])~p_int_median[,1,1])

plot(midpoints,sealevel,type="l")
plot(midpoints,log10(N_est_int[,1,1]),type="l",xlim=rev(range(midpoints)),ylim=c(0,3))
lines(midpoints,2+residuals(fit_seatime)/max(residuals(fit_seatime)))
fit_specest = lm(log10(N_est_int[,1,1])~midpoints)
plot(residuals(fit_specest),2+residuals(fit_seatime)/max(residuals(fit_seatime)))


## Thinking again here.
# The relationship between sea-level and richness.
# There are several possible interpretations.
# Sea-level has direct effect on richness. Could be interpreted as lower sea-level --> more land --> more space for species
# Since there are trends in sea-level this must be detrended

# The relationship between sea-level and sampling.
# Sea-level itself affects sampling. Sea level could affect the amount of areas that could potentially leave fossils.
# 

# Major question: do we need to detrend a time-series with a trend if compares to a time-series without a trend?
# Main point of detrending is to remove spurious correlations due to trend in both series. If there is not trend in one
# and we wish to compare with another with a trend we need not detrend.
plot(sealevel,log10(N_est_int[,1,1]))
# if we detrend, i.e. perform lineare regression across time and use residuals.
par(mfrow=c(2,2))
for (ii in 1:4){
  use = which(!is.na(N_est_int[,ii,1]))
plot(residuals(lm(sealevel[use]~midpoints[use])),residuals(lm(log10(N_est_int[use,ii,1])~midpoints[use])))

}
# Miller sealevels
par(mfrow=c(2,2))
for (ii in 1:4){
  use = which(!is.na(N_est_int[1:18,ii,1]))
  plot(residuals(lm(sealevelM[use]~midpoints[use])),residuals(lm(log10(N_est_int[use,ii,1])~midpoints[use])))
  
}
ii=3
use = which(!is.na(N_est_int[1:18,ii,1]))
fit1<- lm(residuals(lm(sealevelM[use]~midpoints[use]))~residuals(lm(log10(N_est_int[use,ii,1])~midpoints[use])))
fit1<- lm(residuals(lm(log10(N_est_int[use,ii,1])~midpoints[use]))~residuals(lm(sealevelM[use]~midpoints[use])))
summary(fit1)
# So there is a relationship between the detrended Miller data and the detrended log10 richness of Ornithischians.
# Is this something worth reporting?

# Alternative hypotheses are that changes in sea-level affects either sampling or dynamics.
# since arrays are in reverse time use -diff
plot(-diff(sealevel),p_binos[1:26,1,1])
plot(-diff(sealevel),log10(N_est_int[1:26,1,1]))
# Or possibly that changes in sea-levels affect changes in diversity

plot(-diff(sealevel),-diff(log10(N_est_int[,1,1])))

# using millers sealevel
plot(-diff(sealevelM),-diff(log10(N_est_int[1:18,3,1])))
fit<- lm(diff(sealevelM)~diff(log10(N_est_int[1:18,3,1])))
abline(-fit)
summary(fit)



# Can we make a big plot of 4 by 2 with 
# left panels showing detrended sea-level (through linear regression residuals) from Haq
# versus estimated richness in the different groups (detrended log10)
# and right panels use Miller sealevels.
par(mfrow=c(4,2),cex=0.5)
p_values = array(NA,c(4,2))
for (ii in 1:4){

use = which(!is.na(N_est_int[,ii,1]))
plot(residuals(lm(sealevel[use]~midpoints[use])),residuals(lm(log10(N_est_int[use,ii,1])~midpoints[use])))
fit1<- lm(residuals(lm(log10(N_est_int[use,ii,1])~midpoints[use]))~residuals(lm(sealevel[use]~midpoints[use])))
abline(fit1)
p_values[ii,1] = anova(fit1)$Pr[1]
use = which(!is.na(N_est_int[1:18,ii,1]))
plot(residuals(lm(sealevelM[use]~midpoints[use])),residuals(lm(log10(N_est_int[use,ii,1])~midpoints[use])))
fit1<- lm(residuals(lm(log10(N_est_int[use,ii,1])~midpoints[use]))~residuals(lm(sealevelM[use]~midpoints[use])))
abline(fit1)
p_values[ii,2] = anova(fit1)$Pr[1]

}


# All the comparisons of sampling to proxies should be done on the binomial rates, which are directly comparable instead

# Shall we compare binomial rates and estimated richness with
# collections (own data)
# coninental area (from Butler et al 2011)
# sealevel (Miller and Haq)

plot(collections,p_binos[,1,1],log="x")
#proxydata
proxydata = array(NA,dim=c(27,5))
proxydata[,1] = log10(collections);
proxydata[,2] = proxy1;
proxydata[1:26,3] = log10(contarea);
proxydata[,4]   = sealevel;
proxydata[1:18,5] = sealevelM;
rownames(proxydata)<-interval.names
colnames(proxydata)<-c("No collections","Dinosaur bearing formations","Continental area","Sealevel (Haq)","Sealevel (Miller)")
rownames(N_est_int)<-interval.names
colnames(N_est_int)<-Dinogroups
par(mfrow=c(4,5))
for (jj in 1:4){
for (ii in 1:5){
  # plot(proxydata[,ii],log10(N_est_int[,jj,1]))
  # plot(proxydata[,ii],p_binos[,1,1])
  use = which(!is.na(proxydata[,ii]))
  plot(residuals(lm(proxydata[use,ii]~midpoints[use])),residuals(lm(log10(N_est_int[,jj,1])~midpoints))[use])
  
}
}

## Some statistics for the text.
# Collections (often used as proxy)
# detrended log10 collections vs p_binos
fit1 <- lm(residuals(lm(proxydata[,1]~midpoints))~p_binos[,1,1])
par(mfrow=c(1,1))
plot(residuals(lm(proxydata[,1]~midpoints)),p_binos[,4,1])

cor.test(residuals(lm(proxydata[,1]~midpoints)),p_binos[,1,1])
cor.test(residuals(lm(proxydata[,1]~midpoints)),p_binos[,2,1])
cor.test(residuals(lm(proxydata[,1]~midpoints)),p_binos[,3,1])
cor.test(residuals(lm(proxydata[,1]~midpoints)),p_binos[,4,1])


# Seeing how different the 100 sampling rates are
par(mfrow=c(2,2))
for (jj in 1:4){
  if (jj==1){
    p_interval_plot=p_interval_dinos;
  }
  if (jj==2){
    p_interval_plot=p_interval_ornits;
  }
  if (jj==3){
    p_interval_plot=p_interval_sauros;
  }
  if (jj==4){
    p_interval_plot=p_interval_theros;
  }
  
  plot(midpoints,p_interval_plot[1,,1],xlim=rev(range(midpoints)),type="o")
  for (ii in 2:100){
    lines(midpoints,p_interval_plot[ii,,1],type="o")
    
  }
}
