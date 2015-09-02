# Trying to disentangle the 'common cause', i.e. that there is a common driver of both sampling
# and diversity. 

# All the comparisons of sampling to proxies should be done on the binomial rates, which are directly comparable instead

# Shall we compare binomial rates and estimated richness with
# collections (own data)
# coninental area (from Butler et al 2011)
# sealevel (Miller and Haq)

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


# Continental area, sealevel can drive either sampling or diversities.
fit1<-lm(p_binos[,1,1]~midpoints)
summary(fit1) # no trend
fit2<-lm(p_binos[,2,1]~midpoints)
summary(fit2) # no trend
fit3<-lm(p_binos[,3,1]~midpoints)
summary(fit3) # no trend
fit4<-lm(p_binos[,4,1]~midpoints)
summary(fit4) # no trend, i.e.e no need to detrend.

par(mfrow=c(4,2))
useproxy = 5; # index into proxydata
for (ii in 1:4){
  resproxy = residuals(lm(proxydata[,useproxy]~midpoints))
  usprox = which(!is.na(proxydata[,useproxy]))
  use = which(!is.na(N_est_int[usprox,ii,1]))
  
  speclog10 <- log10(N_est_int[use,ii,1]);
  # residuals from a linear detrending of log10 richnesses estimated
  resdiv  = residuals(lm(speclog10~midpoints[use]))
  plot(resproxy[use],resdiv)
  plot(resproxy[use],p_binos[use,ii,1])
       
}

