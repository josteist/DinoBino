

# Some tmp stuff looping over 'sampled' occurrence counts.
p_tmp = array(NA,c(100,3))
p_3_tmp = array(NA,c(100,3));
usenow = dinos;
for (rr in 1:100){

    # Replicating all collection of data.
    J = createDataArrs_v3(usenow);
    # tmpSpec[rr,]=colSums(J$Data>0);
    # tmpOccs[rr,,]=J$Data;
    Occs = J$Data;
    Times = J$Times;
    usen = Occs>0; # which species have occurrences in this stage
    os <- Occs[usen];      # occurrence counts
    dt <- Times[usen]; 
  p_tmp[rr,]=estimatePoiss(dt,os);
  p_3_tmp[rr,] = estimatePoiss(Times[Occs[,3]>0,3],Occs[Occs[,3]>0,3])
}


##
Occs = J$Data;
Times = J$Times;
usen = Occs>0; # which species have occurrences in this stage
os <- Occs[usen];      # occurrence counts
dt <- Times[usen]; 
psss =estimatePoiss(dt,os);

J = createDataArrs_v2(usenow,removedouble=TRUE)
Occs = J$Data;
Times = J$Times;
usen = Occs>0; # which species have occurrences in this stage
os <- Occs[usen];      # occurrence counts
dt <- Times[usen]; 
psss2 =estimatePoiss(dt,os);


## using the tmpOccs now defined, prob theros.
p_tmp2 = array(NA,c(100,27,3));
for (rr in 1:100){
  Occs = tmpOccs[rr,,];
  for (ss in 1:27){
    if (sum((Occs[Occs[,ss]>0,ss])>1)){
      p_tmp2[rr,ss,]  <- estimatePoiss(array(Bins[ss,3],length(Occs[Occs[,ss]>0,ss])),Occs[Occs[,ss]>0,ss])
    }
    
  }
}



## Now ddid the sampling estimate for all 100 replicated occ counts.
p_interval # is using the 'mode' occurrence counts [bin by dino by mle/ci1/2]
p_interval_dinos # sampled; 100 by 27 by 3
# _sauros and _ornits and _theros
par(mfrow=c(2,2))
plot(p_interval[,1,1],colMeans(p_interval_dinos[,,1],na.rm=TRUE))
plot(p_interval[,2,1],colMeans(p_interval_ornits[,,1],na.rm=TRUE))
plot(p_interval[,3,1],colMeans(p_interval_sauros[,,1],na.rm=TRUE))
plot(p_interval[,4,1],colMeans(p_interval_theros[,,1],na.rm=TRUE))


par(mfrow=c(2,2))
hist(p_interval[,1,1]-colMeans(p_interval_dinos[,,1],na.rm=TRUE))
hist(p_interval[,2,1]-colMeans(p_interval_ornits[,,1],na.rm=TRUE))
hist(p_interval[,3,1]-colMeans(p_interval_sauros[,,1],na.rm=TRUE))
hist(p_interval[,4,1]-colMeans(p_interval_theros[,,1],na.rm=TRUE))

## Getting 'median' rates
p_int_median = array(NA,dim=c(27,4,3))
for (ss in 1:27){
  # only using the first 99 to actually get a true median (since with 100 its the mid of two numbers)
  # tix = which(p_interval_dinos[,ss,1]==median(p_interval_dinos[1:99,ss,1],na.rm=TRUE))
  tix = which(p_interval_dinos[,ss,1]==sort(p_interval_dinos[,ss,1])[floor(length(sort(p_interval_dinos[,ss,1]))/2)])
  # if more than one hit, just use the first since they are identical
  p_int_median[ss,1,]=p_interval_dinos[tix[1],ss,]
  
  # tix = which(p_interval_ornits[,ss,1]==median(p_interval_ornits[1:99,ss,1],na.rm=TRUE))
  tix = which(p_interval_ornits[,ss,1]==sort(p_interval_ornits[,ss,1])[floor(length(sort(p_interval_ornits[,ss,1]))/2)])
  # if more than one hit, just use the first since they are identical
  p_int_median[ss,2,]=p_interval_ornits[tix[1],ss,]
  
  # tix = which(p_interval_sauros[,ss,1]==median(p_interval_sauros[1:99,ss,1],na.rm=TRUE))
  tix = which(p_interval_sauros[,ss,1]==sort(p_interval_sauros[,ss,1])[floor(length(sort(p_interval_sauros[,ss,1]))/2)])
  # if more than one hit, just use the first since they are identical
  p_int_median[ss,3,]=p_interval_sauros[tix[1],ss,]
  
  
  # tix = which(p_interval_theros[,ss,1]==median(p_interval_theros[1:99,ss,1],na.rm=TRUE))
  tix = which(p_interval_theros[,ss,1]==sort(p_interval_theros[,ss,1])[floor(length(sort(p_interval_theros[,ss,1]))/2)])
  # if more than one hit, just use the first since they are identical
  p_int_median[ss,4,]=p_interval_theros[tix[1],ss,]
  
}


dimnames(p_int_median)<-list(interval.names,Dinogroups,statnames)
p_binos = array(NA,c(27,4,3));
dimnames(p_binos)<-list(interval.names,Dinogroups,statnames)
for (ss in 1:27){
  for (dd in 1:4){
    p_binos[ss,dd,]=1-exp(-p_int_median[ss,dd,]*Bins[ss,3])
  }
}
SpecRich= array(NA,c(27,4));
dimnames(SpecRich)<-list(interval.names,Dinogroups)
SpecRich[,1] = colSums(dinos_occs>0)
SpecRich[,2] = colSums(ornits_occs>0)
SpecRich[,3] = colSums(sauros_occs>0)
SpecRich[,4] = colSums(theros_occs>0)

N_est_int = array(NA,c(27,4,3)); #bin by group by stat
dimnames(N_est_int)<-list(interval.names,Dinogroups,statnames)
for (ss in 1:27){
  for (dd in 1:4){
    N_est_int[ss,dd,1]=estimatetrue(SpecRich[ss,dd],p_binos[ss,dd,1])[1]
    N_est_int[ss,dd,2]=max(estimatetrue(SpecRich[ss,dd],p_binos[ss,dd,2]))
    N_est_int[ss,dd,3]=min(estimatetrue(SpecRich[ss,dd],p_binos[ss,dd,3]))
    
  }
}


(floor(SpecRich/p_binos[,,1]))
floor(SpecRich[,3]/p_binos[,1,1])

floor(SpecRich[,4]/p_binos[,1,1])


help(diff)
-diff(floor(SpecRich/p_binos[,,1]))/SpecRich[2:27,]

round((-diff(floor(SpecRich/p_binos[,,1]))/SpecRich[2:27,])*100)

# I guess that the previous application a disproportionate amount of dinos spanning
# bin 2 and 3 were placed in bin 2? 

