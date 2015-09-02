# Toying around with collection counts, how many occurrences/entries per collection over time?

unqcid<- unique(dinos$cid)
# so there are 8684 collections
# collecting interval of each collection
Collections <- array(0,27,1);
keepladinian <- array(0,5,1)
tix = 1;
for (cc in 1:8684){
  # Some errors here, since some data have strange ein/lin indexes; e.g. Late Pliocene to Early Pliocene is
  # 96...740, which will span our range. They will only be counted as a collection if the first
  # interval in the range is inside our bins.
dinos[dinos$cid==unqcid[cc],]
  sum(dinos$cid==unqcid[cc])
  range <- seq(min(dinos[dinos$cid==unqcid[cc],]$lin),max(dinos[dinos$cid==unqcid[cc],]$ein,1));
  if (any(range[1]==Bins[,4])){
  for (ss in 1:length(range)){
    Collections[range[ss]-111] <- Collections[range[ss]-111]+1
  }
    
  }
  if (any(range==138)){
    keepladinian[tix] = cc;
    tix = tix+1;
    
  }
}
nocol = cbind(Collections)
dimnames(nocol)<-list(interval.names,"No collections")


# mean occurrence per species. if 1 we cannot find the sampling rate
plot(log(nocol),colSums(J2$Data)/colSums(J2$Data>0))
# I would expect a stronger correlation here. Perhaps this is good, it means that we manage to capture sampling rates 
# irrespective of number of collections?


coloccs = [];
tix = 1;
tmp <-
  for (ii in 1:length(keepladinian){
    tmp<- unqcid[keepladinian[ii]];
    
dinos[tmp,]
## Hmm I suspect the create Dattaarrshave had a bug excluding finds in the ladinian. It should of course not have that.

J1 = createDataArrs(dinos)
J_test = createDataArrs_v2(dinos)