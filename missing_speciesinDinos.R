
# finding the occs in dinos not in any other

use <- setdiff(dinos$oid,union(union(sauros$oid,ornits$oid),theros$oid))

mnas <- dinos[dinos$mra==3,]$oid;
missing <- intersect(use,mnas)
length(missing)

dinos[dinos$oid==missing[2],]

namemism = list(NA,c(191,1));
for (ii in 1:191){
  namemism[ii] = as.character(dinos[dinos$oid==missing[ii],7])
    # dinos[dinos$oid==missing[ii],]$mna
}

##
par(mfrow=c(2,2))
for (dd in 1:4){
  plot(p_int_median[,dd,1],Bins[,3])
  
}