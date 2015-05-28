# Seeing if we can get fossil collection data out.
# cid is collection id.

j <- unique(dinoswobird$cid)

## 8648 collection giving all the 18796 occurrences.
# dinoswobird$mna[dinos$cid==j[2]]

tmp <- dinoswobird$cid==j[2]

j1<-dinoswobird[tmp,]$ein
j2<-dinoswobird[tmp,]$lin
seq(j1[1],j2[1],-1)
seq(dinoswobird[tmp[1],]$ein,dinoswobird[tmp[1]]$lin)


## looping over all unique collections in intervals 112:138
nocol = array(0,27,1)
j <- unique(dinoswobird$cid)
for (ii in 1:length(j)){
  tmp <- dinoswobird$cid==j[ii]
  
  j1<-dinoswobird[tmp,]$ein
  j2<-dinoswobird[tmp,]$lin
  range <-seq(j1[1],j2[1],-1)  
    if (any(range>112 & range<138)){
      nocol[139-range] = nocol[139-range]+1
      # one collection can span several intervals and count in all
      
    }
  
}

nocol = rev(nocol)
# Total diversity
plot(midpoints,(nocol),ylim=c(10,1400),xlim=rev(range(midpoints)),type="l",log="y")

lines(midpoints,MeanSpec[,1,2],log="y")

plot(log10((nocol)),log10(MeanSpec[,1,2]))
plot(log10(sort(nocol[1:26])+1),log10(sort(MeanSpec[1:26,1,2])+1))



##

dinos  <- pbdb_occurrences(limit="all", base_name="Dinosauria",interval="Mesozoic", show=c("phylo", "ident", "time"))
# can use created_befor=2Y 
# will give occurrences created before 2 years ago.

dinos2 <- pbdb_occurrences(limit="all", base_name="Dinosauria",interval="Mesozoic",created_before="5Y", show=c("phylo", "ident", "time"))
dinos1 <- pbdb_occurrences(limit="all", base_name="Dinosauria",interval="Mesozoic",created_after="5Y", show=c("phylo", "ident", "time"))

# This seem to work, and the sum of occurrences before and after 5Y is the sum of all.
# But why does the the ones before 5Y have one more entry


# So are the rates different?
J1 <- createDataArrs_v2(dinos1)
Occs <- J1$Data[J1$Data>0] # Occurrences 
dTs  <- J1$Times[J1$Data>0] # List of durations for each of these occurrences
poisrates1 <- estimatePoiss(dTs,Occs)

J2 <- createDataArrs_v2(dinos2)
Occs <- J2$Data[J2$Data>0] # Occurrences 
dTs  <- J2$Times[J2$Data>0] # List of durations for each of these occurrences
poisrates2 <- estimatePoiss(dTs,Occs)

