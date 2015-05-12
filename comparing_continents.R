

# Trying to also get the continent from PBDB


theros <- pbdb_occurrences(limit="all", base_name="Theropoda", show=c("phylo", "ident", "time","loc"))

# Seems like the continent is not outputted, but the country codes are.
# We can use these. They are listed here;
# https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2

dinoswobird  <- pbdb_occurrences(limit="all", base_name="Dinosauria",exclude_id=36616, show=c("phylo", "ident", "time","loc"))
# 
# sum(theros$cc2=="UK",na.rm=TRUE)
# sum(theros$cc2=="US",na.rm=TRUE)
sum(dinoswobird$cc2=="CA",na.rm=TRUE)
sum(dinoswobird$cc2=="US",na.rm=TRUE)

unique(dinoswobird$cc2)
# Country codes collected for the four main continents. These are not paleocontinents!
NorAm = c("AI","AG","AW","BS","BB","BZ","BM","VG","CA","KY","CR","CU","CW","DM","DO","SV","GL","GD","GP","GT","HT","HN","JM","MQ","MX","PM","MS","CW","KN","NI","PA","PR","KN","LC","PM","VC","TT","TC","VI","US","SX","BQ","SA","SE")
Europ = c("AL",  "AD",  "AT",  "BY",  "BE",  "BA",  "BG",  "HR",  "CY",  "CZ",  "DK",  "EE",  "FO",  "FI",  "FR",  "DE",  "GI",  "GR",  "HU",  "IS",  "IE",  "IT",  "LV",  "LI",  "LT",  "LU",  "MK",  "MT",  "MD",  "MC",  "NL",  "NO",  "PL",  "PT",  "RO",  "RU",  "SM",  "RS",  "SK",  "SI",  "ES",  "SE",  "CH",  "UA",  "GB",  "VA",  "RS",  "IM",  "RS",  "ME")
Asia  = c("AF","AM","AZ","BH","BD","BT","BN","KH","CN","CX","CC","IO","GE","HK","IN","ID","IR","IQ","IL","JP","JO","KZ","KP","KR","KW","KG","LA","LB","MO","MY","MV","MN","MM","NP","OM","PK","PH","QA","SA","SG","LK","SY","TW","TJ","TH","TR","TM","AE","UZ","VN","YE","PS")
SouAm = c("AR","BO","BR","CL","CO","EC","FK","GF","GY","PY","PE","SR","UY","VE");
Ocean = c("AS","AU","NZ","CK","FJ","PF","GU","KI","MP","MH","FM","UM","NR","NC","NZ","NU","NF","PW","PG","MP","SB","TK","TO","TV","VU","UM","WF","WS","TL")
  
CCs = dinoswobird$cc2;

# Oceania
testOC <- function(ii){which(grepl(Ocean[ii],CCs))}
use = array(0,c(1,13726))
for (ii in 1:14){
  j<-testOC(ii)
  use[j] = 1;
}
use = which(use==1);
J_OC = createDataArrs(dinoswobird[use,])

# South America
testSA <- function(ii){which(grepl(SouAm[ii],CCs))}
use = array(0,c(1,13726))
for (ii in 1:14){
  j<-testSA(ii)
  use[j] = 1;
}
use = which(use==1);
J_SA = createDataArrs(dinoswobird[use,])




# Asia
testAs <- function(ii){which(grepl(Asia[ii],CCs))}
use = array(0,c(1,13726));
for (ii in 1:52){
  j<-testAs(ii)
  use[j] = 1;
}
use = which(use==1);
J_AS = createDataArrs(dinoswobird[use,])


# North America
testNA <- function(ii){which(grepl(NorAm[ii],CCs))}
use = array(0,c(1,13726));
for (ii in 1:44){
  j<-testNA(ii)
  use[j] = 1;
}
use = which(use==1);
J_NA = createDataArrs(dinoswobird[use,])






J = createDataArrs(dinoswobird)
J_all = J;
Occs <- J$Data[J$Data>0] # Occurrences 
dTs  <- J$Times[J$Data>0] # List of durations for each of these occurrences
poisrates_all <- estimatePoiss(dTs,Occs)

Occs <- J_NA$Data[J_NA$Data>0] # Occurrences 
dTs  <- J_NA$Times[J_NA$Data>0] # List of durations for each of these occurrences
poisrates_NA <- estimatePoiss(dTs,Occs)

Occs <- J_AS$Data[J_AS$Data>0] # Occurrences 
dTs  <- J_AS$Times[J_AS$Data>0] # List of durations for each of these occurrences
poisrates_AS <- estimatePoiss(dTs,Occs)


Occs <- J_SA$Data[J_SA$Data>0] # Occurrences 
dTs  <- J_SA$Times[J_SA$Data>0] # List of durations for each of these occurrences
poisrates_SA <- estimatePoiss(dTs,Occs)

Occs <- J_OC$Data[J_OC$Data>0] # Occurrences 
dTs  <- J_OC$Times[J_OC$Data>0] # List of durations for each of these occurrences
poisrates_OC <- estimatePoiss(dTs,Occs)


# Hmm, quite large differences here...

poisrates_epoch_all = array(NA,c(6,3,4)) #epoch by [mle,ci1,ci2] by [all Nortamerica asia SouthAmerica 

jj=1;
J = J_all;
for (ii in 1:6){
  tmp = J$Data[,Bins[,6]==ii]
  tmp2 = J$Time[,Bins[,6]==ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  poisrates_epoch_all[ii,,jj] <- estimatePoiss(dTs,Occs)
}

jj=2;
J = J_NA;
for (ii in 1:6){
  tmp = J$Data[,Bins[,6]==ii]
  tmp2 = J$Time[,Bins[,6]==ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  if (sum(Occs>1)>1){
  poisrates_epoch_all[ii,,jj] <- estimatePoiss(dTs,Occs)
  }
}

jj=3;
J = J_AS;
for (ii in 1:6){
  tmp = J$Data[,Bins[,6]==ii]
  tmp2 = J$Time[,Bins[,6]==ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  poisrates_epoch_all[ii,,jj] <- estimatePoiss(dTs,Occs)
}


jj=4;
J = J_SA;
for (ii in 1:6){
  tmp = J$Data[,Bins[,6]==ii]
  tmp2 = J$Time[,Bins[,6]==ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  if (sum(Occs>1)>1){
    poisrates_epoch_all[ii,,jj] <- estimatePoiss(dTs,Occs)
  }
}




# dimnames(poisrates_epoch_all)<-(epoch.names,{"MLE","ci1","ci2"},{"All","NorthAmerica","Asia"}


# Doing some species counting;
# Number of occurrences in each bin
colSums(J_NA$Data)
# NUmber of species in each bin
Spec_cont = array(NA,c(4,27))
Spec_cont[1,] <- colSums(J_all$Data>0)
Spec_cont[2,] <- colSums(J_NA$Data>0)
Spec_cont[3,] <- colSums(J_SA$Data>0)
Spec_cont[4,] <- colSums(J_AS$Data>0)
par(mfrow=c(1,1))
plot(midpoints,Spec_cont[2,],type="l",col=color.list[1],xlim=rev(range(midpoints)))
lines(midpoints,Spec_cont[3,],type="l",col=color.list[2])
lines(midpoints,Spec_cont[3,],type="l",col=color.list[3])
lines(midpoints,Spec_cont[4,],type="l",col=color.list[4])
legend("topleft",c("All","North America","South America","Asia"),pch="l",col=color.list)
# flipping the matrix to use diff
tmp <- t(J_NA$Data[nrow(J_NA$Data):1,])
# diff(tmp) gives the start and end of the different species in [spec by bin]; should add
diff(tmp[,201:210]>0)
