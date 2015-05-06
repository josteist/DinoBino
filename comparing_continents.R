

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



CCs = dinoswobird$cc2;
test <- function(ii){which(grepl(NorAm[ii],CCs))}
testAs <- function(ii){which(grepl(Asia[ii],CCs))}
use = array(0,c(1,13726))
for (ii in 1:44){
  j<-test(ii);
  use[j] = 1;
}

use = which(use==1)
J_NA = createDataArrs(dinoswobird[use,])
use = array(0,c(1,13726));

for (ii in 1:52){
  j<-testAs(ii)
  use[j] = 1;
}
use = which(use==1);
J_AS = createDataArrs(dinoswobird[use,])


J = createDataArrs(dinoswobird)

Occs <- J$Data[J$Data>0] # Occurrences 
dTs  <- J$Times[J$Data>0] # List of durations for each of these occurrences
poisrates_all <- estimatePoiss(dTs,Occs)

Occs <- J_NA$Data[J_NA$Data>0] # Occurrences 
dTs  <- J_NA$Times[J_NA$Data>0] # List of durations for each of these occurrences
poisrates_NA <- estimatePoiss(dTs,Occs)

Occs <- J_AS$Data[J_AS$Data>0] # Occurrences 
dTs  <- J_AS$Times[J_AS$Data>0] # List of durations for each of these occurrences
poisrates_AS <- estimatePoiss(dTs,Occs)
