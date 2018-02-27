########################################################################################################################## 
# Maakt de benodigde data voor alle bestanden waarbij soort is 'histological'.
#
########################################################################################################################## 
benodigdeDataMakenHist <- function(geheleData){
  #De eerste rij wordt verwijderd, omdat hier de histological types staan
  benodigdeData <- geheleData[-1,]
  #De genen (eerste kolom) worden nu de rij namen.
  rownames(benodigdeData) <- benodigdeData[,1]
  #De genen (eerste kolom) worden verwijderd, omdat ze nu al rij namen zijn.
  benodigdeData <- benodigdeData[,-1]
  return(benodigdeData)
}

########################################################################################################################## 
# Maakt de benodigde data voor alle bestanden waarbij er geen histological type is aangegeven.
#
########################################################################################################################## 
benodigdeDataMakenZonderHistType <- function(geheleData){
  #De genen (eerste kolom) worden nu de rij namen.
  rownames(geheleData) <- geheleData[,1]
  #De genen (eerste kolom) worden verwijderd, omdat ze nu al rij namen zijn.
  benodigdeData <- geheleData[,-1]
  return(benodigdeData)
}

