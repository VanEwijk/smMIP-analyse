########################################################################################################################## 
# Een lijst maken met de namen van alle samples.
# 
##########################################################################################################################  
samplesLijstMaken <- function(benodigdeData){
  samples <- as.vector(colnames(benodigdeData))
  return(samples)
}

########################################################################################################################## 
# Deze functie haalt alle bestanden op die dezelfde naam hebben als in de array van de samples.
# Vervolgens worden er bewerkingen gedaan op het file:
#                                         - In de kolom  p..HGVS worden de rijen die gelijk zijn aan p.= en leeg verwijderd
#                   `                     - In de kolom coverage worden de rijen met een percentage lager dan percentageCutOff (10) verwijderd
#                                         - In de kolom Name worden de rijen met de naam ACTB verwijderd
# Nu wordt het bewerkte file opgeslagen.
# Alle mutaties/varianten met hun posities (Name, c. HGVS, p. HGVS) worden in een array gezet.`
# Vervolgens wordt er van die lijst een lijst gemaakt van alleen unieke genen en posities (mutaties).
########################################################################################################################## 
bestandBewerken <- function(samples, percentageCutOff, pathMutatieBestanden, opslaanPath){
  alleMutaties <- array()
  
  for(nummer in 1:length(samples)){
    #Mutatie files ophalen per sample.
    titelBestand <- paste(pathMutatieBestanden, samples[nummer], ".txt", sep = "")
    #De bestanden worden ingelezen.
    mutatieBestand = ldply(titelBestand,function(x) 
      if(grepl(",",readLines(x,n=1))){read.csv(x,sep=",", header = TRUE, na.strings = c("", "NA"))} 
      else if(grepl(";",readLines(x,n=1))){read.csv(x,sep=";", header = TRUE, na.strings = c("", "NA"))} 
      else{read.csv(x,sep="\t", header = TRUE, na.strings = c("", "NA"))})
    soortSample <- samples[nummer]
    
    #De index van de tabs met de volgende namen worden opgehaald.
    p.HGVS <- "p.*HGVS" 
    temp.pHGVS <- regmatches(colnames(mutatieBestand), regexpr(p.HGVS, colnames(mutatieBestand)))
    c.HGVS <- "c.*HGVS" 
    temp.cHGVS <- regmatches(colnames(mutatieBestand), regexpr(c.HGVS, colnames(mutatieBestand)))
    genNaam <- ".*Gene.*" 
    tempGenNaam <- regmatches(colnames(mutatieBestand), regexpr(genNaam, colnames(mutatieBestand)))
    coverage <- ".*Coverage.*" 
    tempCoverage <- regmatches(colnames(mutatieBestand), regexpr(coverage, colnames(mutatieBestand)))
    for(i in 1:ncol(mutatieBestand)){
      if(colnames(mutatieBestand[i]) == temp.pHGVS){
        index.pHGVS <- i
      }
      if(colnames(mutatieBestand[i]) == temp.cHGVS){
        index.cHGVS <- i
      }
      if(colnames(mutatieBestand[i]) == tempGenNaam){
        indexGenNaam <- i
      }
      if(colnames(mutatieBestand[i]) == tempCoverage){
        indexCoverage <- i
      }
    }
    
    
    mutatieBestand <- mutatieBestand[complete.cases(mutatieBestand[, index.pHGVS]), ]
    #Data verwijderen waarbij p.HGVS leeg is of gelijk is aan p.=.
    bestandMinP <- mutatieBestand[!mutatieBestand[,index.pHGVS] == "p.=",]
    bestandMinLeeg <- bestandMinP[!bestandMinP[,index.pHGVS] == "",]
    #Verwijder alle mutaties in het gen ACTB.
    bestandMinACTB <- bestandMinLeeg[!bestandMinLeeg[,indexGenNaam] == "ACTB", ]
    
    #Error opvangen dat sommige bestanden leeg zijn.
    if(nrow(bestandMinACTB) > 0){
      for(rij in 1:nrow(bestandMinACTB)){
        coverage = bestandMinACTB[,indexCoverage][rij]
        #Splitten zodat je inplaats van dit; 92% (3529)   [92% (1756) / 92% (1773)] dit krijgt; 92.
        percentageCharacter <- unlist(strsplit(as.character(coverage), "%"))[[1]][1]
        #Maak van het percentage (nu nog een character) een integer.
        percentageNumeric <- as.numeric(percentageCharacter)
        #Wanneer percentageNumeric lager is dan percentageCutOff wordt er NA neergezet.
        if(percentageNumeric < percentageCutOff){
          bestandMinACTB[,indexCoverage][rij] = NA
        }
      }
    }
    #Verwijder alle rijen die NA bevatten in de kolom coverage.
    bewerkt <- bestandMinACTB[!is.na(bestandMinACTB[,indexCoverage]),]
    #Kolomnamen worden veranderd.
    colnames(bewerkt) <- colnames(mutatieBestand)
    #Titel maken voor het bestand bewerkt.
    titelBewerktBestand <- paste(opslaanPath, "Bewerkt bestand", soortSample, percentageCutOff, sep = " ")
    #Schrijf de bewerkte bestanden weg.
    write.table(bewerkt, file = titelBewerktBestand, row.names=FALSE, na="",col.names=TRUE, sep=";")
    
    #Voeg de naam, c. HGVS en p.HGVS toe aan een lijst.
    naamMutatiePos <- paste(bewerkt[,indexGenNaam], bewerkt[,index.cHGVS], bewerkt[,index.pHGVS], sep = " ; ")
    #Zet alle alle mutaties van 1 bestand in een array.
    lijstMutaties <- unlist(strsplit(as.character(naamMutatiePos), "\t"))
    #Zet alle mutaties van alle bestanden bij elkaar in een array.
    alleMutaties <- unlist(list(alleMutaties, lijstMutaties))
  }
  #Verwijder eventuele NA in de array.
  alleMutaties <- alleMutaties[!is.na(alleMutaties)]
  #Lijst met alle unieke genen.
  uniekeMutaties <- unique(alleMutaties)
  return(uniekeMutaties)
}

