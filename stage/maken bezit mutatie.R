########################################################################################################################## 
# 
#
########################################################################################################################## 
maakBezitMutatieLijst <- function(samples, uniekeMutaties, methodeDistanceMatrix, methodeClusteren,
                                  percentageCutOff, opslaanPath, dataFrameFisherTest, mutatie, mutatiePositie, bezitMutatie){
  
  #Loopt over de samples.
  for(sample in 1:length(samples)){
    naamBestand <- paste(opslaanPath, "Bewerkt bestand", samples[sample], percentageCutOff, sep=" ")
    #Haalt het bestand alleen op als het bestand groter is dan 0.
    if(file.size(naamBestand) > 0){
      bewerktBestand <- ldply(naamBestand,function(x) 
        if(grepl(",",readLines(x,n=1))){read.csv(x,sep=",", header = TRUE, na.strings = c("", "NA"))} 
        else if(grepl(";",readLines(x,n=1))){read.csv(x,sep=";", header = TRUE, na.strings = c("", "NA"))} 
        else{read.csv(x,sep="\t", header = TRUE, na.strings = c("", "NA"))})#read.csv(file = naamBestand, header = TRUE, sep = ";")
      
      #Namen van kolommen uit mutatie bestand.
      #Uit deze kolommen kunnen gegevens worden opgehaald (eventueel)
      p.HGVS <- "p.*HGVS" #15
      temp.pHGVS <- regmatches(colnames(bewerktBestand), regexpr(p.HGVS, colnames(bewerktBestand)))
      c.HGVS <- "c.*HGVS" #14
      temp.cHGVS <- regmatches(colnames(bewerktBestand), regexpr(c.HGVS, colnames(bewerktBestand)))
      genNaam <- ".*Gene.*" #2
      tempGenNaam <- regmatches(colnames(bewerktBestand), regexpr(genNaam, colnames(bewerktBestand)))
      coverage <- ".*Coverage.*" #9
      tempCoverage <- regmatches(colnames(bewerktBestand), regexpr(coverage, colnames(bewerktBestand)))
      for(i in 1:ncol(bewerktBestand)){
        if(colnames(bewerktBestand[i]) == temp.pHGVS){
          index.pHGVS <- i
        }
        if(colnames(bewerktBestand[i]) == temp.cHGVS){
          index.cHGVS <- i
        }
        if(colnames(bewerktBestand[i]) == tempGenNaam){
          indexGenNaam <- i
        }
        if(colnames(bewerktBestand[i]) == tempCoverage){
          indexCoverage <- i
        }
      }
      
      
      #Loopt over de rijen van het bestand van het sample.
      for(rij in 1:nrow(bewerktBestand)){
        #Wanneer de gen naam en de c. HGVS overeenkomen wordt het sample in de vector bezitMutatie geplaatst.
        if(bewerktBestand[rij,indexGenNaam] == mutatiePositie[1] && bewerktBestand[rij, index.cHGVS] == mutatiePositie[2]){
          soortSample <- samples[sample]
          bezitMutatie <- c(bezitMutatie, soortSample)
          
        }
      }
      return(bezitMutatie)
    }else{
      print(paste("LEEG ", samples[sample]))
    }
  }
}

########################################################################################################################## 
# 
#
########################################################################################################################## 
aanroepenFisher <- function(samples, uniekeMutaties, methodeDistanceMatrix, methodeClusteren,
                            percentageCutOff, opslaanPath, dataFrameFisherTest, mutatie, bezitMutatie, groep1Groep2, mutatiePositie){
  print(mutatiePositie[1])
  print("")
  #Wanneer bezitMutatie leeg is hoeft hier ook geen fisher test over uit gevoerd te worden.
  #Alle samples zijn dan namelijk wildtype.
  if(length(bezitMutatie) > 0){
    #Het verschil tussen alle samples en de samples die de mutatie bezitten, hebben de mutatie niet.
    #Deze samples worden opgeslagen in bezitNiet.
    bezitNiet <- setdiff(samples, bezitMutatie)
    
    #Roept fishersExactTest aan.
    #bezitMutaties is unique omdat sommige patienten er twee keer in komen te staan.
    source(paste(wholePath,"fishers exact test.R",sep = ""))
    matrixFisher = fishersExactTest(unique(bezitMutatie), groep1Groep2)
    fTest <- fisher.test(matrixFisher)
    print(fTest)
    #Uit de resultaten van de Fisher test worden de p-values opgehaald.
    pValue <- fTest$p.value
    #Alle gegevens worden aan dataFrameFisherTest toegevoegd.
    dataFrameFisherTest[mutatie, 1] <- mutatiePositie[1] #Naam gen
    dataFrameFisherTest[mutatie, 2] <- mutatiePositie[2] #c. HGVS
    dataFrameFisherTest[mutatie, 3] <- mutatiePositie[3] #p. HGVS
    dataFrameFisherTest[mutatie, 4] <- matrixFisher[1,1] #cluster 1 mutant (aantal)
    dataFrameFisherTest[mutatie, 5] <- matrixFisher[2,1] #cluster 1 wild type (aantal)
    dataFrameFisherTest[mutatie, 6] <- matrixFisher[1,2] #cluster 2 mutant (aantal)
    dataFrameFisherTest[mutatie, 7] <- matrixFisher[2,2] #cluster 2 wild type (aantal)
    dataFrameFisherTest[mutatie, 8] <- pValue #p-value
    return(dataFrameFisherTest)
  }
  
}

