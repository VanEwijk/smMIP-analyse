########################################################################################################################## 
# Een lijst maken met de namen van alle samples.
#
########################################################################################################################## 
lijstPatienten <- function(groep1Groep2){
  samples <- c(groep1Groep2[[1]], groep1Groep2[[2]])
  return(samples)
}

########################################################################################################################## 
# Maakt twee groepen per mutatie: bezitMutatie en bezitNiet.
# Over deze groepen wordt een fisher exact test uitgevoerd.
#
########################################################################################################################## 
KijkenNaarMutaties <- function(samples, groep1Groep2, percentageCutOff, uniekeMutaties, kortDagen, langDagen, naam, opslaanPath, headers, wholePath){
  #Lege dataframe maken voor de informatie voor de fisher testen.
  dataFrameFisherTest <- data.frame()
  
  #Loopt over alle mutaies.
  for(mutatie in 1:length(uniekeMutaties)){
    #Een vector om de samples in te zetten die de mutatie hebben.
    bezitMutatie <- c()
    mutatiePositie <- unlist(strsplit(as.character(uniekeMutaties[mutatie]), " ; "))
    #Loopt over de samples.
    for(sample in 1:length(samples)){
      naamBestand <- paste(opslaanPath ,"Bewerkt bestand", samples[sample], percentageCutOff, sep=" ")
      bewerktBestand <- ldply(naamBestand,function(x)
        if(grepl(",",readLines(x,n=1))){read.csv(x,sep=",", header = TRUE, na.strings = c("", "NA"))}
        else if(grepl(";",readLines(x,n=1))){read.csv(x,sep=";", header = TRUE, na.strings = c("", "NA"))}
        else{read.csv(x,sep="\t", header = TRUE, na.strings = c("", "NA"))})
      #Loopt over de rijen van het bestand van het sample.
      for(rij in 1:nrow(bewerktBestand)){
        #Wanneer de gen naam en de c. HGVS overeenkomen wordt het sample in de vector bezitMutatie geplaatst.
        if(bewerktBestand[rij,2] == mutatiePositie[1] && bewerktBestand[rij, 14] == mutatiePositie[2]){
          soortSample <- samples[sample]
          bezitMutatie <- c(bezitMutatie, soortSample)
        }
      }
    }

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
    }
  }
  #Verwijder alle rijen die NA bevatten in de kolom 
  dataFrameFisherTest <- dataFrameFisherTest[!is.na(dataFrameFisherTest[,1]),]
  #dataFrameFisherTestOrder ordenen op p-value (kolom 8).
  dataFrameFisherTestOrder <- dataFrameFisherTest[order(dataFrameFisherTest[,8]), ]
  #FDR = 0.01
  #Berekenen van multiple Testing Correctie.
  teller = 1
  multipleTestingCorrectie <- c()
  for(i in 1:nrow(dataFrameFisherTestOrder)){
    nieuwewaarde <- (0.01 * teller)/ nrow(dataFrameFisherTestOrder)
    teller = teller + 1
    multipleTestingCorrectie <- c(multipleTestingCorrectie, nieuwewaarde)
  }
  #Multiple Testing Correctie toevoegen aan dataFrameFisherTestOrder
  dataFrameFisherTestOrder$FDR <- multipleTestingCorrectie
  #p.adjust uitvoeren met method = BH (Benjamini & Hochberg).
  dataFrameFisherTestOrder$p.adjust <- p.adjust(dataFrameFisherTestOrder[,8], method = "BH", n = nrow(dataFrameFisherTestOrder))
  #Headers maken voor het bestand dataFrameFisherTestOrder
  dataFrameFisherTestVolledig <- rbind(headers, dataFrameFisherTestOrder)
  #Er wordt een naam voor het file gemaakt en het file wordt weggeschreven.
  rownames(dataFrameFisherTestVolledig) <- NULL
  if(naam == "tweeVSa" || naam == "CA12" || naam == "nierNormaalTumor" || naam == "AO" || naam == "AG" || naam == "GO" || naam == "IDH" || naam == "IDH (zonder IDH1-other en IDH2)"){
    naamdataFrameFisherTestVolledig <- paste(opslaanPath, "000 Fishers Exact test over mutaties", percentageCutOff, "-",  naam, sep = " ")
  }else{
    naamdataFrameFisherTestVolledig <- paste(opslaanPath, "000 Fishers Exact test over mutaties", percentageCutOff, "-",  naam, kortDagen, langDagen, sep = " ")
  }
   write.table(dataFrameFisherTestVolledig, file = naamdataFrameFisherTestVolledig, col.names=FALSE, row.names = FALSE, sep=";")
}

