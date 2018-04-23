########################################################################################################################## 
# Zet per groep de genexpressies van dat groep in een dataframe.
# 
########################################################################################################################## 
alleExpressieWaardenPerGroep <- function(groepEen, groepTwee, benodigdeData, opslaanPath, headerEenTwee, naamWilcoxonBestand, naamMeanEenMinMeanTwee){
  #Allen genen worden nu in een vector/lijst gestopt.
  genen <- rownames(benodigdeData)
  #Maak van het data.frame een nummeric matrix.
  nummericDataFrame <- apply(benodigdeData, 2, as.numeric)
  rownames(nummericDataFrame) <- rownames(benodigdeData)
  
  #Maken van dataframes en counters.
  dataFrameGroepEen <- data.frame(matrix(nrow=nrow(nummericDataFrame),ncol= length(groepEen)))
  countEen = 1
  dataFrameGroepTwee <- data.frame(matrix(nrow=nrow(nummericDataFrame),ncol= length(groepTwee)))
  countTwee = 1
  
  
  #Loopt over alle kolommen van de benodigdeData (over alle samples).
  for(kolom in 1:length(benodigdeData)){
    #Pakt alle expressie waardes van alle samples die in groep 1 zitten en stopt dit in een dataframe.
    #De kolommen van de samples in ddeze groep worden gepakt en in de bovenstaande gemaakte dataframe gestopt.
    for(lijstEen in 1:length(groepEen)){
      if(colnames(benodigdeData[kolom]) == groepEen[lijstEen]){
        dataFrameGroepEen[,countEen] = nummericDataFrame[,kolom]
        countEen <- countEen + 1
      }
    }
    #Pakt alle expressie waardes van alle samples die in groep 2 zitten en stopt dit in een dataframe.
    #De kolommen van de samples in deze groep worden gepakt en in de bovenstaande gemaakte dataframe gestopt.
    for(lijstTwee in 1:length(groepTwee)){
      if(colnames(benodigdeData[kolom]) == groepTwee[lijstTwee]){
        dataFrameGroepTwee[,countTwee] = nummericDataFrame[,kolom]
        countTwee <- countTwee + 1
      }
    }
  }
  
  #Hernoem de rijen en kolommen.
  #De kolommen worden naar de samples vernoemd en de rijen naar de genen.
  colnames(dataFrameGroepEen) <- groepEen
  rownames(dataFrameGroepEen) <- genen
  colnames(dataFrameGroepTwee) <- groepTwee
  rownames(dataFrameGroepTwee) <- genen
  
  #Roept de functie wilcoxonGenen aan.
  wilcoxonGenen(dataFrameGroepEen, dataFrameGroepTwee, nummericDataFrame, benodigdeData, opslaanPath, headerEenTwee, naamWilcoxonBestand, naamMeanEenMinMeanTwee)
}


########################################################################################################################## 
# Wilcoxon test over de genen tussen de groepen
# 
########################################################################################################################## 
wilcoxonGenen <- function(dataFrameGroepEen, dataFrameGroepTwee, nummericDataFrame, benodigdeData, opslaanPath, headers, naamWilcoxonBestand, naamMeanEenMinMeanTwee){
  #Twee lege dataframes worden aangemakat.
  dataFrameWilcoxon <- data.frame()
  dataFrameMeanEenMinMeanTwee <- data.frame()
  #Loopt over alle data heen.
  for(gen in 1:nrow(nummericDataFrame)){
    #Maakt vectors van de expressie waardes per rij per groep.
    expressiesGroepEen <- as.vector(dataFrameGroepEen[gen,])
    expressiesGroepTwee <- as.vector(dataFrameGroepTwee[gen,])
    
    #Genen toevoegen aan dataFrameWilcoxon.
    dataFrameWilcoxon[gen, 1] <- rownames(benodigdeData[gen,])
    
    #GROEP 1
    #Verwijder de gen namen.
    groepEenMinGen <- unlist(lapply(expressiesGroepEen, '[[', 1))
    #Verwijder naam samples.
    groepEenMinSamples <- unname(groepEenMinGen[c(1:length(groepEenMinGen))])
    #Sorteer de genexpressies.
    sortEen <- as.vector(sort(groepEenMinSamples))
    
    #GROEP 2
    #Verwijder de gen namen,
    groepTweeMinGen <- unlist(lapply(expressiesGroepTwee, '[[', 1))
    #Verwijder naam samples.
    groepTweeMinSamples <- unname(groepTweeMinGen[c(1:length(groepTweeMinGen))])
    #Sorteer de genexpressies.
    sortTwee <- as.vector(sort(groepTweeMinSamples))

    
    #Mean berekenen per groep en deze toevoegen aan dataFrameWilcoxon.
    meanEen <- mean(sortEen)
    meanTwee <- mean(sortTwee)
    dataFrameWilcoxon[gen, 2] <- meanEen
    dataFrameWilcoxon[gen, 3] <- meanTwee
    
    #Uitvoeren van Wilcoxon test.
    if(length(sortEen) > 0 && length(sortTwee) > 0){
      wilcoxon <- wilcox.test(sortEen, sortTwee)
      dataFrameWilcoxon[gen, 4] <- wilcoxon$p.value
    }else{
      dataFrameWilcoxon[gen, 4] <- "-"
    }
    
    
    #Voor het sorteren van de Heatmap later worden de genen en means ook toegevoegd aan dataFrameMeanEenMinMeanTwee.
    meanEenMinMeanTwee <- meanEen - meanTwee
    dataFrameMeanEenMinMeanTwee[gen, 1] <- rownames(benodigdeData[gen,])
    dataFrameMeanEenMinMeanTwee[gen, 2] <- meanEenMinMeanTwee
  }
  #dataFrameMeanEenMinMeanTwee wegschrijven, om later op te halen bij het maken van de heatmap .
  write.table(dataFrameMeanEenMinMeanTwee, file = naamMeanEenMinMeanTwee,  row.names=FALSE, na="",col.names=FALSE, sep=";")
  
  #dataFrameWilcoxon sorteren op p-value.
  orderDataFrameWilcoxon <- dataFrameWilcoxon[order(dataFrameWilcoxon[,4]), ]
  
  #FDR = 0.01
  #Berekenen van multiple Testing Correctie.
  teller = 1
  multipleTestingCorrectie <- c()
  for(i in 1:nrow(orderDataFrameWilcoxon)){
    nieuwewaarde <- (0.01 * teller)/ nrow(orderDataFrameWilcoxon)
    teller = teller + 1
    multipleTestingCorrectie <- c(multipleTestingCorrectie, nieuwewaarde)
  }
  #Multiple Testing Correctie toevoegen aan orderDataFrameWilcoxon.
  orderDataFrameWilcoxon$FDR <- multipleTestingCorrectie
  #p.adjust uitvoeren met method = BH (Benjamini & Hochberg).
  orderDataFrameWilcoxon$p.adjust <- p.adjust(orderDataFrameWilcoxon[,4], method = "BH", n = nrow(orderDataFrameWilcoxon))
  
  
  #Headers maken voor het bestand orderDataFrameWilcoxon.
  orderDataFrameWilcoxonVolledig <- rbind(headers, orderDataFrameWilcoxon)
  #Weg schrijven van naamOrderDataFrameWilcoxonVolledig.
  write.table(orderDataFrameWilcoxonVolledig, file = naamWilcoxonBestand, row.names=FALSE, na="",col.names=FALSE, sep=";")
}

