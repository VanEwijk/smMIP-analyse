########################################################################################################################## 
# Wanneer parameter = "overleven";
#     Maakt groepen op dagen. 
#     Wanneer je langer leeft/hebt geleeft dan het ingevoerde aantal dagen behoor je tot de groep2.
#
# Wanneer paramter = "groep";
#     Maakt groepen op cluster.
#
# Wanneer paramter = "CA12";
#     Maakt groepen op CA12 var 1 en CA12 var 2.
#
########################################################################################################################## 
maakGroepen <- function(benodigdeData, groep, kortDagen, langDagen, soortGroep, parameter, opslaanPath){
  #2 vectoren
  groep1 <- c()
  groep2 <- c()
  
  #Wanneer de parameter gelijkt is aan "overleven", dan maakt hij twee groepen; korte en lange overlevers.
  if(parameter == "overleven"){
    for(i in 1:nrow(groep)){
      if(as.numeric(groep$days_survived_from_surgery[i]) <= kortDagen && as.character(groep$course_status[i]) == "Dead"){
        groep1 <- c(groep1, as.character(groep$ID[i]))
      }else if(as.numeric(groep$days_survived_from_surgery[i]) > langDagen){
        groep2 <- c(groep2, as.character(groep$ID[i]))
      }
    }
    #Aanroepen van dataOphalen.
    dataOphalen(benodigdeData, groep1, groep2, soortGroep, kortDagen, langDagen, parameter, opslaanPath)
  }
  #Wanneer de parameter gelijk is aan "groep", dan maakt hij twee groepen; cluster A en cluster B (het zijn er maar twee van cluster B).
  else if(parameter == "groep"){
    for(i in 1:nrow(groep)){
      if(as.character(groep$cluster_A_of_B[i]) == "A"){
        groep1 <- c(groep1, as.character(groep$ID[i]))
      }else if(as.character(groep$cluster_A_of_B[i]) == "B"){
        groep2 <- c(groep2, as.character(groep$ID[i]))
      }
    }
    #Aanroepen van dataOphalen.
    dataOphalen(benodigdeData, groep1, groep2, soortGroep, kortDagen, langDagen, parameter, opslaanPath)
  }
  #Wanneer de parameter gelijk is aan "CA12", dan maakt hij twee groepen; CA12 var 1 en CA12 var 2.
  else if(parameter == "CA12"){
    for(i in 1:nrow(groep)){
      if(as.character(groep$CA12_variant[i]) == "mutant"){
        groep1 <- c(groep1, as.character(groep$ID[i]))
      }else if(as.character(groep$CA12_variant[i]) == "wild_type"){
        groep2 <- c(groep2, as.character(groep$ID[i]))
      }
    }
    #Aanroepen van dataOphalen.
    dataOphalen(benodigdeData, groep1, groep2, soortGroep, kortDagen, langDagen, parameter, opslaanPath)
  }
  
  #Voegt de twee groepen bij elkaar, om in andere bestanden weer mee te werken.
  #Wordt een genestenlijst.
  groep1Groep2 <- list(groep1, groep2)
  return(groep1Groep2)
}



########################################################################################################################## 
# Zet per groep de genexpressies van die groep in een dataframe.
#
########################################################################################################################## 
dataOphalen <- function(benodigdeData, groep1, groep2, soortGroep, kortDagen, langDagen, parameter, opslaanPath){
  #Allen genen worden nu in een vector/lijst gestopt.
  genen <- rownames(benodigdeData)
  
  #Maak van het data.frame een nummeric matrix .
  nummericDataFrame <- apply(benodigdeData, 2, as.numeric) 
  rownames(nummericDataFrame) <- genen
  
  #Maken van dataframes en counters.
  dfGroep1 <- data.frame(matrix(nrow=nrow(nummericDataFrame), ncol = length(groep1)))
  dfGroep2 <- data.frame(matrix(nrow=nrow(nummericDataFrame), ncol = length(groep2)))
  countGroep1 = 1
  countGroep2 = 1
  
  if(length(groep1) > 1 && length(groep2) > 1){
    #Loopt over alle kolommen van de benodigdeData (over alle samples).
    for(kolom in 1:length(benodigdeData)){
      #Pakt alle expressie waardes van alle samples die in groep 1 zitten en stopt dit in een dataframe.
      #De kolommen van de samples in dit cluster worden gepakt en in de bovenstaande gemaakte dataframe gestopt.
      for(k in 1:length(groep1)){
        if(colnames(benodigdeData[kolom]) == groep1[k]){
          dfGroep1[,countGroep1] = nummericDataFrame[,kolom]
          countGroep1 <- countGroep1 + 1
        }
      }
      #Pakt alle expressie waardes van alle samples die in groep 2 zitten en stopt dit in een dataframe.
      #De kolommen van de samples in dit cluster worden gepakt en in de bovenstaande gemaakte dataframe gestopt.
      for(l in 1:length(groep2)){
        if(colnames(benodigdeData[kolom]) == groep2[l]){
          dfGroep2[,countGroep2] = nummericDataFrame[,kolom]
          countGroep2 <- countGroep2 + 1
        }
      }
    }
    
    #Hernoem de rijen en kolommen.
    #De kolommen worden naar de samples vernoemd en de rijen naar de genen.
    colnames(dfGroep1) <- groep1
    rownames(dfGroep1) <- genen
    colnames(dfGroep2) <- groep2
    rownames(dfGroep2) <- genen
    
    #Aanroepen van wilcoxonOverGenen.
    wilcoxonOverGenen(dfGroep1, dfGroep2, nummericDataFrame, benodigdeData, soortGroep, kortDagen, langDagen, parameter, opslaanPath)
    
  }
}


########################################################################################################################## 
# Wilcoxon test over de genen tussen de groepen.
# 
########################################################################################################################## 
wilcoxonOverGenen <- function(dfGroep1, dfGroep2, nummericDataFrame, benodigdeData, soortGroep, kortDagen, langDagen, parameter, opslaanPath){
  #Lege dataframe aanmaken.
  p.valueWilcoxon <- data.frame()
  
  if((parameter == "overleven" || parameter == "CA12") && length(dfGroep1) > 1 && length(dfGroep2) > 1){
    #Loopt over alle data heen.
    for(gen in 1:nrow(nummericDataFrame)){
      #Maakt vectors van de expressie waardes per rij per groep.
      vectorGroep1 <- as.vector(as.numeric(dfGroep1[gen,]))
      vectorGroep2 <- as.vector(as.numeric(dfGroep2[gen,]))
      
      #GROEP 1
      #Verwijder de gen namen.
      groep1MinGenen <- unlist(lapply(vectorGroep1, "[[", 1))
      #Verwijder naam samples.
      groep1MinPatienten <- unname(groep1MinGenen[c(1:length(groep1MinGenen))])
      #Sorteer de genexpressies.
      sortGroep1 <- as.vector(sort(groep1MinPatienten))
      
      #GROEP2
      #Verwijder de gen namen.
      groep2MinGenen <- unlist(lapply(vectorGroep2, "[[", 1))
      #Verwijder naam samples.
      groep2MinPatienten <- unname(groep2MinGenen[c(1:length(groep2MinGenen))])
      #Sorteer de genexpressies.
      sortGroep2 <- as.vector(sort(groep2MinPatienten))
      
      #Genen toevoegen aan dataFrameWilcoxon.
      p.valueWilcoxon[gen, 1] <- rownames(benodigdeData[gen,])
      
      #Mean berekenen per cluster en deze toevoegen aan dataFrameWilcoxon.
      meanGroep1 <- mean(sortGroep1)
      meanGroep2 <- mean(sortGroep2)
      p.valueWilcoxon[gen, 2] <- meanGroep1
      p.valueWilcoxon[gen, 3] <- meanGroep2
      
      #Uitvoeren van Wilcoxon test.
      wilcoxonTest <- wilcox.test(sortGroep1, sortGroep2)
      p.valueWilcoxon[gen, 4] <- wilcoxonTest$p.value
    }
    #order_p.valueWilcoxon sorteren op p-value.
    order_p.valueWilcoxon <- p.valueWilcoxon[order(p.valueWilcoxon[,4]), ]
    #FDR = 0.01
    #Berekenen van multiple Testing Correctie.
    teller = 1
    multipleTestingCorrectie <- c()
    for(i in 1:nrow(order_p.valueWilcoxon)){
      nieuwewaarde <- (0.01 * teller)/ nrow(order_p.valueWilcoxon)
      teller = teller + 1
      multipleTestingCorrectie <- c(multipleTestingCorrectie, nieuwewaarde)
    }
    #Multiple Testing Correctie toevoegen aan order_p.valueWilcoxon.
    order_p.valueWilcoxon$FDR <- multipleTestingCorrectie
    #p.adjust uitvoeren met method = BH (Benjamini & Hochberg).
    order_p.valueWilcoxon$p.adjust <- p.adjust(order_p.valueWilcoxon[,4], method = "BH", n = nrow(order_p.valueWilcoxon))
    
    #Headers maken voor het bestand order_p.valueWilcoxon_MetHeader.
    order_p.valueWilcoxon_MetHeader <- rbind(c("gene", "mean Korte overleving", "mean Lange overleving","p-value", "FDR", "p.adjust"), order_p.valueWilcoxon)
    #Weg schrijven van order_p.valueWilcoxon_MetHeader.
    if(parameter == "CA12"){
      titel <- paste(opslaanPath, " 000 Wilcoxon over overleving ", soortGroep) #NAAM!!!
    }else{
      titel <- paste(opslaanPath, " 000 Wilcoxon over overleving ", soortGroep, " ", kortDagen, " ", langDagen)
    }
    write.table(order_p.valueWilcoxon_MetHeader, file = titel, row.names=FALSE, na="",col.names=FALSE, sep=";")
  }
  
  
  
  
  #Wanneer de parameter gelijk is aan "groep" is er een heel klein verschil tussen wat hij doet.
  #Het verschil zit hem in de mean. Bij "groep" voegt hij ook nog appart de mean van de twee mutanten uit cluster b toe.
  else if(parameter == "groep" && length(dfGroep1) > 1 && length(dfGroep2) > 1){
    #Loopt over alle data heen.
    for(gen in 1:nrow(nummericDataFrame)){
      #Maakt vectors van de expressie waardes per rij per groep.
      vectorGroep1 <- as.vector(as.numeric(dfGroep1[gen,]))
      vectorGroep2 <- as.vector(as.numeric(dfGroep2[gen,]))
      
      #GROEP 1
      #Verwijder de gen namen.
      groep1MinGenen <- unlist(lapply(vectorGroep1, "[[", 1))
      #Verwijder naam samples.
      groep1MinPatienten <- unname(groep1MinGenen[c(1:length(groep1MinGenen))])
      #Sorteer de genexpressies.
      sortGroep1 <- as.vector(sort(groep1MinPatienten))
      
      #GROEP2
      #Verwijder de gen namen.
      groep2MinGenen <- unlist(lapply(vectorGroep2, "[[", 1))
      #Verwijder naam samples.
      groep2MinPatienten <- unname(groep2MinGenen[c(1:length(groep2MinGenen))])
      #Sorteer de genexpressies.
      sortGroep2 <- as.vector(sort(groep2MinPatienten))
      
      #Genen toevoegen aan dataFrameWilcoxon.
      p.valueWilcoxon[gen, 1] <- rownames(benodigdeData[gen,])
      #Mean berekenen per groep en deze toevoegen aan dataFrameWilcoxon.
      #Omdat groep 2 maar twee waardes bevat zijn deze ook appart toegevoegd aan het bestand.
      meanGroep1 <- mean(sortGroep1)
      meanGroep2.1 <- mean(sortGroep2[1])
      meanGroep2.2 <- mean(sortGroep2[2])
      p.valueWilcoxon[gen, 2] <- meanGroep1
      p.valueWilcoxon[gen, 3] <- meanGroep2.1
      p.valueWilcoxon[gen, 4] <- meanGroep2.2
      
      #Uitvoeren van Wilcoxon test.
      wilcoxonTest <- wilcox.test(sortGroep1, sortGroep2)
      p.valueWilcoxon[gen, 5] <- wilcoxonTest$p.value
    }
    #order_p.valueWilcoxon sorteren op p-value.
    order_p.valueWilcoxon <- p.valueWilcoxon[order(p.valueWilcoxon[,5]), ]
    #FDR = 0.01
    #Berekenen van multiple Testing Correctie.
    teller = 1
    multipleTestingCorrectie <- c()
    for(i in 1:nrow(order_p.valueWilcoxon)){
      nieuwewaarde <- (0.01 * teller)/ nrow(order_p.valueWilcoxon)
      teller = teller + 1
      multipleTestingCorrectie <- c(multipleTestingCorrectie, nieuwewaarde)
    }
    #Multiple Testing Correctie toevoegen aan order_p.valueWilcoxon.
    order_p.valueWilcoxon$FDR <- multipleTestingCorrectie
    #p.adjust uitvoeren met method = BH (Benjamini & Hochberg).
    order_p.valueWilcoxon$p.adjust <- p.adjust(order_p.valueWilcoxon[,5], method = "BH", n = nrow(order_p.valueWilcoxon))
    
    #Headers maken voor het bestand order_p.valueWilcoxon_MetHeader.
    order_p.valueWilcoxon_MetHeader <- rbind(c("gene", "mean cluster A", "mean cluster B 1", "mean cluster B 2","p-value", "FDR", "p.adjust"), order_p.valueWilcoxon)
    #Weg schrijven van order_p.valueWilcoxon_MetHeader.
    titel <- paste(opslaanPath, " 000 Wilcoxon over overleving ", soortGroep)
    write.table(order_p.valueWilcoxon_MetHeader, file = titel, row.names = FALSE, na="", col.names = FALSE, sep=";")
    #write.table(order_p.valueWilcoxon_MetHeader, file = titel, row.names=FALSE, na="",col.names=FALSE, sep=";")
  }
  
}

