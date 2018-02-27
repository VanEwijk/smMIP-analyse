########################################################################################################################## 
# Indelen in groepen van histological type.
#
########################################################################################################################## 
scheidenHistType <- function(geheleData, opslaanPath, value1, value2, wholePath){
  #2 vectoren
  groep1 <- c()
  groep2 <- c()
  #De eerste rij van het bestand bestaat uit de soort histological type.
  #Deze worden in histologicalType opgeslagen.
  histologicalType <- geheleData[1, 2:ncol(geheleData)]
  
  for(i in 1:ncol(histologicalType)){
    if(as.character(histologicalType[,i]) == value1){
      groep1 <- c(groep1, as.character(colnames(histologicalType)[i]))
    }else if(as.character(histologicalType[,i]) == value2){
      groep2 <- c(groep2, as.character(colnames(histologicalType)[i]))
    }
  }
  
  source(paste(wholePath, "wilcoxon test over genen.R", sep = ""))
  header <- c("Gene", paste("mean cluster", value1, sep = " "), paste("mean cluster", value2, sep = " "), "p-value", "critical value", "p.adjust")
  naamWilcoxonBestand <- paste(opslaanPath, "000 Wilcoxon over genen -", value1, value2, sep = " ")
  naamMeanEenMinMeanTwee <- paste(opslaanPath, "000 mean cluster een min mean cluster twee -", value1, value2, sep = " ")
  #Roept de functie alleExpressieWaardenPerGroep aan.
  alleExpressieWaardenPerGroep(groep1, groep2, benodigdeData, opslaanPath, header, naamWilcoxonBestand, naamMeanEenMinMeanTwee)
  
  value1Value2 <- list(groep1, groep2)
  return(value1Value2)
}