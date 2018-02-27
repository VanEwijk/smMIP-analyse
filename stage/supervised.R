########################################################################################################################## 
# Maakt twee groepen: mutatie en WT. Pakt per groep alle expressie waarden.
# Returnt uiteindelijk een lijst de mutanten en wildtype : mutatieMutantWildtype.
#
########################################################################################################################## 
alleExpressieWaarden <- function(samples, benodigdeData, opslaanPath, gekozenMutatie, methodeDistanceMatrix, methodeClusteren, percentageCutOff){
  mutatie <- as.character(read.csv(file = paste(opslaanPath, gekozenMutatie, percentageCutOff, sep = " "), check.names = FALSE)[,2])
  WT <- setdiff(samples, mutatie)
  
  #Allen genen worden nu in een vector/lijst gestopt.
  genen <- rownames(benodigdeData)
  #Maak van het data.frame een nummeric matrix.
  nummericDataFrame <- apply(benodigdeData, 2, as.numeric)
  rownames(nummericDataFrame) <- rownames(benodigdeData)
  
  #Maken van dataframes en counters.
  dataFrameMutatie <- data.frame(matrix(nrow=nrow(nummericDataFrame),ncol= length(mutatie)))
  countMutatie = 1
  dataFrameWT <- data.frame(matrix(nrow=nrow(nummericDataFrame),ncol= length(WT)))
  countWT = 1
  
  #De mutanten en wildtype mogen niet leeg zijn.
  if(length(mutatie) > 0 && length(WT) > 0){
    #Loopt over alle kolommen van de benodigdeData (over alle samples).
    for(kolom in 1:length(benodigdeData)){
      #Pakt alle expressie waardes van alle samples die in cluster 1 zitten en stopt dit in een dataframe.
      #De kolommen van de samples in dit cluster worden gepakt en in de bovenstaande gemaakte dataframe gestopt.
      for(lijstMutatie in 1:length(mutatie)){
        if(colnames(benodigdeData[kolom]) == mutatie[lijstMutatie]){
          dataFrameMutatie[,countMutatie] = nummericDataFrame[,kolom]
          countMutatie <- countMutatie + 1
        }
      }
      #Pakt alle expressie waardes van alle samples die in cluster 2 zitten en stopt dit in een dataframe.
      #De kolommen van de samples in dit cluster worden gepakt en in de bovenstaande gemaakte dataframe gestopt.
      for(lijstWT in 1:length(WT)){
        if(colnames(benodigdeData[kolom]) == WT[lijstWT]){
          dataFrameWT[,countWT] = nummericDataFrame[,kolom]
          countWT <- countWT + 1
        }
      }
    }
    
    #Hernoem de rijen en kolommen.
    #De kolommen worden naar de samples vernoemd en de rijen naar de genen.
    colnames(dataFrameMutatie) <- mutatie
    rownames(dataFrameMutatie) <- genen
    colnames(dataFrameWT) <- WT
    rownames(dataFrameWT) <- genen 
    
    #Maakt een nieuw opslaan path
    opslaanP <- paste(opslaanPath, "000", gekozenMutatie, sep = " ")
    source(paste(wholePath, "clusteren van data.R", sep = ""))
    headerMutatie <- c("gene", paste("mean", gekozenMutatie, " mut"), paste("mean" , gekozenMutatie,  "wt"), "p-value wilcoxon test", "FDR", "p.adjust")
    #Voert de wilcoxon test uit over de genen.
    wilcoxonGenen(dataFrameMutatie, dataFrameWT, nummericDataFrame, benodigdeData, opslaanP, headerMutatie)
    
    #Wildtype moet groter dan 1 zijn, omdat er anders niet geclusterd kan worden.
    #Er wordt dan ook geen heatmap gemaakt.
    if(length(WT) > 1){
      # source(paste(wholePath, "clusteren van data.R", sep = ""))
      histogramWT = histogramMaken(dataFrameWT, methodeDistanceMatrix, methodeClusteren)
      dendrogramWT = dendrogramMaken(histogramWT)
      source(paste(wholePath,"heatmap zonder hist_type.R", sep = ""))
      orderOpMeanValuesWT = bewerken(dataFrameWT, opslaanP)
      heatmapMaken(orderOpMeanValuesWT, WT, dendrogramWT)
    }else{
      print(paste("Deze mutatie heeft maar 1 wildtype: ", WT,". Deze kan hierdoor niet geclusterd worden."))
    }
    #Mutant moet groter dan 1 zijn, omdat er anders niet geclusterd kan worden.
    #Er wordt dan ook geen heatmap gemaakt.
    if(length(mutatie) > 1){
      histogramMutatie = histogramMaken(dataFrameMutatie, methodeDistanceMatrix, methodeClusteren)
      dendrogramMutatie = dendrogramMaken(histogramMutatie)
      source(paste(wholePath,"heatmap zonder hist_type.R", sep = ""))
      orderOpMeanValuesMutatie = bewerken(dataFrameMutatie, opslaanP)
      heatmapMaken(orderOpMeanValuesMutatie, mutatie, dendrogramMutatie)
    }else{
      print(paste("Deze mutatie heeft maar 1 patient met deze mutatie: ", mutatie, ". Deze kan hierdoor niet geclusterd worden."))
    }
  }
  #Voegt de mutanten en wildtype samen in een lijst.
  mutatieMutantWildtype <- list(mutatie, WT)
  return(mutatieMutantWildtype)
}