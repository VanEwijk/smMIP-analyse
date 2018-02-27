########################################################################################################################## 
# De gebruiker kan een mutatie kiezen die supervised wordt geclusterd.
# Over de twee groepen (wildtype en mutanten) wordt een wilcoxon test gedaan en fisher exact test.
#
########################################################################################################################## 
gekozenMutatie <- function(uniekeMutaties, samples, benodigdeData, opslaanPath, methodeDistanceMatrix, methodeClusteren, wholePath, percentageCutOff){
  #De mutatie lijst wordt gesorteerd.
  sortUniekeMutatie <- sort(uniekeMutaties)
  #Kies een mutatie waarin u geintereseerd bent.
  indexGekozenMutatie <- menu(gsub(sortUniekeMutatie, pattern = " ; ", replacement = "  -  "), graphics = TRUE, title = "kies een mutatie die u supervised wilt clusteren:")
  splitGekozenMutatie <- unlist(strsplit(as.character(sortUniekeMutatie[indexGekozenMutatie]), " ; "))
  gekozenMutatie <- paste(splitGekozenMutatie[1], gsub(splitGekozenMutatie[2], pattern = ">", replacement = ")"), sep = " ")
  print(indexGekozenMutatie)
  print(splitGekozenMutatie)
  print(gekozenMutatie)
  # #Supervised cluster en testen uitvoeren op CA12 c.875-1_875insTGCAAGTCTGTACTGCGGCAGGACTGAGTCTGG
  # source(paste(wholePath,"supervised.R", sep = ""))
  # mutatieMutantWildtype = alleExpressieWaarden(samples, benodigdeData, opslaanPath, gekozenMutatie, methodeDistanceMatrix, methodeClusteren, percentageCutOff)
  # source(paste(wholePath,"statistiek op groepen (mutaties).R", sep = ""))
  # header <- c("Gene", "c. HGVS", "p. HGVS", paste("mutatie", gekozenMutatie, " mut"), paste("wt", gekozenMutatie, " mut"), paste("mutatie", gekozenMutatie, " wt"), paste("wt", gekozenMutatie, " wt") ,
  #                 "p-value", "FDR", "p.adjust")
  # #De mutanten en wildtype groep moeten beide minstens bestaan uit 1 sample.
  # if(length(mutatieMutantWildtype[[1]]) > 0 && length(mutatieMutantWildtype[[2]]) > 0){
  #   KijkenNaarMutaties(samples, mutatieMutantWildtype, percentageCutOff, uniekeMutaties, "-", "-" ,gekozenMutatie, opslaanPath, header)
  # }else{
  #   #Wanneer dit niet zo is wordt dit weergegeven in het Console.
  #   if(length(mutatieMutantWildtype[[2]]) < 1){
  #     print(paste("Deze mutaties", gekozenMutatie, "heeft alleen mutanten.", sep = " "))
  #   }else{
  #     print(paste("Deze mutaties", gekozenMutatie, "heeft alleen wildtype", sep = " "))
  #   }
  # }
  if(indexGekozenMutatie != 0){
    #Supervised cluster en testen uitvoeren op CA12 c.875-1_875insTGCAAGTCTGTACTGCGGCAGGACTGAGTCTGG
    source(paste(wholePath,"supervised.R", sep = ""))
    mutatieMutantWildtype = alleExpressieWaarden(samples, benodigdeData, opslaanPath, gekozenMutatie, methodeDistanceMatrix, methodeClusteren, percentageCutOff)
    source(paste(wholePath,"statistiek op groepen (mutaties).R", sep = ""))
    header <- c("Gene", "c. HGVS", "p. HGVS", paste("mutatie", gekozenMutatie, " mut"), paste("wt", gekozenMutatie, " mut"), paste("mutatie", gekozenMutatie, " wt"), paste("wt", gekozenMutatie, " wt") ,
                "p-value", "FDR", "p.adjust")
    #De mutanten en wildtype groep moeten beide minstens bestaan uit 1 sample.
    if(length(mutatieMutantWildtype[[1]]) > 0 && length(mutatieMutantWildtype[[2]]) > 0){
      KijkenNaarMutaties(samples, mutatieMutantWildtype, percentageCutOff, uniekeMutaties, "-", "-" ,gekozenMutatie, opslaanPath, header)
    }else{
      #Wanneer dit niet zo is wordt dit weergegeven in het Console.
      if(length(mutatieMutantWildtype[[2]]) < 1){
        print(paste("Deze mutaties", gekozenMutatie, "heeft alleen mutanten.", sep = " "))
      }else{
        print(paste("Deze mutaties", gekozenMutatie, "heeft alleen wildtype", sep = " "))
      }
    }
  }else{
    print("Helaas is er iets fout gegaan.")
    print("kies een nieuwe")
    gekozenMutatie(uniekeMutaties, samples, benodigdeData, opslaanPath, methodeDistanceMatrix, methodeClusteren, wholePath, percentageCutOff)
  }
}