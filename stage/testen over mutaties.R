########################################################################################################################## 
# Loopt over de uniekeMutaties vervolgens over het aantal patienten en vervolgens over de rijen van het bestand.
# Kijkt of elke samples die bepaalde mutatie met dezelfde naam en positie heeft, door in het bestand te kijken of het voorkomt.
# Komt deze in het bestand voor dan wordt het nummer van de samples aan een vector toegevoegd.
# Vervolgens wordt er voor elke gen/positie een dendrogram gemaakt, de takken van de patienten die deze mutatie hebben worden
# rood gekleurd. (deze functie kan je aan en uit)
# Ook wordt er een bestand weggeschreven per gen/positie waarin de patientnummers staan die deze mutatie hebben.
# Er wordt ook een fisher exact test uitgevoerd over de mutaties tussen de cluster.
########################################################################################################################## 
testen <- function(samples, uniekeMutaties, dendrogram, methodeDistanceMatrix, methodeClusteren,
                   percentageCutOff, clusterEenTwee, opslaanPath, plotAanUit, groteData){
  
  #Lege dataframe maken voor de informatie voor de fisher testen.
  dataFrameFisherTest <- data.frame()
  
  lijst0en1 <- c()
  #Loopt over alle mutaies.
  for(mutatie in 1:length(uniekeMutaties)){
    bezitMutatie <- c()
    mutatiePositie <- unlist(strsplit(as.character(uniekeMutaties[mutatie]), " ; "))
    
    
    # source(paste(wholePath,"maken bezit mutatie.R",sep = ""))
    # bezitMutatie <-maakBezitMutatieLijst(samples, uniekeMutaties, methodeDistanceMatrix, methodeClusteren,
    #                                      percentageCutOff, opslaanPath, dataFrameFisherTest, mutatie, mutatiePositie, bezitMutatie)
    # # dataFrameFisherTest <- aanroepenFisher(samples, uniekeMutaties, methodeDistanceMatrix, methodeClusteren,
    #                                        percentageCutOff, opslaanPath, dataFrameFisherTest, mutatie, bezitMutatie, clusterEenTwee, mutatiePositie)
   
    
    
    # Een vector om de samples in te zetten die de mutatie hebben.
    # bezitMutatie <- c()
    # mutatiePositie <- unlist(strsplit(as.character(uniekeMutaties[mutatie]), " ; "))
    # Loopt over de samples.
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
      }else{
        print(paste("LEEG ", samples[sample]))
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
      matrixFisher = fishersExactTest(unique(bezitMutatie), clusterEenTwee)
      fisherTest <- fisher.test(matrixFisher)

      #Uit de resultaten van de Fisher test worden de p-values opgehaald.
      pValue <- fisherTest$p.value
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
    
    #Wanneer de gebruiker heeft gekozen voor AAN worden alle mutaties geplot.
    #Wanneer de gebruiker heeft gekozen voor UIT worden de mutaties niet geplot.
    if(plotAanUit == "AAN"){
      #Wanneer je alle dendrogrammen per mutatie wilt plotten moet je deze functie aanroepen
      plottenVanDendrogram(dendrogram, mutatiePositie, methodeDistanceMatrix, methodeClusteren, percentageCutOff, 
                          pValue, opslaanPath, bezitMutatie)
    }

    #dataFrameBezitMutatie = as.data.frame(bezitMutatie)
    #NulEen wordt aangeroepen
    erInOfNiet = NulEen(bezitMutatie, benodigdeData, uniekeMutaties[mutatie], samples)
    #Alle erInOfNiet lijsten worden toegevoegd aan lijst0en1
    lijst0en1 <- c(lijst0en1, erInOfNiet)
    
    #De / eruit halen anders levert dit errors op.
    mutatiePositie[mutatie] = gsub(mutatiePositie[mutatie], pattern = "/", replacement = "")
    #Bestand wegschrijven.
    naamDataframeBezitMutatie <- paste(opslaanPath, mutatiePositie[1], gsub(mutatiePositie[2], pattern = ">", replacement = ")"), percentageCutOff, sep = " ")
    write.csv(bezitMutatie, file = naamDataframeBezitMutatie)
  }
  #dataFrameFisherTest ordenen op p-value (kolom 8).
  orderDataFrameFisherTest <- dataFrameFisherTest[order(dataFrameFisherTest[,8]), ]
  
  #FDR = 0.01
  #Berekenen van multiple Testing Correctie.
  teller = 1
  multipleTestingCorrectie <- c()
  for(i in 1:nrow(orderDataFrameFisherTest)){
    nieuwewaarde <- (0.01 * teller)/ nrow(orderDataFrameFisherTest)
    teller = teller + 1
    multipleTestingCorrectie <- c(multipleTestingCorrectie, nieuwewaarde)
  }
  #Multiple Testing Correctie toevoegen aan orderDataFrameFisherTest.
  orderDataFrameFisherTest$FDR <- multipleTestingCorrectie
  #p.adjust uitvoeren met method = BH (Benjamini & Hochberg).
  orderDataFrameFisherTest$p.adjust <- p.adjust(orderDataFrameFisherTest[,8], method = "BH", n = nrow(orderDataFrameFisherTest))
  
  #Headers maken voor het bestand orderDataFrameFisherTest.
  orderDataFrameFisherTestMetHeader <- rbind(c("Gene mutatie",  "c. HGVS", "p. HGVS", "mutatie cluster 1",
                                         "wild type cluster 1", "mutatie cluster 2", "wild type cluster 2" ,"p-value Fishers extact test","FDR","p.adjust"), orderDataFrameFisherTest)
  #Er wordt een naam voor het file gemaakt en het file wordt weggeschreven.
  naamorderDataFrameFisherTestMetHeader <- paste(opslaanPath, "000 Fishers Exact test over mutaties", percentageCutOff, sep = " ")
  write.table(orderDataFrameFisherTestMetHeader, file = naamorderDataFrameFisherTestMetHeader,row.names=FALSE, na="",col.names=FALSE, sep=";")

  
  #Maakt dataframe voor de nul en een lijst0en1.
  gesplitteLijst0en1 <- split(lijst0en1, ceiling(seq_along(lijst0en1)/groteData)) #https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
  dataFrameGesplitteLijst0en1 <- as.data.frame(gesplitteLijst0en1)
  dfGesplitteLijst0en1 <- t(dataFrameGesplitteLijst0en1)
  return(dfGesplitteLijst0en1)
}


########################################################################################################################## 
# Maakt een vector met 0 en 1.
# 1 wanneer het de mutatie bevat en 0 wanneer het de mutatie niet bevat.
#
########################################################################################################################## 
NulEen <- function(bezitMutatie, benodigdeData, mutatie, samples){ 
  
  erIn <- c(mutatie)
  nietErIn <- c(mutatie)
  erInOfNiet <- c(mutatie)
  mutatieLijst <- bezitMutatie #https://stackoverflow.com/questions/3823211/convert-a-matrix-to-a-1-dimensional-array
  for(index in 1:length(samples)){
    patient <- samples[index] 
    #Er wordt een 1 op de plek van het sample geplaatst wanneer het sample de mutatie bevat.
    if(any(mutatieLijst == patient)){ #http://r.789695.n4.nabble.com/how-to-test-if-a-vector-contain-a-value-td2326843.html
      erIn <- c(erIn, samples[index])
      erInOfNiet <- c(erInOfNiet, 1)
    }else{
      #Er wordt een 0 op de plek van het sample geplaatst wanneer het sample de mutatie niet bevat.
      nietErIn <- c(nietErIn, samples[index])
      erInOfNiet <- c(erInOfNiet, 0)
    }
  }
  erInOfNietLijst <- as.list(erInOfNiet)
  return(erInOfNiet)
} 


# ########################################################################################################################## 
# # Fisher's Exact Test wordt uitgevoerd.
# #
# ########################################################################################################################## 
# fishersExactTest <- function(bezitMutatie, clusterEenTwee){
#   #Count maken die tellen hoeveel samples in cluster Een en Twee de mutaties bevatten.
#   welInEen = 0
#   welInTwee = 0
#   
#   #Gekeken wordt hoeveel samples per cluster de mutatie bevatten.
#   for(sample in 1:length(bezitMutatie)){
#     #CLUSTER 1
#     for(een in 1:length(clusterEenTwee[[1]])){
#       if(bezitMutatie[sample] == clusterEenTwee[[1]][een]){
#         welInEen = welInEen + 1
#       }
#     }
#     #CLUSTER 2
#     for(twee in 1:length(clusterEenTwee[[2]])){
#       if(bezitMutatie[sample] == clusterEenTwee[[2]][twee]){
#         welInTwee = welInTwee + 1
#       }
#     }
#   }
#   #Het verschil tussen alle samples uit cluster 1 of 2 en het aantal mutaties in dat cluster is het aantal samples dat de mutatie niet heeft.
#   nietInEen <- length(clusterEenTwee[[1]]) - welInEen
#   nietInTwee <- length(clusterEenTwee[[2]]) - welInTwee
#   
#   #Matrix maken om de fisher test later uit te voeren.
#   matrixFisher <- matrix(c(welInEen, nietInEen, welInTwee, nietInTwee), nrow = 2)
#   return(matrixFisher)
# }





########################################################################################################################## 
# Dendrogram per mutatie kleuren en wegschrijven.
# De gekleurde takken zijn de samples die de mutatie bevatten.
# De plot wordt automatisch opgeslagen en wordt niet weergegeven in R zelf.
#
########################################################################################################################## 
plottenVanDendrogram <- function(dendrogram, mutatiePositie, methodeDistanceMatrix, methodeClusteren, percentageCutOff, 
                                 pValue, opslaanPath, bezitMutatie){
  #De naam voor de plot.
  naamPlot <- paste(opslaanPath, mutatiePositie[1], gsub(mutatiePositie[2], pattern = ">", replacement = ")"), methodeClusteren, methodeDistanceMatrix, percentageCutOff, ".png", sep = " ")
  png(filename = naamPlot, width = 800, height = 600)
  par(cex = 0.9)
  #De samples in bezitMutaties worden rood (2) gekleurd.
  dendrogram %>% set("by_labels_branches_col", value = bezitMutatie, TF_values = c(2, Inf)) %>%
    plot()
  axis(2) #x en y as
  par(cex = 1.1)
  #De titel boven de plot.
  titel <- paste(mutatiePositie[1], "\t\t\t", mutatiePositie[2],"\n", mutatiePositie[3], sep =" ")
  cluster <- paste("Soort cluster: ", methodeClusteren)
  pWaarde <- paste("p-value: ", pValue)
  #De subtitel van de plot.
  subtitel <- paste(cluster, pWaarde, sep = "\t\t\t\t")
  title(main = titel)
  par(cex = 0.8)
  #x-as en y-as labeling.
  title(xlab = "sample", sub = subtitel)
  #Legenda wordt toegevoegd.
  legend("topright", legend = c("mutant", "wildtype"), fill = c(2,1))
  dev.off()
}


########################################################################################################################## 
# Maakt het file met 0 en 1 aan.
# En slaat het op als bestand met nul en een.
#
########################################################################################################################## 
maakNulEenFile <- function(benodigdeData, dfGesplitteLijst0en1, opslaanPath){
  #De eerste rij worden de rijnamen.
  rownames(dfGesplitteLijst0en1) <- dfGesplitteLijst0en1[,1]
  datframeNodig <- dfGesplitteLijst0en1[,-1]
  #De kolommen krijgen een naam.
  colnames(datframeNodig) <- colnames(benodigdeData)
  write.csv(datframeNodig, file = paste(opslaanPath, "000 Bestand met nul en een"))
}

