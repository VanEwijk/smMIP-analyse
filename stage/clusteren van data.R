########################################################################################################################## 
# Deze functie zorgt dat er over alle expressiedata een log2 wordt gedaan.
# Ook wordt er een  distance matrix en histogram gemaakt.
# 
########################################################################################################################## 
histogramMaken <- function(benodigdeData, methodeDistanceMatrix, methodeClusteren, soortBestand){
  #De data wordt nummeric gemaakt.
  nummericDataFrame <- apply(benodigdeData, 2, as.numeric)
  #De namen van de genen worden weer de namen van de rijen.
  #rownames(nummericDataFrame) <- rownames(benodigdeData)
  colnames(nummericDataFrame) <- colnames(benodigdeData)
  #Bij alle expressie waarden wordt 0.01 opgeteld. Dit om te voorkomen dat er een log wordt genomen van 0.
  plusGegevens <- nummericDataFrame + 0.01
  #Neem de log2 van alle gen expressies
  logGegevens <- log2(plusGegevens)
  
  # library(ggfortify)
  # logg <- autoplot(prcomp(t(logGegevens)), label = TRUE)
  # #numericc <-autoplot(prcomp(nummericDataFrame))
  # plot(logg)
  # #plot(numericc)
  
  
  
  # omdraaiLogGegevens <- t(logGegevens)
  # na.exclude(omdraaiLogGegevens)
  # View(omdraaiLogGegevens)
  # library(cluster)
  # a <- autoplot(pam(omdraaiLogGegevens, 2), frame = TRUE, frame.type = 'norm', label = TRUE)
  # plot(a)
  # b <- autoplot(fanny(omdraaiLogGegevens, 2), frame = TRUE)
  # plot(b)
  
  # 
  # AA <- autoplot(prcomp(omdraaiLogGegevens), 
  #          data = omdraaiLogGegevens, 
  #          colour = 'blue', 
  #          label = TRUE
  # )
  # plot(AA)
  
  
  #Maak een distance matrix. (Met t() draai je de rijen en kolommen om.)
  distanceMatrix <- dist(t(logGegevens), method = methodeDistanceMatrix)
  #Clusteren van de data
  histogram <- hclust(distanceMatrix, method = methodeClusteren)
  #plotten van histogram
  plot(histogram)
  #title(main = "Histogram")
  return(histogram)
}


########################################################################################################################## 
# Deze functie maakt van het histogram een dendrogram
#
########################################################################################################################## 
dendrogramMaken <- function(histogram, soortBestand){
  #dendrogram maken
  dendrogram <- as.dendrogram(histogram)
  plot(dendrogram)
  #title(main = "Dendrogram")
  return(dendrogram)
}



########################################################################################################################## 
# Deze functies pak de clusters uit de histogram.
# 
########################################################################################################################## 
getClusters <- function(histogram, benodigdeData, opslaanPath, wholePath){
  #Maakt een dataframe met alle clusters.
  #k=2 betekend dat je het histogram in twee clusters verdeeld.
  clusters <- data.frame(Cluster_ID = cutree(histogram, k = 2))
  #Maakt een lijst van alle samples (patientenNummers).
  samples <- unlist(strsplit(rownames(clusters[1]), " "))
  
  #3 vectoren
  clusterEen <- c()
  clusterTwee <- c()
  cluster <- c()
  
  #Wanneer de nummers van het cluster overeenkomen worden ze in een vector gestopt.
  for(rij in 1:nrow(clusters)){
    if(clusters[rij, ] == 1){
      clusterEen <- c(clusterEen, samples[rij])
    }else if(clusters[rij, ] == 2){
        clusterTwee <- c(clusterTwee, samples[rij])
    }else{
        cluster <- c(cluster, samples[rij])
    }
  }
  
  #Cluster 1 en 2 worden geprint
  print("CLUSTER 1:")
  print(clusterEen)
  print("########################################################################################")
  print("CLUSTER 2:")
  print(clusterTwee)
  
  
  headerEenTwee <- c("gene", "mean cluster een", "mean cluster twee", "p-value wilcoxon test", "FDR", "p.adjust")
  naamWilcoxonBestand <- paste(opslaanPath, "000 wilcoxon test over genen - Clusters", sep = " ")
  naamMeanEenMinMeanTwee <- paste(opslaanPath, "000 mean cluster een min mean cluster twee - Clusters", sep = " ")
  
  source(paste(wholePath, "wilcoxon test over genen.R", sep = ""))
  #Roept de functie alleExpressieWaardenPerGroep aan.
  alleExpressieWaardenPerGroep(clusterEen, clusterTwee, benodigdeData, opslaanPath, headerEenTwee, naamWilcoxonBestand, naamMeanEenMinMeanTwee)
  #Voegt de twee cluster bij elkaar, om in andere bestanden weer mee te werken.
  #Wordt een genestenlijst.
  clusterEenTwee <- list(clusterEen, clusterTwee)
  return(clusterEenTwee)
}



