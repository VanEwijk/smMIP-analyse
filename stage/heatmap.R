########################################################################################################################## 
# Maakt het bestand om de heatmap mee te maken. En zorgt dat de genen op de goede volgorde worden gezet.
#
########################################################################################################################## 
maakFileOrder <- function(benodigdeData, opslaanPath){
  #Haalt de meanValues op om hem later op deze volgorde te sorteren.
  naamDataFrameMeanEenMinMeanTwee <- paste(opslaanPath, "000 mean cluster een min mean cluster twee - Clusters", sep = " ")
  meanValues <- read.csv(file = naamDataFrameMeanEenMinMeanTwee, header = FALSE, check.names = FALSE, sep = ";")[,2]
  
  #Maak van het data.frame een nummeric matrix .
  nummericBenodigdeData <- apply(benodigdeData, 2, as.numeric)
  
  #Bij alle expressie waarden wordt 0.01 bij opgeteld. Dit om te voorkomen dat er een log wordt genomen van 0.
  dataPlus <- nummericBenodigdeData + 0.01
  #Neem de log2 van alle gen expressies
  dataLog <- log2(dataPlus)
  
  rownames(dataLog) <- rownames(benodigdeData)
  
  #Data omdraaien (kolommen, rijen) en de meanValues er aan plakken
  draaiomDataLog <- t(dataLog)
  matrixMetMeanValues <- rbind(draaiomDataLog, meanValues)
  #omdraaien
  draaiomMatrixMetMeanValues <- t(matrixMetMeanValues)
  orderDraaiomMatrixMetMeanValues <- draaiomMatrixMetMeanValues[order(draaiomMatrixMetMeanValues[,ncol(draaiomMatrixMetMeanValues)]), ]
  #Haal de mean Values nu weer weg (verwijder laatste kolom)
  orderOpMeanValues <- orderDraaiomMatrixMetMeanValues[,-(ncol(draaiomMatrixMetMeanValues))]
  return(orderOpMeanValues)
}

########################################################################################################################## 
# Maakt verschillende heatmaps die je vervolgens kan samenvoegen tot 1 in paint of een andere programma.
#
########################################################################################################################## 
heatmapMaken <- function(orderOpMeanValues, samples, dendrogram, kleurenTabel, kleurenWtMt, soortBestand, statusHistological){
  dend <- dendrogram %>% set("branches_lwd", 1)
  
  #Wanneer soortBestand gelijk is aan "N 66" dan moet hij in cluster B de twee grootste cluster omdraaien.
  if(soortBestand == "N 66"){
    dend <- dend
    tmp <- dend[[2]][[2]]  
    dend[[2]][[2]] <- dend[[2]][[1]] 
    dend[[2]][[1]] <- tmp; 
  }
  
  if(statusHistological == "histological"){
    #Wanneer soortBestand gelijk is aan "N 66" of "N 75" moeten er twee balen komen.
    #Een voor IDH type en de andere voor histological type.
    if(soortBestand == "N 66" || soortBestand == "N 75"){
      kleurenVoorBalkjes = cbind(kleurenTabel, kleurenWtMt)
      colnames(kleurenVoorBalkjes) = c("histological type", "IDH status")
      rownames(kleurenVoorBalkjes) = samples
      
      #Heatmap voor de balkjes
      heatmapBalkjes <-heatmap.plus(orderOpMeanValues, Rowv = FALSE, Colv = dend ,
                                    margins = c(5.5,5), col = bluered(200),
                                    cexRow = 0.6, cexCol = 0.8, scale = "row",
                                    xlab = "patienten", ylab = "ROI",
                                    ColSideColors = kleurenVoorBalkjes)
      
      
      #Heatmap voor de legenda
      legendaVoorAfbeelding <-heatmap.2(orderOpMeanValues, Rowv = FALSE, Colv = dend ,col = bluered,
                                        trace = "none", density = "none", margins = c(5.5,5), sepcolor="grey", cexRow = 0.6,
                                        cexCol = 0.8, scale = "row", dendrogram = "column", xlab = "sample", ylab = "gene",
                                        colsep=1:ncol(orderOpMeanValues), rowsep=1:nrow(orderOpMeanValues),
                                        keysize = 1 ,ColSideColors = kleurenWtMt)
      # legend("topright", legend = c("IDH1-R132H", "IDH1-Y183C", "IDH1-V178I", "IDH1-R132H / IDH1-V178I", "IDH2",
      #                               "IDH-WT", "astrocytoma", "oligodendroglioma",
      #                               "glioblastoma", "metastasis", "ependymoma", "Variant glioma", "Ganglioglioma",
      #                               "Lymphoproliferative disease", "Dysembryoplastic neuroepethelial tumor", "other histologicial type"),
      #        fill = c("#FF0000", "#00994C", "#FF9933", "#B2FF66", "#000000",
      #                 "#3333FF", "#EFE93E","#C0C0C0", "#666600", "#CC6600"
      #                 , "#B5EF89", "#9933FF", "#99FFFF", "#FF99CC", "#660000", "#FFFFFF"))
      legend("topright", legend = c("IDH1-R132H", "IDH1-other", "IDH2-MT", 
                                    "IDH-WT", "astrocytoma", "oligodendroglioma",
                                    "glioblastoma", "metastasis", "ependymoma", "Variant glioma", "Ganglioglioma",
                                    "Lymphoproliferative disease", "Dysembryoplastic neuroepethelial tumor", "other histologicial type"),
             fill = c("#FF0000", "#FFFFFF", "#000000", 
                      "#3333FF", "#EFE93E","#C0C0C0", "#666600", "#CC6600"
                      , "#B5EF89", "#9933FF", "#99FFFF", "#FF99CC", "#660000", "#FFFFFF"))
      
      
    }else{
      #HeatMap met legenda.
      legendaHist <- heatmap.2(orderOpMeanValues, Rowv = FALSE, Colv = dend ,col = bluered,
                               trace = "none", density = "none", margins = c(5.5,5), sepcolor="grey", cexRow = 0.6,
                               cexCol = 0.8, scale = "row", dendrogram = "column", xlab = "sample", ylab = "gene",
                               colsep=1:ncol(orderOpMeanValues), rowsep=1:nrow(orderOpMeanValues),
                               keysize = 1 ,ColSideColors = kleurenTabel)
      legend("topright", legend = c("astrocytoma", "oligodendroglioma","glioblastoma", "metastasis", "ependymoma", 
                                    "Variant glioma", "Ganglioglioma",  "Lymphoproliferative disease", 
                                    "Dysembryoplastic neuroepethelial tumor", "other histologicial type"),
             fill = c("#EFE93E","#C0C0C0", "#666600", "#CC6600", "#B5EF89", "#9933FF", "#99FFFF", "#FF99CC", "#660000", "#FFFFFF"))
      
    }
    #Heatmap die gebruikt gaat worden
    heatmapKleuren <-heatmap.2(orderOpMeanValues, Rowv = FALSE, Colv = dend ,col = bluered,
                               trace = "none", density = "none", margins = c(5.5,5), sepcolor="grey", cexRow = 0.6,
                               cexCol = 0.8, scale = "row", dendrogram = "column", xlab = "sample", ylab = "gene",
                               colsep=1:ncol(orderOpMeanValues), rowsep=1:nrow(orderOpMeanValues),
                               keysize = 1 ,ColSideColors = kleurenTabel)
  }else{
    #Heatmap die gebruikt gaat worden
    heatmapKleuren <-heatmap.2(orderOpMeanValues, Rowv = FALSE, Colv = dend ,col = bluered,
                               trace = "none", density = "none", margins = c(7.5,5.5), sepcolor="grey", cexRow = 0.6,
                               cexCol = 0.7, srtCol = 50, scale = "row", dendrogram = "column", xlab = "sample", ylab = "gene",
                               colsep=1:ncol(orderOpMeanValues), rowsep=1:nrow(orderOpMeanValues),
                               keysize = 1)
  }

  
  
  
}


# ########################################################################################################################## 
# # Maakt een heatmap.
# #
# ########################################################################################################################## 
# heatmapMaken<- function(orderOpMeanValues, samples, dendrogram){
#   dend <- dendrogram %>% set("branches_lwd", 1)
#   
#   heatmap <- heatmap.2(orderOpMeanValues, Rowv = FALSE, Colv = dend, col = bluered,
#                        trace = "none", density = "none", margins = c(7,5), sepcolor = "grey",
#                        cexRow = 0.4, cexCol = 0.5, scale = "row", dendrogram = "column", 
#                        xlab = "sample", ylab = "gene", colsep = 1:ncol(orderOpMeanValues),
#                        sepwidth = c(0.0005,0.005), keysize = 1)
# }