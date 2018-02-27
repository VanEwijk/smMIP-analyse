########################################################################################################################## 
# 
#
########################################################################################################################## 

scheidenNieren <- function(geheleData, opslaanPath, wholePath){
  print("doet het")
  
  #2 vectoren
  groep1 <- c()
  groep2 <- c()
  #De eerste rij van het bestand bestaat uit de soort histological type.
  #Deze worden in histologicalType opgeslagen.
  histologicalType <- geheleData[1, 2:ncol(geheleData)]
  
  
 
  #Patroon maken.
  patroonTumor <- "_Tumor"
  #Pak de index wanneer de kolomnaam hetzelfde is als het patroonTumor
  colIndexTagTumor <- grep(patroonTumor, colnames(histologicalType))
  #Pak de naam van de kolom wanneer de kolomnaam hetzelfde is als het patroonTumor
  colNaamTumor <- grep(patroonTumor, colnames(histologicalType), value = TRUE)
  
  #Patroon maken.
  patroonNormal <- "_normal"
  #Pak de index wanneer de kolomnaam hetzelfde is als het patroonNormal
  colIndexTagNormal <- grep(patroonNormal, colnames(histologicalType))
  #Pak de naam van de kolom wanneer de kolomnaam hetzelfde is als het patroonNormal
  colNaamNormal <- grep(patroonNormal, colnames(histologicalType), value = TRUE)
  
  print(colNaamNormal)
  print(colNaamTumor)
  
  source(paste(wholePath, "wilcoxon test over genen.R", sep = ""))
  header <- c("gene", "mean normaal", "mean tumor", "p-value wilcoxon test", "FDR", "p.adjust")
  naamWilcoxonBestand <- paste(opslaanPath, "000 wilcoxon test over genen - Normaal Tumor")
  naamMeanEenMinMeanTwee <- paste(opslaanPath, "000 mean cluster een min mean cluster twee - Normaal Tumor", sep = " ")
  #Roept de functie alleExpressieWaardenPerGroep aan.
  alleExpressieWaardenPerGroep(colNaamNormal, colNaamTumor, benodigdeData, opslaanPath, header, naamWilcoxonBestand, naamMeanEenMinMeanTwee)
  
  
  
  normalTumor <- list(colNaamNormal, colNaamTumor)
  return(normalTumor)
}