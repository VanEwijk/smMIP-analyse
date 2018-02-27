########################################################################################################################## 
# Indelen in groepen van nieren.
#
########################################################################################################################## 

scheidenNieren <- function(geheleData, opslaanPath, wholePath){
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
  
  
  source(paste(wholePath, "wilcoxon test over genen.R", sep = ""))
  header <- c("Gene", "mean normaal", "mean tumor", "p-value", "critical value", "p.adjust")
  naamWilcoxonBestand <- paste(opslaanPath, "000 Wilcoxon over genen - Normaalweefsel Tumorweefsel")
  naamMeanEenMinMeanTwee <- paste(opslaanPath, "000 mean cluster een min mean cluster twee - Normaalweefsel Tumorweefsel", sep = " ")
  #Roept de functie alleExpressieWaardenPerGroep aan.
  alleExpressieWaardenPerGroep(colNaamNormal, colNaamTumor, benodigdeData, opslaanPath, header, naamWilcoxonBestand, naamMeanEenMinMeanTwee)
  
  normalTumor <- list(colNaamNormal, colNaamTumor)
  return(normalTumor)
}