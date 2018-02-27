########################################################################################################################## 
# Indelen in groepen van IDH. Met de rest (ook IDH1-other en IDH2) als wildtype.
#
########################################################################################################################## 

scheidenIDH <- function(benodigdeData, opslaanPath, wholePath){
    #Haalt de volgende bestanden op.
    IDH1.395GA <- as.character(read.csv(file = paste(opslaanPath, "IDH1 c.395G)A 10", sep = " "), check.names = FALSE)[,2])
    #Beschouwd al de bovenstaande bestanden als mutanten.
    MT <- unique(as.character(unlist(list(IDH1.395GA))))
    #Het verschil tussen de mutanten en alle sampels zijn de wild types.
    WT <- setdiff(samples, MT)
    
    print(WT)
    print(length(WT))
    print(MT)
    print(length(MT))
    
    source(paste(wholePath, "wilcoxon test over genen.R", sep = ""))
    header <- c("gene", "mean IDH-WT", "mean IDH1-R132H", "p-value wilcoxon test", "FDR", "p.adjust")
    naamWilcoxonBestand <- paste(opslaanPath, "000 wilcoxon test over genen - IDH1", sep = " ")
    naamMeanEenMinMeanTwee <- paste(opslaanPath, "000 mean cluster een min mean cluster twee - IDH1", sep = " ")
    #Roept de functie alleExpressieWaardenPerGroep aan.
    alleExpressieWaardenPerGroep(WT, MT, benodigdeData, opslaanPath, header, naamWilcoxonBestand, naamMeanEenMinMeanTwee)
    
    wildtypeMutatieIDH <- list(WT, MT)
    return(wildtypeMutatieIDH)
  }
  
########################################################################################################################## 
# Indelen in groepen van IDH. Met de rest (ook IDH1-other en IDH2) niet als wildtype.
#
##########################################################################################################################
scheidenIDHwt <- function(benodigdeData, opslaanPath, wholePath){
  #Haalt de volgende bestanden op.
  IDH1.395GA <- as.character(read.csv(file = paste(opslaanPath, "IDH1 c.395G)A 10", sep = " "), check.names = FALSE)[,2])
  IDH1.532GA <- as.character(read.csv(file = paste(opslaanPath, "IDH1 c.532G)A 10", sep = " "), check.names = FALSE)[,2])
  IDH1.548AG <- as.character(read.csv(file = paste(opslaanPath, "IDH1 c.548A)G 10", sep = " "), check.names = FALSE)[,2])
  IDH2.515GA <- as.character(read.csv(file = paste(opslaanPath, "IDH2 c.515G)A 10", sep = " "), check.names = FALSE)[,2])
  IDH2.514AT <- as.character(read.csv(file = paste(opslaanPath, "IDH2 c.514A)T 10", sep = " "), check.names = FALSE)[,2])
  IDH2.515GT <- as.character(read.csv(file = paste(opslaanPath, "IDH2 c.515G)T 10", sep = " "), check.names = FALSE)[,2])
  #Alleen IDH1-R132H
  onlyIDH1R132H <- setdiff(IDH1.395GA, IDH1.532GA)
  
  
  #Beschouwd al de bovenstaande bestanden als mutanten.
  MT <- unique(as.character(unlist(list(IDH1.395GA, IDH1.532GA, IDH1.548AG, IDH2.515GA, IDH2.514AT, IDH2.515GT))))
  #Het verschil tussen de mutanten en alle sampels zijn de wild types.
  WT <- setdiff(samples, MT)
  
  source(paste(wholePath, "wilcoxon test over genen.R", sep = ""))
  header <- c("gene", "mean IDH-WT", "mean IDH1-R132H", "p-value wilcoxon test", "FDR", "p.adjust")
  naamWilcoxonBestand <- paste(opslaanPath, "000 wilcoxon test over genen - IDH1 wt (zonder IDH1-other en IDH2)", sep = " ")
  naamMeanEenMinMeanTwee <- paste(opslaanPath, "000 mean cluster een min mean cluster twee - IDH1 wt (zonder IDH1-other en IDH2)", sep = " ")
  #Roept de functie alleExpressieWaardenPerGroep aan.
  alleExpressieWaardenPerGroep(WT, onlyIDH1R132H, benodigdeData, opslaanPath, header, naamWilcoxonBestand, naamMeanEenMinMeanTwee)
  
  wildtypeZonderIDHotherMutatieIDH <- list(WT, onlyIDH1R132H)
  print(length(WT))
  print(length(onlyIDH1R132H))
  return(wildtypeZonderIDHotherMutatieIDH)
}
  
  
 