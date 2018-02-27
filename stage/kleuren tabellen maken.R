########################################################################################################################## 
# Groepen maken op histological type.
#
########################################################################################################################## 
histologicalType <- function(geheleData, samples){
  #De eerste rij van het bestand bestaat uit de soort histological type.
  #Deze worden in histologicalType opgeslagen.
  histologicalType <- geheleData[1, 2:ncol(geheleData)]
  
  #Vectoren per histological type.
  astrocyten <- c()
  oligodendrocyten <- c()
  GB <- c()
  ependymoma <- c()
  variantGlioma <- c()
  ganglioglioma <- c()
  metastasis <- c()
  dysembryoplasticNeuroepethelialTumor <- c()
  lymphoproliferativeDisease <- c()
  #Wanneer er een niet gedifineert histological type voorkomt wordt hij beschouwd als onbekend.
  nogNietBekend <- c()
  #Verdeelt de samples in groepen.
  for(kolom in 1:length(histologicalType)){
    if(as.character(histologicalType[,kolom]) == "A"){
      astrocyten <- c(astrocyten, colnames(histologicalType[kolom]))
    }else if (as.character(histologicalType[,kolom]) == "O"){
      oligodendrocyten <- c(oligodendrocyten, colnames(histologicalType[kolom]))
    }else if(as.character(histologicalType[,kolom]) == "G"){
      GB <- c(GB, colnames(histologicalType[kolom]))
    }else if(as.character(histologicalType[,kolom]) == "M"){
      metastasis <- c(metastasis, colnames(histologicalType[kolom]))
    }else if(as.character(histologicalType[,kolom]) == "E"){
      ependymoma <- c(ependymoma, colnames(histologicalType[kolom]))
    }else if(as.character(histologicalType[,kolom]) == "V"){
      variantGlioma <- c(variantGlioma, colnames(histologicalType[kolom]))
    }else if(as.character(histologicalType[,kolom]) == "GG"){
      ganglioglioma <- c(ganglioglioma, colnames(histologicalType[kolom]))
    }else if(as.character(histologicalType[,kolom]) == "L"){
      lymphoproliferativeDisease <- c(lymphoproliferativeDisease, colnames(histologicalType[kolom]))
    }else if(as.character(histologicalType[,kolom]) == "DNET"){
      dysembryoplasticNeuroepethelialTumor <- c(dysembryoplasticNeuroepethelialTumor, colnames(histologicalType[kolom]))
    }else{
      nogNietBekend <- c(nogNietBekend, colnames(histologicalType[kolom]))
    }
  }
  #Roept de functie kleurenTabelMaken aan.
  kleurenTabelHistologicalType = kleurenTabelMakenHisto(astrocyten, oligodendrocyten, GB, metastasis,ependymoma, 
                                   variantGlioma, ganglioglioma,  dysembryoplasticNeuroepethelialTumor,
                                   lymphoproliferativeDisease, samples, nogNietBekend)
  return(kleurenTabelHistologicalType)
}

########################################################################################################################## 
# Maakt de kleurenTabel voor histological types.
#
########################################################################################################################## 
kleurenTabelMakenHisto <- function(astro, oligo, gliobla, metastasis, ependymoma, 
                              variantGlioma, ganglioglioma, dysembryoplasticNeuroepethelialTumor,
                              lymphoproliferativeDisease, samples, nogNietBekend){
  kleurenTabelHistologicalType <- samples
  #Zet op de plek van het sample een bepaalde kleurcode.
  #Kleurcodes van R kan je vinden op: https://www.rapidtables.com/web/color/RGB_Color.html
  for(pat in 1:length(kleurenTabelHistologicalType)){
    astrocyten <- grep(kleurenTabelHistologicalType[pat], astro, value = TRUE)
    oligodendrocyten <- grep(kleurenTabelHistologicalType[pat], oligo, value = TRUE)
    glioblastoma <- grep(kleurenTabelHistologicalType[pat], gliobla, value = TRUE)
    metasta <- grep(kleurenTabelHistologicalType[pat], metastasis, value = TRUE)
    epebdymo <- grep(kleurenTabelHistologicalType[pat], ependymoma, value = TRUE)
    varGlio <- grep(kleurenTabelHistologicalType[pat], variantGlioma, value = TRUE)
    ganggliogli<- grep(kleurenTabelHistologicalType[pat], ganglioglioma, value = TRUE)
    lymphoDisease <- grep(kleurenTabelHistologicalType[pat], lymphoproliferativeDisease, value = TRUE)
    DyNET<- grep(kleurenTabelHistologicalType[pat], dysembryoplasticNeuroepethelialTumor, value = TRUE)
    onbekend <- grep(kleurenTabelHistologicalType[pat], nogNietBekend, value = TRUE)
    if(length(astrocyten) > 0){
      kleurenTabelHistologicalType = gsub(kleurenTabelHistologicalType, pattern = astrocyten, replacement = "#EFE93E")
    }else if(length(oligodendrocyten) > 0){
      kleurenTabelHistologicalType = gsub(kleurenTabelHistologicalType, pattern = oligodendrocyten, replacement = "#C0C0C0")
    }else if(length(glioblastoma) > 0){
      kleurenTabelHistologicalType = gsub(kleurenTabelHistologicalType, pattern = glioblastoma, replacement = "#666600")
    }    else if(length(metasta) > 0){
      kleurenTabelHistologicalType = gsub(kleurenTabelHistologicalType, pattern = metasta, replacement = "#CC6600")
    }else if(length(epebdymo) > 0){
      kleurenTabelHistologicalType = gsub(kleurenTabelHistologicalType, pattern = epebdymo, replacement = "#B5EF89")
    }else if(length(varGlio) > 0){
      kleurenTabelHistologicalType = gsub(kleurenTabelHistologicalType, pattern = varGlio, replacement = "#9933FF")
    }    else if(length(ganggliogli) > 0){
      kleurenTabelHistologicalType = gsub(kleurenTabelHistologicalType, pattern = ganggliogli, replacement = "#99FFFF")
    }else if(length(lymphoDisease) > 0){
      kleurenTabelHistologicalType = gsub(kleurenTabelHistologicalType, pattern = lymphoDisease, replacement = "#FF99CC")
    }else if(length(DyNET) > 0){
      kleurenTabelHistologicalType = gsub(kleurenTabelHistologicalType, pattern = DyNET, replacement = "#660000")
    }else if(length(onbekend) > 0){
      kleurenTabelHistologicalType = gsub(kleurenTabelHistologicalType, pattern = onbekend, replacement = "#FFFFFF")
    }
  }
  return(kleurenTabelHistologicalType)
}

########################################################################################################################## 
# Maakt groepen op IDHwt en IDHmt.
#
########################################################################################################################## 
bestandenOpenenIDH <- function(samples, opslaanPath){
  #Haalt de volgende bestanden op.
  IDH1.395GA <- as.character(read.csv(file = paste(opslaanPath, "IDH1 c.395G)A 10", sep = " "), check.names = FALSE)[,2])
  IDH1.532GA <- as.character(read.csv(file = paste(opslaanPath, "IDH1 c.532G)A 10", sep = " "), check.names = FALSE)[,2])
  IDH1.548AG <- as.character(read.csv(file = paste(opslaanPath, "IDH1 c.548A)G 10", sep = " "), check.names = FALSE)[,2])
  IDH2.515GA <- as.character(read.csv(file = paste(opslaanPath, "IDH2 c.515G)A 10", sep = " "), check.names = FALSE)[,2])
  IDH2.514AT <- as.character(read.csv(file = paste(opslaanPath, "IDH2 c.514A)T 10", sep = " "), check.names = FALSE)[,2])
  IDH2.515GT <- as.character(read.csv(file = paste(opslaanPath, "IDH2 c.515G)T 10", sep = " "), check.names = FALSE)[,2])
  #Beschouwd al de bovenstaande bestanden als mutanten.
  MT <- unique(as.character(unlist(list(IDH1.395GA, IDH1.532GA, IDH1.548AG, IDH2.515GA, IDH2.514AT, IDH2.515GT))))
  #Het verschil tussen de mutanten en alle sampels zijn de wild types.
  WT <- setdiff(samples, MT)
  
  #HARDCODE
  IDH1.395GA.IDH1.532GA <- intersect(IDH1.395GA, IDH1.532GA)
  IDH1.395GA <- setdiff(IDH1.395GA, IDH1.395GA.IDH1.532GA)
  IDH1.532GA <- setdiff(IDH1.532GA, IDH1.395GA.IDH1.532GA)
  
  #Roept de functie kleurenLijstWtMt aan.
  kleurenTabelStatusIDH = kleurenTabelMakenIDH(IDH1.395GA, IDH1.532GA, IDH1.548AG, IDH1.395GA.IDH1.532GA, IDH2.515GA,
                                 IDH2.514AT, IDH2.515GT, WT, samples)
  return(kleurenTabelStatusIDH)
}

########################################################################################################################## 
# Maakt de kleurenTabel voor IDHwt en IDHmt.
#
########################################################################################################################## 
kleurenTabelMakenIDH <- function(IDH1.395GA, IDH1.532GA, IDH1.548AG, IDH1.395GA.IDH1.532GA, IDH2.515GA,
                             IDH2.514AT, IDH2.515GT, WT, samples){
  kleurenTabelStatusIDH <- samples
  #Zet op de plek van het sample een bepaalde kleurcode.
  for(pat in 1:length(kleurenTabelStatusIDH)){
    IDH1.R132H <- grep(kleurenTabelStatusIDH[pat], IDH1.395GA, value = TRUE)
    IDH1.V178I <- grep(kleurenTabelStatusIDH[pat], IDH1.532GA, value = TRUE)
    IDH1.Y183C <- grep(kleurenTabelStatusIDH[pat], IDH1.548AG, value = TRUE)
    IDH1.R132H.IDH1.V178I <- grep(kleurenTabelStatusIDH[pat], IDH1.395GA.IDH1.532GA, value = TRUE)
    IDH2.R172K <- grep(kleurenTabelStatusIDH[pat], IDH2.515GA, value = TRUE)
    IDH2.R172W <- grep(kleurenTabelStatusIDH[pat], IDH2.514AT, value = TRUE)
    IDH2.R172M <- grep(kleurenTabelStatusIDH[pat], IDH2.515GT, value = TRUE)
    wildType <- grep(kleurenTabelStatusIDH[pat], WT, value = TRUE)
    if(length(IDH1.R132H) > 0){
      kleurenTabelStatusIDH = gsub(kleurenTabelStatusIDH, pattern = IDH1.R132H, replacement = "#FF0000")
    }else if(length(IDH1.V178I) > 0){
      kleurenTabelStatusIDH = gsub(kleurenTabelStatusIDH, pattern = IDH1.V178I, replacement = "#FFFFFF") 
    }else if(length(IDH1.Y183C) > 0){
      kleurenTabelStatusIDH = gsub(kleurenTabelStatusIDH, pattern = IDH1.Y183C, replacement = "#FFFFFF") 
    }else if(length(IDH1.R132H.IDH1.V178I) > 0){
      kleurenTabelStatusIDH = gsub(kleurenTabelStatusIDH, pattern = IDH1.R132H.IDH1.V178I, replacement = "#FF0000")
    }else if(length(IDH2.R172K) > 0){
      kleurenTabelStatusIDH = gsub(kleurenTabelStatusIDH, pattern = IDH2.R172K, replacement = "#000000") 
    }else if(length(IDH2.R172W) > 0){
      kleurenTabelStatusIDH = gsub(kleurenTabelStatusIDH, pattern = IDH2.R172W, replacement = "#000000") 
    }else if(length(IDH2.R172M) > 0){
      kleurenTabelStatusIDH = gsub(kleurenTabelStatusIDH, pattern = IDH2.R172M, replacement = "#000000") 
    }else if(length(wildType) > 0){
      kleurenTabelStatusIDH = gsub(kleurenTabelStatusIDH, pattern = wildType, replacement = "#3333FF")
    }
  }
  return(kleurenTabelStatusIDH)
}
