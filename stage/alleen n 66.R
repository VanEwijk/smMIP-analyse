########################################################################################################################## 
# Dit script wordt alleen aangroepen wanneer "soortBestand" gelijk is aan N 66.
# Wanneer je het bestand van n=66 hebt ingeladen, vergelijkt hij ook de lange en korte overlevers van cluster A, B en AB.
# Ook vergelijkt hij cluster A tegen de twee IDH1-R132H uit cluster B.
# Ook vergelijkt hij CA12 var1 tegen CA12 var2 (in IDHwt).
#
########################################################################################################################## 
alleenN66 <- function(benodigdeData, opslaanPath, percentageCutOff, uniekeMutaties, wholePath, geheleData, statusHistological){
  #Aantal dagen korte en lange overlevers wordt opgehaald.
  kortClusterA <- as.numeric(dlgInput("Korte overlevers cluster A aantal dagen:", Sys.info()["kortClusterA"])$res)
  langClusterA <- as.numeric(dlgInput("Lange overlevers cluster A aantal dagen:", Sys.info()["langClusterA"])$res)
  kortClusterB <- as.numeric(dlgInput("Korte overlevers cluster B aantal dagen:", Sys.info()["kortClusterB"])$res)
  langClusterB <- as.numeric(dlgInput("Lange overlevers cluster B aantal dagen:", Sys.info()["langClusterB"])$res)
  kortClusterAB <- as.numeric(dlgInput("Korte overlevers cluster A en B aantal dagen:", Sys.info()["kortClusterAB"])$res)
  langClusterAB <- as.numeric(dlgInput("Lange overlevers cluster A en B aantal dagen:", Sys.info()["langClusterAB"])$res)

  #Bestanden inlezen. Dit zijn de bestanden over de overleving en andere informatie van de patienten (klinische informatie).
  clusterA <- read.csv(file = paste(wholePath,"Kopie van Kopie van smMIP_Glioma_excel_export_20171113124934 met idh bewerkt cluster A.csv", sep = ""),
                       header = TRUE, sep = ",", na.strings=c("","NA"))
  clusterB <- read.csv(file = paste(wholePath,"Kopie van Kopie van smMIP_Glioma_excel_export_20171113124934 met idh bewerkt cluster B.csv", sep = ""),
                       header = TRUE, sep = ",", na.strings=c("","NA"))
  clusterAB <- read.csv(file = paste(wholePath,"Kopie van Kopie van smMIP_Glioma_excel_export_20171113124934 met idh bewerkt cluster A en B.csv", sep = ""),
                        header = TRUE, sep = ",", na.strings=c("","NA"))
  CA12IDHwt <- read.csv(file = paste(wholePath,"Kopie van Kopie van smMIP_Glioma_excel_export_20171113124934 met idh bewerkt IDHWT.csv", sep = ""),
                        header = TRUE, sep = ",", na.strings=c("","NA"))
  tweeMutantTegenA <- read.csv(file = paste(wholePath,"Kopie van smMIP_Glioma_excel_export_20171113124934 met idh bewerkt 2 mutant tegen A.csv", sep = ""),
                               header = TRUE, sep = ",")
  #Headers voor de bestanden met lange en korte overlevers (A, B en AB)
  headers <- c("Gene", "c. HGVS", "p. HGVS", "MT kort overleven","WT kort overleven", "MT lang overleven", "WT lang overleven" ,
                "p-value", "critical value", "p.adjust")

  source(paste(wholePath,"indeling in groepen.R", sep = ""))
  source(paste(wholePath,"statistiek op groepen (mutaties).R", sep = ""))
  #Lange vs korte overleving cluster A
  kortLangOverlevingA = maakGroepen(benodigdeData, clusterA, kortClusterA,  langClusterA, "clusterA", "overleven", opslaanPath)
  samplesA = lijstPatienten(kortLangOverlevingA)
  if(length(kortLangOverlevingA[[1]]) > 1 && length(kortLangOverlevingA[[2]]) > 1){
    KijkenNaarMutaties(samplesA, kortLangOverlevingA, percentageCutOff, uniekeMutaties, kortClusterA, langClusterA ,"clusterA", opslaanPath, headers, wholePath)
  }
  na.exclude(clusterB)
  #Lange vs korte overleving cluster B
  kortLangOverlevingB = maakGroepen(benodigdeData, clusterB, kortClusterB, langClusterB,"clusterB", "overleven", opslaanPath)
  samplesB = lijstPatienten(kortLangOverlevingB)
  if(length(kortLangOverlevingB[[1]]) > 1 && length(kortLangOverlevingB[[2]]) > 1){
    KijkenNaarMutaties(samplesB, kortLangOverlevingB, percentageCutOff, uniekeMutaties, kortClusterB, langClusterB ,"clusterB", opslaanPath, headers, wholePath)
  }
  na.exclude(clusterAB)
  #Lange vs korte overleving cluster AB
  kortLangOverlevingAB = maakGroepen(benodigdeData, clusterAB, kortClusterAB, langClusterAB,"clusterAB", "overleven", opslaanPath)
  samplesAB = lijstPatienten(kortLangOverlevingAB)
  if(length(kortLangOverlevingAB[[1]]) > 1 && length(kortLangOverlevingAB[[2]]) > 1){
    KijkenNaarMutaties(samplesAB, kortLangOverlevingAB, percentageCutOff, uniekeMutaties, kortClusterAB, langClusterAB ,"clusterAB", opslaanPath, headers, wholePath)
  }

  dagenGroep <- " "
  na.exclude(CA12IDHwt)
  headersCA12 <- c("Gene (in IDH-WT groep)", "c. HGVS", "p. HGVS", "MT CA12-insert","WT CA12-insert", "MT CA12-del", "WT CA12-del" ,
               "p-value", "critical value", "p.adjust")
  #CA12 var1 tegen CA12 var2 (in IDHwt)
  groep1Groep2CA12 = maakGroepen(benodigdeData, CA12IDHwt, dagenGroep, dagenGroep,"CA12 IDHwt", "CA12", opslaanPath)
  samplesCA12 = lijstPatienten(groep1Groep2CA12)
  if(length(groep1Groep2CA12[[1]]) > 1 && length(groep1Groep2CA12[[2]]) > 1){
    KijkenNaarMutaties(samplesCA12, groep1Groep2CA12, percentageCutOff, uniekeMutaties, dagenGroep, dagenGroep ,"CA12", opslaanPath, headersCA12, wholePath)
  }


  #cluster A tegen twee mutanten cluster B
  groep1Groep2 = maakGroepen(benodigdeData, tweeMutantTegenA, dagenGroep, dagenGroep,"tweeVSa", "groep", opslaanPath)
  samplesTweeMutantA = lijstPatienten(groep1Groep2)
  headersTweeMutantA <- c("Gene", "c. HGVS", "p. HGVS", "MT cluster A","WT cluster A", "MT cluster B", "WT cluster B" ,
                          "p-value", "critical value", "p.adjust")
  if(length(groep1Groep2[[1]]) > 1 && length(groep1Groep2[[2]]) > 1){
    KijkenNaarMutaties(samplesTweeMutantA, groep1Groep2, percentageCutOff, uniekeMutaties, dagenGroep, dagenGroep ,"tweeVSa", opslaanPath, headersTweeMutantA, wholePath)
  }

  if(statusHistological == "histological"){
    #Indelen in groepen op basis van histological type (A, O, G)
    #En op deze groepen ook een wilcoxon test doen.
    source(paste(wholePath,"indeling in groepen hist_Type.R", sep = ""))
    source(paste(wholePath,"statistiek op groepen (mutaties).R", sep = ""))
    valueAvalueO = scheidenHistType(geheleData, opslaanPath, "A", "O", wholePath)
    samplesAO = lijstPatienten(valueAvalueO)
    headerAO <- c("Gene", "c. HGVS", "p. HGVS", "MT A","WT A", "MT O", "WT O" ,
                 "p-value", "critical value", "p.adjust")
    if(length(valueAvalueO[[1]]) > 1 && length(valueAvalueO[[2]]) > 1){
      KijkenNaarMutaties(samplesAO, valueAvalueO, percentageCutOff, uniekeMutaties, "", "" ,"AO", opslaanPath, headerAO, wholePath)
    }

    valueAvalueG = scheidenHistType(geheleData, opslaanPath, "A", "G", wholePath)
    samplesAG = lijstPatienten(valueAvalueG)
    headerAG <- c("Gene", "c. HGVS", "p. HGVS", "MT A","WT A", "MT G", "WT G" ,
                 "p-value", "critical value", "p.adjust")
    if(length(valueAvalueG[[1]]) > 1 && length(valueAvalueG[[2]]) > 1){
      KijkenNaarMutaties(samplesAG, valueAvalueG, percentageCutOff, uniekeMutaties, "", "" ,"AG", opslaanPath, headerAG, wholePath)
    }

    valueGvalueO = scheidenHistType(geheleData, opslaanPath, "G", "O", wholePath)
    samplesGO = lijstPatienten(valueGvalueO)
    headerGO <- c("Gene", "c. HGVS", "p. HGVS", "MT G","WT G", "MT O", "WT O" ,
                 "p-value", "critical value", "p.adjust")
    if(length(valueGvalueO[[1]]) > 1 && length(valueGvalueO[[2]]) > 1){
      KijkenNaarMutaties(samplesGO, valueGvalueO, percentageCutOff, uniekeMutaties, "", "" ,"GO", opslaanPath, headerGO, wholePath)
    }
  }
  # Indelen in groepen IDH-R132H en IDHwt (met IDH2 en IDH1-other als WT)
  # En op deze groepen ook een wilcoxon test doen
  source(paste(wholePath,"indeling in groepen IDH1 R132H.R", sep = ""))
  source(paste(wholePath,"statistiek op groepen (mutaties).R", sep = ""))
  wildtypeMutatieIDH = scheidenIDH_Met(benodigdeData, opslaanPath, wholePath)
  samplesIDH = lijstPatienten(wildtypeMutatieIDH)
  headerIDH <- c("Gene", "c. HGVS", "p. HGVS", "MT IDH-WT (met IDH2 en IDH1-other)","WT IDH-WT (met IDH2 en IDH1-other)", "MT IDH-R132H", "WT IDH-R132H" ,
                "p-value", "critical value", "p.adjust")
  if(length(wildtypeMutatieIDH[[1]]) > 1 && length(wildtypeMutatieIDH[[2]]) > 1){
    KijkenNaarMutaties(samplesIDH, wildtypeMutatieIDH, percentageCutOff, uniekeMutaties, "", "" ,"IDH", opslaanPath, headerIDH, wholePath)
  }
  # Indelen in groepen IDH-R132H en IDHwt (zonder IDH2 en IDH1-other als WT)
  # En op deze groepen ook een wilcoxon test doen
  # source(paste(wholePath,"indeling in groepen IDH1 R132H.R", sep = ""))
  # source(paste(wholePath,"statistiek op groepen (mutaties).R", sep = ""))
  wildtypeZonderIDHotherMutatieIDH = scheidenIDH_Zonder(benodigdeData, opslaanPath, wholePath)
  samplesIDHzonder = lijstPatienten(wildtypeZonderIDHotherMutatieIDH)
  headerIDHzonder <- c("Gene", "c. HGVS", "p. HGVS", "MT IDH-WT (zonder IDH2 en IDH1-other)","WT IDH-WT (zonder IDH2 en IDH1-other)", "MT IDH-R132H", "WT IDH-R132H" ,
                 "p-value", "critical value", "p.adjust")
  if(length(wildtypeZonderIDHotherMutatieIDH[[1]]) > 1 && length(wildtypeZonderIDHotherMutatieIDH[[2]]) > 1){
    KijkenNaarMutaties(samplesIDHzonder, wildtypeZonderIDHotherMutatieIDH, percentageCutOff, uniekeMutaties, "", "" ,"IDH (zonder IDH1-other en IDH2)", opslaanPath, headerIDHzonder, wholePath)
  }
}

