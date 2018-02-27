#Verschillende librarys die eerst geimporteerd moeten worden.
library(dendextend)
library("gplots")
library(heatmap.plus)
library(svDialogs) #library voor het ophalen van inputGebruiker
library(plyr) #Bestand lezen ;, , ,\t


#Bestand ophalen FPM bestand (moeten .csv bestanden zijn).
path <- file.choose()
#Soort slash wordt opgehaald, dit kan namelijk per computer verschillen.
#Het ingeladen bestand moet in een mapje staan met de naam: stage.
#In de map stage moet een map zitten met de mutatie bestanden met de naam: SeqNext raw data files.
#In de map stage moet ook een map zitten met de naam: Resultaten.
soortSlash <- substr(path, regexec(".stage", path)[[1]][1], regexec(".stage", path)[[1]][1])
#Maakt een pad, om verschllende bestanden op te halen.
wholePath <- paste(strsplit(path, "stage")[[1]][1], "stage", soortSlash , sep = "") #Om de andere bestanden (met functies) op te halen.
#Haalt het bestand dat de gebruiker heeft gekozen op.
geheleData <- ldply(path,function(x) 
  if(grepl(",",readLines(x,n=1))){read.csv(x,sep=",", header = TRUE, na.strings = c("", "NA"))} 
  else if(grepl(";",readLines(x,n=1))){read.csv(x,sep=";", header = TRUE, na.strings = c("", "NA"))} 
  else{read.csv(x,sep="\t", header = TRUE, na.strings = c("", "NA"))})
#Aangeven wat voor een soort bestand het is.
soortBestand <- switch(menu(c("n = 66","n = 75", "Nieren", "ADV", "anders"), graphics = TRUE, title = "soort bestand") + 1,
       cat("Nothing done\n"), "N 66", "N 75", "Nieren", "ADV", "anders")

if(soortBestand == "anders"){
  soortBestand <- as.character(dlgInput("Wat voor een data staat er in het bestand?:", Sys.info()["soortBestand"])$res)
}


#Het path voor het opslaan van de resultaten.
opslaanPath <- paste(paste(strsplit(path, "stage")[[1]][1], "stage" , soortSlash, "Resultaten", soortSlash, sep = ""), soortBestand, sep = "") #Pad om op te slaan, met ervoor wat voor een bestand het is.
#Wanneer je in de heatmap histological types hebt staan moet soort 'histological' zijn, zo niet dan moet soort 'niet histological' zijn.
statusHistological <- switch(menu(c("histological", "niet histological"), graphics = TRUE, title = "histological type") + 1,
                       cat("Nothing done\n"), "histological", "")
#Plotjes maken per mutatie aan of uit.
plotAanUit <- switch(menu(c("plot AAN", "plot UIT"), graphics = TRUE, title = "plotten?") + 1,
                cat("Nothing done\n"), "AAN", "UIT")
#Hier moet het path staan die naar alle mutatiefiles leidt.
pathMutatieBestanden <- paste(wholePath, "SeqNext raw data files", soortSlash, sep = "")

# #Parameters die de gebruiker kan aanpassen
# methodeDistanceMatrix = "manhattan" #"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
# methodeClusteren = "ward.D2" #"ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid".
# #Cutoff voor percentage in mutatiefile.
# percentageCutOff = 10

methodeDistanceMatrix <- switch(menu(c("manhattan","euclidean", "maximum", "canberra", "binary", "minkowski"), graphics = TRUE, title = "Distance matrix") + 1,
                       cat("Nothing done\n"), "manhattan","euclidean", "maximum", "canberra", "binary", "minkowski")

methodeClusteren <- switch(menu(c("ward.D2" ,"ward.D", "single", "complete", "average", "mcquitty", "median", "centroid"), graphics = TRUE, title = "Cluster Methode") + 1,
                                cat("Nothing done\n"), "ward.D2" ,"ward.D", "single", "complete", "average", "mcquitty", "median", "centroid")
#Cutoff voor percentage in mutatiefile.
percentageCutOff <- as.numeric(dlgInput("Percentage cutoff voor mutatiefile:", Sys.info()["percentageCutOff"])$res)



########################################################################################################################## 
# Begin van aanroepen andere files.
# Door source worden andere bestanden aangeroepen.
##########################################################################################################################
#Het ligt aan het soort bestand (histological of niet) welke functie er wordt gebruikt.
source(paste(wholePath, "benodigde data maken.R", sep = ""))
if(statusHistological == "histological"){
  benodigdeData = benodigdeDataMakenHist(geheleData)
}else{
  benodigdeData = benodigdeDataMakenZonderHistType(geheleData)
}

#Het histogram, dendrogram en clusters worden hier gemaakt.
source(paste(wholePath, "clusteren van data.R", sep = ""))
histogram = histogramMaken(benodigdeData, methodeDistanceMatrix, methodeClusteren, soortBestand)
dendrogram = dendrogramMaken(histogram, soortBestand)
clusterEenTwee = getClusters(histogram, benodigdeData, opslaanPath, wholePath)

#Er wordt een lijst met sample nummers gemaakt -> samples.
#En er wordt een lijst gemaakt met unieke mutaties.
source(paste(wholePath,"bewerken bestand1.R", sep = ""))
samples = samplesLijstMaken(benodigdeData)
uniekeMutaties = bestandBewerken(samples, percentageCutOff, pathMutatieBestanden, opslaanPath)

#Er wordt een fisher test op de mutaties over het dendrogram uitgevoerd.
#Je kunt ook dendrogrammen laten maken met de takken van de patienten gekleurd met die mutatie.
#Deze functie kan je in dit bestand aan en uit zetten <- voor deze keuze komt bij het runnen van het programma een popup.
source(paste(wholePath,"testen over mutaties.R",sep = ""))
dfGesplitteLijst0en1 = testen(samples, uniekeMutaties, dendrogram, methodeDistanceMatrix, methodeClusteren,
       percentageCutOff, clusterEenTwee, opslaanPath, plotAanUit, length(benodigdeData)+1)
maakNulEenFile(benodigdeData, dfGesplitteLijst0en1, opslaanPath)


#Wanneer je het bestand van n=66 hebt ingeladen, vergelijkt hij ook de lange en korte overlevers van cluster A, B en AB.
#Ook vergelijkt hij cluster A tegen de twee IDH1-R132H uit cluster B.
#Ook vergelijkt hij CA12 var1 tegen CA12 var2 (in IDHwt)
#De genen worden ook vergelijken tussen A - O, A - G en G - O.
#De genen worden ook vergelijken tussen IDH1-R132H en wildtype.
if(soortBestand == "N 66" ){
  source(paste(wholePath,"alleen n 66.R",sep = ""))
  alleenN66(benodigdeData, opslaanPath, percentageCutOff, uniekeMutaties, wholePath, geheleData, statusHistological)
}

#Wanneer soort is 'histological' dan moeten er voor de IDHwt en histological types kleuren tabellen worden gemaakt voor de kleurenbalken.
if(statusHistological == "histological"){
  source(paste(wholePath,"kleuren tabellen maken.R", sep = ""))
  kleurenTabelHistologicalType = histologicalType(geheleData, samples)
  #Alleen voor n=66 en n=75 wordt er ook een tabel gemaakt voor IDHwt en IDH mutant.
  #Voor de andere is deze tabel 'kleurenTabelStatusIDH' leeg.
  if(soortBestand == "N 66" || soortBestand == "N 75"){
    kleurenTabelStatusIDH = bestandenOpenenIDH(samples, opslaanPath)
  }else{
    kleurenTabelStatusIDH = ""
  }
  # #De heatmaps worden gemaakt.
  # source(paste(wholePath,"heatmap met hist_type.R", sep = ""))
  # orderOpMeanValues = bewerken(benodigdeData, opslaanPath)
  # heatmapMaken(orderOpMeanValues, samples, dendrogram, kleurenTabelHistologicalType, kleurenTabelStatusIDH, soortBestand)
}else{
  kleurenTabelStatusIDH = ""
  kleurenTabelHistologicalType = ""
}
# else{
#   #Hierbij wordt er alleen een heatmap gemaakt en geen gekleurde balken.
#   source(paste(wholePath,"heatmap zonder hist_type.R", sep = ""))
#   orderOpMeanValues = bewerken(benodigdeData, opslaanPath)
#   heatmapMaken(orderOpMeanValues, samples, dendrogram)
# }
source(paste(wholePath,"heatmap.R", sep = ""))
orderOpMeanValues = maakFileOrder(benodigdeData, opslaanPath)
heatmapMaken(orderOpMeanValues, samples, dendrogram, kleurenTabelHistologicalType, kleurenTabelStatusIDH, soortBestand, statusHistological)


# #Supervised clusteren van de mt/var tegen wt.
# #Werkt op sommige punten nog niet
# source(paste(wholePath,"mutatie die is gekozen.R", sep = ""))
# gekozenMutatie(uniekeMutaties, samples, benodigdeData, opslaanPath, methodeDistanceMatrix, methodeClusteren, wholePath, percentageCutOff)

#Correlatie tussen verschillende nier samples wordt berekend.
#Wanneer het soortBestand niet gelijk is aan Nieren wordt de correlatie tussen alle samples berekend uit de benodigdeData.
source(paste(wholePath,"correlatie.R", sep = ""))
if(soortBestand == "Nieren"){
  berekenenCorraltieNieren(benodigdeData, opslaanPath)
  #Vergelijkt de normale weefsels van de nieren met de tumor weefsels.
  source(paste(wholePath,"indeling in groepen nieren.R", sep = ""))
  normalTumor = scheidenNieren(geheleData, opslaanPath, wholePath)
  source(paste(wholePath,"statistiek op groepen (mutaties).R", sep = ""))
  headers <- c("Gene", "c. HGVS", "p. HGVS", "mutatie normaal weefsel","wt normaal weefsel", "mutatie tumor weefsel", "wt tumor weefsel" ,
               "p-value", "FDR", "p.adjust")
  if(length(normalTumor[[1]]) > 1 && length(normalTumor[[2]]) > 1){
    KijkenNaarMutaties(colnames(benodigdeData), normalTumor, percentageCutOff, uniekeMutaties, "", "" ,"nierNormaalTumor", opslaanPath, headers, wholePath)
  }
}else{
  berekenenCorraltie(benodigdeData, opslaanPath)
}

