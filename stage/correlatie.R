########################################################################################################################## 
# Correlaties worden berekend wanneer soortBestand is Nieren.
#
##########################################################################################################################
berekenenCorraltieNieren <- function(benodigdeData, opslaanPath){
  sample <- colnames(benodigdeData)
  
  #Door de "?" pakt hij de eerste "_" die hij tegenkomt, dus niet die na Tumor.
  pattern_tag <- "_.*?_"
  
  alleTags  <- c()
  
  #name is elke kolomnaam
  for(name in sample){
    #temp dat zijn de A, B, etc die hij allemaal vind.
    temp <- regmatches(name, regexpr(pattern_tag, name))
    #Hierdoor hou je alleen A, B, etc over en niet de "_" meer.
    tag <- gsub("_","", temp)
    #Alle tags toevoegen aan een lijst 
    alleTags <- c(alleTags, tag)
  }
  
  #Alleen de unieke overhouden 
  uniekeTags <- unique(alleTags)
  
  #Correlatie berekenen over gehele data.
  correlatieGeheleData <- cor(benodigdeData)
  naamBestand <- paste(opslaanPath, "00 correlatie over gehele data", sep = " ")
  write.table(correlatieGeheleData, file = naamBestand, sep = ",")
  
  #Je hebt _Tumor en _normal weefsel
  soortWeefsel <- c("", "_Tumor", "_normal")
  #Er zijn 4 verschillende weefsels wanneer je alle weefsels pakt.
  #Er zijn 3 verschillende tumor weefsels.
  #Er is 1 normaal weefsel.
  lengteRijen <- c(4, 3, 1)
  normaalVsTumor <- data.frame(matrix(nrow = nrow(benodigdeData)))
  
  
  for(index in 1:length(soortWeefsel)){
    #Haalt het soort weefsel op.
    weefsel <- soortWeefsel[index]
    allePatienten <- data.frame(matrix(nrow = (nrow(benodigdeData) * lengteRijen[index])))
    namen <- c("")
    tumorNaam <- c("")
    tumoren <- data.frame(matrix(nrow = nrow(benodigdeData)))
    #Naam weefsel veranderen.
    if(weefsel == ""){
      weefsel <- "Alle weefsels"
    }
    #Loop over de unieke Tags.
    for(tags in uniekeTags){
      #Patroon maken.
      patroon <- paste0("ccRCC_", tags, soortWeefsel[index])
      #Pak de index wanneer de kolomnaam hetzelfde is als het patroon.
      colIndexTag <- grep(patroon, colnames(benodigdeData))
      #Pak de naam van de kolom wanneer de kolomnaam hetzelfde is als het patroon.
      colNaam <- grep(patroon, colnames(benodigdeData), value = TRUE)
      
      dataFrame<-NULL
      
      vectorPatient <- c()
      #Loop over aantal kolomnamen met het patroon.
      for(ind in 1:length(colIndexTag)){
        dataFrame <- cbind(dataFrame,benodigdeData[,colIndexTag[ind]])
        vectorPatient <- c(vectorPatient, benodigdeData[,colIndexTag[ind]])
        #Wanneer het weefsel gelijk is aan "_Tumor" voegt hij de kolommen samen.
        #Ook voegt hij alle kolomnamen toe aan een vector.
        if(soortWeefsel[index] == "_Tumor"){
          tumoren <- cbind(tumoren, benodigdeData[,colIndexTag[ind]])
          tumorNaam <- c(tumorNaam, colNaam[ind])
        }
      }
      #Correlatie kan alleen berekend worden wanneer je minimaal 2 groepen met elkaar vergelijkt.
      if(length(colIndexTag) > 1){
        colnames(dataFrame) <- colNaam
        correlatie <- cor(dataFrame)
        naamBestand <- paste(opslaanPath, "00 correlatie", tags, weefsel, sep = " ")
        write.table(correlatie, file = naamBestand, sep = ",")
      }
      allePatienten <- cbind(allePatienten, vectorPatient)
      namen <- c(namen, patroon)
      
      #Wanneer het soortWeefsel gelijk is aan "_Tumor" wordt er een gemiddelde berekend van de 3 tumoren van 1 patient.
      if(soortWeefsel[index] == "_Tumor"){
        gemiddeldeTumoren <- as.data.frame(rowMeans(dataFrame))
        colnames(gemiddeldeTumoren) <- patroon
        normaalVsTumor <- cbind(normaalVsTumor, gemiddeldeTumoren)
      }
      #Wanneer het soortWeefsel gelijk is aan "_normal" wordt er een gemiddelde berekend van de 3 tumoren van 1 patient.
      if(soortWeefsel[index] == "_normal"){
        gemiddeldeTumoren <- as.data.frame(rowMeans(dataFrame))
        colnames(gemiddeldeTumoren) <- patroon
        normaalVsTumor <- cbind(normaalVsTumor, gemiddeldeTumoren)
      }
    }
    #De kolomnamen worden veranderd.
    colnames(allePatienten) <- namen
    allePatienten <- allePatienten[,-1]
    #Correlatie wordt berekend over alle patienten, waarbij alle weefsels onder elkaar zijn geplakt.
    correlatieGroot <- cor(allePatienten)
    naamBestand <- paste(opslaanPath, "00 correlatie", weefsel, "Alles bij elkaar gepakt", sep = " ")
    write.table(correlatieGroot, file = naamBestand, sep = ",")
    #Wanneer soortWeefsel gelijk is aan "_Tumor" wordt de correlatie berekend over alleen de tumoren.
    if(soortWeefsel[index] == "_Tumor"){
      colnames(tumoren) <- tumorNaam
      tumoren <- tumoren[,-1]
      correlatieTumoren <- cor(tumoren)
      naamBestand <- paste(opslaanPath, "00 correlatie", weefsel, sep = "")
      write.table(correlatieTumoren, file = naamBestand, sep = ",")
    }
  }
  #Correlatie wordt berekend over het gemiddelde van de tumorweefsels tegen de normale weefsels.
  normaalVsTumor <- normaalVsTumor[,-1]
  correlatieNormaalVsTumor <- cor(normaalVsTumor)
  naamBestand <- paste(opslaanPath, "00 correlatie met gemiddelde tumoren vs normaal weefsel", sep = " ")
  write.table(correlatieNormaalVsTumor, file = naamBestand, sep = ",")
  
}

########################################################################################################################## 
# Correlaties worden berekend over de gehele data.
#
##########################################################################################################################
berekenenCorraltie <- function(benodigdeData, opslaanPath){
  #De data wordt nummeric gemaakt.
  nummericDataFrame <- apply(benodigdeData, 2, as.numeric)
  #De namen van de genen worden weer de namen van de rijen.
  rownames(nummericDataFrame) <- rownames(benodigdeData)
  correlatieGeheleData <- cor(nummericDataFrame)
  naamBestand <- paste(opslaanPath, "00 correlatie over gehele data", sep = " ")
  write.table(correlatieGeheleData, file = naamBestand, sep = ",")
}

