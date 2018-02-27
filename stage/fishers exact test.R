########################################################################################################################## 
# Fisher's Exact Test wordt uitgevoerd.
#
########################################################################################################################## 
fishersExactTest <- function(bezitMutatie, groep1Groep2){
  #Count maken die tellen hoeveel samples in cluster Een en Twee de mutaties bevatten.
  welInGroep1 = 0
  welInGroep2 = 0
  
  #Gekeken wordt hoeveel samples per cluster de mutatie bevatten.
  for(sample in 1:length(bezitMutatie)){
    #GROEP 1
    for(gr1 in 1:length(groep1Groep2[[1]])){
      if(bezitMutatie[sample] == groep1Groep2[[1]][gr1]){
        welInGroep1 = welInGroep1 + 1
      }
    }
    #GROEP 2
    for(gr2 in 1:length(groep1Groep2[[2]])){
      if(bezitMutatie[sample] == groep1Groep2[[2]][gr2]){
        welInGroep2 = welInGroep2 + 1
      }
    }
  }
  #Het verschil tussen alle samples uit cluster 1 of 2 en het aantal mutaties in dat cluster is het aantal samples dat de mutatie niet heeft.
  nietInGroep1 <- length(groep1Groep2[[1]]) - welInGroep1
  nietInGroep2 <- length(groep1Groep2[[2]]) - welInGroep2
  
  #Matrix maken om de fisher test later uit te voeren.
  matrixFisher <- matrix(c(welInGroep1, nietInGroep1, welInGroep2, nietInGroep2), nrow = 2)
  return(matrixFisher)
}
