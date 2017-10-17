library("rEDM")
library("ggplot2")
library("hexbin")
library("RColorBrewer")
library("MASS")

# ubuntu
# setwd("/home/bo/Desktop/deeplearning/project")
# mac
setwd("/Users/alexpb/Desktop/Course/820/deeplearning/Project")


standardDFFor4VariableSystem <- function(strScenarioName, strWPathAndFile, strXPathAndFile, strYPathAndFile, strZPathAndFile)
{
  wSimpleCorrelation <- read.csv(strWPathAndFile, header=FALSE)
  xSimpleCorrelation <- read.csv(strXPathAndFile, header=FALSE)
  ySimpleCorrelation <- read.csv(strYPathAndFile, header=FALSE)
  zSimpleCorrelation <- read.csv(strZPathAndFile, header=FALSE)
  
  dfSimpleCorrelation <- data.frame(w=wSimpleCorrelation$V1,
                                    x=xSimpleCorrelation$V1,
                                    y=ySimpleCorrelation$V1, 
                                    z=zSimpleCorrelation$V1)
  
  dfSimpleCorrelation$deltaW = c(0, diff(dfSimpleCorrelation$w))
  dfSimpleCorrelation$stdDeltaW = scale(dfSimpleCorrelation$deltaW)
  
  dfSimpleCorrelation$deltaX = c(0, diff(dfSimpleCorrelation$x))
  dfSimpleCorrelation$stdDeltaX = scale(dfSimpleCorrelation$deltaX)
  
  dfSimpleCorrelation$deltaY = c(0, diff(dfSimpleCorrelation$y))
  dfSimpleCorrelation$stdDeltaY = scale(dfSimpleCorrelation$deltaY)
  
  dfSimpleCorrelation$deltaZ = c(0, diff(dfSimpleCorrelation$z))
  dfSimpleCorrelation$stdDeltaZ = scale(dfSimpleCorrelation$deltaZ)
  
  dfSimpleCorrelation$scenarioName = strScenarioName
  
  return(dfSimpleCorrelation);
}





createPerRealizationPlot <- function(ccmXYIncludingNonnegColumn, strXDescription=ccmXYIncludingNonnegColumn$lib_column[1], strYDescription=ccmXYIncludingNonnegColumn$target_column[1], isTruncatePerRealizationRho=FALSE, isDisplayMAE=FALSE)
{
  EEmbeddingDimension = ccmXYIncludingNonnegColumn$E[1] 
  countSamplesPerL = ccmXYIncludingNonnegColumn$num_pred[1]
  
  strTitle = paste(strXDescription, " xmap ", strYDescription, "(E", EEmbeddingDimension, ",", countSamplesPerL, " samples)")
  
  gPlot <- createGGPlotInternal(ccmXYIncludingNonnegColumn, isDisplayMAE, isTruncatePerRealizationRho)
  gPlot + geom_point() +  geom_smooth(se=TRUE) + stat_quantile(quantiles=c(0.05, 0.25,0.5,0.75,0.95)) + ggtitle(strTitle)
}





performStandardCCM <- function(df, strXVariableName, strYVariableName, strXDescription=strXVariableName, strYDescription=strYVariableName, EEmbeddingDimension=6, tauDelay=1, countSamplesPerL=1000, minL=10, stepL = 10, stepLBelow250 = max(1.0,stepL/10.0), isShowPerRealizationPlot=TRUE, isTruncatePerRealizationRho=FALSE)
{
  maxL <- length(df[[strXVariableName]])
  if (maxL < 250)
    vecLibrarySizes <- seq(minL, maxL, stepLBelow250)
  else if (minL < 250)
    vecLibrarySizes <- c(seq(minL, 250-1, stepLBelow250), seq(250, maxL, stepL))
  else
    vecLibrarySizes <- seq(minL, maxL, stepL)
  
  ccmXY <- ccm( df, 
                E=EEmbeddingDimension, 
                tau=tauDelay,
                lib_column=strXVariableName, 
                target_column=strYVariableName, 
                lib_sizes=vecLibrarySizes, 
                num_samples=countSamplesPerL)
  
  ccmXY$nonnegRho = pmax(0, ccmXY$rho)
  
  if (isShowPerRealizationPlot)
    createPerRealizationPlot(ccmXY, strXDescription, strYDescription, isTruncatePerRealizationRho)
  
  return(ccmXY)
}




performSymmetricAllPairsAnalysisFor4VariableSystem <- function(df, 
                                                               EEmbeddingDimension, 
                                                               strWVariableName="w", 
                                                               strXVariableName="x",
                                                               strYVariableName="y",
                                                               strZVariableName="z", 
                                                               strWDescription = strWVariableName, 
                                                               strXDescription = strXVariableName, 
                                                               strYDescription = strYVariableName, 
                                                               strZDescription = strZVariableName, 
                                                               stepL=10, 
                                                               stepLBelow250 = max(1.0,stepL/10.0),
                                                               countSamplesPerL=1000, 
                                                               tauDelay = 1)
{
  # to clean up the code and allow us to focus on the essential detail, we define some closures up front
  performCCM <- function(strVariableName1, strVariableName2, strDescription1, strDescription2 ) { performStandardCCM(df, strVariableName1, strVariableName2, strDescription1, strDescription2, stepL=stepL, stepLBelow250=stepLBelow250, countSamplesPerL=countSamplesPerL, EEmbeddingDimension = EEmbeddingDimension, tauDelay = tauDelay) }
  
  #  print("#3")
  #  print(paste(strXVariableName, ",strYVariableName=", strYVariableName))
  # because dfSimpleCorrelation is a dataframe, dfSimpleCorrelation$scenarioName$scenarioName is of the same length as
  # other members
  # first, with no standardization or first-differencing
  
  # we proceed down the upper-diagonal matrix of possible matches
  # pairwise with w
  # w and x
  resultsWX = performCCM(strWVariableName, strXVariableName, strWDescription, strXDescription)
  resultsXW = performCCM(strXVariableName, strWVariableName, strXDescription, strWDescription)
  
  # w and y
  resultsWY = performCCM(strWVariableName, strYVariableName, strWDescription, strYDescription)
  resultsYW = performCCM(strYVariableName, strWVariableName, strYDescription, strWDescription)
  
  # w and z
  resultsWZ = performCCM(strWVariableName, strZVariableName, strWDescription, strZDescription)
  resultsZW = performCCM(strZVariableName, strWVariableName, strZDescription, strWDescription)
  
  
  # now on to pairwise with x  
  # x and y
  resultsXY = performCCM(strXVariableName, strYVariableName, strXDescription, strYDescription)
  resultsYX = performCCM(strYVariableName, strXVariableName, strYDescription, strXDescription)
  
  # x and z
  resultsXZ = performCCM(strXVariableName, strZVariableName, strXDescription, strZDescription)
  resultsZX = performCCM(strZVariableName, strXVariableName, strZDescription, strXDescription)
  
  # finally, now for pairwise with y
  # y and z
  resultsYZ = performCCM(strYVariableName, strZVariableName, strYDescription, strZDescription)
  resultsZY = performCCM(strZVariableName, strYVariableName, strZDescription, strYDescription)
  
  return(list(  resultsWX=resultsWX,
                resultsXW=resultsXW,
                resultsWY=resultsWY,
                resultsYW=resultsYW,
                resultsWZ=resultsWZ,
                resultsZW=resultsZW,
                resultsXY=resultsXY,
                resultsYX=resultsYX,
                resultsXZ=resultsXZ,
                resultsZX=resultsZX,
                resultsYZ=resultsYZ,
                resultsZY=resultsZY))
}


# now perform the analysis for all of
#   * the raw variables
#   * the first-differenced versions of the variables
#   * when just the flow variable is first differenced
performStandardNonSweepAnalysisFor4VariableSystem <- function(df, EEmbeddingDimension, stepL=10, stepLBelow250 = max(1.0,stepL/10.0), countSamplesPerL = 1000, tauDelay = 1)
{
  
  # ok, first conduct analysis without first-differencing
  sensitivityResultsForRawVariables <- performSymmetricAllPairsAnalysisFor4VariableSystem(df, EEmbeddingDimension, "w", "x", "y", "z", stepL=stepL, stepLBelow250=stepLBelow250, countSamplesPerL=countSamplesPerL, tauDelay = tauDelay)
  print("Finished analysis for non-differenced variables, now starting analysis for uniformly differenced variables...") 
  
  # now conduct analysis with symmetric first-differencing
  sensitivityResultsForAllVariablesFirstDifferenced <- performSymmetricAllPairsAnalysisFor4VariableSystem(df, EEmbeddingDimension, "stdDeltaW", "stdDeltaX", "stdDeltaY", "stdDeltaZ", "STD(DELTA w)",   "STD(DELTA x)",  "STD(DELTA y)",  "STD(DELTA z)", stepL=stepL, stepLBelow250=stepLBelow250, countSamplesPerL=countSamplesPerL, tauDelay = tauDelay)
  print("Finished analysis for uniformly differenced variables, starting analysis for stock-flow based differenced values....") 
  
  # now conduct analysis with first-differencing only for the flows
  sensitivityResultsForFlowVariablesFirstDifferenced <- performStockFlowAllPairsAnalysisFor4VariableSystem(df, EEmbeddingDimension, stepL=stepL, stepLBelow250=stepLBelow250, countSamplesPerL=countSamplesPerL, tauDelay = tauDelay)
  
  print("Finished non-sweep analysis....") 
  
  return(list(sensitivityResultsForRawVariables=sensitivityResultsForRawVariables,
              sensitivityResultsForAllVariablesFirstDifferenced=sensitivityResultsForAllVariablesFirstDifferenced,
              sensitivityResultsForFlowVariablesFirstDifferenced=sensitivityResultsForFlowVariablesFirstDifferenced))  
}







performStockFlowAllPairsAnalysisFor4VariableSystem <- function(df, 
                                                               EEmbeddingDimension, 
                                                               strWVariableName="w",
                                                               strDeltaWVariableName="deltaW",
                                                               strXVariableName="x",
                                                               strDeltaXVariableName="deltaX",
                                                               strYVariableName="y",
                                                               strDeltaYVariableName="deltaY",
                                                               strZVariableName="z",  
                                                               strDeltaZVariableName="deltaZ",
                                                               strWDescription = strWVariableName, 
                                                               strXDescription = strXVariableName, 
                                                               strYDescription = strYVariableName, 
                                                               strZDescription = strZVariableName, 
                                                               stepL=10, 
                                                               stepLBelow250 = max(1.0,stepL/10.0),
                                                               countSamplesPerL=1000, 
                                                               tauDelay=1,
                                                               minL=10)
{
  # to clean up the code and allow us to focus on the essential detail, we define some closures up front
  performCCM <- function(strVariableName1, strVariableName2, strDescription1, strDescription2 ) { performStandardCCM(df, strVariableName1, strVariableName2, paste("DELTA ", strDescription1), strDescription2, stepL=stepL, stepLBelow250=stepLBelow250, countSamplesPerL=countSamplesPerL, EEmbeddingDimension = EEmbeddingDimension, tauDelay = tauDelay, minL = minL) }
  
  #  print("#3")
  #  print(paste(strXVariableName, ",strYVariableName=", strYVariableName))
  # because dfSimpleCorrelation is a dataframe, dfSimpleCorrelation$scenarioName$scenarioName is of the same length as
  # other members
  # first, with no standardization or first-differencing
  
  # we proceed down the upper-diagonal matrix of possible matches
  # pairwise with w
  # w and x
  resultsWX = performCCM(strDeltaWVariableName, strXVariableName, strWDescription, strXDescription)
  resultsXW = performCCM(strDeltaXVariableName, strWVariableName, strXDescription, strWDescription)
  
  # w and y
  resultsWY = performCCM(strDeltaWVariableName, strYVariableName, strWDescription, strYDescription)
  resultsYW = performCCM(strDeltaYVariableName, strWVariableName, strYDescription, strWDescription)
  
  # w and z
  resultsWZ = performCCM(strDeltaWVariableName, strZVariableName, strWDescription, strZDescription)
  resultsZW = performCCM(strDeltaZVariableName, strWVariableName, strZDescription, strWDescription)
  
  
  # pairwise with x
  # x and y
  resultsXY = performCCM(strDeltaXVariableName, strYVariableName, strXDescription, strYDescription)
  resultsYX = performCCM(strDeltaYVariableName, strXVariableName, strYDescription, strXDescription)
  
  # x and z
  resultsXZ = performCCM(strDeltaXVariableName, strZVariableName, strXDescription, strZDescription)
  resultsZX = performCCM(strDeltaZVariableName, strXVariableName, strZDescription, strXDescription)
  
  
  # pairwise with y
  # y and z
  resultsYZ = performCCM(strDeltaYVariableName, strZVariableName, strYDescription, strZDescription)
  resultsZY = performCCM(strDeltaZVariableName, strYVariableName, strZDescription, strYDescription)
  
  return(list(resultsWX=resultsWX,
              resultsXW=resultsXW,
              resultsWY=resultsWY,
              resultsYW=resultsYW,
              resultsWZ=resultsWZ,
              resultsZW=resultsZW,
              resultsXY=resultsXY,
              resultsYX=resultsYX,
              resultsXZ=resultsXZ,
              resultsZX=resultsZX,
              resultsYZ=resultsYZ,
              resultsZY=resultsZY))
}




# The below is useful when all of the following apply
#   1) We know that we wish to perform a stock-flow first-differencing scheme (Rather than 2 different schemeas, as supported by performStandardSweepAnalysisFor4VariableSystem)
#       the standard "symmetric" mechanism prevents use of a stock-flow first-differencing scheme, as it assumes that the variable on the left side of the xmap is the same as that on the right in the other element of the pair 
#   2) We want to performa sweep through the different Es
#       Note that all 3 schemes (stock-flow diff, diff on both sides, raw variables on both sides) are supported for NON-sweep analysis by performStandardNonSweepAnalysisFor4VariableSystem
#   3) We don't want to perform an analysis by tau
performStandardStockFlowCCMByEFor4VariableSystem <-  function(df, 
                                                              vecEmbeddingDimensionE=c(2,4,6,8,10,12,14),  
                                                              tauDelay=1,
                                                              strWVariableName="w",
                                                              strDeltaWVariableName="deltaW",
                                                              strXVariableName="x",
                                                              strDeltaXVariableName="deltaX",
                                                              strYVariableName="y",
                                                              strDeltaYVariableName="deltaY",
                                                              strZVariableName="z",  
                                                              strDeltaZVariableName="deltaZ",
                                                              strWDescription = strWVariableName, 
                                                              strXDescription = strXVariableName, 
                                                              strYDescription = strYVariableName, 
                                                              strZDescription = strZVariableName, 
                                                              stepL=10, 
                                                              stepLBelow250 = max(1.0,stepL/10.0),
                                                              countSamplesPerL=1000, 
                                                              minL=10)
{
  #	vecPairs <-  sapply(vecEmbeddingDimensionE, function (EEmbeddingDimension) {  performStandardCCMXYPair(df, strXVariableName, strYVariableName, strXDescription, strYDescription, EEmbeddingDimension, tau, countSamplesPerL, minL, stepL, isShowPerRealizationPlot, isTruncatePerRealizationRho, isShowCCMMeansPlot, strTitleForCCMPairMeanPlot) })
  countDimensions <- length(vecEmbeddingDimensionE)
  
  # After a little bit of exploration was unable to find much more elegant (e.g., one-line) way to capture the below.  I hope that this
  # can be refactored into something less ugly down the road...
  vecCCM_resultsByE = list()
  
  for (i in seq(1, length(vecEmbeddingDimensionE)))
  {
    EEmbeddingDimension <- vecEmbeddingDimensionE[i]
    vecCCM_resultsByE[[i]] <- performStockFlowAllPairsAnalysisFor4VariableSystem(df, 
                                                                                 EEmbeddingDimension, 
                                                                                 strWVariableName = strWVariableName,
                                                                                 strDeltaWVariableName = strDeltaWVariableName,
                                                                                 strXVariableName = strXVariableName,
                                                                                 strDeltaXVariableName = strDeltaXVariableName,
                                                                                 strYVariableName = strYVariableName,
                                                                                 strDeltaYVariableName = strDeltaYVariableName,
                                                                                 strZVariableName = strZVariableName,  
                                                                                 strDeltaZVariableName = strDeltaZVariableName,
                                                                                 strWDescription = strWDescription, 
                                                                                 strXDescription = strXDescription, 
                                                                                 strYDescription = strYDescription, 
                                                                                 strZDescription = strZDescription, 
                                                                                 minL=minL, 
                                                                                 stepL=stepL, 
                                                                                 stepLBelow250=stepLBelow250,
                                                                                 countSamplesPerL=countSamplesPerL, 
                                                                                 tauDelay = tauDelay)
  } 
  
  return(list(vecCCM_resultsBySweptVariable=vecCCM_resultsByE, 
              strParameterBeingVaried="E", 
              valuesForVariedParameter=vecEmbeddingDimensionE,
              countSamplesPerL=countSamplesPerL, 
              stepL=stepL, 
              minL=minL))
}








createGGPlotInternal <- function(ccmXYIncludingNonnegColumn, isDisplayMAE, isTruncatePerRealizationRho)
{
  if (isDisplayMAE)
    gPlot <- ggplot(ccmXYIncludingNonnegColumn, aes(x=lib_size, y=mae))
  else
  {
    if (isTruncatePerRealizationRho)
      gPlot <- ggplot(ccmXYIncludingNonnegColumn, aes(x=lib_size, y=nonnegRho))
    else
      gPlot <- ggplot(ccmXYIncludingNonnegColumn, aes(x=lib_size, y=rho)) 
  }
  
  return(gPlot)  
}



createPerRealizationDensityPlot <- function(ccmXYIncludingNonnegColumn, strXDescription=ccmXYIncludingNonnegColumn$lib_column[1], strYDescription=ccmXYIncludingNonnegColumn$target_column[1], isTruncatePerRealizationRho=FALSE, isDisplayMAE=FALSE)
{
  EEmbeddingDimension = ccmXYIncludingNonnegColumn$E[1] 
  countSamplesPerL = ccmXYIncludingNonnegColumn$num_pred[1]
  tau = ccmXYIncludingNonnegColumn$tau[1] 
  
  vecColours <- colorRampPalette(rev(brewer.pal(10,'Spectral')))(50)
  strTitle = paste(strXDescription, " xmap ", strYDescription, "(E", EEmbeddingDimension, ",", "tau", tau, ",", countSamplesPerL, " samples)")
  
  uniqueLibrarySizes <- unique(ccmXYIncludingNonnegColumn$lib_size)
  maxLibrarySize <- max(uniqueLibrarySizes)
  # because there can be an extra small number in the transition between from the short step size to the longer one, to get the smallest
  # standard step size, get the SECOND smallest one
  smallestStandardStepL <- diff(uniqueLibrarySizes)[2]
  
  countDistinctLibrarySizes <- maxLibrarySize/smallestStandardStepL
  countRhoBins <- 100
  
  gPlot <- createGGPlotInternal(ccmXYIncludingNonnegColumn, isDisplayMAE, isTruncatePerRealizationRho)
  
  #  print(countDistinctLibrarySizes)
  result <- gPlot + stat_bin2d(bins=c(countDistinctLibrarySizes, countRhoBins)) + scale_fill_gradientn(colours=vecColours) + ggtitle(strTitle);
  print(result);
  #  countDistinctLibrarySizes <- length(unique(ccmXYIncludingNonnegColumn$lib_size))
}



createDensityPlotsFromAllPairs4DAnalysis <- function(resultsAllPairsAnalysis, isDisplayMAE=FALSE)
{
  # we proceed down the upper-diagonal matrix of possible matches
  # pairwise with w
  # w and x
  createPerRealizationDensityPlot(resultsAllPairsAnalysis$resultsWX, isDisplayMAE=isDisplayMAE);
  createPerRealizationDensityPlot(resultsAllPairsAnalysis$resultsXW, isDisplayMAE=isDisplayMAE);
  
  # w and y
  createPerRealizationDensityPlot(resultsAllPairsAnalysis$resultsWY, isDisplayMAE=isDisplayMAE);
  createPerRealizationDensityPlot(resultsAllPairsAnalysis$resultsYW, isDisplayMAE=isDisplayMAE);
  
  # w and z
  createPerRealizationDensityPlot(resultsAllPairsAnalysis$resultsWZ, isDisplayMAE=isDisplayMAE)
  createPerRealizationDensityPlot(resultsAllPairsAnalysis$resultsZW, isDisplayMAE=isDisplayMAE)
  
  
  # pairwise with x
  # x and y
  createPerRealizationDensityPlot(resultsAllPairsAnalysis$resultsXY, isDisplayMAE=isDisplayMAE);
  createPerRealizationDensityPlot(resultsAllPairsAnalysis$resultsYX, isDisplayMAE=isDisplayMAE);
  
  # x and z
  createPerRealizationDensityPlot(resultsAllPairsAnalysis$resultsXZ, isDisplayMAE=isDisplayMAE)
  createPerRealizationDensityPlot(resultsAllPairsAnalysis$resultsZX, isDisplayMAE=isDisplayMAE)
  
  # finally, pairwise with y
  # y and z
  createPerRealizationDensityPlot(resultsAllPairsAnalysis$resultsYZ, isDisplayMAE=isDisplayMAE)
  createPerRealizationDensityPlot(resultsAllPairsAnalysis$resultsZY, isDisplayMAE=isDisplayMAE)
}



# The below can be used to plot results generated by performStandardStockFlowCCMByEFor4VariableSystem (and similar functions in the future using a PARTICULAR set of variables, rather than exploring many such sets of variables)
createDensityPlotsFromStandardSweepSingleAllPairsAnalysisFor4VariableSystem <- function(resultsStandardSweepSingleAllPairsAnalysisFor4VariableSystem, isDisplayMAE=FALSE)
{
  vecCCM_resultsBySweptVariable = resultsStandardSweepSingleAllPairsAnalysisFor4VariableSystem$vecCCM_resultsBySweptVariable
  for (i in seq(1, length(vecCCM_resultsBySweptVariable)))
  {
    createDensityPlotsFromAllPairs4DAnalysis(vecCCM_resultsBySweptVariable[[i]], isDisplayMAE=isDisplayMAE)
  } 
}


createDensityPlotsFromStandardNonSweepAnalysisFor4VariableSystem <- function(resultsStandardNonSweepAnalysisFor4VariableSystem, isDisplayMAE=FALSE)
{
  createDensityPlotsFromAllPairs4DAnalysis(resultsStandardNonSweepAnalysisFor4VariableSystem$sensitivityResultsForRawVariables, isDisplayMAE=isDisplayMAE)
  createDensityPlotsFromAllPairs4DAnalysis(resultsStandardNonSweepAnalysisFor4VariableSystem$sensitivityResultsForAllVariablesFirstDifferenced, isDisplayMAE=isDisplayMAE)
  createDensityPlotsFromAllPairs4DAnalysis(resultsStandardNonSweepAnalysisFor4VariableSystem$sensitivityResultsForFlowVariablesFirstDifferenced, isDisplayMAE=isDisplayMAE)
}




#   now the code to perform the analysis
#   the below loads in the CSV files generated by the scala code that I already shared with you

preStr <- "./data/dualPredpreyVariant"
variantStrW <- "_prey1.csv"
variantStrX <- "_predator1.csv"
variantStrY <- "_prey2.csv"
variantStrZ <- "_predator2.csv"
totalNum <- 1


for(item in c(1 : totalNum)){
  # for every set -> causal category
  print(paste("This is model", toString(item)))
  start_time <- Sys.time()

  #read file
  wStr <- paste(preStr, toString(item), variantStrW, sep="")
  xStr <- paste(preStr, toString(item), variantStrX, sep="")
  yStr <- paste(preStr, toString(item), variantStrY, sep="")
  zStr <- paste(preStr, toString(item), variantStrZ, sep="")
  description <- paste("Dual Pred Prey Variant ", toString(item), sep="")



  #   below performs the actual analysis for a SINGLE value of E (E=4)
  #resultDualPredPreyVariant6E4 <- performStandardNonSweepAnalysisFor4VariableSystem(dfDualPredPreyVariant6, stepL=200, EEmbeddingDimension = 4, countSamplesPerL=1000)

  #   the below creatses a density plot for the actual analysis for the SINGLE value of E (E=4) computed above
  #createDensityPlotsFromStandardNonSweepAnalysisFor4VariableSystem(resultDualPredPreyVariant6E4)

  dfDualPredPreyVariant <- standardDFFor4VariableSystem(description, wStr, xStr, yStr, zStr)

  #   the below performs the analysis for multiple values of E (here, E={2,4,8,16}):
  resultsSymmetricSweepAnalysisPredPreyVariant <- performStandardStockFlowCCMByEFor4VariableSystem(dfDualPredPreyVariant, 
                                                                                                  vecEmbeddingDimensionE=c(2,4,8,16),  
                                                                                                  stepL = 200, 
                                                                                                  countSamplesPerL = 1000)

  # the below creates the density plots from this
  createDensityPlotsFromStandardSweepSingleAllPairsAnalysisFor4VariableSystem(resultsSymmetricSweepAnalysisPredPreyVariant)

  end_time <- Sys.time()

  print(paste("Duration: ", toString(end_time - start_time)))
}
