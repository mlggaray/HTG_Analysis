
saveFileAsExcel <- function(table, filename, nameOfSheet){
  options(java.parameters = "-Xmx128G")
 library(r2excel)
  
  wb <- createWorkbook(type="xlsx")
  sheet <- createSheet(wb, sheetName = nameOfSheet)
  xlsx.addTable(wb, sheet, table, startCol=1)
  saveWorkbook(wb, filename)
}

findMinimumSample <- function(sampleTable, types)
{
  sc1 <- length(which(sampleTable$diagnosis == types[1]))
  sc2 <- length(which(sampleTable$diagnosis == types[2]))
  minNumberOfSamples <- 0
  if(sc1 > sc2) {
    minNumberOfSamples <- sc2
  } else {
    minNumberOfSamples <- sc1
  }
  return(minNumberOfSamples)
}


generateHeatMapSimple <- function(excelFile, pngFile){
  suppressMessages(library(grid))
  suppressMessages(library(circlize))
  suppressMessages(library(RColorBrewer))
  suppressMessages(library(pheatmap))
  suppressMessages(library(readxl))
  
  transTable <- read_excel(excelFile)
  diagnosis <- as.character(transTable[1,-1])
  transTable <- transTable[-1,]
  samplesToUse <- as.character(colnames(transTable[,-1]))
  genesNamesForHM <- transTable$geneSymbols
  
  highly_variable_lcpm <- as.data.frame(transTable)
  rownames(highly_variable_lcpm) <- highly_variable_lcpm$geneSymbols
  columns <- colnames(highly_variable_lcpm)
  columnsIdToKeep <- which(columns %in% samplesToUse)
  highly_variable_lcpm <- highly_variable_lcpm[,columnsIdToKeep]
  for(i in 1:length(columnsIdToKeep)){
    highly_variable_lcpm[,i] <- as.numeric(highly_variable_lcpm[,i])
  }
  
  diagnosisToUse <- as.factor(diagnosis)
  types <- levels(diagnosisToUse)
  
  samplesPertype <- c()
  for(i in 1:length(types)){
    samplesPertype <- c(samplesPertype, length(which(diagnosisToUse %in% types[i])))
  }
  
  Diagnosis <- diagnosisToUse
  gapsToUse <- samplesPertype[length(samplesPertype) - 1]
  
  
  smallestValue  <- min(apply(highly_variable_lcpm,2,min))
  biggestValue  <- max(apply(highly_variable_lcpm,2,max))
  range <- abs(smallestValue) + abs(biggestValue)
  divi <- range / 3
  first <- smallestValue + divi
  second <- first + 0.001
  third <- second + divi
  fourth <- third + 0.001
  
  
  #  mypalette <- colorRampPalette(c("#25D811", "#FFFFFF", "#CE0C0C"))(n = 299)
  mypalette <- colorRampPalette(c("#1142d8", "#FFFFFF", "#CE0C0C"))(n = 299)
  col_breaks = c(seq(smallestValue,first,length = 100), 
                 seq(second,third,length = 100),
                 seq(fourth,biggestValue,length = 100))
  
  
  
  someMatrix <- as.matrix(highly_variable_lcpm)
  col.diagnosis <- data.frame(Diagnosis)
  
  rownames(col.diagnosis) <- colnames(someMatrix)
  
  
  pheatmap(someMatrix,
           scale="none",
           color = mypalette,
           #           breaks = col_breaks,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           clustering_distance_rows = "correlation",#Pearson's
           clustering_method = "average",
           gaps_col = gapsToUse,
           cellheight = 6,
           cellwidth = 20,
           border_color=NA,
           fontsize_row = 6,
           filename = pngFile,
           annotation_col = col.diagnosis
  )

}



generateHeatMapSimpleNoColClustering <- function(excelFile, pngFile){
  suppressMessages(library(grid))
  suppressMessages(library(circlize))
  suppressMessages(library(RColorBrewer))
  suppressMessages(library(pheatmap))
  suppressMessages(library(readxl))
  
  transTable <- read_excel(excelFile)
  diagnosis <- as.character(transTable[1,-1])
  transTable <- transTable[-1,]
  samplesToUse <- as.character(colnames(transTable[,-1]))
  genesNamesForHM <- transTable$geneSymbols
  
  highly_variable_lcpm <- as.data.frame(transTable)
  rownames(highly_variable_lcpm) <- highly_variable_lcpm$geneSymbols
  columns <- colnames(highly_variable_lcpm)
  columnsIdToKeep <- which(columns %in% samplesToUse)
  highly_variable_lcpm <- highly_variable_lcpm[,columnsIdToKeep]
  for(i in 1:length(columnsIdToKeep)){
    highly_variable_lcpm[,i] <- as.numeric(highly_variable_lcpm[,i])
  }
  
  diagnosisToUse <- as.factor(diagnosis)
  types <- levels(diagnosisToUse)
  
  samplesPertype <- c()
  for(i in 1:length(types)){
    samplesPertype <- c(samplesPertype, length(which(diagnosisToUse %in% types[i])))
  }
  
  Diagnosis <- diagnosisToUse
  gapsToUse <- samplesPertype[length(samplesPertype) - 1]
  
  
  smallestValue  <- min(apply(highly_variable_lcpm,2,min))
  biggestValue  <- max(apply(highly_variable_lcpm,2,max))
  range <- abs(smallestValue) + abs(biggestValue)
  divi <- range / 3
  first <- smallestValue + divi
  second <- first + 0.001
  third <- second + divi
  fourth <- third + 0.001
  
  
  #  mypalette <- colorRampPalette(c("#25D811", "#FFFFFF", "#CE0C0C"))(n = 299)
  mypalette <- colorRampPalette(c("#1142d8", "#FFFFFF", "#CE0C0C"))(n = 299)
  col_breaks = c(seq(smallestValue,first,length = 100), 
                 seq(second,third,length = 100),
                 seq(fourth,biggestValue,length = 100))
  
  
  
  someMatrix <- as.matrix(highly_variable_lcpm)
  col.diagnosis <- data.frame(Diagnosis)
  
  rownames(col.diagnosis) <- colnames(someMatrix)
  
  
  pheatmap(someMatrix,
           scale="none",
           color = mypalette,
           #           breaks = col_breaks,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           clustering_distance_rows = "correlation",#Pearson's
           clustering_method = "average",
           gaps_col = gapsToUse,
           cellheight = 6,
           cellwidth = 20,
           border_color=NA,
           fontsize_row = 6,
           filename = pngFile,
           annotation_col = col.diagnosis
  )
  
}





generateHeatMapSimpleAlternative <- function(excelFile, pngFileAlt, showRowNames=TRUE){
  library(markdown)
  options(markdown.HTML.options = c(options('markdown.HTML.options')[[1]], "toc"))
  
  library(knitr)
  knitr::opts_chunk$set(
    error = FALSE,
    tidy  = FALSE,
    message = FALSE,
    fig.align = "center",
    fig.width = 5,
    fig.height = 5)
  options(markdown.HTML.stylesheet = "custom.css")
  
  options(width = 100)
  
  suppressMessages(library(ComplexHeatmap))
  suppressMessages(library(circlize))
  suppressMessages(library(readxl))
  
  transTable <- read_excel(excelFile)
  diagnosis <- as.character(transTable[1,-1])
  transTable <- transTable[-1,]
  samplesToUse <- as.character(colnames(transTable[,-1]))
  genesNamesForHM <- transTable$geneSymbols
  
  
  hitsToMap <- as.data.frame(transTable)
  rownames(hitsToMap) <- hitsToMap$geneSymbols
  columns <- colnames(hitsToMap)
  columnsIdToKeep <- which(columns %in% samplesToUse)
  hitsToMap <- hitsToMap[,columnsIdToKeep]
  for(i in 1:length(columnsIdToKeep)){
    hitsToMap[,i] <- as.numeric(hitsToMap[,i])
  }
  
  
  ha = HeatmapAnnotation(df = data.frame(group = diagnosis))
  
  numberOfRows <- dim(hitsToMap)[1]

  
  plotHeight <- 1200 
  if(showRowNames){
    hm <- Heatmap(hitsToMap, top_annotation = ha, row_title_gp = gpar(fontsize = 4))
    plotHeight <- (length(genesNamesForHM) * 10) + 400 
  } else {
    hm <- Heatmap(hitsToMap, top_annotation = ha, show_row_names = FALSE, row_title_gp = gpar(fontsize = 4))
    plotHeight <- (length(genesNamesForHM) * 5) + 400 
  }
  
  
  png(file=pngFileAlt, width = 900, height = plotHeight)
  print(hm)
  dev.off()
}


transfromTSVIntoXLSX <- function(tsv, xlsx, nameOfSheet, nameOfFirstColumn)
{
  celFilesTab <- read.csv(tsv, header = FALSE, stringsAsFactors = FALSE,  sep="\t")
  numberOfColumns <- dim(celFilesTab)[2] - 1
  colnames(celFilesTab) <- c("geneSymbols",as.character(celFilesTab[1,1:numberOfColumns]))
  celFilesTab <- celFilesTab[-1,]
  saveFileAsExcel(celFilesTab,xlsx, nameOfSheet)
}


generateVolcanoPlot <- function(genes, title, pValue, lfc, significantTranscripts){
  library(ggplot2)
  library(ggrepel)
  
  
  genes$regulation <- rep("2. noSignificant", dim(genes)[1])
  transcripts <- rownames(genes)
  significantIds <- which(transcripts %in% significantTranscripts)
  temp <- genes[significantIds,]
  posIds <- which(temp$pvalue < pValue & temp$log2FoldChange >= (lfc) & temp$Gene != "")
  posTranscripts <- rownames(temp[posIds,])
  pos <- which(transcripts %in% posTranscripts)
  negIds <- which(temp$pvalue < pValue & temp$log2FoldChange <= -(lfc) & temp$Gene != "")
  negTranscripts <- rownames(temp[negIds,])
  neg <- which(transcripts %in% negTranscripts)

  padValue <- 0.3
  colNumber <- which(colnames(genes) == "regulation")
  genes[pos,colNumber] = "3. upRegulated"
  genes[neg,colNumber] = "1. downRegulated"
  minAdjPval<- -log10(min(genes$pvalue)) + 1
  lfcPlus_label = paste0("+ ", round(2^lfc, digits = 2), " fold")
  lfcPlus_df=data.frame(name=lfcPlus_label, x = lfc + padValue, y = minAdjPval)
  lfcMinus_label = paste0("- ", round(2^lfc, digits = 2), " fold")
  lfcMinus_df=data.frame(name=lfcMinus_label, x = -lfc + padValue, y = minAdjPval)
  pval_label = paste0("pvalue = ",pValue )
  yintercept <- (-log10(pValue))
  minlog2 <- min(genes$log2FoldChange)
  sizeOfX0 <- round(minlog2 - 0.5, digits = 1)
  pval_df=data.frame(name=pval_label, x = sizeOfX0 + 0.65, y = (yintercept + 0.2))
  addText <- subset(genes, pvalue < (pValue * 0.01) & abs(log2FoldChange) >= (lfc + 1.0) & Gene != "")
  thePlot <- ggplot(genes, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = regulation)) +  scale_color_manual(values = c("green", "grey", "red")) +
    labs(x="log2 Fold Change", 
         y="-log10(p-value)", 
         title=title, 
         caption = "MLGG") +
    geom_hline(aes(yintercept = yintercept), colour = "#FFA500", linetype = "dashed", show.legend = TRUE) +  
    geom_vline(aes(xintercept = lfc), colour = "blue", linetype="dashed") + 
    geom_vline(aes(xintercept = -lfc), colour = "blue", linetype="dashed") + 
    geom_label(data=pval_df, mapping=aes(x=x, y=y, label=name), colour="#FFA500", size=4) + 
    geom_label(data=lfcPlus_df, mapping=aes(x=x, y=y, label=name), colour="blue", size=4) +
    geom_label(data=lfcMinus_df, mapping=aes(x=x, y=y, label=name), colour="blue", size=4) +
    geom_text_repel(
      data = addText,
      aes(label = Gene),
      size = 3,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")
    ) + 
    theme_bw() +
    theme(legend.position="top", legend.direction="horizontal", 
          legend.justification = "right", 
          legend.background = element_rect(colour = "black", size = 0.5), 
          plot.title = element_text(size=16, colour = "BLACK", vjust= -10)) 
  return(thePlot)
  
}


getSampleInfoInCaseControlFile <- function(y, ids){
  controlCounts <- y[, ids]
  sampleInfo <- data.frame(controlCounts$samples)
  sampleInfo$sampleName <- rownames(sampleInfo)
  sampleInfo <- sampleInfo[,c(4,1)]
  sampleInfo$group <- factor(sampleInfo$group)
  typeOfSamples <- levels(sampleInfo$group)
  sampleInfo$type <- "empty"
  typeI <- which(sampleInfo$group %in% typeOfSamples[1])
  sampleInfo[typeI,3] <- "control"
  sampleInfo[-typeI,3] <- "case"
  return(sampleInfo)
}

coefOfVar <- function(values){
  stdev <- sd(values)
  meanValue <- mean(values)
  t <- (   (stdev/meanValue) * 100 )
  notAnumbers <- which(is.nan(t))
  t[notAnumbers] <- 0
  return(t)
}

transcriptsToRemoveVector <- function(efiControlForScatterPlot)
{
  
  remove <- which(round(efiControlForScatterPlot$caseAverage) < 3 & round(efiControlForScatterPlot$controlAverage) < 3 )  
  removeTable <- efiControlForScatterPlot[remove,]
  
  remove <- which(removeTable$controlCoefOfVar > 100 | removeTable$caseCoefOfVar > 100)
  removeTable <- removeTable[remove,]
  
  transcriptsToRemove <- rownames(removeTable)
  
  return(transcriptsToRemove)
}

multipleVectorByVector <- function(myVector, myFactorVector){
  newVector <- c()
  sizeV <- length(myVector)
  sizeF <- length(myFactorVector)
  if(sizeV == sizeF){
    for(i in 1:sizeV){
      v <- myVector[i]
      f <- myFactorVector[i]
      newVector <- c(newVector, (v * f))
    }
  }else{
    print("erros sizeV is different from size F")
  }
  return(newVector)
}

generateFactors <- function(myVector, maxValue){
  newVector <- c()
  for(v in myVector) newVector <- c(newVector, (v / maxValue))
  return(newVector)
}

multipleVectorByFactor <- function(myVector, myFactor){
  newVector <- c()
  for(v in myVector) newVector <- c(newVector, (v * myFactor))
  return(newVector)
}


multipleDataFrameByVectorFactors <- function(myDataFrame, myVectorFactors){
  numberOfRows = dim(myDataFrame)[1]
  numberOfColumns =  dim(dataSheet)[2]
  tempMatrix <- matrix(data = NA, nrow = numberOfRows, ncol = numberOfColumns)
  for(i in 1:numberOfColumns){
    tempMatrix[,i] <- multipleVectorByFactor(myDataFrame[,i], myVectorFactors[i])
  }
  tempDataFrame <- as.data.frame(tempMatrix)
  colnames(tempDataFrame) <- colnames(myDataFrame)
  rownames(tempDataFrame) <- rownames(myDataFrame)
  return(tempDataFrame)
}

generateVolcanoPlotTable <- function(fit){
  
  fitNames <- rownames(fit)
  lengthfit <- length(fitNames)
  fitTableForVolcanoPlot <- topTable(fit, genelist=fit$geneSymbols,number = lengthfit , adjust.method= "BH", sort.by="p")
  if(is.null(fit$geneSymbols)){
    fitTableForVolcanoPlot$ID <- rownames(fitTableForVolcanoPlot)
    fitTableForVolcanoPlot <- fitTableForVolcanoPlot[,c(7,1:6)]
  }
  colnames(fitTableForVolcanoPlot) <- c("Gene","log2FoldChange", "aveExpr", "tTest", "pvalue", "padj", "logOdds")
  return(fitTableForVolcanoPlot)
}

generateCaseControlFile <- function(yOriginalCounts, idsToUse, fitName){

  fitTableForVolcanoPlot <- generateVolcanoPlotTable(fitName)

  library(dplyr)
  controlCounts <- yOriginalCounts[, idsToUse]
  rowCounts <- data.frame(controlCounts$counts)
  sampleInfo <- data.frame(controlCounts$samples)
  allColRowCounts <- rownames(controlCounts$samples)
  numberOfSamples <- dim(sampleInfo)[1]
  sampleInfo$group <- factor(sampleInfo$group)
  typeOfSamples <- levels(sampleInfo$group)
  normalizationFactors <- as.numeric(sampleInfo$norm.factors)
  for(i in 1:numberOfSamples) rowCounts[,i] <- rowCounts[,i] *  normalizationFactors[i]

  numberOfRows <- dim(rowCounts)[1]
  type1 <- which(sampleInfo$group == typeOfSamples[1])
  type1names <- rownames(sampleInfo[type1,])
  type1Ids <- which(allColRowCounts %in% type1names)
  type1AdjustedCounts <- rowCounts[,type1Ids]
  type1Length = length(type1)
  if(type1Length > 1){
    meanValue <- apply(type1AdjustedCounts,1, mean, na.rm = TRUE)
    stdev <-  apply(type1AdjustedCounts,1, sd, na.rm = TRUE)
    coefVar <- apply(type1AdjustedCounts,1, coefOfVar)
  } else if(type1Length == 1){
    meanValue <- type1AdjustedCounts
    stdev <- rep(0,numberOfRows)
    coefVar <- rep(0,numberOfRows)
  } else{
    print("ERROR at generateCaseControlFile() the number of record in type1 is null")
  }
  type1AdjustedCounts$mean <- meanValue
  type1AdjustedCounts$coefOfVar <- coefVar
  type1AdjustedCounts$log2 <- log2(round(type1AdjustedCounts$mean))
  infValues <- which(type1AdjustedCounts$log2 == -Inf)
  log2ColumNumber <- which(colnames(type1AdjustedCounts) == "log2")
  type1AdjustedCounts[infValues,log2ColumNumber] = 0
  type1AdjustedCounts$transcripts <- rownames(type1AdjustedCounts)
  type1AdjustedCounts <- arrange(type1AdjustedCounts, transcripts)
  rownames(type1AdjustedCounts) <- type1AdjustedCounts$transcripts
  
  type2 <- which(sampleInfo$group == typeOfSamples[2])
  type2names <- rownames(sampleInfo[type2,])
  type2Ids <- which(allColRowCounts %in% type2names)
  type2AdjustedCounts <- rowCounts[,type2Ids]
  type2Length = length(type2)
  if(type2Length > 1){
    meanValue <- apply(type2AdjustedCounts,1, mean, na.rm = TRUE)
    stdev <-  apply(type2AdjustedCounts,1, sd, na.rm = TRUE)
    coefVar <- apply(type2AdjustedCounts,1, coefOfVar)
  } else if(type2Length == 1){
    meanValue <- type2AdjustedCounts
    stdev <- rep(0,numberOfRows)
    coefVar <- rep(0,numberOfRows)
  } else{
    print("ERROR at generateCaseControlFile() the number of record in type2 is null")
  }
  type2AdjustedCounts$mean <- meanValue
  type2AdjustedCounts$coefOfVar <- coefVar
  type2AdjustedCounts$log2 <- log2(round(type2AdjustedCounts$mean))
  infValues <- which(type2AdjustedCounts$log2 == -Inf)
  log2ColumNumber <- which(colnames(type2AdjustedCounts) == "log2")
  type2AdjustedCounts[infValues,log2ColumNumber] = 0
  type2AdjustedCounts$transcripts <- rownames(type2AdjustedCounts)
  type2AdjustedCounts <- arrange(type2AdjustedCounts, transcripts)
  rownames(type2AdjustedCounts) <- type2AdjustedCounts$transcripts

  
  fitTableForVolcanoPlot$transcripts <-  rownames(fitTableForVolcanoPlot)
  fitTableForVolcanoPlot <- arrange(fitTableForVolcanoPlot, transcripts)
  rownames(fitTableForVolcanoPlot) <- fitTableForVolcanoPlot$transcripts
# appears that case and control are reversed check TODO
  caseControlFile <- data.frame(Gene = fitTableForVolcanoPlot$Gene, 
                                log2FoldChange = fitTableForVolcanoPlot$log2FoldChange,
                                aveExpr = fitTableForVolcanoPlot$aveExpr, tTest = fitTableForVolcanoPlot$tTest,
                                pvalue = fitTableForVolcanoPlot$pvalue, padj = fitTableForVolcanoPlot$padj,
                                logOdds = fitTableForVolcanoPlot$logOdds, 
                                caseAverage = type2AdjustedCounts$log2, 
                                controlAverage = type1AdjustedCounts$log2, 
                                caseCoefOfVar = type2AdjustedCounts$coefOfVar, 
                                controlCoefOfVar = type1AdjustedCounts$coefOfVar)
  rownames(caseControlFile) <- fitTableForVolcanoPlot$transcripts

   return(caseControlFile)
  
}

tableWithDiffs <- function(theList, fileName, logValue=TRUE, topValue=6, lowValue=4){
  
  myDGEList <- calcNormFactors(theList)
  normalizedCounts <- as.data.frame(cpm(myDGEList, normalized.lib.sizes=TRUE, log = logValue))
  controlsrows <- c("NEG_CTRL_ANT1","NEG_CTRL_ANT2", "NEG_CTRL_ANT3","NEG_CTRL_ANT4",
                    "POS_CTRL_POS1", "POS_CTRL_POS2","POS_CTRL_POS3", "POS_CTRL_POS4")
  
  nameOfRows <- rownames(normalizedCounts)
  
  
  #sapply(normalizedCounts, max, na.rm = T)
  
  samplesInColumn <- colnames(normalizedCounts)
  
  controlOnlyIds <- grep("Control", samplesInColumn, ignore.case = TRUE)
  controlsOnlyNames <- samplesInColumn[controlOnlyIds]
  noControlNames <- samplesInColumn[-controlOnlyIds]
  noControlIds <- which(samplesInColumn %in% noControlNames)
  
  table <- normalizedCounts
  
  for(i in noControlIds){
    ids <- which(table[,i] > topValue)
    table <- table[ids,]
  }
  
  for(i in controlOnlyIds){
    ids <- which(table[,i] < lowValue)
    table <- table[ids,]
  }
  
  passFile <- paste0(fileName,"ControlPassed.tsv")
  write.table(table, passFile, sep="\t")
  
  return(table)
}


tableWithCounts <- function(group, diagnosis){
  theList <- group
  x <-  as.matrix(diagnosis)
  colnames(x) <- "diagnosis"
  x <- as.data.frame(t(x))
  colnames(x) <- colnames(theList)
  newdf <- rbind(theList, x)
  return(newdf)
}


writeCountTables <- function(theList, fileName){
  
  noLogFile <- paste0(fileName,"Counts.tsv")
  logFile <- paste0(fileName,"LogCounts.tsv")
  normalizedCounts <- as.data.frame(cpm(theList, normalized.lib.sizes=TRUE, log = TRUE))
  transcripts <- rownames(normalizedCounts)
  geneSymbols <- returnListOfGenes(transcripts)
  normalizedCounts$geneSymbols <- geneSymbols
  write.table(normalizedCounts, logFile, sep="\t")
  normalizedCounts <- as.data.frame(cpm(theList, normalized.lib.sizes=TRUE, log = FALSE))
  transcripts <- rownames(normalizedCounts)
  geneSymbols <- returnListOfGenes(transcripts)
  normalizedCounts$geneSymbols <- geneSymbols
  write.table(normalizedCounts, noLogFile, sep="\t")
}

# writeDiagnosis <- function(diag, colnameList, fileName){
#   diagFileName <- paste0(fileName,"_diagnosis.tsv")
#   diagnosis <- as.character(diag)
#   x <-  as.matrix(diagnosis)
#   colnames(x) <- "diagnosis"
#   x <- as.data.frame(t(x))
#   colnames(x) <- colnameList
#   write.table(x, diagFileName, sep="\t")
# }



returnListOfGenes <- function(listOfTranscripts){
  load("/home/mlggaray/fitTranscripts2Gene.Rdata")
  genesSymbols <- c()
  identicalFits <- FALSE
  
  sameLength <- (length(listOfTranscripts) == length(fitTranscripts2Gene$transcripts))
  if(sameLength){
    identicalFits <- all(listOfTranscripts == fitTranscripts2Gene$transcripts)
  } 
  
  if(identicalFits){
    genesSymbols <- fitTranscripts2Gene$geneSymbols
  }else{
    load("/home/mlggaray/tx2geneHash.Rdata")
    for(v in listOfTranscripts) genesSymbols <- c(genesSymbols, tx2geneHash[[v]])
  }
  return(genesSymbols)
}

reduceFile <- function(controlCaseFile, transcriptsToRemove){
  
  splotTrans <- rownames(controlCaseFile)
  removeScatterPlotIds <- which(splotTrans %in% transcriptsToRemove)
  if(length(removeScatterPlotIds) > 0){
    smallerCaseControlFile <- controlCaseFile[-removeScatterPlotIds,]
  }else{
    smallerCaseControlFile <- controlCaseFile
  }
  return(smallerCaseControlFile)
}

generateFilesForProject <- function(yOriginalCounts, idsToUse, projectDir, projectFiles,nameOfProject, fitName, lfc, pval, significantTranscripts){

  caseControlFileName <- paste0(projectDir, projectFiles, ".tsv")
  removedTranscriptFileName <- paste0(projectDir, projectFiles, "_removed.tsv")
  volcanoFile <-  paste0(projectDir, "volcanoPlot-", projectFiles, ".png")
  volcanoPlotTitle <- paste0("Volcano Plot of ", nameOfProject)


  typeTable <- getSampleInfoInCaseControlFile(yOriginalCounts, idsToUse)
  if(length(levels(typeTable$group)) != 2) print("ERROR")
  caseControlFile <- generateCaseControlFile(yOriginalCounts, idsToUse,fitName)
#  write.table(caseControlFile, caseControlFileName, sep="\t")

  transcriptsToRemove <- transcriptsToRemoveVector(caseControlFile)
  smallerCaseControlFile <- reduceFile(caseControlFile, transcriptsToRemove)
  removedTable <- caseControlFile[transcriptsToRemove,]

#  write.table(removedTable, removedTranscriptFileName, sep="\t")

  myVolcanoPlot <- generateVolcanoPlot(smallerCaseControlFile, volcanoPlotTitle, pval, lfc, significantTranscripts)
  ggsave(volcanoFile, plot = myVolcanoPlot, width = 300, height = 250, units = "mm")

  
  return(transcriptsToRemove)
}


generateHeatMap <- function(transcripts, sampleIdToUse, samplesToUse, lcpm, types, 
                            fileName, transcriptsToRemove, fileNameWithData, theHashFile, useGenes = TRUE){
  suppressMessages(library(grid))
  suppressMessages(library(circlize))
  suppressMessages(library(RColorBrewer))
  suppressMessages(library(pheatmap))
  suppressMessages(library(hash))
  
  load(theHashFile)
  
  lastGroup <- unique(transcripts)

  samplesToUse$sampleNames <- rownames(samplesToUse)
  samplesToUse <- samplesToUse[sampleIdToUse,]
  samplesToUse$diagnosis <- factor(samplesToUse$diagnosis)

  areTypesCorrectIds <- which(levels(samplesToUse$diagnosis) %in% types)
  areCorrect <- length(areTypesCorrectIds) == length(types)
  if(!areCorrect){
    print(paste0("the types are not correct! ", types, levels(samplesToUse$diagnosis)))
    return(-9)
  }
  samplesToUse$diagnosis <- factor(samplesToUse$diagnosis, levels = types)
  samplesPertype <- c()
  for(i in 1:length(types)){
    samplesPertype <- c(samplesPertype, length(which(samplesToUse$diagnosis %in% types[i])))
  }
  diagnosisToUse <- samplesToUse$diagnosis
  samplesToUseNames <- rownames(samplesToUse)[sampleIdToUse]
  logcounts <- lcpm[,sampleIdToUse]
  bignames <- rownames(logcounts)
  found  <- bignames[which(bignames %in% lastGroup)]
  lastGroup <- found
  
  #useGenes = TRUE
  remTran = TRUE
  samplesToUse <- samplesToUse
  Diagnosis <- diagnosisToUse
  gapsToUse <- samplesPertype[length(samplesPertype) - 1]
  
  fileToUse <- fileName
  
  transcriptsToUse <- lastGroup
  if(length(transcriptsToRemove) > 0){
   removeTranscripts <- transcriptsToRemove
    
    removeTranscriptsIndex <- which(transcriptsToUse %in% removeTranscripts)
    if(length(removeTranscriptsIndex) > 0){
      transcriptsToUse <- transcriptsToUse[-removeTranscriptsIndex]
    }
      }
  
  

  highly_variable_lcpm <- logcounts[transcriptsToUse,]
  
  list <- rownames(highly_variable_lcpm)
  if(is.null(list)) return(0)
  genesNamesForHM <- c()
  for(v in list){
#    transf <- as.character(ids2geneHash[[v]])
    transf <- as.character(tx2geneHash[[v]])    #error here is hardcoded
    if(is.null(transf) || identical(transf,character(0))){ 
      transf <- v
    }else if( length(which(genesNamesForHM %in% transf) ) > 0) {
      transf <- paste0(transf, "_", v)
    }
    genesNamesForHM <- c(genesNamesForHM, transf)
  }
  
  if(is.null(genesNamesForHM)){
    genesNamesForHM <- list
  }
  
  if(useGenes){
    rownames(highly_variable_lcpm) <- genesNamesForHM
  }

  if(!is.null(fileNameWithData)){
    write.table(highly_variable_lcpm, fileNameWithData, sep="\t")
  }
  
  
  smallestValue  <- min(apply(highly_variable_lcpm,2,min))
  biggestValue  <- max(apply(highly_variable_lcpm,2,max))
  range <- abs(smallestValue) + abs(biggestValue)
  divi <- range / 3
  first <- smallestValue + divi
  second <- first + 0.001
  third <- second + divi
  fourth <- third + 0.001
  
  
#  mypalette <- colorRampPalette(c("#25D811", "#FFFFFF", "#CE0C0C"))(n = 299)
  mypalette <- colorRampPalette(c("#1142d8", "#FFFFFF", "#CE0C0C"))(n = 299)
  col_breaks = c(seq(smallestValue,first,length = 100), 
                 seq(second,third,length = 100),
                 seq(fourth,biggestValue,length = 100))
  
  
  
  someMatrix <- as.matrix(highly_variable_lcpm)
  col.diagnosis <- data.frame(Diagnosis)
  
  rownames(col.diagnosis) <- colnames(someMatrix)
  
  
  pheatmap(someMatrix,
           scale="none",
           color = mypalette,
#           breaks = col_breaks,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           clustering_distance_rows = "correlation",#Pearson's
           clustering_method = "average",
           gaps_col = gapsToUse,
           cellheight = 6,
           cellwidth = 20,
           border_color=NA,
           fontsize_row = 6,
           filename = fileToUse,
           annotation_col = col.diagnosis
  )
  
  
}


getTableWithSignificantTranscripts <- function(dte, fit,index, sign){

  all <- which(dte[,index] != 0)
  table <- topTable(fit, genelist=fit$geneSymbols, number = length(all), adjust.method= "BH", sort.by="p")
  # fixed The table only have a subset of the names but the dte has all the elements as rows for example 2436 genes
  allNames <- rownames(dte)
  # fixed the table names only have the name of the significant genes 33 genes
  tableNames <- rownames(table)
    
  if(sign == "plus") {
    ids <- which(dte[,index] == 1)
  }else if(sign == "minus"){
    ids <- which(dte[,index] == -1)
  }else{
    ids <- c(which(dte[,index] == -1), which(dte[,index] == 1)) 
  }
  # fixed the ids is a numeric vector from which no names associated the ids are from the dte for example 15 down regulated
  namesOfIds  <- allNames[ids]
  # fixed extracted the names
  idsToExtract <- which(tableNames %in% namesOfIds)
  # fixed now we extract the names from the table

  newTable <- table[idsToExtract,]
  
  return(newTable)
}

reorderIdsUsingTypes <- function(sampleTable, types, sampleIds){
  newSampleTable <- sampleTable
  newSampleTable$ids <- c(1:dim(sampleTable)[1])
  newSampleTable <- newSampleTable[sampleIds,]
  first <- newSampleTable[which(newSampleTable$diagnosis == types[1]),]$ids
  second <- newSampleTable[which(newSampleTable$diagnosis == types[2]),]$ids
  reorderedIds <- c(first, second)
  

  theIds <- newSampleTable[which(reorderedIds %in% sampleIds),]$ids
  pass <- all(newSampleTable$ids == reorderedIds)
  
  if(!pass) print("An ERROR in reorderIds")

  return(reorderedIds)
}

removeFile <- function(fileName){
  #Define the file name that will be deleted
  fn <- fileName
  #Check its existence
  if (file.exists(fn)) 
    #Delete file if it exists
    file.remove(fn)
  
}


processProject <- function(yOriginalCounts, originalIds, projectFiles, nameOfProject, fileName,fileWithListOfDataUsed,
         significantTableName, types, lfc, pval,indexNumber, fitName, nameDte, sampleTable, 
         lcpm, projectDir, theHashFile="fastqFiles/tx2geneHTGHash.Rdata", showRowNames=TRUE){

  debug <- FALSE
  
  if(debug){
    yOriginalCounts <- yGene
    originalIds <- sampleIds
    theHashFile="fastqFiles/tx2geneHTGHash.Rdata"
    showRowNames=TRUE
  }
  
  fileWithListOfDataUsedTSV <- paste0(fileWithListOfDataUsed, ".tsv")
  fileWithListOfDataUsedExcel <- paste0(fileWithListOfDataUsed, ".xlsx")
  idsToUse <- reorderIdsUsingTypes(sampleTable, types, originalIds)
  downControlSignificant <- getTableWithSignificantTranscripts(nameDte, fitName, indexNumber, "minus" )
  upControlSignificant <- getTableWithSignificantTranscripts(nameDte, fitName, indexNumber, "plus" )
  topTableSigResults <- getTableWithSignificantTranscripts(nameDte, fitName, indexNumber, "all" )
  significantTranscripts <-rownames(topTableSigResults)

  saveFileAsExcel(topTableSigResults, significantTableName, "topSignificantHits")
  transcriptsToRemove <- generateFilesForProject(yOriginalCounts, originalIds , projectDir, 
                                                 projectFiles,nameOfProject, fitName, lfc, pval, significantTranscripts)
  inListToRemove <- which(significantTranscripts %in% transcriptsToRemove)
  transInListToRemoveTable <- topTableSigResults[inListToRemove,]
  transInListToRemove <- significantTranscripts[inListToRemove]
  generateHeatMap(significantTranscripts, originalIds, sampleTable, lcpm, types, 
                  fileName, transInListToRemove, fileWithListOfDataUsedTSV,theHashFile, TRUE)
  if(file.exists(fileWithListOfDataUsedTSV)){
    transfromTSVIntoXLSX(fileWithListOfDataUsedTSV, fileWithListOfDataUsedExcel, "heatMap", "geneSymbols")
  }

  removeFile(fileWithListOfDataUsedTSV)

}