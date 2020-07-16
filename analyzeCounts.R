suppressMessages(library(dplyr))
suppressMessages(library(hash))
suppressMessages(library(limma))
suppressMessages(library(edgeR))

library(optparse)
options(java.parameters = "-Xmx4024m")
suppressMessages(library(r2excel))



stopQuietly <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "))
  stop(simpleError(blankMsg))
} # stopQuietly()




  
  
  option_list = list(
    make_option(c("-p", "--projectHome"), type="character", default=NULL,
                help="the full path of the projec directory. REQUIRED if you are in the directory use -p $( pwd -LP . )",
                metavar="character"),
    make_option(c("-a", "--analisisDir"), type="character", default=NULL,
                help="analisisDir"),
    make_option(c("-c", "--counts"), type="character", default="countsFromHTG.Rdata",
                help="counts file"),
    make_option(c("-t", "--tableAndDesign"), type="character", default="tableAndDesign.Rdata",
                help="tableAndDesign file"),
    make_option(c("-g", "--hashTable"), type="character", default="fastqFiles/tx2geneHTGHash.Rdata",
                help="hashTable"),
    make_option(c("-v", "--pvalue"), type="double", default=0.01,
                help="pvalue"),
    make_option(c("-l", "--log2FoldChange"), type="double", default=1.0,
                help="log2FoldChange")
  )
  
  
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  
  
  ### Rscript --vanilla htgDataProcessing.R -p "projectHome" -a "analysisDir" -c "countsFromHTG.Rdata" -t "tableAndDesign.Rdata" -g "fastqFiles/tx2geneHTGHash.Rdata" -v  0.01 -l 1.0
  ### Rscript --vanilla htgDataProcessing.R --help 
  
  
  if ( is.null(opt$projectHome) || is.null(opt$analisisDir) ){
    print_help(opt_parser)
    stop("Arguments must be supplied, use full path", call.= FALSE)
  }
  
  workingDir <- opt$projectHome
  setwd(workingDir)

  source("scripts/minRNASeqFuncLib.R")

  tableAndDesign <- opt$tableAndDesign
  load(tableAndDesign)
  countsFromHTG <- opt$counts
  load(countsFromHTG)

  projectDir <- paste0(opt$analisisDir, "/")
  if(!dir.exists(projectDir))
	  dir.create(projectDir)

  samplenames <- dataSheet$sampleName
  lcpm <- cpm(yGene, log=TRUE)  # just the log2 counts Counts per Million
  allRowNames <- samplenames

  hashTable <- opt$hashTable
  pval <- opt$pvalue
  lfc <- opt$log2FoldChange

  sampleTable <- data.frame(samplenames=samplenames, diagnosis=mydiag)
  rownames(sampleTable) <- sampleTable$samplenames

  lungSarcNumbers <- which(sampleTable$diagnosis %in% as.factor("lung_sarc"))
  lungCocciNumbers <- which(sampleTable$diagnosis %in% as.factor("lung_cocci"))
  lungControlNumbers <- which(sampleTable$diagnosis %in% as.factor("lung_control"))
  lymphSarcNumbers <- which(sampleTable$diagnosis %in% as.factor("lymph_sarc"))
  lymphControlNumbers <- which(sampleTable$diagnosis %in% as.factor("lymph_control"))
  lymphTBNumbers <- which(sampleTable$diagnosis %in% as.factor("lymph_TB"))

lungCoccivslungControlNumbers <- c(lungCocciNumbers, lungControlNumbers)
lungCoccivslungControl <- lcpm[, lungCoccivslungControlNumbers]

lungCoccivslungControlSampleTable <- data.frame(samplenames=samplenames[lungCoccivslungControlNumbers], diagnosis=as.character(mydiag[lungCoccivslungControlNumbers]))
rownames(lungCoccivslungControlSampleTable) <- lungCoccivslungControlSampleTable$samplenames
diag <- as.factor(lungCoccivslungControlSampleTable$diagnosis)
designlungCoccivslungControl <- model.matrix(~0+diag)
colnames(designlungCoccivslungControl) <- gsub("diag", "", colnames(designlungCoccivslungControl))


contr.matrixControl <- makeContrasts(
  lungCoccivslungControl = lung_cocci - lung_control,
  levels = colnames(designlungCoccivslungControl))

lungCoccivslungControlGroup <- lcpm[, lungCoccivslungControlNumbers]
yLungCoccivslungControl <- yGene[, lungCoccivslungControlNumbers]


vVoomLungCoccivslungControl <- voom(yLungCoccivslungControl, designlungCoccivslungControl, plot=FALSE)
# vVoom is now ready for lmFit() see limma User's Guide

vVoomLungCoccivslungControl <- lmFit(vVoomLungCoccivslungControl, designlungCoccivslungControl)
vVoomLungCoccivslungControl <- contrasts.fit(vVoomLungCoccivslungControl, contrasts=contr.matrixControl)
efitLungCoccivslungControl <- eBayes(vVoomLungCoccivslungControl)
efitLungCoccivslungControl$geneSymbols <- names(efitLungCoccivslungControl$genes)
dteLungCoccivslungControl <- decideTests(efitLungCoccivslungControl, method = "separate", adjust.method = "none", p.value = pval, lfc=lfc)


# a simplify method to get the transcripts
sampleIds <- lungCoccivslungControlNumbers
projectFiles="lungCocciVsLungControl"
nameOfProject <- "Lung Cocci Vs Lung Control"
fileName <- paste0(projectDir, projectFiles, "_heatmap.png")
fileWithListOfDataUsed <- paste0(projectDir, projectFiles, "_heatmap_data")
significantTableName <- paste0(projectDir, projectFiles, "_significantTable.xlsx")

types <- c("lung_cocci", "lung_control")
lfc = lfc
pval = pval
indexNumber <- 1
fitName <- efitLungCoccivslungControl[,indexNumber]
nameDte <- dteLungCoccivslungControl
summaryDte <- summary(nameDte)
summaryDte 
found <- summaryDte[1,1] + summaryDte[3,1]


if(found > 1 && findMinimumSample(sampleTable, types) > 1){
sigTranscripts <- processProject(yGene,sampleIds, projectFiles, nameOfProject, fileName,fileWithListOfDataUsed,
                           significantTableName, types, lfc, pval,indexNumber, fitName, 
               nameDte,sampleTable, lcpm, projectDir, hashTable)
}

###----------------------  lung control vs lung sarco
allRowNames <- samplenames

sampleTable <- data.frame(samplenames=samplenames, diagnosis=as.character(mydiag))
rownames(sampleTable) <- sampleTable$samplenames

lungSarcovslungControlNumbers <- c(lungSarcNumbers, lungControlNumbers)
lungSarcovslungControl <- lcpm[, lungSarcovslungControlNumbers]

lungSarcovslungControlSampleTable <- data.frame(samplenames=samplenames[lungSarcovslungControlNumbers], diagnosis=as.character(mydiag[lungSarcovslungControlNumbers]))
rownames(lungSarcovslungControlSampleTable) <- lungSarcovslungControlSampleTable$samplenames
diag <- as.factor(lungSarcovslungControlSampleTable$diagnosis)
designlungSarcovslungControl <- model.matrix(~0+diag)
colnames(designlungSarcovslungControl) <- gsub("diag", "", colnames(designlungSarcovslungControl))


contr.matrixControl <- makeContrasts(
  lungSarcovslungControl = lung_sarc - lung_control,
  levels = colnames(designlungSarcovslungControl))

lungSarcovslungControlGroup <- lcpm[, lungSarcovslungControlNumbers]
yLungSarcovslungControl <- yGene[, lungSarcovslungControlNumbers]


vVoomLungSarcovslungControl <- voom(yLungSarcovslungControl, designlungSarcovslungControl, plot=FALSE)
# vVoom is now ready for lmFit() see limma User's Guide

vVoomLungSarcovslungControl <- lmFit(vVoomLungSarcovslungControl, designlungSarcovslungControl)
vVoomLungSarcovslungControl <- contrasts.fit(vVoomLungSarcovslungControl, contrasts=contr.matrixControl)
efitLungSarcovslungControl <- eBayes(vVoomLungSarcovslungControl)
efitLungSarcovslungControl$geneSymbols <- names(efitLungSarcovslungControl$genes)
dteLungSarcovslungControl <- decideTests(efitLungSarcovslungControl, method = "separate", adjust.method = "none", p.value = pval, lfc=lfc)


# a simplify method to get the transcripts
sampleIds <- lungSarcovslungControlNumbers
projectFiles <- "lungSarcoVsLungControl"
nameOfProject <- "Lung Sarco Vs Lung Control"
fileName <- paste0(projectDir, projectFiles, "_heatmap.png")
fileWithListOfDataUsed <- paste0(projectDir, projectFiles, "_heatmap_data")
significantTableName <- paste0(projectDir, projectFiles, "_significantTable.xlsx")
types <- c("lung_sarc", "lung_control")
lfc = lfc
pval = pval
indexNumber <- 1
fitName <- efitLungSarcovslungControl[,indexNumber]
nameDte <- dteLungSarcovslungControl
summaryDte <- summary(nameDte)
summaryDte
found <- summaryDte[1,1] + summaryDte[3,1]


if(found > 1 && findMinimumSample(sampleTable, types) > 1){
  sigTranscripts <- processProject(yGene,sampleIds, projectFiles, nameOfProject, fileName,fileWithListOfDataUsed,
                                   significantTableName, types, lfc, pval,indexNumber, fitName, 
                                   nameDte,sampleTable, lcpm, projectDir, hashTable)
}

###----------------------  lymphSarc vs lymphControl

allRowNames <- samplenames

sampleTable <- data.frame(samplenames=samplenames, diagnosis=as.character(mydiag))
rownames(sampleTable) <- sampleTable$samplenames
lymphSarcovslymphControlNumbers <- c(lymphSarcNumbers, lymphControlNumbers)
lymphSarcovslymphControl <- lcpm[, lymphSarcovslymphControlNumbers]

lymphSarcovslymphControlSampleTable <- data.frame(samplenames=samplenames[lymphSarcovslymphControlNumbers], diagnosis=as.character(mydiag[lymphSarcovslymphControlNumbers]))
rownames(lymphSarcovslymphControlSampleTable) <- lymphSarcovslymphControlSampleTable$samplenames
diag <- as.factor(lymphSarcovslymphControlSampleTable$diagnosis)
designlymphSarcovslymphControl <- model.matrix(~0+diag)
colnames(designlymphSarcovslymphControl) <- gsub("diag", "", colnames(designlymphSarcovslymphControl))


contr.matrixControl <- makeContrasts(
  lymphSarcovslymphControl = lymph_sarc - lymph_control,
  levels = colnames(designlymphSarcovslymphControl))

lymphSarcovslymphControlGroup <- lcpm[, lymphSarcovslymphControlNumbers]
ylymphSarcovslymphControl <- yGene[, lymphSarcovslymphControlNumbers]


vVoomlymphSarcovslymphControl <- voom(ylymphSarcovslymphControl, designlymphSarcovslymphControl, plot=FALSE)
# vVoom is now ready for lmFit() see limma User's Guide

vVoomlymphSarcovslymphControl <- lmFit(vVoomlymphSarcovslymphControl, designlymphSarcovslymphControl)
vVoomlymphSarcovslymphControl <- contrasts.fit(vVoomlymphSarcovslymphControl, contrasts=contr.matrixControl)
efitlymphSarcovslymphControl <- eBayes(vVoomlymphSarcovslymphControl)
efitlymphSarcovslymphControl$geneSymbols <- names(efitlymphSarcovslymphControl$genes)
dtelymphSarcovslymphControl <- decideTests(efitlymphSarcovslymphControl, method = "separate", adjust.method = "none", p.value = pval, lfc=lfc)


# a simplify method to get the transcripts
sampleIds <- lymphSarcovslymphControlNumbers
projectFiles <- "lymphSarcovslymphControl"
nameOfProject <- "Lymph Sarco Vs Lymph Control"
fileName <- paste0(projectDir, projectFiles, "_heatmap.png")
fileWithListOfDataUsed <- paste0(projectDir, projectFiles, "_heatmap_data")
significantTableName <- paste0(projectDir, projectFiles, "_significantTable.xlsx")
types <- c("lymph_sarc", "lymph_control")
lfc = lfc
pval = pval

indexNumber <- 1
fitName <- efitlymphSarcovslymphControl[,indexNumber]
nameDte <- dtelymphSarcovslymphControl
summaryDte <- summary(nameDte)
summaryDte
found <- summaryDte[1,1] + summaryDte[3,1]


if(found > 1 && findMinimumSample(sampleTable, types) > 1){
  sigTranscripts <- processProject(yGene,sampleIds, projectFiles, nameOfProject, fileName,fileWithListOfDataUsed,
                                   significantTableName, types, lfc, pval,indexNumber, fitName, 
                                   nameDte,sampleTable, lcpm, projectDir, hashTable)
}


# lymphTB vs lymphControl
allRowNames <- samplenames

sampleTable <- data.frame(samplenames=samplenames, diagnosis=as.character(mydiag))
rownames(sampleTable) <- sampleTable$samplenames
lymphTBvslymphControlNumbers <- c(lymphTBNumbers, lymphControlNumbers)
lymphTBvslymphControl <- lcpm[, lymphTBvslymphControlNumbers]

lymphTBvslymphControlSampleTable <- data.frame(samplenames=samplenames[lymphTBvslymphControlNumbers], diagnosis=as.character(mydiag[lymphTBvslymphControlNumbers]))
rownames(lymphTBvslymphControlSampleTable) <- lymphTBvslymphControlSampleTable$samplenames
diag <- as.factor(lymphTBvslymphControlSampleTable$diagnosis)
designlymphTBvslymphControl <- model.matrix(~0+diag)
colnames(designlymphTBvslymphControl) <- gsub("diag", "", colnames(designlymphTBvslymphControl))


contr.matrixControl <- makeContrasts(
  lymphTBvslymphControl = lymph_TB - lymph_control,
  levels = colnames(designlymphTBvslymphControl))

lymphTBvslymphControlGroup <- lcpm[, lymphTBvslymphControlNumbers]
ylymphTBvslymphControl <- yGene[, lymphTBvslymphControlNumbers]


vVoomlymphTBvslymphControl <- voom(ylymphTBvslymphControl, designlymphTBvslymphControl, plot=FALSE)
# vVoom is now ready for lmFit() see limma User's Guide

vVoomlymphTBvslymphControl <- lmFit(vVoomlymphTBvslymphControl, designlymphTBvslymphControl)
vVoomlymphTBvslymphControl <- contrasts.fit(vVoomlymphTBvslymphControl, contrasts=contr.matrixControl)
efitlymphTBvslymphControl <- eBayes(vVoomlymphTBvslymphControl)
efitlymphTBvslymphControl$geneSymbols <- names(efitlymphTBvslymphControl$genes)
dtelymphTBvslymphControl <- decideTests(efitlymphTBvslymphControl, method = "separate", adjust.method = "none", p.value = pval, lfc=lfc)


# a simplify method to get the transcripts
sampleIds <- lymphTBvslymphControlNumbers
projectFiles <- "lymphTBvslymphControl"
nameOfProject <- "Lymph TB Vs Lymph Control"
fileName <- paste0(projectDir, projectFiles, "_heatmap.png")
fileWithListOfDataUsed <- paste0(projectDir, projectFiles, "_heatmap_data")
significantTableName <- paste0(projectDir, projectFiles, "_significantTable.xlsx")
types <- c("lymph_TB", "lymph_control")
lfc = lfc
pval = pval
indexNumber <- 1
fitName <- efitlymphTBvslymphControl[,indexNumber]
nameDte <- dtelymphTBvslymphControl
summaryDte <- summary(nameDte)
summaryDte
found <- summaryDte[1,1] + summaryDte[3,1]


if(found > 1 && findMinimumSample(sampleTable, types) > 1){
  sigTranscripts <- processProject(yGene,sampleIds, projectFiles, nameOfProject, fileName,fileWithListOfDataUsed,
                                   significantTableName, types, lfc, pval,indexNumber, fitName, 
                                   nameDte,sampleTable, lcpm, projectDir, hashTable)
}


