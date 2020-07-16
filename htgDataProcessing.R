suppressMessages(library(tximport))
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
                help="the full path of the projec directory. REQUIRED "),
    make_option(c("-f", "--fileInfo"), type="character", default=NULL,
                help="fileInfo"),
    make_option(c("-g", "--geneModels"), type="character", default=NULL,
                help="geneModels"),
    make_option(c("-s", "--serializeData"), type="logical", 
                action="store_true", default=TRUE, 
                dest="serializeData", help="an optional flag to save serialized data")
  )
  
  
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  
  
  ### Rscript --vanilla htgDataProcessing.R -p "projectHome" -f "fileInfo.tsv" -g "fastqFiles/probesToGenes.csv" 
  ### Rscript --vanilla htgDataProcessing.R --help 
  
  
  if ( is.null(opt$projectHome) || is.null(opt$fileInfo) || is.null(opt$geneModels) ){
    print_help(opt_parser)
    stop("Arguments must be supplied, use full path", call.= FALSE)
  }
  
  project_home <- opt$projectHome
  print(project_home)
  fileInfo <- opt$fileInfo
  print(fileInfo)
  fileInfo <- paste0(project_home, "/", fileInfo) 
  print(fileInfo)
  geneModel <- paste0(project_home, "/", opt$geneModels) 
  print(geneModel)
  serializeData <- opt$serialize


# check that the project_home exists
if(dir.exists(project_home)){
  cat("The HTG project exists!", sep="\n", file = stderr())
} else{
  cat("Bad or missing HTG project path", sep="\n", file = stderr())
  stopQuietly()
}



dataSheet <- read.csv(fileInfo, stringsAsFactors = FALSE, sep = "\t")
files <- dataSheet$location
samplenames <- dataSheet$sampleName
tissue <- dataSheet$tissue
diagnosis <- dataSheet$diagnosis
class <- paste0(dataSheet$tissue, "_", dataSheet$diagnosis)
mydiag <- as.factor(class)
all(file.exists(files))
design <- model.matrix(~0+mydiag)
colnames(design) <- gsub("mydiag", "", colnames(design))
tx2gene <- read.csv(geneModel, stringsAsFactors = FALSE,  sep="\t")
colnames(tx2gene) <- c("ensembl_transcript_id", "GENEID")
txiGene <- tximport(files, type = "salmon", tx2gene = tx2gene)
colnames(txiGene$counts) <- samplenames
colnames(txiGene$abundance) <- samplenames
colnames(txiGene$length) <- samplenames

cts <- txiGene$counts
normMat <- txiGene$length
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
yGene <- DGEList(cts)
yGene$offset <- t(t(log(normMat)) + o)
rownames(yGene$samples)

yGene$samples$group <- mydiag
yGene <- calcNormFactors(yGene)
yGene$genes <- apply(txiGene$length,1, mean, na.rm = TRUE)
dim(yGene)



if(serializeData){
	save(yGene, file="countsFromHTG.Rdata")	
	save(dataSheet, mydiag, design, file="tableAndDesign.Rdata")
}

