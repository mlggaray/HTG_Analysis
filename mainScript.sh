#!/bin/bash

if [ $# -ne 3 ]; then
    echo "illegal number of parameters"
    echo "Script name: "$0
    echo "expectedArguments: salmon-full-path to binary, test that the binary is working"
    echo "BASEDIR(parent directory where the decompressed fastqFiles.tar.bz2 are located as provided by Nancy Casanova) "
    echo "analysis directory name for example <analysisResults>"
    exit 1
fi


SALMON=$1
currentDir=$2
ANALISISDIR=$3


cd $currentDir



scripts/generateTemplate.sh $SALMON $currentDir > scripts/generateScriptToRunSalmon.sh
chmod 755 scripts/generateScriptToRunSalmon.sh
scripts/generateScriptToRunSalmon.sh > scripts/runSalmon.sh
chmod 755 scripts/runSalmon.sh
scripts/runSalmon.sh
scripts/generateRFiles.sh $currentDir
Rscript --vanilla scripts/htgDataProcessing.R -p $currentDir -f fileInfo.tsv -g fastqFiles/probesToGenes.csv
Rscript --vanilla scripts/analyzeCounts.R -p $currentDir -a $ANALISISDIR -c "countsFromHTG.Rdata" -t "tableAndDesign.Rdata" -g "fastqFiles/tx2geneHTGHash.Rdata" -v  0.01 -l 1.0
