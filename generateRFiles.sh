#!/bin/bash

if [ $# -ne 1 ]; then
    echo "illegal number of parameters"
    echo "Script name: "$0
    echo "expectedArguments: BASEDIR(directory where the decompressed htg fastq files are located as provided by Nancy Casanova) "
    exit 1
fi


BASEDIR=$1"/"

cd $BASEDIR

echo "#psudoID sampleName diagnosis tissue location"|sed -e 's/ /\t/g' > fileInfoHeader.txt
find $BASEDIR  -name 'quant.sf' -exec echo {} {} \;|sed -e 's/ \/.*quantification\//\t/' -e 's/\tlung\/sarc\//\tlung\tsarc\t/' -e 's/\tlymph\/sarc\//\tlymph\tsarc\t/' -e 's/\tlung\/cocci\//\tlung\tcocci\t/' -e 's/\tlung\/control\//\tlung\tcontrol\t/' -e 's/\tlymph\/TB\//\tlymph\tTB\t/' -e 's/\tlymph\/control\//\tlymph\tcontrol\t/' -e 's/_qual\/quant.sf$//' |awk -F '\t' '{print $4"\t"$3"\t"$2"\t"$1}' > fileInfoData.txt

cat fileInfoHeader.txt fileInfoData.txt > fileInfoTempFile.txt
scripts/reNumberFile.rb fileInfoTempFile.txt fileInfo.tsv

rm -rf fileInfoHeader.txt fileInfoData.txt fileInfoTempFile.txt
