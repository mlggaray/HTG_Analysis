#!/bin/bash

if [ $# -ne 2 ]; then
    echo "illegal number of parameters"
    echo "Script name: "$0
    echo "expectedArguments: salmon-full-path BASEDIR(parent directory where the decompressed fastqFiles.tar.bz2 are located as provided by Nancy Casanova) "
    echo " example: scripts/$0 /home/userName/bin/salmon /home/userName/HTGData > scripts/generateScriptToRunSalmon.sh"
    echo "chmod 755 scripts/generateScriptToRunSalmon.sh"
    echo "scripts/generateScriptToRunSalmon.sh > scripts/runSalmon.sh"
    echo "chmod 755 scripts/runSalmon.sh"
    exit 1
fi


SALMON=$1
BASEDIR=$2"/"
HTGPROBES=$BASEDIR"fastqFiles/probes/goodHTGProbes.fa"
SALMONLIB=$BASEDIR"fastqFiles/probes/goodHTGProbesLib"

cd $BASEDIR

echo "#!/bin/bash"
echo "SALMON=$SALMON"
echo "BASEDIR=$BASEDIR"
echo "HTGPROBES=$HTGPROBES"
echo "SALMONLIB=$SALMONLIB"
echo ""
echo "cd $BASEDIR"
echo ""
echo "rm -rf $BASEDIR"quantification""
echo ""
echo "mkdir $BASEDIR"quantification""
echo "mkdir $BASEDIR"quantification/lung""
echo "mkdir $BASEDIR"quantification/lung/cocci""
echo "mkdir $BASEDIR"quantification/lung/control""
echo "mkdir $BASEDIR"quantification/lung/sarc""
echo ""
echo "mkdir $BASEDIR"quantification/lymph""

echo "mkdir $BASEDIR"quantification/lymph/control""
echo "mkdir $BASEDIR"quantification/lymph/sarc""
echo "mkdir $BASEDIR"quantification/lymph/TB""
echo ""
echo 'echo $SALMON index -t $HTGPROBES -i $SALMONLIB'

find $BASEDIR -name '*.fastq.gz' -exec echo 'echo $SALMON quant -i $SALMONLIB -l SR -r ' {} {} \; |sed -e 's/.fastq.gz \/.*fastqFiles/.fastq.gz --output $BASEDIR"quantification/' -e 's/.fastq.gz$/_qual/' -e 's/\/.*fastqFiles/$BASEDIR"fastqFiles/'
