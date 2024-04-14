#!/bin/sh
echo $#
if [ $# -le 2 ]; then
    echo "usage: bash sim.sh [datadir] [program] [prefix] [run] ([maxloop] [skip])"
    exit
fi

DATADIR=$1
ROOTDIR=$DATADIR
ANADIR=$ANALYZERDIR/e73

RESOL=0.005
RESOL=0.018

echo $DATADIR
PROG=$2
RUN=$(printf %03d $4)

BIN=$ANADIR/bin/$PROG
#CONF=$ANADIR/param/conf/analyzer_e73.conf
#CONF=$ANADIR/param/conf/analyzer_noresol.conf
CONF=$ANADIR/param/conf/analyzer_e80.conf

tmp=$(echo "$RESOL * 1000" | bc)
tmp=$(printf %03d $tmp)

INPUT=$DATADIR/$3_$RUN.root
#OUTPUT=$ROOTDIR/$3_${PROG}_$RUN\_resol$tmp.root
OUTPUT=$ROOTDIR/$3_${PROG}_$RUN.root
CDCTREE=$ROOTDIR/$3_ExCDCTree_$RUN.root

MAXLOOP=-1
if [ $# -eq 5 ]; then
MAXLOOP=$5
fi
if [ $# -eq 6 ]; then
MAXLOOP=$5
SKIP=$6
fi

if [ -r $INPUT ]; then
    $BIN --conf=$CONF --infile=$INPUT --outfile=$OUTPUT --maxloop=$MAXLOOP --skip=$SKIP  --cdctree=$CDCTREE --run=$RUN --cdcresol=$RESOL
elif [ -r $INPUT.gz ]; then
    $BIN --conf=$CONF --infile=$INPUT.gz --outfile=$OUTPUT --maxloop=$MAXLOOP --skip=$SKIP  --cdctree=$CDCTREE --run=$RUN --cdcresol=$RESOL
else
    echo $INPUT "does not exist"
fi
