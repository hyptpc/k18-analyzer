#!/bin/sh
DATADIR=$E73DATADIR/Run85
ROOTDIR=$E73ROOTDIR
ANADIR=$ANALYZERDIR/e73
echo $DATADIR
PROG=$1
RUN=$(printf %05d $3)

BIN=$ANADIR/bin/$PROG
CONF=$ANADIR/param/conf/analyzer.conf
INPUT=$DATADIR/run$RUN.dat
#OUTPUT=$ROOTDIR/$PROG\_$RUN.root
OUTPUT=$ROOTDIR/$2\_$RUN.root
#CDCTREE=$ROOTDIR/CDCTree_$RUN.root
CDCTREE=$ROOTDIR/ExCDCTree_$RUN.root

MAXLOOP=-1
echo $#
if [ $# -eq 3 ]; then
MAXLOOP=$4
fi
if [ $# -eq 4 ]; then
MAXLOOP=$4
SKIP=$5
fi

if [ -r $INPUT ]; then
    $BIN --conf=$CONF --infile=$INPUT --outfile=$OUTPUT --maxloop=$MAXLOOP --skip=$SKIP --cdctree=$CDCTREE --run=$RUN
elif [ -r $INPUT.gz ]; then
    $BIN --conf=$CONF --infile=$INPUT.gz --outfile=$OUTPUT --maxloop=$MAXLOOP --skip=$SKIP --cdctree=$CDCTREE --run=$RUN
else
    echo $INPUT "does not exist"
fi
