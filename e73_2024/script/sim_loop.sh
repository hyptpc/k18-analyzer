#!/bin/sh
echo $#
if [ $# -le 2 ]; then
    echo "usage: bash ana.sh [program] [output prefix] [run] ([maxloop] [skip])"
    exit
fi

DATADIR=$KNUCLROOTDIR/geant/e73_quasifree/
ROOTDIR=$DATADIR
ANADIR=$ANALYZERDIR/e73
echo $DATADIR
LOGDIR=$ROOTDIR/log
mkdir -p $LOGDIR

PROG=$1
BIN=$ANADIR/bin/$PROG
CONF=$ANADIR/param/conf/analyzer_e73.conf
Q=s

MAXLOOP=-1
if [ $# -eq 5 ]; then
MAXLOOP=$5
fi
if [ $# -eq 6 ]; then
MAXLOOP=$5
SKIP=$6
fi

for run in $(seq $3 $4)
do
RUN=$(printf %03d $run)
INPUT=$DATADIR/$2_$RUN.root
OUTPUT=$ROOTDIR/$2_${PROG}_$RUN.root
CDCTREE=$ROOTDIR/$2_ExCDCTree_$RUN.root
LOGNAME=$LOGDIR/$PROG\_$RUN.log
if [ -r $INPUT ]; then
    bsub -q $Q -o $LOGNAME $BIN --conf=$CONF --infile=$INPUT --outfile=$OUTPUT --maxloop=$MAXLOOP --skip=$SKIP  --cdctree=$CDCTREE --run=$RUN
elif [ -r $INPUT.gz ]; then
    bsub -q $Q -o $LOGNAME $BIN --conf=$CONF --infile=$INPUT.gz --outfile=$OUTPUT --maxloop=$MAXLOOP --skip=$SKIP  --cdctree=$CDCTREE --run=$RUN
else
    echo $INPUT "does not exist"
fi
done
