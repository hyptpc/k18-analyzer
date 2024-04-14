#!/bin/sh
echo $#
if [ $# -le 2 ]; then
    echo "usage: bash ana.sh [program] [output prefix] [run] (--maxloop [maxloop] --skip [skip] --mc --cosmic)"
    exit
fi

DATADIR=$E73DATADIR/Run91
ROOTDIR=$HOME/data/e73_2024
ANADIR=$ANALYZERDIR/e73
echo $DATADIR

OPT=""
PROG=$1
OUTPRE=$2
RUN=$(printf %05d $3)
shift 3
while (( $# > 0 ))
do
    case $1 in
	--cosmic | --mc)
	    OPT+="$1 "
	    ;;
	--maxloop | --skip)
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "$1 requires an argument." 1>&2
		exit 1
	    else
		OPT+="$1=$2 "
		shift
	    fi
	    ;;
	*)
	    echo "invalid argument $1"
	    exit 1
	    ;;
    esac
    shift
done
echo $OPT

BIN=$ANADIR/bin/$PROG
CONF=$ANADIR/param/conf/analyzer_e73_2024.conf
#CONF=$ANADIR/param/conf/analyzer.online
INPUT=$DATADIR/run$RUN.dat
OUTPUT=$ROOTDIR/${OUTPRE}_$RUN.root
CDCTREE=$ROOTDIR/ExCDCTree_$RUN.root

if [ -r $INPUT ]; then
    $BIN --conf=$CONF --infile=$INPUT --outfile=$OUTPUT --cdctree=$CDCTREE $OPT --run=$RUN
elif [ -r $INPUT.gz ]; then
    echo $BIN --conf=$CONF --infile=$INPUT.gz --outfile=$OUTPUT --cdctree=$CDCTREE $OPT --run=$RUN
    $BIN --conf=$CONF --infile=$INPUT.gz --outfile=$OUTPUT --cdctree=$CDCTREE $OPT --run=$RUN
else
    echo $INPUT "does not exist"
fi
