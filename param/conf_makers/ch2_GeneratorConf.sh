#!/bin/bash
#rm -rf ../conf_ch2target/*
set -o noclobber
for run in {5641..5666}
do
    file="../conf_ch2target/analyzer_0$run.conf"
    #file="../conf_beamstudy/analyzer_0$run.conf"
    echo "$file"
    rm -rf $file
    pK18=`grep pK18 /group/had/sks/E42/JPARC2021May/K18MagnetParam/K18MagParam_0$run | cut -c17-28`
    Kurama_F=`grep KURAMA_F /group/had/sks/E42/JPARC2021May/K18MagnetParam/K18MagParam_0$run | cut -c18-28`
    SHS_F=`grep SHS_F /group/had/sks/E42/JPARC2021May/K18MagnetParam/K18MagParam_0$run | cut -c18-28`
    echo "UNPACK:		param/UNPACK/unpacker.xml">>$file
    echo "DIGIT:		param/DIGIT/digit_data_e42_20210423.xml	">>$file
    echo "CMAP:		param/CMAP/channel_map_e42_20210423.xml">>$file
    echo "DCGEO:		param/DCGEO/DCGeomParam_E42_20231216">>$file
    #echo "DCGEO:		param/DCGEO/DCGeomParam_E42_beamstudy">>$file
    echo "DCDRFT:		param/DCDRFT/DCDriftParam_E42_20230614">>$file
    echo "DCTDC:		param/DCTDC/DCTdcCalib_E42_20230614">>$file
    echo "HDPRM:		param/HDPRM/HodoParam_E42_20240228">>$file
    echo "HDPHC:		param/HDPHC/HodoPHCParam_htof">>$file
    echo "USER:		param/USER/UserParam_E42_htof">>$file
    echo "TPCPRM:		param/TPCPRM/TPCParam_20230825">>$file
    echo "TPCPOS:		param/TPCPOS/TPCPositionCorrectionMap_thetacalib">>$file
    echo "BH2FLT:		param/BH2FLT/BH2Filter.param">>$file
    echo "BH1MTH:		param/BH1MTH/BH1Matching.param">>$file
    echo "">>$file
    #echo "#FLDCALC:	0.749779 #KURAMA">>$file
    echo "FLDCALC:	0.744601 #KURAMA+HS">>$file
    echo "FLDNMR:		$Kurama_F">>$file
    echo "FLDMAP:		param/KURAMA/KuramaFieldMap_20210526">>$file
    echo "">>$file
    echo "#HSFLDCALIB:	0.993623">>$file
    echo "HSFLDCALIB:	1.00049">>$file
    echo "HSFLDCALC:	1.0">>$file
    #echo "HSFLDHALL:	$SHS_F">>$file
    echo "HSFLDHALL:	1.0">>$file
    echo "#HSFLDMAP:	param/KURAMA/ShsFieldMap_20210526">>$file
    echo "#HSFLDMAP:	param/KURAMA/ShsFieldMap_20210526_Extrapolated">>$file
    echo "HSFLDMAP:	param/KURAMA/ShsFieldMap_20210526_Adjusted">>$file
    echo "">>$file
    echo "#PK18:		0.5970048">>$file
    echo "#K18TM:		param/K18TM/K18MatrixParam_0.6_p_05311">>$file
    echo "PK18:		$pK18">>$file
    echo "K18TM:		param/K18TM/K18MatrixParam_2019Feb_plus">>$file
    echo "#TPCGDML:		param/geometry/hyptpcGeo_targetDiamond.gdml">>$file
    echo "#TPCGDML:		param/geometry/hyptpcGeo_targetEmpty.gdml">>$file
    echo "TPCGDML:		param/geometry/hyptpcGeo_targetCH2.gdml">>$file
done
