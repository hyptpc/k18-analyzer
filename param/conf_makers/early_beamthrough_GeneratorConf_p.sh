#!/bin/bash
set -o noclobber
for run in {5197,5205,5210,5215} #p Vgem 316V
################
do
    file="../conf_beamthrough_highgain/analyzer_0$run.conf"
    echo "$file"
    pK18=`grep pK18 /group/had/sks/E42/JPARC2021May/K18MagnetParam/K18MagParam_0$run | cut -c17-28`
    Kurama_F=`grep KURAMA_F /group/had/sks/E42/JPARC2021May/K18MagnetParam/K18MagParam_0$run | cut -c18-28`
    SHS_F=`grep SHS_F /group/had/sks/E42/JPARC2021May/K18MagnetParam/K18MagParam_0$run | cut -c18-28`
    echo "UNPACK:		param/UNPACK/unpacker.xml">>$file
    echo "DIGIT:		param/DIGIT/digit_data_e42_20210423.xml	">>$file
    echo "CMAP:		param/CMAP/channel_map_e42_20210423.xml">>$file
    #echo "DCGEO:		param/DCGEO/DCGeomParam_E42_20230614">>$file
    #echo "DCGEO:		param/DCGEO/DCGeomParam_E42_20230725">>$file
    echo "DCGEO:		param/DCGEO/DCGeomParam_E42_20230910">>$file
    echo "DCDRFT:		param/DCDRFT/DCDriftParam_E42_20230614">>$file
    echo "DCTDC:		param/DCTDC/DCTdcCalib_E42_20230614">>$file
    #echo "HDPRM:		param/HDPRM/HodoParam_E42_20220316">>$file
    echo "HDPRM:		param/HDPRM/HodoParam_E42_ch2target">>$file
    echo "HDPHC:		param/HDPHC/HodoPHCParam_htof">>$file
    #echo "USER:		param/USER/UserParam_E42_early_beamthrough_pi">>$file
    echo "USER:		param/USER/UserParam_E42_early_beamthrough_p">>$file
    echo "TPCPRM:		param/TPCPRM/TPCParam_20230825">>$file
    #echo "TPCPRM:		param/TPCPRM/TPCParam_20230825_NoPadOff">>$file
    echo "TPCPOS:		param/TPCPOS/TPCPositionCorrectionMap_thetacalib">>$file
    echo "#TPCPOS:		param/TPCPOS/TPCPositionCorrectionMap_0">>$file
    echo "BH2FLT:		param/BH2FLT/BH2Filter.param">>$file
    echo "BH1MTH:		param/BH1MTH/BH1Matching.param">>$file
    echo "">>$file
    echo "FLDCALC:	0.749779">>$file
    echo "FLDNMR:		$Kurama_F">>$file
    echo "FLDMAP:		param/KURAMA/KuramaFieldMap_20210526">>$file
    echo "">>$file
    echo "HSFLDCALIB:	0.94838">>$file
    #echo "HSFLDCALIB:	1.0">>$file
    echo "HSFLDCALC:	0.90607">>$file
    echo "HSFLDHALL:	$SHS_F">>$file
    echo "HSFLDMAP:	param/KURAMA/ShsFieldMap_20210526">>$file
    echo "">>$file
    echo "#PK18:		0.5970048">>$file
    echo "#K18TM:		param/K18TM/K18MatrixParam_0.6_p_05311">>$file
    echo "PK18:		$pK18">>$file
    echo "K18TM:		param/K18TM/K18MatrixParam_2019Feb_plus">>$file
    echo "#TPCGDML:		param/geometry/hyptpcGeo_targetDiamond.gdml">>$file
    echo "TPCGDML:		param/geometry/hyptpcGeo_targetEmpty.gdml">>$file
    echo "#TPCGDML:		param/geometry/hyptpcGeo_targetCH2.gdml">>$file
done
