#!/bin/bash
#rm -rf ../conf_defocus_hsoff/*
#rm -rf ../conf_defocus_hson/*
#rm -rf ../conf_beamthrough/*
#rm -rf ../conf_beamthrough_highgain/*
set -o noclobber
#for run in {5808,5811,5814,5818,5821,5828} #pi- after BH2 modification
#for run in {5831,5837,5838,5840,5844,5851,5855,5858,5862,5864} #pi+ after BH2 modification
#for run in {5835,5842,5846,5847,5856,5860,5866} #proton after BH2 modification
#for run in {5217,5219,5220,5221} #pi- Vgem 316V
#for run in {5197,5205,5210,5215} #p Vgem 316V
#for run in {5191,5200,5207,5212} #pi+ Vgem 316V
#for run in {5191,5200,5207,5212,5197,5205,5210,5215,5217,5219,5220,5221} #highgian all
#for run in {5835,5842,5846,5847,5856,5860,5866,5831,5837,5838,5840,5844,5851,5855,5858,5862,5864,5808,5811,5814,5818,5821,5828} #all
################
for run in {5744..5762} #defocus hsoff pi-
#for run in {5723..5742} #defocus hson pi-
################
do
    file="../conf_defocus_hsoff/analyzer_0$run.conf"
    #file="../conf_defocus_hson/analyzer_0$run.conf"
    #file="../conf_beamthrough/analyzer_0$run.conf"
    #file="../conf_beamthrough_highgain/analyzer_0$run.conf"
    rm $file
    echo "$file"
    pK18=`grep pK18 /group/had/sks/E42/JPARC2021May/K18MagnetParam/K18MagParam_0$run | cut -c17-28`
    Kurama_F=`grep KURAMA_F /group/had/sks/E42/JPARC2021May/K18MagnetParam/K18MagParam_0$run | cut -c18-28`
    SHS_F=`grep SHS_F /group/had/sks/E42/JPARC2021May/K18MagnetParam/K18MagParam_0$run | cut -c18-28`
    echo "UNPACK:		param/UNPACK/unpacker.xml">>$file
    echo "DIGIT:		param/DIGIT/digit_data_e42_20210423.xml	">>$file
    echo "CMAP:		param/CMAP/channel_map_e42_20210423.xml">>$file
    echo "DCGEO:		param/DCGEO/DCGeomParam_E42_20231216">>$file
    echo "DCDRFT:		param/DCDRFT/DCDriftParam_E42_20230614">>$file
    echo "DCTDC:		param/DCTDC/DCTdcCalib_E42_20230614">>$file
    echo "HDPRM:		param/HDPRM/HodoParam_E42_20240228">>$file
    echo "HDPHC:		param/HDPHC/HodoPHCParam_htof">>$file
    echo "USER:		param/USER/UserParam_E42_beamthrough_pi">>$file
    echo "TPCPRM:		param/TPCPRM/TPCParam_20230825">>$file
    echo "TPCPOS:		param/TPCPOS/TPCPositionCorrectionMap_thetacalib">>$file
    echo "#TPCPOS:		param/TPCPOS/TPCPositionCorrectionMap_0">>$file
    echo "BH2FLT:		param/BH2FLT/BH2Filter.param">>$file
    echo "BH1MTH:		param/BH1MTH/BH1Matching.param">>$file
    echo "">>$file
    #echo "FLDCALC:	0.749779">>$file
    echo "FLDCALC:	0.744601">>$file
    echo "FLDNMR:		$Kurama_F">>$file
    echo "FLDMAP:		param/KURAMA/KuramaFieldMap_20210526">>$file
    echo "">>$file
    echo "HSFLDCALIB:	1.00049">>$file
    echo "HSFLDCALC:	1.0">>$file
    echo "HSFLDHALL:	$SHS_F">>$file
    #echo "HSFLDHALL:	1.0">>$file
    #echo "HSFLDMAP:	param/KURAMA/ShsFieldMap_20210526">>$file
    echo "HSFLDMAP:	param/KURAMA/ShsFieldMap_20210526_Adjusted">>$file
    echo "">>$file
    echo "#K18TM:		param/K18TM/K18MatrixParam_0.6_p_05311">>$file
    echo "PK18:		$pK18">>$file
    echo "K18TM:		param/K18TM/K18MatrixParam_2019Feb_plus">>$file
    echo "#TPCGDML:		param/geometry/hyptpcGeo_targetDiamond.gdml">>$file
    echo "TPCGDML:		param/geometry/hyptpcGeo_targetEmpty.gdml">>$file
    echo "#TPCGDML:		param/geometry/hyptpcGeo_targetCH2.gdml">>$file
done
