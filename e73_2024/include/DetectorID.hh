#ifndef DETECTOR_ID_HH
#define DETECTOR_ID_HH

#include <iostream>
#include <Rtypes.h>

// Counters ___________________________________________________________
const int DetIdCDH    =   0;
const int DetIdBHD    =   1;
const int DetIdBHT    =   1;
const int DetIdT0     =   2;
const int DetIdAC     =   3;
const int DetIdT0new  =   4;
const int DetIdT1     =   4;
const int DetIdE0     =   5;
const int DetIdDEF    =   6;
const int DetIdStart  =  11;
const int DetIdStop   =  12;
const int DetIdCVC    =  14;
const int DetIdSDD    =  20;
const int DetIdSDDVeto=  25;
const int DetIdNC     =  32;
const int DetIdBVC    =  33;
const int DetIdPC     =  35;
const int DetIdIH     =  42;
const int DetIdVFT    =  53;

const Int_t NumOfSegAC = 1;

const int DetIdPbG    =  29;
const int DetIdPbF2   =  30;
const int DetIdVeto   =  31;

const int DetIdBTC    =  33;
const int DetIdCHCbarrel  =  61;
const int DetIdCHCcapF    =  62;
const int DetIdCHCcapB    =  63;
const int DetIdNCbarrel  =  66;
#if E73
const int DetIdVeto1  =  31;
const int DetIdVeto0  =  32;
const int DetIdFinger =  32;
const int DetIdLeak   =  34;
#endif
#if T98
const int DetIdRC     =  98;
#endif
#if E15
const int DetIdNC     =  32;
const int DetIdLB     =  36;
const int DetIdWVC    =  37;
const int DetIdBPD    =  41;
const int DetIdBD     =  90;
#endif
// const int NumOfSegBHD    =  11;
// const int NumOfSegT0     =  5;
// const int NumOfSegT0new  =  1;
const int DetIdVmeRm     =  81;
const int DetIdScaler    =  91;
const int DetIdTrigFlag  =  99;
const Int_t NumOfSegTrigFlag = 32;

// Chambers
const int DetIdCDC    = 100;
const int DetIdBLC1a  = 101;
const int DetIdBLC1b  = 102;
const int DetIdBLC2a  = 103;
const int DetIdBLC2b  = 104;
const int DetIdBPC    = 105;
const int DetIdBPC2   = 105;
const int DetIdBPC1   = 106;
const int DetIdBPCmini= 106;
const int DetIdSDC    = 106;
const int DetIdFDC    = 107;
const int DetIdBLC1   = 111;
const int DetIdBLC2   = 112;
const int DetIdBPC0   = 113;

// Hall etc.
const int DetIdHall       	= 200;
const int DetIdFloor		= 201;
// const int DetIdBeamDump		= 110;
// const int DetIdSideDump		= 111;
// const int DetIdNShield		= 112;
// const int DetIdSideCon		= 113;
// const int DetIdDoorCon		= 114;
// const int DetIdDumpCon		= 115;
const int DetIdDoraemon		= 220;
// const int DetIdUSWK		= 121;
// const int DetIdCDSBobbin	= 125;
// const int DetIdCDSCoil		= 126;
const int DetIdCDCCFRP		= 230;
const int DetIdCDCMylar		= 231;
// const int DetIdCDCEndCap	= 132;
// Target etc.
const int DetIdTarSys		= 240;
const int DetIdRadS		= 241;
const int DetIdTarCFRP		= 242;
const int DetIdTarCap		= 243;
const int DetIdTarRing		= 244;
// const int DetIdTarChm		= 145;
// const int DetIdBeamWindow	= 146;
// const int DetIdScatWindow	= 147;
// const int DetIdMagShield	= 148;

const int DetIdTarCell		= 250;
const int DetIdTarget		= 251;
const int DetIdCellTube		= 252;
const int DetIdCellFlange	= 253;
// const int DetIdBShield		= 154;
// const int DetIdCellWindow	= 154;
// const int DetIdBFrange		= 155;
// const int DetIdCellRing		= 156;
const int DetIdFiducial		= 260;

const int NumOfPlaneVmeRm=3;

// enum eTriggerFlag
//   {
//     kSpillStart =  1,
//     kSpillEnd   =  2,
//     kBeam1      =  3,
//     kBeam2      =  4,
//     kBeam_f     =  5,
//     kPion       =  6,
//     kPion_f     =  7,
//     kKaon1      =  8,
//     kKaon2      =  9,
//     kKaon3      =  10,
//     kKaon1_f       = 11,
//     kKaon2_f       = 12,
//     kKaonStart     = 13,
//     kKaonStartStop = 14,
//     kStartStop     = 15,
//     kMisc          = 16
//   };

#endif
