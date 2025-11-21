//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Nov  6 16:58:31 2025 by ROOT version 6.32.04
// from TTree tpc/tree of GenfitCarbon
// found on file: ../../k18-analyzer/rootfiles/Prod/KFCor_th50/run05667_GenfitCarbonVer10.root
//////////////////////////////////////////////////////////

#ifndef tpc_h
#define tpc_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class tpc {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           status;
   Int_t           runnum;
   Int_t           evnum;
   vector<int>     *trigpat;
   vector<int>     *trigflag;
   Int_t           nhHtof;
   vector<double>  *HtofSeg;
   vector<double>  *tHtof;
   vector<double>  *dtHtof;
   vector<double>  *deHtof;
   vector<double>  *posHtof;
   Int_t           nclTpc;
   Int_t           remain_nclTpc;
   vector<double>  *remain_cluster_x;
   vector<double>  *remain_cluster_y;
   vector<double>  *remain_cluster_z;
   vector<double>  *remain_cluster_de;
   vector<int>     *remain_cluster_size;
   vector<int>     *remain_cluster_layer;
   vector<int>     *remain_cluster_row_center;
   vector<double>  *remain_cluster_mrow;
   vector<double>  *remain_cluster_de_center;
   vector<double>  *remain_cluster_x_center;
   vector<double>  *remain_cluster_y_center;
   vector<double>  *remain_cluster_z_center;
   vector<int>     *remain_cluster_houghflag;
   Int_t           ntTpc;
   vector<int>     *isInTarget;
   vector<int>     *nhtrack;
   vector<int>     *trackid;
   vector<int>     *isXi;
   vector<int>     *isBeam;
   vector<int>     *isKurama;
   vector<int>     *isK18;
   vector<int>     *isAccidental;
   vector<int>     *isMultiloop;
   vector<int>     *charge;
   vector<int>     *pid;
   vector<double>  *chisqr;
   vector<double>  *pval;
   vector<double>  *helix_cx;
   vector<double>  *helix_cy;
   vector<double>  *helix_z0;
   vector<double>  *helix_r;
   vector<double>  *helix_dz;
   vector<double>  *dE;
   vector<double>  *dEdx;
   vector<double>  *mom0;
   vector<double>  *path;
   vector<int>     *isElectron;
   vector<double>  *nsigma_triton;
   vector<double>  *nsigma_deutron;
   vector<double>  *nsigma_proton;
   vector<double>  *nsigma_kaon;
   vector<double>  *nsigma_pion;
   vector<double>  *nsigma_electron;
   vector<vector<double> > *hitlayer;
   vector<vector<double> > *hitpos_x;
   vector<vector<double> > *hitpos_y;
   vector<vector<double> > *hitpos_z;
   vector<vector<double> > *calpos_x;
   vector<vector<double> > *calpos_y;
   vector<vector<double> > *calpos_z;
   vector<vector<double> > *residual;
   vector<vector<double> > *residual_x;
   vector<vector<double> > *residual_y;
   vector<vector<double> > *residual_z;
   vector<vector<double> > *resolution_x;
   vector<vector<double> > *resolution_y;
   vector<vector<double> > *resolution_z;
   vector<vector<double> > *pathhit;
   vector<vector<double> > *alpha;
   vector<vector<double> > *track_cluster_de;
   vector<vector<double> > *track_cluster_size;
   vector<vector<double> > *track_cluster_mrow;
   vector<vector<double> > *track_cluster_de_center;
   vector<vector<double> > *track_cluster_x_center;
   vector<vector<double> > *track_cluster_y_center;
   vector<vector<double> > *track_cluster_z_center;
   vector<vector<double> > *track_cluster_row_center;
   vector<vector<double> > *helix_t;
   vector<int>     *chargeIndistinguishable;
   vector<double>  *chisqr_inverted;
   vector<double>  *pval_inverted;
   vector<double>  *helix_cx_inverted;
   vector<double>  *helix_cy_inverted;
   vector<double>  *helix_z0_inverted;
   vector<double>  *helix_r_inverted;
   vector<double>  *helix_dz_inverted;
   vector<double>  *mom0_inverted;
   vector<int>     *pid_inverted;
   Int_t           nvtxTpc;
   vector<double>  *vtx_x;
   vector<double>  *vtx_y;
   vector<double>  *vtx_z;
   vector<double>  *vtx_dist;
   vector<double>  *vtx_angle;
   vector<vector<double> > *vtxid;
   vector<vector<double> > *vtxmom_theta;
   vector<vector<double> > *vtxpos_x;
   vector<vector<double> > *vtxpos_y;
   vector<vector<double> > *vtxpos_z;
   vector<vector<double> > *vtxmom_x;
   vector<vector<double> > *vtxmom_y;
   vector<vector<double> > *vtxmom_z;
   vector<double>  *GFcharge;
   vector<double>  *GFchisqr;
   vector<double>  *GFtof;
   vector<double>  *GFpval;
   vector<int>     *GFfitstatus;
   vector<int>     *GFpdgcode;
   vector<int>     *GFnhtrack;
   vector<vector<double> > *GFlayer;
   vector<vector<double> > *GFpos_x;
   vector<vector<double> > *GFpos_y;
   vector<vector<double> > *GFpos_z;
   vector<vector<double> > *GFmom;
   vector<vector<double> > *GFmom_x;
   vector<vector<double> > *GFmom_y;
   vector<vector<double> > *GFmom_z;
   vector<vector<double> > *GFresidual_x;
   vector<vector<double> > *GFresidual_y;
   vector<vector<double> > *GFresidual_z;
   vector<vector<double> > *GFresidual_p;
   vector<vector<double> > *GFresidual_px;
   vector<vector<double> > *GFresidual_py;
   vector<vector<double> > *GFresidual_pz;
   Int_t           GFntTpc_target;
   Double_t        GFprodvtx_x;
   Double_t        GFprodvtx_y;
   Double_t        GFprodvtx_z;
   vector<int>     *GFinside;
   vector<double>  *GFtracklen;
   vector<double>  *GFtrack2vtxdist;
   vector<double>  *GFcalctof;
   vector<double>  *GFsegHtof;
   vector<double>  *GFtofHtof;
   vector<double>  *GFtdiffHtof;
   vector<double>  *GFposHtof;
   vector<double>  *GFposx;
   vector<double>  *GFposy;
   vector<double>  *GFposz;
   vector<double>  *GFinvbeta;
   vector<double>  *GFm2;
   vector<double>  *nsigma_tritonHtof;
   vector<double>  *nsigma_deutronHtof;
   vector<double>  *nsigma_protonHtof;
   vector<double>  *nsigma_kaonHtof;
   vector<double>  *nsigma_pionHtof;
   vector<double>  *nsigma_electronHtof;
   vector<int>     *isLambda;
   vector<int>     *ncombiLambda;
   vector<double>  *distLambda;
   vector<double>  *angleLambda;
   vector<double>  *bestmassLambda;
   vector<vector<double> > *massLambda;
   vector<vector<double> > *vtxLambda_x;
   vector<vector<double> > *vtxLambda_y;
   vector<vector<double> > *vtxLambda_z;
   vector<vector<double> > *momLambda;
   vector<vector<double> > *momLambda_x;
   vector<vector<double> > *momLambda_y;
   vector<vector<double> > *momLambda_z;
   vector<vector<double> > *decaysidLambda;
   vector<vector<double> > *decaysmomLambda;
   vector<vector<double> > *decaysmomLambda_x;
   vector<vector<double> > *decaysmomLambda_y;
   vector<vector<double> > *decaysmomLambda_z;
   Int_t           nvtxTpcClustered;
   vector<double>  *clusteredVtx_x;
   vector<double>  *clusteredVtx_y;
   vector<double>  *clusteredVtx_z;
   vector<vector<double> > *clusteredVtxid;
   Int_t           ncombiReconFailedLambda;
   vector<int>     *ReconFailedLambdaPId;
   vector<int>     *ReconFailedLambdaPiId;
   vector<double>  *ReconFailedLambdaMass;
   vector<double>  *ReconFailedLambdaDecayVtx_x;
   vector<double>  *ReconFailedLambdaDecayVtx_y;
   vector<double>  *ReconFailedLambdaDecayVtx_z;
   vector<double>  *ReconFailedLambdaMom;
   vector<double>  *ReconFailedLambdaMom_x;
   vector<double>  *ReconFailedLambdaMom_y;
   vector<double>  *ReconFailedLambdaMom_z;
   vector<double>  *ReconFailedLambdaPMom;
   vector<double>  *ReconFailedLambdaPMom_x;
   vector<double>  *ReconFailedLambdaPMom_y;
   vector<double>  *ReconFailedLambdaPMom_z;
   vector<double>  *ReconFailedLambdaPiMom;
   vector<double>  *ReconFailedLambdaPiMom_x;
   vector<double>  *ReconFailedLambdaPiMom_y;
   vector<double>  *ReconFailedLambdaPiMom_z;
   vector<double>  *ReconFailedLambdavtxCloseDist;
   Int_t           ncombiPiPair;
   vector<int>     *PiPairPipId;
   vector<int>     *PiPairPimId;
   vector<double>  *PiPairPipMom;
   vector<double>  *PiPairPipMom_x;
   vector<double>  *PiPairPipMom_y;
   vector<double>  *PiPairPipMom_z;
   vector<double>  *PiPairPimMom;
   vector<double>  *PiPairPimMom_x;
   vector<double>  *PiPairPimMom_y;
   vector<double>  *PiPairPimMom_z;
   vector<double>  *PiPairMom;
   vector<double>  *PiPairMom_x;
   vector<double>  *PiPairMom_y;
   vector<double>  *PiPairMom_z;
   vector<double>  *PiPairReconLambdaMass;
   vector<double>  *PiPairReconMass;
   vector<double>  *PiPairCloseDist;
   Int_t           ntK18;
   vector<double>  *pK18;
   vector<double>  *chisqrK18;
   vector<double>  *xtgtK18;
   vector<double>  *ytgtK18;
   vector<double>  *utgtK18;
   vector<double>  *vtgtK18;
   vector<int>     *isgoodTPCK18;
   vector<double>  *chisqrTPCK18;
   vector<double>  *qTPCK18;
   vector<double>  *pTPCK18;
   vector<double>  *xtgtTPCK18;
   vector<double>  *ytgtTPCK18;
   vector<double>  *utgtTPCK18;
   vector<double>  *vtgtTPCK18;
   vector<double>  *thetaTPCK18;
   vector<double>  *lhtofTPCK18;
   vector<double>  *xhtofTPCK18;
   vector<double>  *yhtofTPCK18;
   vector<vector<double> > *lvpTPCK18;
   vector<vector<double> > *xvpTPCK18;
   vector<vector<double> > *yvpTPCK18;
   Int_t           ntKurama;
   vector<double>  *pKurama;
   vector<double>  *qKurama;
   vector<double>  *chisqrKurama;
   vector<double>  *m2Kurama;
   vector<double>  *xtgtKurama;
   vector<double>  *ytgtKurama;
   vector<double>  *utgtKurama;
   vector<double>  *vtgtKurama;
   vector<double>  *thetaKurama;
   vector<double>  *pathwcKurama;
   vector<double>  *xin;
   vector<double>  *yin;
   vector<double>  *zin;
   vector<double>  *pxin;
   vector<double>  *pyin;
   vector<double>  *pzin;
   vector<int>     *isgoodTPCKurama;
   vector<int>     *kflagTPCKurama;
   vector<double>  *chisqrTPCKurama;
   vector<double>  *pTPCKurama;
   vector<double>  *qTPCKurama;
   vector<double>  *m2TPCKurama;
   vector<double>  *xtgtTPCKurama;
   vector<double>  *ytgtTPCKurama;
   vector<double>  *utgtTPCKurama;
   vector<double>  *vtgtTPCKurama;
   vector<double>  *thetaTPCKurama;
   vector<double>  *pathTPCKurama;
   vector<double>  *lhtofTPCKurama;
   vector<double>  *xhtofTPCKurama;
   vector<double>  *yhtofTPCKurama;
   vector<vector<double> > *lvpTPCKurama;
   vector<vector<double> > *xvpTPCKurama;
   vector<vector<double> > *yvpTPCKurama;
   vector<int>     *isgoodTPC;
   vector<int>     *insideTPC;
   vector<double>  *vtxTPC;
   vector<double>  *vtyTPC;
   vector<double>  *vtzTPC;
   vector<double>  *closeDistTPC;
   vector<double>  *MissMassTPC;
   vector<double>  *MissMassCorrTPC;
   vector<double>  *MissMassCorrDETPC;
   vector<double>  *pOrgTPC;
   vector<double>  *pCalcTPC;
   vector<double>  *pCorrTPC;
   vector<double>  *pCorrDETPC;
   vector<double>  *xbTPC;
   vector<double>  *ybTPC;
   vector<double>  *ubTPC;
   vector<double>  *vbTPC;
   vector<double>  *xsTPC;
   vector<double>  *ysTPC;
   vector<double>  *usTPC;
   vector<double>  *vsTPC;
   Int_t           nKK;
   vector<int>     *Kflag;
   vector<double>  *MissMass;
   vector<double>  *MissMassCorr;
   vector<double>  *MissMassCorrDE;
   vector<double>  *vtx;
   vector<double>  *vty;
   vector<double>  *vtz;
   vector<double>  *pOrg;
   vector<double>  *pCalc;
   vector<double>  *pCorr;
   vector<double>  *pCorrDE;
   vector<double>  *xb;
   vector<double>  *yb;
   vector<double>  *ub;
   vector<double>  *vb;
   vector<double>  *xs;
   vector<double>  *ys;
   vector<double>  *us;
   vector<double>  *vs;
   vector<double>  *KmMom_x;
   vector<double>  *KmMom_y;
   vector<double>  *KmMom_z;
   vector<double>  *KpMom_x;
   vector<double>  *KpMom_y;
   vector<double>  *KpMom_z;
   vector<double>  *BE;
   vector<double>  *BETPC;
   vector<double>  *BE_LL;
   vector<double>  *BETPC_LL;
   Double_t        GFKKXiProductionVtx_x;
   Double_t        GFKKXiProductionVtx_y;
   Double_t        GFKKXiProductionVtx_z;
   Double_t        GFKKLLProductionVtx_x;
   Double_t        GFKKLLProductionVtx_y;
   Double_t        GFKKLLProductionVtx_z;
   Double_t        GFKKLProductionVtx_x1;
   Double_t        GFKKLProductionVtx_y1;
   Double_t        GFKKLProductionVtx_z1;
   Double_t        GFKKLProductionVtx_x2;
   Double_t        GFKKLProductionVtx_y2;
   Double_t        GFKKLProductionVtx_z2;
   Double_t        GFKKLProductionVtx_x;
   Double_t        GFKKLProductionVtx_y;
   Double_t        GFKKLProductionVtx_z;
   Bool_t          Xiflag;
   Bool_t          XiPflag;
   Double_t        XiMass;
   Double_t        XiDecayVtx_x;
   Double_t        XiDecayVtx_y;
   Double_t        XiDecayVtx_z;
   Double_t        XiMom;
   Double_t        XiMom_x;
   Double_t        XiMom_y;
   Double_t        XiMom_z;
   Double_t        XiVtxCloseDist;
   Double_t        LambdaMass;
   Double_t        LambdaDecayVtx_x;
   Double_t        LambdaDecayVtx_y;
   Double_t        LambdaDecayVtx_z;
   Double_t        LambdaMom;
   Double_t        LambdaMom_x;
   Double_t        LambdaMom_y;
   Double_t        LambdaMom_z;
   Double_t        LambdaVtxCloseDist;
   Double_t        XiTarget_x;
   Double_t        XiTarget_y;
   Double_t        XiTarget_z;
   Double_t        XiTargetMom;
   Double_t        XiTargetMom_x;
   Double_t        XiTargetMom_y;
   Double_t        XiTargetMom_z;
   Double_t        XiTargetCloseDist;
   Double_t        GFXiMass;
   Double_t        GFXiDecayVtx_x;
   Double_t        GFXiDecayVtx_y;
   Double_t        GFXiDecayVtx_z;
   Double_t        GFXiMom;
   Double_t        GFXiMom_x;
   Double_t        GFXiMom_y;
   Double_t        GFXiMom_z;
   Double_t        GFXiVtxCloseDist;
   Double_t        GFXiKKVtx_x;
   Double_t        GFXiKKVtx_y;
   Double_t        GFXiKKVtx_z;
   Double_t        GFXiKKVtxMom;
   Double_t        GFXiKKVtxMom_x;
   Double_t        GFXiKKVtxMom_y;
   Double_t        GFXiKKVtxMom_z;
   Double_t        GFXiKKVtxCloseDist;
   Double_t        GFXiProductionVtx_x;
   Double_t        GFXiProductionVtx_y;
   Double_t        GFXiProductionVtx_z;
   Double_t        GFXiProductionVtxMom;
   Double_t        GFXiProductionVtxMom_x;
   Double_t        GFXiProductionVtxMom_y;
   Double_t        GFXiProductionVtxMom_z;
   Double_t        GFXiProductionVtxCloseDist;
   Double_t        GFXiTrackLen;
   Double_t        GFXiTof;
   Double_t        GFXiMomLoss;
   Double_t        GFXiExcitation;
   Double_t        GFXiTarget_x;
   Double_t        GFXiTarget_y;
   Double_t        GFXiTarget_z;
   Double_t        GFXiTargetMom;
   Double_t        GFXiTargetMom_x;
   Double_t        GFXiTargetMom_y;
   Double_t        GFXiTargetMom_z;
   Double_t        GFXiTargetCloseDist;
   Double_t        GFXiTargetCenter_x;
   Double_t        GFXiTargetCenter_y;
   Double_t        GFXiTargetCenter_z;
   Double_t        GFXiTargetCenterMom;
   Double_t        GFXiTargetCenterMom_x;
   Double_t        GFXiTargetCenterMom_y;
   Double_t        GFXiTargetCenterMom_z;
   Double_t        GFXiTargetCenterCloseDist;
   Bool_t          LLflag;
   Double_t        LambdaLambdaVtx_x;
   Double_t        LambdaLambdaVtx_y;
   Double_t        LambdaLambdaVtx_z;
   Double_t        LambdaLambdaCloseDist;
   Double_t        LambdaMass1;
   Double_t        LambdaDecayVtx_x1;
   Double_t        LambdaDecayVtx_y1;
   Double_t        LambdaDecayVtx_z1;
   Double_t        LambdaMom1;
   Double_t        LambdaMom_x1;
   Double_t        LambdaMom_y1;
   Double_t        LambdaMom_z1;
   Double_t        LambdaVtxCloseDist1;
   Double_t        LambdaTargetCloseDist1;
   Double_t        LambdaTargetCloseVtx_x1;
   Double_t        LambdaTargetCloseVtx_y1;
   Double_t        LambdaTargetCloseVtx_z1;
   Double_t        LambdaMass2;
   Double_t        LambdaDecayVtx_x2;
   Double_t        LambdaDecayVtx_y2;
   Double_t        LambdaDecayVtx_z2;
   Double_t        LambdaMom2;
   Double_t        LambdaMom_x2;
   Double_t        LambdaMom_y2;
   Double_t        LambdaMom_z2;
   Double_t        LambdaVtxCloseDist2;
   Double_t        LambdaTargetCloseDist2;
   Double_t        LambdaTargetCloseVtx_x2;
   Double_t        LambdaTargetCloseVtx_y2;
   Double_t        LambdaTargetCloseVtx_z2;
   Double_t        GFLLExcitation;
   Double_t        GFLambdaMass1;
   Double_t        GFLambdaDecayVtx_x1;
   Double_t        GFLambdaDecayVtx_y1;
   Double_t        GFLambdaDecayVtx_z1;
   Double_t        GFLambdaMom1;
   Double_t        GFLambdaMom_x1;
   Double_t        GFLambdaMom_y1;
   Double_t        GFLambdaMom_z1;
   Double_t        GFLambdaVtxCloseDist1;
   Double_t        GFLambdaMass2;
   Double_t        GFLambdaDecayVtx_x2;
   Double_t        GFLambdaDecayVtx_y2;
   Double_t        GFLambdaDecayVtx_z2;
   Double_t        GFLambdaMom2;
   Double_t        GFLambdaMom_x2;
   Double_t        GFLambdaMom_y2;
   Double_t        GFLambdaMom_z2;
   Double_t        GFLambdaVtxCloseDist2;
   Double_t        GFLambdaMass1_Alt;
   Double_t        GFLambdaDecayVtx_x1_Alt;
   Double_t        GFLambdaDecayVtx_y1_Alt;
   Double_t        GFLambdaDecayVtx_z1_Alt;
   Double_t        GFLambdaMom1_Alt;
   Double_t        GFLambdaMom_x1_Alt;
   Double_t        GFLambdaMom_y1_Alt;
   Double_t        GFLambdaMom_z1_Alt;
   Double_t        GFLambdaVtxCloseDist1_Alt;
   Double_t        GFLambdaTargetCloseDist1_Alt;
   Double_t        GFLambdaTarget_x1_Alt;
   Double_t        GFLambdaTarget_y1_Alt;
   Double_t        GFLambdaTarget_z1_Alt;
   Double_t        GFLambdaTargetCenterCloseDist1_Alt;
   Double_t        GFLambdaTargetCenter_x1_Alt;
   Double_t        GFLambdaTargetCenter_y1_Alt;
   Double_t        GFLambdaTargetCenter_z1_Alt;
   Double_t        GFLambdaMass2_Alt;
   Double_t        GFLambdaDecayVtx_x2_Alt;
   Double_t        GFLambdaDecayVtx_y2_Alt;
   Double_t        GFLambdaDecayVtx_z2_Alt;
   Double_t        GFLambdaMom2_Alt;
   Double_t        GFLambdaMom_x2_Alt;
   Double_t        GFLambdaMom_y2_Alt;
   Double_t        GFLambdaMom_z2_Alt;
   Double_t        GFLambdaVtxCloseDist2_Alt;
   Double_t        GFLambdaTargetCloseDist2_Alt;
   Double_t        GFLambdaTarget_x2_Alt;
   Double_t        GFLambdaTarget_y2_Alt;
   Double_t        GFLambdaTarget_z2_Alt;
   Double_t        GFLambdaTargetCenterCloseDist2_Alt;
   Double_t        GFLambdaTargetCenter_x2_Alt;
   Double_t        GFLambdaTargetCenter_y2_Alt;
   Double_t        GFLambdaTargetCenter_z2_Alt;
   Double_t        GFLambdaLambdaVtx_x;
   Double_t        GFLambdaLambdaVtx_y;
   Double_t        GFLambdaLambdaVtx_z;
   Double_t        GFLambdaLambdaCloseDist;
   Double_t        GFLambdaTarget_x1;
   Double_t        GFLambdaTarget_y1;
   Double_t        GFLambdaTarget_z1;
   Double_t        GFLambdaTargetCloseDist1;
   Double_t        GFLambdaTarget_x2;
   Double_t        GFLambdaTarget_y2;
   Double_t        GFLambdaTarget_z2;
   Double_t        GFLambdaTargetCloseDist2;
   Double_t        GFLambdaTargetCenter_x1;
   Double_t        GFLambdaTargetCenter_y1;
   Double_t        GFLambdaTargetCenter_z1;
   Double_t        GFLambdaTargetCenterCloseDist1;
   Double_t        GFLambdaTargetCenter_x2;
   Double_t        GFLambdaTargetCenter_y2;
   Double_t        GFLambdaTargetCenter_z2;
   Double_t        GFLambdaTargetCenterCloseDist2;
   Double_t        GFLambdaProductionVtx_x1;
   Double_t        GFLambdaProductionVtx_y1;
   Double_t        GFLambdaProductionVtx_z1;
   Double_t        GFLambdaProductionVtxCloseDist1;
   Double_t        GFLambdaTrackLen1;
   Double_t        GFLambdaTof1;
   Double_t        GFLambdaProductionVtx_x2;
   Double_t        GFLambdaProductionVtx_y2;
   Double_t        GFLambdaProductionVtx_z2;
   Double_t        GFLambdaProductionVtxCloseDist2;
   Double_t        GFLambdaTrackLen2;
   Double_t        GFLambdaTof2;
   Double_t        KFLLExcitation;
   Double_t        KFKKLLProductionVtxChisqr;
   Double_t        KFKKLLProductionVtx_x;
   Double_t        KFKKLLProductionVtx_y;
   Double_t        KFKKLLProductionVtx_z;
   Double_t        KFKKLProductionVtx_x1;
   Double_t        KFKKLProductionVtx_y1;
   Double_t        KFKKLProductionVtx_z1;
   Double_t        KFKKLProductionVtx_x2;
   Double_t        KFKKLProductionVtx_y2;
   Double_t        KFKKLProductionVtx_z2;
   Double_t        KFLambdaLambdaVtx_x;
   Double_t        KFLambdaLambdaVtx_y;
   Double_t        KFLambdaLambdaVtx_z;
   Double_t        KFLambdaLambdaCloseDist;
   Double_t        KFLambdaProductionVtx_x1;
   Double_t        KFLambdaProductionVtx_y1;
   Double_t        KFLambdaProductionVtx_z1;
   Double_t        KFLambdaProductionVtxCloseDist1;
   Double_t        KFLambdaTrackLen1;
   Double_t        KFLambdaTof1;
   Double_t        KFLambdaProductionVtx_x2;
   Double_t        KFLambdaProductionVtx_y2;
   Double_t        KFLambdaProductionVtx_z2;
   Double_t        KFLambdaProductionVtxCloseDist2;
   Double_t        KFLambdaTrackLen2;
   Double_t        KFLambdaTof2;
   Double_t        KFLambdaPval1;
   vector<double>  *KFLambdaPull1;
   Double_t        KFLambdaChisqr1;
   Double_t        KFLambdaMom1;
   Double_t        KFLambdaMom_x1;
   Double_t        KFLambdaMom_y1;
   Double_t        KFLambdaMom_z1;
   Double_t        KFLambdaPval2;
   vector<double>  *KFLambdaPull2;
   Double_t        KFLambdaChisqr2;
   Double_t        KFLambdaMom2;
   Double_t        KFLambdaMom_x2;
   Double_t        KFLambdaMom_y2;
   Double_t        KFLambdaMom_z2;
   Bool_t          Lflag;
   Double_t        LambdaTargetCloseDist;
   Double_t        LambdaTargetCloseVtx_x;
   Double_t        LambdaTargetCloseVtx_y;
   Double_t        LambdaTargetCloseVtx_z;
   Double_t        GFLambdaMass;
   Double_t        GFLambdaDecayVtx_x;
   Double_t        GFLambdaDecayVtx_y;
   Double_t        GFLambdaDecayVtx_z;
   Double_t        GFLambdaMom;
   Double_t        GFLambdaMom_x;
   Double_t        GFLambdaMom_y;
   Double_t        GFLambdaMom_z;
   Double_t        GFLambdaVtxCloseDist;
   Double_t        GFLambdaTargetCloseDist;
   Double_t        GFLambdaTarget_x;
   Double_t        GFLambdaTarget_y;
   Double_t        GFLambdaTarget_z;
   Double_t        GFLambdaTargetCenter_x;
   Double_t        GFLambdaTargetCenter_y;
   Double_t        GFLambdaTargetCenter_z;
   Double_t        GFLambdaTargetCenterCloseDist;
   Double_t        GFLambdaProductionVtx_x;
   Double_t        GFLambdaProductionVtx_y;
   Double_t        GFLambdaProductionVtx_z;
   Double_t        GFLambdaProductionVtxCloseDist;
   Double_t        GFLambdaTrackLen;
   Double_t        GFLambdaTof;
   Double_t        KFKKLProductionVtx_x;
   Double_t        KFKKLProductionVtx_y;
   Double_t        KFKKLProductionVtx_z;
   Double_t        KFLambdaProductionVtx_x;
   Double_t        KFLambdaProductionVtx_y;
   Double_t        KFLambdaProductionVtx_z;
   Double_t        KFLambdaProductionVtxCloseDist;
   Double_t        KFLambdaTrackLen;
   Double_t        KFLambdaTof;
   Bool_t          LPhiflag;
   Double_t        PhiMass;
   Double_t        PhiDecayVtx_x;
   Double_t        PhiDecayVtx_y;
   Double_t        PhiDecayVtx_z;
   Double_t        PhiCosKK;
   Double_t        PhiMom;
   Double_t        PhiMom_x;
   Double_t        PhiMom_y;
   Double_t        PhiMom_z;
   Double_t        PhiVtxCloseDist;
   vector<int>     *PhiDecaysTrackId;
   vector<double>  *PhiDecaysMom;
   vector<double>  *PhiDecaysMom_x;
   vector<double>  *PhiDecaysMom_y;
   vector<double>  *PhiDecaysMom_z;
   Double_t        GFPhiMass;
   Double_t        GFPhiDecayVtx_x;
   Double_t        GFPhiDecayVtx_y;
   Double_t        GFPhiDecayVtx_z;
   Double_t        GFPhiCosKK;
   Double_t        GFPhiMom;
   Double_t        GFPhiMom_x;
   Double_t        GFPhiMom_y;
   Double_t        GFPhiMom_z;
   Double_t        GFPhiVtxCloseDist;
   Double_t        GFPhiProductionVtxCloseDist;
   Double_t        GFPhiKmMassSquare;
   Double_t        GFPhiKmInvbeta;
   Double_t        GFPhiKpMassSquare;
   Double_t        GFPhiKpInvbeta;
   vector<double>  *GFPhiDecaysMom;
   vector<double>  *GFPhiDecaysMom_x;
   vector<double>  *GFPhiDecaysMom_y;
   vector<double>  *GFPhiDecaysMom_z;
   Double_t        GFPhiMassKurama;
   Double_t        GFPhiMomKurama;
   Double_t        GFPhiMomKurama_x;
   Double_t        GFPhiMomKurama_y;
   Double_t        GFPhiMomKurama_z;
   Bool_t          LPiflag;
   Bool_t          PiPiflag;
   Bool_t          Pimflag;
   Bool_t          Emptyflag;
   vector<int>     *GFXiDecaysNhit;
   vector<double>  *GFXiDecaysChisqr;
   vector<double>  *GFXiDecaysCharge;
   vector<double>  *GFXiDecaysTof;
   vector<double>  *GFXiDecaysPval;
   vector<int>     *GFXiDecaysPdgcode;
   vector<int>     *GFXiDecaysHtofId;
   vector<double>  *GFXiDecaysTrackLen;
   vector<double>  *GFXiDecaysTrackTof;
   vector<double>  *GFXiDecaysMassSquare;
   vector<double>  *GFXiDecaysInvbeta;
   vector<double>  *GFXiDecaysMom;
   vector<double>  *GFXiDecaysMom_x;
   vector<double>  *GFXiDecaysMom_y;
   vector<double>  *GFXiDecaysMom_z;
   vector<double>  *GFXiDecaysMomCM;
   vector<double>  *GFXiDecaysMomCM_x;
   vector<double>  *GFXiDecaysMomCM_y;
   vector<double>  *GFXiDecaysMomCM_z;
   vector<double>  *GFXiDecaysMomLoss;
   vector<double>  *GFXiDecaysELoss;
   vector<int>     *XiDecaysTrackId;
   vector<double>  *XiDecaysMom;
   vector<double>  *XiDecaysMom_x;
   vector<double>  *XiDecaysMom_y;
   vector<double>  *XiDecaysMom_z;
   vector<double>  *XiDecaysMomCM;
   vector<double>  *XiDecaysMomCM_x;
   vector<double>  *XiDecaysMomCM_y;
   vector<double>  *XiDecaysMomCM_z;
   vector<int>     *GFLLDecaysNhit;
   vector<double>  *GFLLDecaysChisqr;
   vector<double>  *GFLLDecaysCharge;
   vector<double>  *GFLLDecaysTof;
   vector<double>  *GFLLDecaysPval;
   vector<int>     *GFLLDecaysPdgcode;
   vector<int>     *GFLLDecaysHtofId;
   vector<double>  *GFLLDecaysTrackLen;
   vector<double>  *GFLLDecaysTrackTof;
   vector<double>  *GFLLDecaysMassSquare;
   vector<double>  *GFLLDecaysInvbeta;
   vector<double>  *GFLLDecaysMom;
   vector<double>  *GFLLDecaysMom_x;
   vector<double>  *GFLLDecaysMom_y;
   vector<double>  *GFLLDecaysMom_z;
   vector<double>  *GFLLDecaysMomCM;
   vector<double>  *GFLLDecaysMomCM_x;
   vector<double>  *GFLLDecaysMomCM_y;
   vector<double>  *GFLLDecaysMomCM_z;
   vector<double>  *GFLLDecaysMomLoss;
   vector<double>  *GFLLDecaysELoss;
   vector<int>     *LLDecaysTrackId;
   vector<double>  *LLDecaysMom;
   vector<double>  *LLDecaysMom_x;
   vector<double>  *LLDecaysMom_y;
   vector<double>  *LLDecaysMom_z;
   vector<double>  *LLDecaysMomCM;
   vector<double>  *LLDecaysMomCM_x;
   vector<double>  *LLDecaysMomCM_y;
   vector<double>  *LLDecaysMomCM_z;
   vector<int>     *GFDecaysNhit;
   vector<double>  *GFDecaysChisqr;
   vector<double>  *GFDecaysCharge;
   vector<double>  *GFDecaysTof;
   vector<double>  *GFDecaysPval;
   vector<int>     *GFDecaysPdgcode;
   vector<int>     *GFDecaysHtofId;
   vector<double>  *GFDecaysTrackLen;
   vector<double>  *GFDecaysTrackTof;
   vector<double>  *GFDecaysMassSquare;
   vector<double>  *GFDecaysInvbeta;
   vector<double>  *GFDecaysMom;
   vector<double>  *GFDecaysMom_x;
   vector<double>  *GFDecaysMom_y;
   vector<double>  *GFDecaysMom_z;
   vector<double>  *GFDecaysMomCM;
   vector<double>  *GFDecaysMomCM_x;
   vector<double>  *GFDecaysMomCM_y;
   vector<double>  *GFDecaysMomCM_z;
   vector<double>  *GFDecaysMomLoss;
   vector<double>  *GFDecaysELoss;
   vector<int>     *DecaysTrackId;
   vector<double>  *DecaysMom;
   vector<double>  *DecaysMom_x;
   vector<double>  *DecaysMom_y;
   vector<double>  *DecaysMom_z;
   vector<double>  *DecaysMomCM;
   vector<double>  *DecaysMomCM_x;
   vector<double>  *DecaysMomCM_y;
   vector<double>  *DecaysMomCM_z;
   Int_t           AccidentalMultiplicity;
   vector<int>     *AccidentalTrackId;
   Int_t           LLResidualsMultiplicity;
   Int_t           LLpMultiplicity;
   Int_t           LLpipMultiplicity;
   Int_t           LLpimMultiplicity;
   Int_t           LLepMultiplicity;
   Int_t           LLemMultiplicity;
   Int_t           LLppipMultiplicity;
   vector<int>     *LLResidualsTrackId;
   vector<double>  *LLResidualsMassSquare;
   vector<double>  *LLResidualsInvbeta;
   vector<double>  *LLResidualsCloseDistTgt;
   vector<double>  *LLResidualsGFProductionVtxCloseDist;
   vector<double>  *LLResidualsKFProductionVtxCloseDist;
   vector<double>  *LLResidualsMom;
   vector<double>  *LLResidualsMom_x;
   vector<double>  *LLResidualsMom_y;
   vector<double>  *LLResidualsMom_z;
   vector<double>  *LLResidualsCharge;
   Int_t           XiResidualsMultiplicity;
   Int_t           XipMultiplicity;
   Int_t           XipipMultiplicity;
   Int_t           XipimMultiplicity;
   Int_t           XiepMultiplicity;
   Int_t           XiemMultiplicity;
   Int_t           XippipMultiplicity;
   vector<int>     *XiResidualsTrackId;
   vector<double>  *XiResidualsMassSquare;
   vector<double>  *XiResidualsInvbeta;
   vector<double>  *XiResidualsCloseDistTgt;
   vector<double>  *XiResidualsGFProductionVtxCloseDist;
   vector<double>  *XiResidualsKFProductionVtxCloseDist;
   vector<double>  *XiResidualsMom;
   vector<double>  *XiResidualsMom_x;
   vector<double>  *XiResidualsMom_y;
   vector<double>  *XiResidualsMom_z;
   vector<double>  *XiResidualsCharge;
   Int_t           ResidualsMultiplicity;
   Int_t           pMultiplicity;
   Int_t           pipMultiplicity;
   Int_t           pimMultiplicity;
   Int_t           epMultiplicity;
   Int_t           emMultiplicity;
   Int_t           ppipMultiplicity;
   vector<int>     *ResidualsTrackId;
   vector<double>  *ResidualsMassSquare;
   vector<double>  *ResidualsInvbeta;
   vector<double>  *ResidualsCloseDistTgt;
   vector<double>  *ResidualsGFProductionVtxCloseDist;
   vector<double>  *ResidualsKFProductionVtxCloseDist;
   vector<double>  *ResidualsMom;
   vector<double>  *ResidualsMom_x;
   vector<double>  *ResidualsMom_y;
   vector<double>  *ResidualsMom_z;
   vector<double>  *ResidualsCharge;
   Double_t        KFLambdaMomPpi;
   Double_t        KFLambdaMomPpi_x;
   Double_t        KFLambdaMomPpi_y0;
   Double_t        KFLambdaMomPpi_z0;
   Double_t        KFLambdaMom;
   Double_t        KFLambdaMom_x;
   Double_t        KFLambdaMom_y;
   Double_t        KFLambdaMom_z;
   Double_t        KFLambdaChisqr;
   Double_t        KFLambdaPval;
   vector<double>  *KFLambdaPull;
   vector<vector<double> > *KFLambdaCovMatrix;
   Double_t        KFXiVtxCloseDist;
   Double_t        KFXiMom;
   Double_t        KFXiMom_x;
   Double_t        KFXiMom_y;
   Double_t        KFXiMom_z;
   Double_t        KFXiChisqr;
   Double_t        KFXiPval;
   Double_t        KFXiMass;
   Double_t        KFXiDecayVtx_x;
   Double_t        KFXiDecayVtx_y;
   Double_t        KFXiDecayVtx_z;
   vector<double>  *KFXiPull;
   vector<vector<double> > *KFXiCovMatrix;
   Double_t        KFKKXiProductionVtxChisqr;
   Double_t        KFKKXiProductionVtx_x;
   Double_t        KFKKXiProductionVtx_y;
   Double_t        KFKKXiProductionVtx_z;
   Double_t        KFXiProductionVtx_x;
   Double_t        KFXiProductionVtx_y;
   Double_t        KFXiProductionVtx_z;
   Double_t        KFXiProductionVtxMom;
   Double_t        KFXiProductionVtxMom_x;
   Double_t        KFXiProductionVtxMom_y;
   Double_t        KFXiProductionVtxMom_z;
   Double_t        KFXiProductionVtxCloseDist;
   Double_t        KFXiTrackLen;
   Double_t        KFXiTof;
   Double_t        KFXiMomLoss;
   Double_t        KFXiExcitation;
   Double_t        KFKpXiProductionVtx_x;
   Double_t        KFKpXiProductionVtx_y;
   Double_t        KFKpXiProductionVtx_z;
   Double_t        KFXiProductionVtx_x_KpXi;
   Double_t        KFXiProductionVtx_y_KpXi;
   Double_t        KFXiProductionVtx_z_KpXi;
   Double_t        KFXiProductionVtxMom_KpXi;
   Double_t        KFXiProductionVtxMom_x_KpXi;
   Double_t        KFXiProductionVtxMom_y_KpXi;
   Double_t        KFXiProductionVtxMom_z_KpXi;
   Double_t        KFXiProductionVtxCloseDist_KpXi;
   Double_t        KFXiProductionVtx_x_KK;
   Double_t        KFXiProductionVtx_y_KK;
   Double_t        KFXiProductionVtx_z_KK;
   Double_t        KFXiProductionVtxMom_KK;
   Double_t        KFXiProductionVtxMom_x_KK;
   Double_t        KFXiProductionVtxMom_y_KK;
   Double_t        KFXiProductionVtxMom_z_KK;
   Double_t        KFXiProductionVtxCloseDist_KK;
   vector<double>  *KFLLDecaysMom;
   vector<double>  *KFLLDecaysMom_x;
   vector<double>  *KFLLDecaysMom_y;
   vector<double>  *KFLLDecaysMom_z;
   vector<double>  *KFLLDecaysMomCM;
   vector<double>  *KFLLDecaysMomCM_x;
   vector<double>  *KFLLDecaysMomCM_y;
   vector<double>  *KFLLDecaysMomCM_z;
   vector<double>  *KFXiDecaysMom;
   vector<double>  *KFXiDecaysMom_x;
   vector<double>  *KFXiDecaysMom_y;
   vector<double>  *KFXiDecaysMom_z;
   vector<double>  *KFXiDecaysMomCM;
   vector<double>  *KFXiDecaysMomCM_x;
   vector<double>  *KFXiDecaysMomCM_y;
   vector<double>  *KFXiDecaysMomCM_z;
   vector<double>  *KFDecaysMom;
   vector<double>  *KFDecaysMom_x;
   vector<double>  *KFDecaysMom_y;
   vector<double>  *KFDecaysMom_z;
   vector<double>  *KFDecaysMomCM;
   vector<double>  *KFDecaysMomCM_x;
   vector<double>  *KFDecaysMomCM_y;
   vector<double>  *KFDecaysMomCM_z;
   Int_t           LLnGamma;
   vector<int>     *LLGammaEpTrackId;
   vector<int>     *LLGammaEmTrackId;
   vector<double>  *LLGammaMomId;
   vector<double>  *LLGammaEpMom_x;
   vector<double>  *LLGammaEpMom_y;
   vector<double>  *LLGammaEpMom_z;
   vector<double>  *LLGammaEmMom;
   vector<double>  *LLGammaEmMom_x;
   vector<double>  *LLGammaEmMom_y;
   vector<double>  *LLGammaEmMom_z;
   vector<double>  *LLGammaMom;
   vector<double>  *LLGammaMom_x;
   vector<double>  *LLGammaMom_y;
   vector<double>  *LLGammaMom_z;
   vector<double>  *LLGammaVtxCloseDist;
   vector<double>  *LLGammaDecayVtx_x;
   vector<double>  *LLGammaDecayVtx_y;
   vector<double>  *LLGammaDecayVtx_z;
   Int_t           XinGamma;
   vector<int>     *XiGammaEpTrackId;
   vector<int>     *XiGammaEmTrackId;
   vector<double>  *XiGammaMomId;
   vector<double>  *XiGammaEpMom_x;
   vector<double>  *XiGammaEpMom_y;
   vector<double>  *XiGammaEpMom_z;
   vector<double>  *XiGammaEmMom;
   vector<double>  *XiGammaEmMom_x;
   vector<double>  *XiGammaEmMom_y;
   vector<double>  *XiGammaEmMom_z;
   vector<double>  *XiGammaMom;
   vector<double>  *XiGammaMom_x;
   vector<double>  *XiGammaMom_y;
   vector<double>  *XiGammaMom_z;
   vector<double>  *XiGammaVtxCloseDist;
   vector<double>  *XiGammaDecayVtx_x;
   vector<double>  *XiGammaDecayVtx_y;
   vector<double>  *XiGammaDecayVtx_z;
   Int_t           nGamma;
   vector<int>     *GammaEpTrackId;
   vector<int>     *GammaEmTrackId;
   vector<double>  *GammaMomId;
   vector<double>  *GammaEpMom_x;
   vector<double>  *GammaEpMom_y;
   vector<double>  *GammaEpMom_z;
   vector<double>  *GammaEmMom;
   vector<double>  *GammaEmMom_x;
   vector<double>  *GammaEmMom_y;
   vector<double>  *GammaEmMom_z;
   vector<double>  *GammaMom;
   vector<double>  *GammaMom_x;
   vector<double>  *GammaMom_y;
   vector<double>  *GammaMom_z;
   vector<double>  *GammaVtxCloseDist;
   vector<double>  *GammaDecayVtx_x;
   vector<double>  *GammaDecayVtx_y;
   vector<double>  *GammaDecayVtx_z;

   // List of branches
   TBranch        *b_status;   //!
   TBranch        *b_runnum;   //!
   TBranch        *b_evnum;   //!
   TBranch        *b_trigpat;   //!
   TBranch        *b_trigflag;   //!
   TBranch        *b_nhHtof;   //!
   TBranch        *b_HtofSeg;   //!
   TBranch        *b_tHtof;   //!
   TBranch        *b_dtHtof;   //!
   TBranch        *b_deHtof;   //!
   TBranch        *b_posHtof;   //!
   TBranch        *b_nclTpc;   //!
   TBranch        *b_remain_nclTpc;   //!
   TBranch        *b_remain_cluster_x;   //!
   TBranch        *b_remain_cluster_y;   //!
   TBranch        *b_remain_cluster_z;   //!
   TBranch        *b_remain_cluster_de;   //!
   TBranch        *b_remain_cluster_size;   //!
   TBranch        *b_remain_cluster_layer;   //!
   TBranch        *b_remain_cluster_row_center;   //!
   TBranch        *b_remain_cluster_mrow;   //!
   TBranch        *b_remain_cluster_de_center;   //!
   TBranch        *b_remain_cluster_x_center;   //!
   TBranch        *b_remain_cluster_y_center;   //!
   TBranch        *b_remain_cluster_z_center;   //!
   TBranch        *b_remain_cluster_houghflag;   //!
   TBranch        *b_ntTpc;   //!
   TBranch        *b_isInTarget;   //!
   TBranch        *b_nhtrack;   //!
   TBranch        *b_trackid;   //!
   TBranch        *b_isXi;   //!
   TBranch        *b_isBeam;   //!
   TBranch        *b_isKurama;   //!
   TBranch        *b_isK18;   //!
   TBranch        *b_isAccidental;   //!
   TBranch        *b_isMultiloop;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_pid;   //!
   TBranch        *b_chisqr;   //!
   TBranch        *b_pval;   //!
   TBranch        *b_helix_cx;   //!
   TBranch        *b_helix_cy;   //!
   TBranch        *b_helix_z0;   //!
   TBranch        *b_helix_r;   //!
   TBranch        *b_helix_dz;   //!
   TBranch        *b_dE;   //!
   TBranch        *b_dEdx;   //!
   TBranch        *b_mom0;   //!
   TBranch        *b_path;   //!
   TBranch        *b_isElectron;   //!
   TBranch        *b_nsigma_triton;   //!
   TBranch        *b_nsigma_deutron;   //!
   TBranch        *b_nsigma_proton;   //!
   TBranch        *b_nsigma_kaon;   //!
   TBranch        *b_nsigma_pion;   //!
   TBranch        *b_nsigma_electron;   //!
   TBranch        *b_hitlayer;   //!
   TBranch        *b_hitpos_x;   //!
   TBranch        *b_hitpos_y;   //!
   TBranch        *b_hitpos_z;   //!
   TBranch        *b_calpos_x;   //!
   TBranch        *b_calpos_y;   //!
   TBranch        *b_calpos_z;   //!
   TBranch        *b_residual;   //!
   TBranch        *b_residual_x;   //!
   TBranch        *b_residual_y;   //!
   TBranch        *b_residual_z;   //!
   TBranch        *b_resolution_x;   //!
   TBranch        *b_resolution_y;   //!
   TBranch        *b_resolution_z;   //!
   TBranch        *b_pathhit;   //!
   TBranch        *b_alpha;   //!
   TBranch        *b_track_cluster_de;   //!
   TBranch        *b_track_cluster_size;   //!
   TBranch        *b_track_cluster_mrow;   //!
   TBranch        *b_track_cluster_de_center;   //!
   TBranch        *b_track_cluster_x_center;   //!
   TBranch        *b_track_cluster_y_center;   //!
   TBranch        *b_track_cluster_z_center;   //!
   TBranch        *b_track_cluster_row_center;   //!
   TBranch        *b_helix_t;   //!
   TBranch        *b_chargeIndistinguishable;   //!
   TBranch        *b_chisqr_inverted;   //!
   TBranch        *b_pval_inverted;   //!
   TBranch        *b_helix_cx_inverted;   //!
   TBranch        *b_helix_cy_inverted;   //!
   TBranch        *b_helix_z0_inverted;   //!
   TBranch        *b_helix_r_inverted;   //!
   TBranch        *b_helix_dz_inverted;   //!
   TBranch        *b_mom0_inverted;   //!
   TBranch        *b_pid_inverted;   //!
   TBranch        *b_nvtxTpc;   //!
   TBranch        *b_vtx_x;   //!
   TBranch        *b_vtx_y;   //!
   TBranch        *b_vtx_z;   //!
   TBranch        *b_vtx_dist;   //!
   TBranch        *b_vtx_angle;   //!
   TBranch        *b_vtxid;   //!
   TBranch        *b_vtxmom_theta;   //!
   TBranch        *b_vtxpos_x;   //!
   TBranch        *b_vtxpos_y;   //!
   TBranch        *b_vtxpos_z;   //!
   TBranch        *b_vtxmom_x;   //!
   TBranch        *b_vtxmom_y;   //!
   TBranch        *b_vtxmom_z;   //!
   TBranch        *b_GFcharge;   //!
   TBranch        *b_GFchisqr;   //!
   TBranch        *b_GFtof;   //!
   TBranch        *b_GFpval;   //!
   TBranch        *b_GFfitstatus;   //!
   TBranch        *b_GFpdgcode;   //!
   TBranch        *b_GFnhtrack;   //!
   TBranch        *b_GFlayer;   //!
   TBranch        *b_GFpos_x;   //!
   TBranch        *b_GFpos_y;   //!
   TBranch        *b_GFpos_z;   //!
   TBranch        *b_GFmom;   //!
   TBranch        *b_GFmom_x;   //!
   TBranch        *b_GFmom_y;   //!
   TBranch        *b_GFmom_z;   //!
   TBranch        *b_GFresidual_x;   //!
   TBranch        *b_GFresidual_y;   //!
   TBranch        *b_GFresidual_z;   //!
   TBranch        *b_GFresidual_p;   //!
   TBranch        *b_GFresidual_px;   //!
   TBranch        *b_GFresidual_py;   //!
   TBranch        *b_GFresidual_pz;   //!
   TBranch        *b_GFntTpc_target;   //!
   TBranch        *b_GFprodvtx_x;   //!
   TBranch        *b_GFprodvtx_y;   //!
   TBranch        *b_GFprodvtx_z;   //!
   TBranch        *b_GFinside;   //!
   TBranch        *b_GFtracklen;   //!
   TBranch        *b_GFtrack2vtxdist;   //!
   TBranch        *b_GFcalctof;   //!
   TBranch        *b_GFsegHtof;   //!
   TBranch        *b_GFtofHtof;   //!
   TBranch        *b_GFtdiffHtof;   //!
   TBranch        *b_GFposHtof;   //!
   TBranch        *b_GFposx;   //!
   TBranch        *b_GFposy;   //!
   TBranch        *b_GFposz;   //!
   TBranch        *b_GFinvbeta;   //!
   TBranch        *b_GFm2;   //!
   TBranch        *b_nsigma_tritonHtof;   //!
   TBranch        *b_nsigma_deutronHtof;   //!
   TBranch        *b_nsigma_protonHtof;   //!
   TBranch        *b_nsigma_kaonHtof;   //!
   TBranch        *b_nsigma_pionHtof;   //!
   TBranch        *b_nsigma_electronHtof;   //!
   TBranch        *b_isLambda;   //!
   TBranch        *b_ncombiLambda;   //!
   TBranch        *b_distLambda;   //!
   TBranch        *b_angleLambda;   //!
   TBranch        *b_bestmassLambda;   //!
   TBranch        *b_massLambda;   //!
   TBranch        *b_vtxLambda_x;   //!
   TBranch        *b_vtxLambda_y;   //!
   TBranch        *b_vtxLambda_z;   //!
   TBranch        *b_momLambda;   //!
   TBranch        *b_momLambda_x;   //!
   TBranch        *b_momLambda_y;   //!
   TBranch        *b_momLambda_z;   //!
   TBranch        *b_decaysidLambda;   //!
   TBranch        *b_decaysmomLambda;   //!
   TBranch        *b_decaysmomLambda_x;   //!
   TBranch        *b_decaysmomLambda_y;   //!
   TBranch        *b_decaysmomLambda_z;   //!
   TBranch        *b_nvtxTpcClustered;   //!
   TBranch        *b_clusteredVtx_x;   //!
   TBranch        *b_clusteredVtx_y;   //!
   TBranch        *b_clusteredVtx_z;   //!
   TBranch        *b_clusteredVtxid;   //!
   TBranch        *b_ncombiReconFailedLambda;   //!
   TBranch        *b_ReconFailedLambdaPId;   //!
   TBranch        *b_ReconFailedLambdaPiId;   //!
   TBranch        *b_ReconFailedLambdaMass;   //!
   TBranch        *b_ReconFailedLambdaDecayVtx_x;   //!
   TBranch        *b_ReconFailedLambdaDecayVtx_y;   //!
   TBranch        *b_ReconFailedLambdaDecayVtx_z;   //!
   TBranch        *b_ReconFailedLambdaMom;   //!
   TBranch        *b_ReconFailedLambdaMom_x;   //!
   TBranch        *b_ReconFailedLambdaMom_y;   //!
   TBranch        *b_ReconFailedLambdaMom_z;   //!
   TBranch        *b_ReconFailedLambdaPMom;   //!
   TBranch        *b_ReconFailedLambdaPMom_x;   //!
   TBranch        *b_ReconFailedLambdaPMom_y;   //!
   TBranch        *b_ReconFailedLambdaPMom_z;   //!
   TBranch        *b_ReconFailedLambdaPiMom;   //!
   TBranch        *b_ReconFailedLambdaPiMom_x;   //!
   TBranch        *b_ReconFailedLambdaPiMom_y;   //!
   TBranch        *b_ReconFailedLambdaPiMom_z;   //!
   TBranch        *b_ReconFailedLambdavtxCloseDist;   //!
   TBranch        *b_ncombiPiPair;   //!
   TBranch        *b_PiPairPipId;   //!
   TBranch        *b_PiPairPimId;   //!
   TBranch        *b_PiPairPipMom;   //!
   TBranch        *b_PiPairPipMom_x;   //!
   TBranch        *b_PiPairPipMom_y;   //!
   TBranch        *b_PiPairPipMom_z;   //!
   TBranch        *b_PiPairPimMom;   //!
   TBranch        *b_PiPairPimMom_x;   //!
   TBranch        *b_PiPairPimMom_y;   //!
   TBranch        *b_PiPairPimMom_z;   //!
   TBranch        *b_PiPairMom;   //!
   TBranch        *b_PiPairMom_x;   //!
   TBranch        *b_PiPairMom_y;   //!
   TBranch        *b_PiPairMom_z;   //!
   TBranch        *b_PiPairReconLambdaMass;   //!
   TBranch        *b_PiPairReconMass;   //!
   TBranch        *b_PiPairCloseDist;   //!
   TBranch        *b_ntK18;   //!
   TBranch        *b_pK18;   //!
   TBranch        *b_chisqrK18;   //!
   TBranch        *b_xtgtK18;   //!
   TBranch        *b_ytgtK18;   //!
   TBranch        *b_utgtK18;   //!
   TBranch        *b_vtgtK18;   //!
   TBranch        *b_isgoodTPCK18;   //!
   TBranch        *b_chisqrTPCK18;   //!
   TBranch        *b_qTPCK18;   //!
   TBranch        *b_pTPCK18;   //!
   TBranch        *b_xtgtTPCK18;   //!
   TBranch        *b_ytgtTPCK18;   //!
   TBranch        *b_utgtTPCK18;   //!
   TBranch        *b_vtgtTPCK18;   //!
   TBranch        *b_thetaTPCK18;   //!
   TBranch        *b_lhtofTPCK18;   //!
   TBranch        *b_xhtofTPCK18;   //!
   TBranch        *b_yhtofTPCK18;   //!
   TBranch        *b_lvpTPCK18;   //!
   TBranch        *b_xvpTPCK18;   //!
   TBranch        *b_yvpTPCK18;   //!
   TBranch        *b_ntKurama;   //!
   TBranch        *b_pKurama;   //!
   TBranch        *b_qKurama;   //!
   TBranch        *b_chisqrKurama;   //!
   TBranch        *b_m2Kurama;   //!
   TBranch        *b_xtgtKurama;   //!
   TBranch        *b_ytgtKurama;   //!
   TBranch        *b_utgtKurama;   //!
   TBranch        *b_vtgtKurama;   //!
   TBranch        *b_thetaKurama;   //!
   TBranch        *b_pathwcKurama;   //!
   TBranch        *b_xin;   //!
   TBranch        *b_yin;   //!
   TBranch        *b_zin;   //!
   TBranch        *b_pxin;   //!
   TBranch        *b_pyin;   //!
   TBranch        *b_pzin;   //!
   TBranch        *b_isgoodTPCKurama;   //!
   TBranch        *b_kflagTPCKurama;   //!
   TBranch        *b_chisqrTPCKurama;   //!
   TBranch        *b_pTPCKurama;   //!
   TBranch        *b_qTPCKurama;   //!
   TBranch        *b_m2TPCKurama;   //!
   TBranch        *b_xtgtTPCKurama;   //!
   TBranch        *b_ytgtTPCKurama;   //!
   TBranch        *b_utgtTPCKurama;   //!
   TBranch        *b_vtgtTPCKurama;   //!
   TBranch        *b_thetaTPCKurama;   //!
   TBranch        *b_pathTPCKurama;   //!
   TBranch        *b_lhtofTPCKurama;   //!
   TBranch        *b_xhtofTPCKurama;   //!
   TBranch        *b_yhtofTPCKurama;   //!
   TBranch        *b_lvpTPCKurama;   //!
   TBranch        *b_xvpTPCKurama;   //!
   TBranch        *b_yvpTPCKurama;   //!
   TBranch        *b_isgoodTPC;   //!
   TBranch        *b_insideTPC;   //!
   TBranch        *b_vtxTPC;   //!
   TBranch        *b_vtyTPC;   //!
   TBranch        *b_vtzTPC;   //!
   TBranch        *b_closeDistTPC;   //!
   TBranch        *b_MissMassTPC;   //!
   TBranch        *b_MissMassCorrTPC;   //!
   TBranch        *b_MissMassCorrDETPC;   //!
   TBranch        *b_pOrgTPC;   //!
   TBranch        *b_pCalcTPC;   //!
   TBranch        *b_pCorrTPC;   //!
   TBranch        *b_pCorrDETPC;   //!
   TBranch        *b_xbTPC;   //!
   TBranch        *b_ybTPC;   //!
   TBranch        *b_ubTPC;   //!
   TBranch        *b_vbTPC;   //!
   TBranch        *b_xsTPC;   //!
   TBranch        *b_ysTPC;   //!
   TBranch        *b_usTPC;   //!
   TBranch        *b_vsTPC;   //!
   TBranch        *b_nKK;   //!
   TBranch        *b_Kflag;   //!
   TBranch        *b_MissMass;   //!
   TBranch        *b_MissMassCorr;   //!
   TBranch        *b_MissMassCorrDE;   //!
   TBranch        *b_vtx;   //!
   TBranch        *b_vty;   //!
   TBranch        *b_vtz;   //!
   TBranch        *b_pOrg;   //!
   TBranch        *b_pCalc;   //!
   TBranch        *b_pCorr;   //!
   TBranch        *b_pCorrDE;   //!
   TBranch        *b_xb;   //!
   TBranch        *b_yb;   //!
   TBranch        *b_ub;   //!
   TBranch        *b_vb;   //!
   TBranch        *b_xs;   //!
   TBranch        *b_ys;   //!
   TBranch        *b_us;   //!
   TBranch        *b_vs;   //!
   TBranch        *b_KmMom_x;   //!
   TBranch        *b_KmMom_y;   //!
   TBranch        *b_KmMom_z;   //!
   TBranch        *b_KpMom_x;   //!
   TBranch        *b_KpMom_y;   //!
   TBranch        *b_KpMom_z;   //!
   TBranch        *b_BE;   //!
   TBranch        *b_BETPC;   //!
   TBranch        *b_BE_LL;   //!
   TBranch        *b_BETPC_LL;   //!
   TBranch        *b_GFKKXiProductionVtx_x;   //!
   TBranch        *b_GFKKXiProductionVtx_y;   //!
   TBranch        *b_GFKKXiProductionVtx_z;   //!
   TBranch        *b_GFKKLLProductionVtx_x;   //!
   TBranch        *b_GFKKLLProductionVtx_y;   //!
   TBranch        *b_GFKKLLProductionVtx_z;   //!
   TBranch        *b_GFKKLProductionVtx_x1;   //!
   TBranch        *b_GFKKLProductionVtx_y1;   //!
   TBranch        *b_GFKKLProductionVtx_z1;   //!
   TBranch        *b_GFKKLProductionVtx_x2;   //!
   TBranch        *b_GFKKLProductionVtx_y2;   //!
   TBranch        *b_GFKKLProductionVtx_z2;   //!
   TBranch        *b_GFKKLProductionVtx_x;   //!
   TBranch        *b_GFKKLProductionVtx_y;   //!
   TBranch        *b_GFKKLProductionVtx_z;   //!
   TBranch        *b_Xiflag;   //!
   TBranch        *b_XiPflag;   //!
   TBranch        *b_XiMass;   //!
   TBranch        *b_XiDecayVtx_x;   //!
   TBranch        *b_XiDecayVtx_y;   //!
   TBranch        *b_XiDecayVtx_z;   //!
   TBranch        *b_XiMom;   //!
   TBranch        *b_XiMom_x;   //!
   TBranch        *b_XiMom_y;   //!
   TBranch        *b_XiMom_z;   //!
   TBranch        *b_XiVtxCloseDist;   //!
   TBranch        *b_LambdaMass;   //!
   TBranch        *b_LambdaDecayVtx_x;   //!
   TBranch        *b_LambdaDecayVtx_y;   //!
   TBranch        *b_LambdaDecayVtx_z;   //!
   TBranch        *b_LambdaMom;   //!
   TBranch        *b_LambdaMom_x;   //!
   TBranch        *b_LambdaMom_y;   //!
   TBranch        *b_LambdaMom_z;   //!
   TBranch        *b_LambdaVtxCloseDist;   //!
   TBranch        *b_XiTarget_x;   //!
   TBranch        *b_XiTarget_y;   //!
   TBranch        *b_XiTarget_z;   //!
   TBranch        *b_XiTargetMom;   //!
   TBranch        *b_XiTargetMom_x;   //!
   TBranch        *b_XiTargetMom_y;   //!
   TBranch        *b_XiTargetMom_z;   //!
   TBranch        *b_XiTargetCloseDist;   //!
   TBranch        *b_GFXiMass;   //!
   TBranch        *b_GFXiDecayVtx_x;   //!
   TBranch        *b_GFXiDecayVtx_y;   //!
   TBranch        *b_GFXiDecayVtx_z;   //!
   TBranch        *b_GFXiMom;   //!
   TBranch        *b_GFXiMom_x;   //!
   TBranch        *b_GFXiMom_y;   //!
   TBranch        *b_GFXiMom_z;   //!
   TBranch        *b_GFXiVtxCloseDist;   //!
   TBranch        *b_GFXiKKVtx_x;   //!
   TBranch        *b_GFXiKKVtx_y;   //!
   TBranch        *b_GFXiKKVtx_z;   //!
   TBranch        *b_GFXiKKVtxMom;   //!
   TBranch        *b_GFXiKKVtxMom_x;   //!
   TBranch        *b_GFXiKKVtxMom_y;   //!
   TBranch        *b_GFXiKKVtxMom_z;   //!
   TBranch        *b_GFXiKKVtxCloseDist;   //!
   TBranch        *b_GFXiProductionVtx_x;   //!
   TBranch        *b_GFXiProductionVtx_y;   //!
   TBranch        *b_GFXiProductionVtx_z;   //!
   TBranch        *b_GFXiProductionVtxMom;   //!
   TBranch        *b_GFXiProductionVtxMom_x;   //!
   TBranch        *b_GFXiProductionVtxMom_y;   //!
   TBranch        *b_GFXiProductionVtxMom_z;   //!
   TBranch        *b_GFXiProductionVtxCloseDist;   //!
   TBranch        *b_GFXiTrackLen;   //!
   TBranch        *b_GFXiTof;   //!
   TBranch        *b_GFXiMomLoss;   //!
   TBranch        *b_GFXiExcitation;   //!
   TBranch        *b_GFXiTarget_x;   //!
   TBranch        *b_GFXiTarget_y;   //!
   TBranch        *b_GFXiTarget_z;   //!
   TBranch        *b_GFXiTargetMom;   //!
   TBranch        *b_GFXiTargetMom_x;   //!
   TBranch        *b_GFXiTargetMom_y;   //!
   TBranch        *b_GFXiTargetMom_z;   //!
   TBranch        *b_GFXiTargetCloseDist;   //!
   TBranch        *b_GFXiTargetCenter_x;   //!
   TBranch        *b_GFXiTargetCenter_y;   //!
   TBranch        *b_GFXiTargetCenter_z;   //!
   TBranch        *b_GFXiTargetCenterMom;   //!
   TBranch        *b_GFXiTargetCenterMom_x;   //!
   TBranch        *b_GFXiTargetCenterMom_y;   //!
   TBranch        *b_GFXiTargetCenterMom_z;   //!
   TBranch        *b_GFXiTargetCenterCloseDist;   //!
   TBranch        *b_LLflag;   //!
   TBranch        *b_LambdaLambdaVtx_x;   //!
   TBranch        *b_LambdaLambdaVtx_y;   //!
   TBranch        *b_LambdaLambdaVtx_z;   //!
   TBranch        *b_LambdaLambdaCloseDist;   //!
   TBranch        *b_LambdaMass1;   //!
   TBranch        *b_LambdaDecayVtx_x1;   //!
   TBranch        *b_LambdaDecayVtx_y1;   //!
   TBranch        *b_LambdaDecayVtx_z1;   //!
   TBranch        *b_LambdaMom1;   //!
   TBranch        *b_LambdaMom_x1;   //!
   TBranch        *b_LambdaMom_y1;   //!
   TBranch        *b_LambdaMom_z1;   //!
   TBranch        *b_LambdaVtxCloseDist1;   //!
   TBranch        *b_LambdaTargetCloseDist1;   //!
   TBranch        *b_LambdaTargetCloseVtx_x1;   //!
   TBranch        *b_LambdaTargetCloseVtx_y1;   //!
   TBranch        *b_LambdaTargetCloseVtx_z1;   //!
   TBranch        *b_LambdaMass2;   //!
   TBranch        *b_LambdaDecayVtx_x2;   //!
   TBranch        *b_LambdaDecayVtx_y2;   //!
   TBranch        *b_LambdaDecayVtx_z2;   //!
   TBranch        *b_LambdaMom2;   //!
   TBranch        *b_LambdaMom_x2;   //!
   TBranch        *b_LambdaMom_y2;   //!
   TBranch        *b_LambdaMom_z2;   //!
   TBranch        *b_LambdaVtxCloseDist2;   //!
   TBranch        *b_LambdaTargetCloseDist2;   //!
   TBranch        *b_LambdaTargetCloseVtx_x2;   //!
   TBranch        *b_LambdaTargetCloseVtx_y2;   //!
   TBranch        *b_LambdaTargetCloseVtx_z2;   //!
   TBranch        *b_GFLLExcitation;   //!
   TBranch        *b_GFLambdaMass1;   //!
   TBranch        *b_GFLambdaDecayVtx_x1;   //!
   TBranch        *b_GFLambdaDecayVtx_y1;   //!
   TBranch        *b_GFLambdaDecayVtx_z1;   //!
   TBranch        *b_GFLambdaMom1;   //!
   TBranch        *b_GFLambdaMom_x1;   //!
   TBranch        *b_GFLambdaMom_y1;   //!
   TBranch        *b_GFLambdaMom_z1;   //!
   TBranch        *b_GFLambdaVtxCloseDist1;   //!
   TBranch        *b_GFLambdaMass2;   //!
   TBranch        *b_GFLambdaDecayVtx_x2;   //!
   TBranch        *b_GFLambdaDecayVtx_y2;   //!
   TBranch        *b_GFLambdaDecayVtx_z2;   //!
   TBranch        *b_GFLambdaMom2;   //!
   TBranch        *b_GFLambdaMom_x2;   //!
   TBranch        *b_GFLambdaMom_y2;   //!
   TBranch        *b_GFLambdaMom_z2;   //!
   TBranch        *b_GFLambdaVtxCloseDist2;   //!
   TBranch        *b_GFLambdaMass1_Alt;   //!
   TBranch        *b_GFLambdaDecayVtx_x1_Alt;   //!
   TBranch        *b_GFLambdaDecayVtx_y1_Alt;   //!
   TBranch        *b_GFLambdaDecayVtx_z1_Alt;   //!
   TBranch        *b_GFLambdaMom1_Alt;   //!
   TBranch        *b_GFLambdaMom_x1_Alt;   //!
   TBranch        *b_GFLambdaMom_y1_Alt;   //!
   TBranch        *b_GFLambdaMom_z1_Alt;   //!
   TBranch        *b_GFLambdaVtxCloseDist1_Alt;   //!
   TBranch        *b_GFLambdaTargetCloseDist1_Alt;   //!
   TBranch        *b_GFLambdaTarget_x1_Alt;   //!
   TBranch        *b_GFLambdaTarget_y1_Alt;   //!
   TBranch        *b_GFLambdaTarget_z1_Alt;   //!
   TBranch        *b_GFLambdaTargetCenterCloseDist1_Alt;   //!
   TBranch        *b_GFLambdaTargetCenter_x1_Alt;   //!
   TBranch        *b_GFLambdaTargetCenter_y1_Alt;   //!
   TBranch        *b_GFLambdaTargetCenter_z1_Alt;   //!
   TBranch        *b_GFLambdaMass2_Alt;   //!
   TBranch        *b_GFLambdaDecayVtx_x2_Alt;   //!
   TBranch        *b_GFLambdaDecayVtx_y2_Alt;   //!
   TBranch        *b_GFLambdaDecayVtx_z2_Alt;   //!
   TBranch        *b_GFLambdaMom2_Alt;   //!
   TBranch        *b_GFLambdaMom_x2_Alt;   //!
   TBranch        *b_GFLambdaMom_y2_Alt;   //!
   TBranch        *b_GFLambdaMom_z2_Alt;   //!
   TBranch        *b_GFLambdaVtxCloseDist2_Alt;   //!
   TBranch        *b_GFLambdaTargetCloseDist2_Alt;   //!
   TBranch        *b_GFLambdaTarget_x2_Alt;   //!
   TBranch        *b_GFLambdaTarget_y2_Alt;   //!
   TBranch        *b_GFLambdaTarget_z2_Alt;   //!
   TBranch        *b_GFLambdaTargetCenterCloseDist2_Alt;   //!
   TBranch        *b_GFLambdaTargetCenter_x2_Alt;   //!
   TBranch        *b_GFLambdaTargetCenter_y2_Alt;   //!
   TBranch        *b_GFLambdaTargetCenter_z2_Alt;   //!
   TBranch        *b_GFLambdaLambdaVtx_x;   //!
   TBranch        *b_GFLambdaLambdaVtx_y;   //!
   TBranch        *b_GFLambdaLambdaVtx_z;   //!
   TBranch        *b_GFLambdaLambdaCloseDist;   //!
   TBranch        *b_GFLambdaTarget_x1;   //!
   TBranch        *b_GFLambdaTarget_y1;   //!
   TBranch        *b_GFLambdaTarget_z1;   //!
   TBranch        *b_GFLambdaTargetCloseDist1;   //!
   TBranch        *b_GFLambdaTarget_x2;   //!
   TBranch        *b_GFLambdaTarget_y2;   //!
   TBranch        *b_GFLambdaTarget_z2;   //!
   TBranch        *b_GFLambdaTargetCloseDist2;   //!
   TBranch        *b_GFLambdaTargetCenter_x1;   //!
   TBranch        *b_GFLambdaTargetCenter_y1;   //!
   TBranch        *b_GFLambdaTargetCenter_z1;   //!
   TBranch        *b_GFLambdaTargetCenterCloseDist1;   //!
   TBranch        *b_GFLambdaTargetCenter_x2;   //!
   TBranch        *b_GFLambdaTargetCenter_y2;   //!
   TBranch        *b_GFLambdaTargetCenter_z2;   //!
   TBranch        *b_GFLambdaTargetCenterCloseDist2;   //!
   TBranch        *b_GFLambdaProductionVtx_x1;   //!
   TBranch        *b_GFLambdaProductionVtx_y1;   //!
   TBranch        *b_GFLambdaProductionVtx_z1;   //!
   TBranch        *b_GFLambdaProductionVtxCloseDist1;   //!
   TBranch        *b_GFLambdaTrackLen1;   //!
   TBranch        *b_GFLambdaTof1;   //!
   TBranch        *b_GFLambdaProductionVtx_x2;   //!
   TBranch        *b_GFLambdaProductionVtx_y2;   //!
   TBranch        *b_GFLambdaProductionVtx_z2;   //!
   TBranch        *b_GFLambdaProductionVtxCloseDist2;   //!
   TBranch        *b_GFLambdaTrackLen2;   //!
   TBranch        *b_GFLambdaTof2;   //!
   TBranch        *b_KFLLExcitation;   //!
   TBranch        *b_KFKKLLProductionVtxChisqr;   //!
   TBranch        *b_KFKKLLProductionVtx_x;   //!
   TBranch        *b_KFKKLLProductionVtx_y;   //!
   TBranch        *b_KFKKLLProductionVtx_z;   //!
   TBranch        *b_KFKKLProductionVtx_x1;   //!
   TBranch        *b_KFKKLProductionVtx_y1;   //!
   TBranch        *b_KFKKLProductionVtx_z1;   //!
   TBranch        *b_KFKKLProductionVtx_x2;   //!
   TBranch        *b_KFKKLProductionVtx_y2;   //!
   TBranch        *b_KFKKLProductionVtx_z2;   //!
   TBranch        *b_KFLambdaLambdaVtx_x;   //!
   TBranch        *b_KFLambdaLambdaVtx_y;   //!
   TBranch        *b_KFLambdaLambdaVtx_z;   //!
   TBranch        *b_KFLambdaLambdaCloseDist;   //!
   TBranch        *b_KFLambdaProductionVtx_x1;   //!
   TBranch        *b_KFLambdaProductionVtx_y1;   //!
   TBranch        *b_KFLambdaProductionVtx_z1;   //!
   TBranch        *b_KFLambdaProductionVtxCloseDist1;   //!
   TBranch        *b_KFLambdaTrackLen1;   //!
   TBranch        *b_KFLambdaTof1;   //!
   TBranch        *b_KFLambdaProductionVtx_x2;   //!
   TBranch        *b_KFLambdaProductionVtx_y2;   //!
   TBranch        *b_KFLambdaProductionVtx_z2;   //!
   TBranch        *b_KFLambdaProductionVtxCloseDist2;   //!
   TBranch        *b_KFLambdaTrackLen2;   //!
   TBranch        *b_KFLambdaTof2;   //!
   TBranch        *b_KFLambdaPval1;   //!
   TBranch        *b_KFLambdaPull1;   //!
   TBranch        *b_KFLambdaChisqr1;   //!
   TBranch        *b_KFLambdaMom1;   //!
   TBranch        *b_KFLambdaMom_x1;   //!
   TBranch        *b_KFLambdaMom_y1;   //!
   TBranch        *b_KFLambdaMom_z1;   //!
   TBranch        *b_KFLambdaPval2;   //!
   TBranch        *b_KFLambdaPull2;   //!
   TBranch        *b_KFLambdaChisqr2;   //!
   TBranch        *b_KFLambdaMom2;   //!
   TBranch        *b_KFLambdaMom_x2;   //!
   TBranch        *b_KFLambdaMom_y2;   //!
   TBranch        *b_KFLambdaMom_z2;   //!
   TBranch        *b_Lflag;   //!
   TBranch        *b_LambdaTargetCloseDist;   //!
   TBranch        *b_LambdaTargetCloseVtx_x;   //!
   TBranch        *b_LambdaTargetCloseVtx_y;   //!
   TBranch        *b_LambdaTargetCloseVtx_z;   //!
   TBranch        *b_GFLambdaMass;   //!
   TBranch        *b_GFLambdaDecayVtx_x;   //!
   TBranch        *b_GFLambdaDecayVtx_y;   //!
   TBranch        *b_GFLambdaDecayVtx_z;   //!
   TBranch        *b_GFLambdaMom;   //!
   TBranch        *b_GFLambdaMom_x;   //!
   TBranch        *b_GFLambdaMom_y;   //!
   TBranch        *b_GFLambdaMom_z;   //!
   TBranch        *b_GFLambdaVtxCloseDist;   //!
   TBranch        *b_GFLambdaTargetCloseDist;   //!
   TBranch        *b_GFLambdaTarget_x;   //!
   TBranch        *b_GFLambdaTarget_y;   //!
   TBranch        *b_GFLambdaTarget_z;   //!
   TBranch        *b_GFLambdaTargetCenter_x;   //!
   TBranch        *b_GFLambdaTargetCenter_y;   //!
   TBranch        *b_GFLambdaTargetCenter_z;   //!
   TBranch        *b_GFLambdaTargetCenterCloseDist;   //!
   TBranch        *b_GFLambdaProductionVtx_x;   //!
   TBranch        *b_GFLambdaProductionVtx_y;   //!
   TBranch        *b_GFLambdaProductionVtx_z;   //!
   TBranch        *b_GFLambdaProductionVtxCloseDist;   //!
   TBranch        *b_GFLambdaTrackLen;   //!
   TBranch        *b_GFLambdaTof;   //!
   TBranch        *b_KFKKLProductionVtx_x;   //!
   TBranch        *b_KFKKLProductionVtx_y;   //!
   TBranch        *b_KFKKLProductionVtx_z;   //!
   TBranch        *b_KFLambdaProductionVtx_x;   //!
   TBranch        *b_KFLambdaProductionVtx_y;   //!
   TBranch        *b_KFLambdaProductionVtx_z;   //!
   TBranch        *b_KFLambdaProductionVtxCloseDist;   //!
   TBranch        *b_KFLambdaTrackLen;   //!
   TBranch        *b_KFLambdaTof;   //!
   TBranch        *b_LPhiflag;   //!
   TBranch        *b_PhiMass;   //!
   TBranch        *b_PhiDecayVtx_x;   //!
   TBranch        *b_PhiDecayVtx_y;   //!
   TBranch        *b_PhiDecayVtx_z;   //!
   TBranch        *b_PhiCosKK;   //!
   TBranch        *b_PhiMom;   //!
   TBranch        *b_PhiMom_x;   //!
   TBranch        *b_PhiMom_y;   //!
   TBranch        *b_PhiMom_z;   //!
   TBranch        *b_PhiVtxCloseDist;   //!
   TBranch        *b_PhiDecaysTrackId;   //!
   TBranch        *b_PhiDecaysMom;   //!
   TBranch        *b_PhiDecaysMom_x;   //!
   TBranch        *b_PhiDecaysMom_y;   //!
   TBranch        *b_PhiDecaysMom_z;   //!
   TBranch        *b_GFPhiMass;   //!
   TBranch        *b_GFPhiDecayVtx_x;   //!
   TBranch        *b_GFPhiDecayVtx_y;   //!
   TBranch        *b_GFPhiDecayVtx_z;   //!
   TBranch        *b_GFPhiCosKK;   //!
   TBranch        *b_GFPhiMom;   //!
   TBranch        *b_GFPhiMom_x;   //!
   TBranch        *b_GFPhiMom_y;   //!
   TBranch        *b_GFPhiMom_z;   //!
   TBranch        *b_GFPhiVtxCloseDist;   //!
   TBranch        *b_GFPhiProductionVtxCloseDist;   //!
   TBranch        *b_GFPhiKmMassSquare;   //!
   TBranch        *b_GFPhiKmInvbeta;   //!
   TBranch        *b_GFPhiKpMassSquare;   //!
   TBranch        *b_GFPhiKpInvbeta;   //!
   TBranch        *b_GFPhiDecaysMom;   //!
   TBranch        *b_GFPhiDecaysMom_x;   //!
   TBranch        *b_GFPhiDecaysMom_y;   //!
   TBranch        *b_GFPhiDecaysMom_z;   //!
   TBranch        *b_GFPhiMassKurama;   //!
   TBranch        *b_GFPhiMomKurama;   //!
   TBranch        *b_GFPhiMomKurama_x;   //!
   TBranch        *b_GFPhiMomKurama_y;   //!
   TBranch        *b_GFPhiMomKurama_z;   //!
   TBranch        *b_LPiflag;   //!
   TBranch        *b_PiPiflag;   //!
   TBranch        *b_Pimflag;   //!
   TBranch        *b_Emptyflag;   //!
   TBranch        *b_GFXiDecaysNhit;   //!
   TBranch        *b_GFXiDecaysChisqr;   //!
   TBranch        *b_GFXiDecaysCharge;   //!
   TBranch        *b_GFXiDecaysTof;   //!
   TBranch        *b_GFXiDecaysPval;   //!
   TBranch        *b_GFXiDecaysPdgcode;   //!
   TBranch        *b_GFXiDecaysHtofId;   //!
   TBranch        *b_GFXiDecaysTrackLen;   //!
   TBranch        *b_GFXiDecaysTrackTof;   //!
   TBranch        *b_GFXiDecaysMassSquare;   //!
   TBranch        *b_GFXiDecaysInvbeta;   //!
   TBranch        *b_GFXiDecaysMom;   //!
   TBranch        *b_GFXiDecaysMom_x;   //!
   TBranch        *b_GFXiDecaysMom_y;   //!
   TBranch        *b_GFXiDecaysMom_z;   //!
   TBranch        *b_GFXiDecaysMomCM;   //!
   TBranch        *b_GFXiDecaysMomCM_x;   //!
   TBranch        *b_GFXiDecaysMomCM_y;   //!
   TBranch        *b_GFXiDecaysMomCM_z;   //!
   TBranch        *b_GFXiDecaysMomLoss;   //!
   TBranch        *b_GFXiDecaysELoss;   //!
   TBranch        *b_XiDecaysTrackId;   //!
   TBranch        *b_XiDecaysMom;   //!
   TBranch        *b_XiDecaysMom_x;   //!
   TBranch        *b_XiDecaysMom_y;   //!
   TBranch        *b_XiDecaysMom_z;   //!
   TBranch        *b_XiDecaysMomCM;   //!
   TBranch        *b_XiDecaysMomCM_x;   //!
   TBranch        *b_XiDecaysMomCM_y;   //!
   TBranch        *b_XiDecaysMomCM_z;   //!
   TBranch        *b_GFLLDecaysNhit;   //!
   TBranch        *b_GFLLDecaysChisqr;   //!
   TBranch        *b_GFLLDecaysCharge;   //!
   TBranch        *b_GFLLDecaysTof;   //!
   TBranch        *b_GFLLDecaysPval;   //!
   TBranch        *b_GFLLDecaysPdgcode;   //!
   TBranch        *b_GFLLDecaysHtofId;   //!
   TBranch        *b_GFLLDecaysTrackLen;   //!
   TBranch        *b_GFLLDecaysTrackTof;   //!
   TBranch        *b_GFLLDecaysMassSquare;   //!
   TBranch        *b_GFLLDecaysInvbeta;   //!
   TBranch        *b_GFLLDecaysMom;   //!
   TBranch        *b_GFLLDecaysMom_x;   //!
   TBranch        *b_GFLLDecaysMom_y;   //!
   TBranch        *b_GFLLDecaysMom_z;   //!
   TBranch        *b_GFLLDecaysMomCM;   //!
   TBranch        *b_GFLLDecaysMomCM_x;   //!
   TBranch        *b_GFLLDecaysMomCM_y;   //!
   TBranch        *b_GFLLDecaysMomCM_z;   //!
   TBranch        *b_GFLLDecaysMomLoss;   //!
   TBranch        *b_GFLLDecaysELoss;   //!
   TBranch        *b_LLDecaysTrackId;   //!
   TBranch        *b_LLDecaysMom;   //!
   TBranch        *b_LLDecaysMom_x;   //!
   TBranch        *b_LLDecaysMom_y;   //!
   TBranch        *b_LLDecaysMom_z;   //!
   TBranch        *b_LLDecaysMomCM;   //!
   TBranch        *b_LLDecaysMomCM_x;   //!
   TBranch        *b_LLDecaysMomCM_y;   //!
   TBranch        *b_LLDecaysMomCM_z;   //!
   TBranch        *b_GFDecaysNhit;   //!
   TBranch        *b_GFDecaysChisqr;   //!
   TBranch        *b_GFDecaysCharge;   //!
   TBranch        *b_GFDecaysTof;   //!
   TBranch        *b_GFDecaysPval;   //!
   TBranch        *b_GFDecaysPdgcode;   //!
   TBranch        *b_GFDecaysHtofId;   //!
   TBranch        *b_GFDecaysTrackLen;   //!
   TBranch        *b_GFDecaysTrackTof;   //!
   TBranch        *b_GFDecaysMassSquare;   //!
   TBranch        *b_GFDecaysInvbeta;   //!
   TBranch        *b_GFDecaysMom;   //!
   TBranch        *b_GFDecaysMom_x;   //!
   TBranch        *b_GFDecaysMom_y;   //!
   TBranch        *b_GFDecaysMom_z;   //!
   TBranch        *b_GFDecaysMomCM;   //!
   TBranch        *b_GFDecaysMomCM_x;   //!
   TBranch        *b_GFDecaysMomCM_y;   //!
   TBranch        *b_GFDecaysMomCM_z;   //!
   TBranch        *b_GFDecaysMomLoss;   //!
   TBranch        *b_GFDecaysELoss;   //!
   TBranch        *b_DecaysTrackId;   //!
   TBranch        *b_DecaysMom;   //!
   TBranch        *b_DecaysMom_x;   //!
   TBranch        *b_DecaysMom_y;   //!
   TBranch        *b_DecaysMom_z;   //!
   TBranch        *b_DecaysMomCM;   //!
   TBranch        *b_DecaysMomCM_x;   //!
   TBranch        *b_DecaysMomCM_y;   //!
   TBranch        *b_DecaysMomCM_z;   //!
   TBranch        *b_AccidentalMultiplicity;   //!
   TBranch        *b_AccidentalTrackId;   //!
   TBranch        *b_LLResidualsMultiplicity;   //!
   TBranch        *b_LLpMultiplicity;   //!
   TBranch        *b_LLpipMultiplicity;   //!
   TBranch        *b_LLpimMultiplicity;   //!
   TBranch        *b_LLepMultiplicity;   //!
   TBranch        *b_LLemMultiplicity;   //!
   TBranch        *b_LLppipMultiplicity;   //!
   TBranch        *b_LLResidualsTrackId;   //!
   TBranch        *b_LLResidualsMassSquare;   //!
   TBranch        *b_LLResidualsInvbeta;   //!
   TBranch        *b_LLResidualsCloseDistTgt;   //!
   TBranch        *b_LLResidualsGFProductionVtxCloseDist;   //!
   TBranch        *b_LLResidualsKFProductionVtxCloseDist;   //!
   TBranch        *b_LLResidualsMom;   //!
   TBranch        *b_LLResidualsMom_x;   //!
   TBranch        *b_LLResidualsMom_y;   //!
   TBranch        *b_LLResidualsMom_z;   //!
   TBranch        *b_LLResidualsCharge;   //!
   TBranch        *b_XiResidualsMultiplicity;   //!
   TBranch        *b_XipMultiplicity;   //!
   TBranch        *b_XipipMultiplicity;   //!
   TBranch        *b_XipimMultiplicity;   //!
   TBranch        *b_XiepMultiplicity;   //!
   TBranch        *b_XiemMultiplicity;   //!
   TBranch        *b_XippipMultiplicity;   //!
   TBranch        *b_XiResidualsTrackId;   //!
   TBranch        *b_XiResidualsMassSquare;   //!
   TBranch        *b_XiResidualsInvbeta;   //!
   TBranch        *b_XiResidualsCloseDistTgt;   //!
   TBranch        *b_XiResidualsGFProductionVtxCloseDist;   //!
   TBranch        *b_XiResidualsKFProductionVtxCloseDist;   //!
   TBranch        *b_XiResidualsMom;   //!
   TBranch        *b_XiResidualsMom_x;   //!
   TBranch        *b_XiResidualsMom_y;   //!
   TBranch        *b_XiResidualsMom_z;   //!
   TBranch        *b_XiResidualsCharge;   //!
   TBranch        *b_ResidualsMultiplicity;   //!
   TBranch        *b_pMultiplicity;   //!
   TBranch        *b_pipMultiplicity;   //!
   TBranch        *b_pimMultiplicity;   //!
   TBranch        *b_epMultiplicity;   //!
   TBranch        *b_emMultiplicity;   //!
   TBranch        *b_ppipMultiplicity;   //!
   TBranch        *b_ResidualsTrackId;   //!
   TBranch        *b_ResidualsMassSquare;   //!
   TBranch        *b_ResidualsInvbeta;   //!
   TBranch        *b_ResidualsCloseDistTgt;   //!
   TBranch        *b_ResidualsGFProductionVtxCloseDist;   //!
   TBranch        *b_ResidualsKFProductionVtxCloseDist;   //!
   TBranch        *b_ResidualsMom;   //!
   TBranch        *b_ResidualsMom_x;   //!
   TBranch        *b_ResidualsMom_y;   //!
   TBranch        *b_ResidualsMom_z;   //!
   TBranch        *b_ResidualsCharge;   //!
   TBranch        *b_KFLambdaMomPpi;   //!
   TBranch        *b_KFLambdaMomPpi_x;   //!
   TBranch        *b_KFLambdaMomPpi_y0;   //!
   TBranch        *b_KFLambdaMomPpi_z0;   //!
   TBranch        *b_KFLambdaMom;   //!
   TBranch        *b_KFLambdaMom_x;   //!
   TBranch        *b_KFLambdaMom_y;   //!
   TBranch        *b_KFLambdaMom_z;   //!
   TBranch        *b_KFLambdaChisqr;   //!
   TBranch        *b_KFLambdaPval;   //!
   TBranch        *b_KFLambdaPull;   //!
   TBranch        *b_KFLambdaCovMatrix;   //!
   TBranch        *b_KFXiVtxCloseDist;   //!
   TBranch        *b_KFXiMom;   //!
   TBranch        *b_KFXiMom_x;   //!
   TBranch        *b_KFXiMom_y;   //!
   TBranch        *b_KFXiMom_z;   //!
   TBranch        *b_KFXiChisqr;   //!
   TBranch        *b_KFXiPval;   //!
   TBranch        *b_KFXiMass;   //!
   TBranch        *b_KFXiDecayVtx_x;   //!
   TBranch        *b_KFXiDecayVtx_y;   //!
   TBranch        *b_KFXiDecayVtx_z;   //!
   TBranch        *b_KFXiPull;   //!
   TBranch        *b_KFXiCovMatrix;   //!
   TBranch        *b_KFKKXiProductionVtxChisqr;   //!
   TBranch        *b_KFKKXiProductionVtx_x;   //!
   TBranch        *b_KFKKXiProductionVtx_y;   //!
   TBranch        *b_KFKKXiProductionVtx_z;   //!
   TBranch        *b_KFXiProductionVtx_x;   //!
   TBranch        *b_KFXiProductionVtx_y;   //!
   TBranch        *b_KFXiProductionVtx_z;   //!
   TBranch        *b_KFXiProductionVtxMom;   //!
   TBranch        *b_KFXiProductionVtxMom_x;   //!
   TBranch        *b_KFXiProductionVtxMom_y;   //!
   TBranch        *b_KFXiProductionVtxMom_z;   //!
   TBranch        *b_KFXiProductionVtxCloseDist;   //!
   TBranch        *b_KFXiTrackLen;   //!
   TBranch        *b_KFXiTof;   //!
   TBranch        *b_KFXiMomLoss;   //!
   TBranch        *b_KFXiExcitation;   //!
   TBranch        *b_KFKpXiProductionVtx_x;   //!
   TBranch        *b_KFKpXiProductionVtx_y;   //!
   TBranch        *b_KFKpXiProductionVtx_z;   //!
   TBranch        *b_KFXiProductionVtx_x_KpXi;   //!
   TBranch        *b_KFXiProductionVtx_y_KpXi;   //!
   TBranch        *b_KFXiProductionVtx_z_KpXi;   //!
   TBranch        *b_KFXiProductionVtxMom_KpXi;   //!
   TBranch        *b_KFXiProductionVtxMom_x_KpXi;   //!
   TBranch        *b_KFXiProductionVtxMom_y_KpXi;   //!
   TBranch        *b_KFXiProductionVtxMom_z_KpXi;   //!
   TBranch        *b_KFXiProductionVtxCloseDist_KpXi;   //!
   TBranch        *b_KFXiProductionVtx_x_KK;   //!
   TBranch        *b_KFXiProductionVtx_y_KK;   //!
   TBranch        *b_KFXiProductionVtx_z_KK;   //!
   TBranch        *b_KFXiProductionVtxMom_KK;   //!
   TBranch        *b_KFXiProductionVtxMom_x_KK;   //!
   TBranch        *b_KFXiProductionVtxMom_y_KK;   //!
   TBranch        *b_KFXiProductionVtxMom_z_KK;   //!
   TBranch        *b_KFXiProductionVtxCloseDist_KK;   //!
   TBranch        *b_KFLLDecaysMom;   //!
   TBranch        *b_KFLLDecaysMom_x;   //!
   TBranch        *b_KFLLDecaysMom_y;   //!
   TBranch        *b_KFLLDecaysMom_z;   //!
   TBranch        *b_KFLLDecaysMomCM;   //!
   TBranch        *b_KFLLDecaysMomCM_x;   //!
   TBranch        *b_KFLLDecaysMomCM_y;   //!
   TBranch        *b_KFLLDecaysMomCM_z;   //!
   TBranch        *b_KFXiDecaysMom;   //!
   TBranch        *b_KFXiDecaysMom_x;   //!
   TBranch        *b_KFXiDecaysMom_y;   //!
   TBranch        *b_KFXiDecaysMom_z;   //!
   TBranch        *b_KFXiDecaysMomCM;   //!
   TBranch        *b_KFXiDecaysMomCM_x;   //!
   TBranch        *b_KFXiDecaysMomCM_y;   //!
   TBranch        *b_KFXiDecaysMomCM_z;   //!
   TBranch        *b_KFDecaysMom;   //!
   TBranch        *b_KFDecaysMom_x;   //!
   TBranch        *b_KFDecaysMom_y;   //!
   TBranch        *b_KFDecaysMom_z;   //!
   TBranch        *b_KFDecaysMomCM;   //!
   TBranch        *b_KFDecaysMomCM_x;   //!
   TBranch        *b_KFDecaysMomCM_y;   //!
   TBranch        *b_KFDecaysMomCM_z;   //!
   TBranch        *b_LLnGamma;   //!
   TBranch        *b_LLGammaEpTrackId;   //!
   TBranch        *b_LLGammaEmTrackId;   //!
   TBranch        *b_LLGammaMomId;   //!
   TBranch        *b_LLGammaEpMom_x;   //!
   TBranch        *b_LLGammaEpMom_y;   //!
   TBranch        *b_LLGammaEpMom_z;   //!
   TBranch        *b_LLGammaEmMom;   //!
   TBranch        *b_LLGammaEmMom_x;   //!
   TBranch        *b_LLGammaEmMom_y;   //!
   TBranch        *b_LLGammaEmMom_z;   //!
   TBranch        *b_LLGammaMom;   //!
   TBranch        *b_LLGammaMom_x;   //!
   TBranch        *b_LLGammaMom_y;   //!
   TBranch        *b_LLGammaMom_z;   //!
   TBranch        *b_LLGammaVtxCloseDist;   //!
   TBranch        *b_LLGammaDecayVtx_x;   //!
   TBranch        *b_LLGammaDecayVtx_y;   //!
   TBranch        *b_LLGammaDecayVtx_z;   //!
   TBranch        *b_XinGamma;   //!
   TBranch        *b_XiGammaEpTrackId;   //!
   TBranch        *b_XiGammaEmTrackId;   //!
   TBranch        *b_XiGammaMomId;   //!
   TBranch        *b_XiGammaEpMom_x;   //!
   TBranch        *b_XiGammaEpMom_y;   //!
   TBranch        *b_XiGammaEpMom_z;   //!
   TBranch        *b_XiGammaEmMom;   //!
   TBranch        *b_XiGammaEmMom_x;   //!
   TBranch        *b_XiGammaEmMom_y;   //!
   TBranch        *b_XiGammaEmMom_z;   //!
   TBranch        *b_XiGammaMom;   //!
   TBranch        *b_XiGammaMom_x;   //!
   TBranch        *b_XiGammaMom_y;   //!
   TBranch        *b_XiGammaMom_z;   //!
   TBranch        *b_XiGammaVtxCloseDist;   //!
   TBranch        *b_XiGammaDecayVtx_x;   //!
   TBranch        *b_XiGammaDecayVtx_y;   //!
   TBranch        *b_XiGammaDecayVtx_z;   //!
   TBranch        *b_nGamma;   //!
   TBranch        *b_GammaEpTrackId;   //!
   TBranch        *b_GammaEmTrackId;   //!
   TBranch        *b_GammaMomId;   //!
   TBranch        *b_GammaEpMom_x;   //!
   TBranch        *b_GammaEpMom_y;   //!
   TBranch        *b_GammaEpMom_z;   //!
   TBranch        *b_GammaEmMom;   //!
   TBranch        *b_GammaEmMom_x;   //!
   TBranch        *b_GammaEmMom_y;   //!
   TBranch        *b_GammaEmMom_z;   //!
   TBranch        *b_GammaMom;   //!
   TBranch        *b_GammaMom_x;   //!
   TBranch        *b_GammaMom_y;   //!
   TBranch        *b_GammaMom_z;   //!
   TBranch        *b_GammaVtxCloseDist;   //!
   TBranch        *b_GammaDecayVtx_x;   //!
   TBranch        *b_GammaDecayVtx_y;   //!
   TBranch        *b_GammaDecayVtx_z;   //!

   tpc(TTree *tree=0);
   virtual ~tpc();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef tpc_cxx
tpc::tpc(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../k18-analyzer/rootfiles/Prod/KFCor_th50/run05667_GenfitCarbonVer10.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../../k18-analyzer/rootfiles/Prod/KFCor_th50/run05667_GenfitCarbonVer10.root");
      }
      f->GetObject("tpc",tree);

   }
   Init(tree);
}

tpc::~tpc()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tpc::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tpc::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void tpc::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   trigpat = 0;
   trigflag = 0;
   HtofSeg = 0;
   tHtof = 0;
   dtHtof = 0;
   deHtof = 0;
   posHtof = 0;
   remain_cluster_x = 0;
   remain_cluster_y = 0;
   remain_cluster_z = 0;
   remain_cluster_de = 0;
   remain_cluster_size = 0;
   remain_cluster_layer = 0;
   remain_cluster_row_center = 0;
   remain_cluster_mrow = 0;
   remain_cluster_de_center = 0;
   remain_cluster_x_center = 0;
   remain_cluster_y_center = 0;
   remain_cluster_z_center = 0;
   remain_cluster_houghflag = 0;
   isInTarget = 0;
   nhtrack = 0;
   trackid = 0;
   isXi = 0;
   isBeam = 0;
   isKurama = 0;
   isK18 = 0;
   isAccidental = 0;
   isMultiloop = 0;
   charge = 0;
   pid = 0;
   chisqr = 0;
   pval = 0;
   helix_cx = 0;
   helix_cy = 0;
   helix_z0 = 0;
   helix_r = 0;
   helix_dz = 0;
   dE = 0;
   dEdx = 0;
   mom0 = 0;
   path = 0;
   isElectron = 0;
   nsigma_triton = 0;
   nsigma_deutron = 0;
   nsigma_proton = 0;
   nsigma_kaon = 0;
   nsigma_pion = 0;
   nsigma_electron = 0;
   hitlayer = 0;
   hitpos_x = 0;
   hitpos_y = 0;
   hitpos_z = 0;
   calpos_x = 0;
   calpos_y = 0;
   calpos_z = 0;
   residual = 0;
   residual_x = 0;
   residual_y = 0;
   residual_z = 0;
   resolution_x = 0;
   resolution_y = 0;
   resolution_z = 0;
   pathhit = 0;
   alpha = 0;
   track_cluster_de = 0;
   track_cluster_size = 0;
   track_cluster_mrow = 0;
   track_cluster_de_center = 0;
   track_cluster_x_center = 0;
   track_cluster_y_center = 0;
   track_cluster_z_center = 0;
   track_cluster_row_center = 0;
   helix_t = 0;
   chargeIndistinguishable = 0;
   chisqr_inverted = 0;
   pval_inverted = 0;
   helix_cx_inverted = 0;
   helix_cy_inverted = 0;
   helix_z0_inverted = 0;
   helix_r_inverted = 0;
   helix_dz_inverted = 0;
   mom0_inverted = 0;
   pid_inverted = 0;
   vtx_x = 0;
   vtx_y = 0;
   vtx_z = 0;
   vtx_dist = 0;
   vtx_angle = 0;
   vtxid = 0;
   vtxmom_theta = 0;
   vtxpos_x = 0;
   vtxpos_y = 0;
   vtxpos_z = 0;
   vtxmom_x = 0;
   vtxmom_y = 0;
   vtxmom_z = 0;
   GFcharge = 0;
   GFchisqr = 0;
   GFtof = 0;
   GFpval = 0;
   GFfitstatus = 0;
   GFpdgcode = 0;
   GFnhtrack = 0;
   GFlayer = 0;
   GFpos_x = 0;
   GFpos_y = 0;
   GFpos_z = 0;
   GFmom = 0;
   GFmom_x = 0;
   GFmom_y = 0;
   GFmom_z = 0;
   GFresidual_x = 0;
   GFresidual_y = 0;
   GFresidual_z = 0;
   GFresidual_p = 0;
   GFresidual_px = 0;
   GFresidual_py = 0;
   GFresidual_pz = 0;
   GFinside = 0;
   GFtracklen = 0;
   GFtrack2vtxdist = 0;
   GFcalctof = 0;
   GFsegHtof = 0;
   GFtofHtof = 0;
   GFtdiffHtof = 0;
   GFposHtof = 0;
   GFposx = 0;
   GFposy = 0;
   GFposz = 0;
   GFinvbeta = 0;
   GFm2 = 0;
   nsigma_tritonHtof = 0;
   nsigma_deutronHtof = 0;
   nsigma_protonHtof = 0;
   nsigma_kaonHtof = 0;
   nsigma_pionHtof = 0;
   nsigma_electronHtof = 0;
   isLambda = 0;
   ncombiLambda = 0;
   distLambda = 0;
   angleLambda = 0;
   bestmassLambda = 0;
   massLambda = 0;
   vtxLambda_x = 0;
   vtxLambda_y = 0;
   vtxLambda_z = 0;
   momLambda = 0;
   momLambda_x = 0;
   momLambda_y = 0;
   momLambda_z = 0;
   decaysidLambda = 0;
   decaysmomLambda = 0;
   decaysmomLambda_x = 0;
   decaysmomLambda_y = 0;
   decaysmomLambda_z = 0;
   clusteredVtx_x = 0;
   clusteredVtx_y = 0;
   clusteredVtx_z = 0;
   clusteredVtxid = 0;
   ReconFailedLambdaPId = 0;
   ReconFailedLambdaPiId = 0;
   ReconFailedLambdaMass = 0;
   ReconFailedLambdaDecayVtx_x = 0;
   ReconFailedLambdaDecayVtx_y = 0;
   ReconFailedLambdaDecayVtx_z = 0;
   ReconFailedLambdaMom = 0;
   ReconFailedLambdaMom_x = 0;
   ReconFailedLambdaMom_y = 0;
   ReconFailedLambdaMom_z = 0;
   ReconFailedLambdaPMom = 0;
   ReconFailedLambdaPMom_x = 0;
   ReconFailedLambdaPMom_y = 0;
   ReconFailedLambdaPMom_z = 0;
   ReconFailedLambdaPiMom = 0;
   ReconFailedLambdaPiMom_x = 0;
   ReconFailedLambdaPiMom_y = 0;
   ReconFailedLambdaPiMom_z = 0;
   ReconFailedLambdavtxCloseDist = 0;
   PiPairPipId = 0;
   PiPairPimId = 0;
   PiPairPipMom = 0;
   PiPairPipMom_x = 0;
   PiPairPipMom_y = 0;
   PiPairPipMom_z = 0;
   PiPairPimMom = 0;
   PiPairPimMom_x = 0;
   PiPairPimMom_y = 0;
   PiPairPimMom_z = 0;
   PiPairMom = 0;
   PiPairMom_x = 0;
   PiPairMom_y = 0;
   PiPairMom_z = 0;
   PiPairReconLambdaMass = 0;
   PiPairReconMass = 0;
   PiPairCloseDist = 0;
   pK18 = 0;
   chisqrK18 = 0;
   xtgtK18 = 0;
   ytgtK18 = 0;
   utgtK18 = 0;
   vtgtK18 = 0;
   isgoodTPCK18 = 0;
   chisqrTPCK18 = 0;
   qTPCK18 = 0;
   pTPCK18 = 0;
   xtgtTPCK18 = 0;
   ytgtTPCK18 = 0;
   utgtTPCK18 = 0;
   vtgtTPCK18 = 0;
   thetaTPCK18 = 0;
   lhtofTPCK18 = 0;
   xhtofTPCK18 = 0;
   yhtofTPCK18 = 0;
   lvpTPCK18 = 0;
   xvpTPCK18 = 0;
   yvpTPCK18 = 0;
   pKurama = 0;
   qKurama = 0;
   chisqrKurama = 0;
   m2Kurama = 0;
   xtgtKurama = 0;
   ytgtKurama = 0;
   utgtKurama = 0;
   vtgtKurama = 0;
   thetaKurama = 0;
   pathwcKurama = 0;
   xin = 0;
   yin = 0;
   zin = 0;
   pxin = 0;
   pyin = 0;
   pzin = 0;
   isgoodTPCKurama = 0;
   kflagTPCKurama = 0;
   chisqrTPCKurama = 0;
   pTPCKurama = 0;
   qTPCKurama = 0;
   m2TPCKurama = 0;
   xtgtTPCKurama = 0;
   ytgtTPCKurama = 0;
   utgtTPCKurama = 0;
   vtgtTPCKurama = 0;
   thetaTPCKurama = 0;
   pathTPCKurama = 0;
   lhtofTPCKurama = 0;
   xhtofTPCKurama = 0;
   yhtofTPCKurama = 0;
   lvpTPCKurama = 0;
   xvpTPCKurama = 0;
   yvpTPCKurama = 0;
   isgoodTPC = 0;
   insideTPC = 0;
   vtxTPC = 0;
   vtyTPC = 0;
   vtzTPC = 0;
   closeDistTPC = 0;
   MissMassTPC = 0;
   MissMassCorrTPC = 0;
   MissMassCorrDETPC = 0;
   pOrgTPC = 0;
   pCalcTPC = 0;
   pCorrTPC = 0;
   pCorrDETPC = 0;
   xbTPC = 0;
   ybTPC = 0;
   ubTPC = 0;
   vbTPC = 0;
   xsTPC = 0;
   ysTPC = 0;
   usTPC = 0;
   vsTPC = 0;
   Kflag = 0;
   MissMass = 0;
   MissMassCorr = 0;
   MissMassCorrDE = 0;
   vtx = 0;
   vty = 0;
   vtz = 0;
   pOrg = 0;
   pCalc = 0;
   pCorr = 0;
   pCorrDE = 0;
   xb = 0;
   yb = 0;
   ub = 0;
   vb = 0;
   xs = 0;
   ys = 0;
   us = 0;
   vs = 0;
   KmMom_x = 0;
   KmMom_y = 0;
   KmMom_z = 0;
   KpMom_x = 0;
   KpMom_y = 0;
   KpMom_z = 0;
   BE = 0;
   BETPC = 0;
   BE_LL = 0;
   BETPC_LL = 0;
   KFLambdaPull1 = 0;
   KFLambdaPull2 = 0;
   PhiDecaysTrackId = 0;
   PhiDecaysMom = 0;
   PhiDecaysMom_x = 0;
   PhiDecaysMom_y = 0;
   PhiDecaysMom_z = 0;
   GFPhiDecaysMom = 0;
   GFPhiDecaysMom_x = 0;
   GFPhiDecaysMom_y = 0;
   GFPhiDecaysMom_z = 0;
   GFXiDecaysNhit = 0;
   GFXiDecaysChisqr = 0;
   GFXiDecaysCharge = 0;
   GFXiDecaysTof = 0;
   GFXiDecaysPval = 0;
   GFXiDecaysPdgcode = 0;
   GFXiDecaysHtofId = 0;
   GFXiDecaysTrackLen = 0;
   GFXiDecaysTrackTof = 0;
   GFXiDecaysMassSquare = 0;
   GFXiDecaysInvbeta = 0;
   GFXiDecaysMom = 0;
   GFXiDecaysMom_x = 0;
   GFXiDecaysMom_y = 0;
   GFXiDecaysMom_z = 0;
   GFXiDecaysMomCM = 0;
   GFXiDecaysMomCM_x = 0;
   GFXiDecaysMomCM_y = 0;
   GFXiDecaysMomCM_z = 0;
   GFXiDecaysMomLoss = 0;
   GFXiDecaysELoss = 0;
   XiDecaysTrackId = 0;
   XiDecaysMom = 0;
   XiDecaysMom_x = 0;
   XiDecaysMom_y = 0;
   XiDecaysMom_z = 0;
   XiDecaysMomCM = 0;
   XiDecaysMomCM_x = 0;
   XiDecaysMomCM_y = 0;
   XiDecaysMomCM_z = 0;
   GFLLDecaysNhit = 0;
   GFLLDecaysChisqr = 0;
   GFLLDecaysCharge = 0;
   GFLLDecaysTof = 0;
   GFLLDecaysPval = 0;
   GFLLDecaysPdgcode = 0;
   GFLLDecaysHtofId = 0;
   GFLLDecaysTrackLen = 0;
   GFLLDecaysTrackTof = 0;
   GFLLDecaysMassSquare = 0;
   GFLLDecaysInvbeta = 0;
   GFLLDecaysMom = 0;
   GFLLDecaysMom_x = 0;
   GFLLDecaysMom_y = 0;
   GFLLDecaysMom_z = 0;
   GFLLDecaysMomCM = 0;
   GFLLDecaysMomCM_x = 0;
   GFLLDecaysMomCM_y = 0;
   GFLLDecaysMomCM_z = 0;
   GFLLDecaysMomLoss = 0;
   GFLLDecaysELoss = 0;
   LLDecaysTrackId = 0;
   LLDecaysMom = 0;
   LLDecaysMom_x = 0;
   LLDecaysMom_y = 0;
   LLDecaysMom_z = 0;
   LLDecaysMomCM = 0;
   LLDecaysMomCM_x = 0;
   LLDecaysMomCM_y = 0;
   LLDecaysMomCM_z = 0;
   GFDecaysNhit = 0;
   GFDecaysChisqr = 0;
   GFDecaysCharge = 0;
   GFDecaysTof = 0;
   GFDecaysPval = 0;
   GFDecaysPdgcode = 0;
   GFDecaysHtofId = 0;
   GFDecaysTrackLen = 0;
   GFDecaysTrackTof = 0;
   GFDecaysMassSquare = 0;
   GFDecaysInvbeta = 0;
   GFDecaysMom = 0;
   GFDecaysMom_x = 0;
   GFDecaysMom_y = 0;
   GFDecaysMom_z = 0;
   GFDecaysMomCM = 0;
   GFDecaysMomCM_x = 0;
   GFDecaysMomCM_y = 0;
   GFDecaysMomCM_z = 0;
   GFDecaysMomLoss = 0;
   GFDecaysELoss = 0;
   DecaysTrackId = 0;
   DecaysMom = 0;
   DecaysMom_x = 0;
   DecaysMom_y = 0;
   DecaysMom_z = 0;
   DecaysMomCM = 0;
   DecaysMomCM_x = 0;
   DecaysMomCM_y = 0;
   DecaysMomCM_z = 0;
   AccidentalTrackId = 0;
   LLResidualsTrackId = 0;
   LLResidualsMassSquare = 0;
   LLResidualsInvbeta = 0;
   LLResidualsCloseDistTgt = 0;
   LLResidualsGFProductionVtxCloseDist = 0;
   LLResidualsKFProductionVtxCloseDist = 0;
   LLResidualsMom = 0;
   LLResidualsMom_x = 0;
   LLResidualsMom_y = 0;
   LLResidualsMom_z = 0;
   LLResidualsCharge = 0;
   XiResidualsTrackId = 0;
   XiResidualsMassSquare = 0;
   XiResidualsInvbeta = 0;
   XiResidualsCloseDistTgt = 0;
   XiResidualsGFProductionVtxCloseDist = 0;
   XiResidualsKFProductionVtxCloseDist = 0;
   XiResidualsMom = 0;
   XiResidualsMom_x = 0;
   XiResidualsMom_y = 0;
   XiResidualsMom_z = 0;
   XiResidualsCharge = 0;
   ResidualsTrackId = 0;
   ResidualsMassSquare = 0;
   ResidualsInvbeta = 0;
   ResidualsCloseDistTgt = 0;
   ResidualsGFProductionVtxCloseDist = 0;
   ResidualsKFProductionVtxCloseDist = 0;
   ResidualsMom = 0;
   ResidualsMom_x = 0;
   ResidualsMom_y = 0;
   ResidualsMom_z = 0;
   ResidualsCharge = 0;
   KFLambdaPull = 0;
   KFLambdaCovMatrix = 0;
   KFXiPull = 0;
   KFXiCovMatrix = 0;
   KFLLDecaysMom = 0;
   KFLLDecaysMom_x = 0;
   KFLLDecaysMom_y = 0;
   KFLLDecaysMom_z = 0;
   KFLLDecaysMomCM = 0;
   KFLLDecaysMomCM_x = 0;
   KFLLDecaysMomCM_y = 0;
   KFLLDecaysMomCM_z = 0;
   KFXiDecaysMom = 0;
   KFXiDecaysMom_x = 0;
   KFXiDecaysMom_y = 0;
   KFXiDecaysMom_z = 0;
   KFXiDecaysMomCM = 0;
   KFXiDecaysMomCM_x = 0;
   KFXiDecaysMomCM_y = 0;
   KFXiDecaysMomCM_z = 0;
   KFDecaysMom = 0;
   KFDecaysMom_x = 0;
   KFDecaysMom_y = 0;
   KFDecaysMom_z = 0;
   KFDecaysMomCM = 0;
   KFDecaysMomCM_x = 0;
   KFDecaysMomCM_y = 0;
   KFDecaysMomCM_z = 0;
   LLGammaEpTrackId = 0;
   LLGammaEmTrackId = 0;
   LLGammaMomId = 0;
   LLGammaEpMom_x = 0;
   LLGammaEpMom_y = 0;
   LLGammaEpMom_z = 0;
   LLGammaEmMom = 0;
   LLGammaEmMom_x = 0;
   LLGammaEmMom_y = 0;
   LLGammaEmMom_z = 0;
   LLGammaMom = 0;
   LLGammaMom_x = 0;
   LLGammaMom_y = 0;
   LLGammaMom_z = 0;
   LLGammaVtxCloseDist = 0;
   LLGammaDecayVtx_x = 0;
   LLGammaDecayVtx_y = 0;
   LLGammaDecayVtx_z = 0;
   XiGammaEpTrackId = 0;
   XiGammaEmTrackId = 0;
   XiGammaMomId = 0;
   XiGammaEpMom_x = 0;
   XiGammaEpMom_y = 0;
   XiGammaEpMom_z = 0;
   XiGammaEmMom = 0;
   XiGammaEmMom_x = 0;
   XiGammaEmMom_y = 0;
   XiGammaEmMom_z = 0;
   XiGammaMom = 0;
   XiGammaMom_x = 0;
   XiGammaMom_y = 0;
   XiGammaMom_z = 0;
   XiGammaVtxCloseDist = 0;
   XiGammaDecayVtx_x = 0;
   XiGammaDecayVtx_y = 0;
   XiGammaDecayVtx_z = 0;
   GammaEpTrackId = 0;
   GammaEmTrackId = 0;
   GammaMomId = 0;
   GammaEpMom_x = 0;
   GammaEpMom_y = 0;
   GammaEpMom_z = 0;
   GammaEmMom = 0;
   GammaEmMom_x = 0;
   GammaEmMom_y = 0;
   GammaEmMom_z = 0;
   GammaMom = 0;
   GammaMom_x = 0;
   GammaMom_y = 0;
   GammaMom_z = 0;
   GammaVtxCloseDist = 0;
   GammaDecayVtx_x = 0;
   GammaDecayVtx_y = 0;
   GammaDecayVtx_z = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("status", &status, &b_status);
   fChain->SetBranchAddress("runnum", &runnum, &b_runnum);
   fChain->SetBranchAddress("evnum", &evnum, &b_evnum);
   fChain->SetBranchAddress("trigpat", &trigpat, &b_trigpat);
   fChain->SetBranchAddress("trigflag", &trigflag, &b_trigflag);
   fChain->SetBranchAddress("nhHtof", &nhHtof, &b_nhHtof);
   fChain->SetBranchAddress("HtofSeg", &HtofSeg, &b_HtofSeg);
   fChain->SetBranchAddress("tHtof", &tHtof, &b_tHtof);
   fChain->SetBranchAddress("dtHtof", &dtHtof, &b_dtHtof);
   fChain->SetBranchAddress("deHtof", &deHtof, &b_deHtof);
   fChain->SetBranchAddress("posHtof", &posHtof, &b_posHtof);
   fChain->SetBranchAddress("nclTpc", &nclTpc, &b_nclTpc);
   fChain->SetBranchAddress("remain_nclTpc", &remain_nclTpc, &b_remain_nclTpc);
   fChain->SetBranchAddress("remain_cluster_x", &remain_cluster_x, &b_remain_cluster_x);
   fChain->SetBranchAddress("remain_cluster_y", &remain_cluster_y, &b_remain_cluster_y);
   fChain->SetBranchAddress("remain_cluster_z", &remain_cluster_z, &b_remain_cluster_z);
   fChain->SetBranchAddress("remain_cluster_de", &remain_cluster_de, &b_remain_cluster_de);
   fChain->SetBranchAddress("remain_cluster_size", &remain_cluster_size, &b_remain_cluster_size);
   fChain->SetBranchAddress("remain_cluster_layer", &remain_cluster_layer, &b_remain_cluster_layer);
   fChain->SetBranchAddress("remain_cluster_row_center", &remain_cluster_row_center, &b_remain_cluster_row_center);
   fChain->SetBranchAddress("remain_cluster_mrow", &remain_cluster_mrow, &b_remain_cluster_mrow);
   fChain->SetBranchAddress("remain_cluster_de_center", &remain_cluster_de_center, &b_remain_cluster_de_center);
   fChain->SetBranchAddress("remain_cluster_x_center", &remain_cluster_x_center, &b_remain_cluster_x_center);
   fChain->SetBranchAddress("remain_cluster_y_center", &remain_cluster_y_center, &b_remain_cluster_y_center);
   fChain->SetBranchAddress("remain_cluster_z_center", &remain_cluster_z_center, &b_remain_cluster_z_center);
   fChain->SetBranchAddress("remain_cluster_houghflag", &remain_cluster_houghflag, &b_remain_cluster_houghflag);
   fChain->SetBranchAddress("ntTpc", &ntTpc, &b_ntTpc);
   fChain->SetBranchAddress("isInTarget", &isInTarget, &b_isInTarget);
   fChain->SetBranchAddress("nhtrack", &nhtrack, &b_nhtrack);
   fChain->SetBranchAddress("trackid", &trackid, &b_trackid);
   fChain->SetBranchAddress("isXi", &isXi, &b_isXi);
   fChain->SetBranchAddress("isBeam", &isBeam, &b_isBeam);
   fChain->SetBranchAddress("isKurama", &isKurama, &b_isKurama);
   fChain->SetBranchAddress("isK18", &isK18, &b_isK18);
   fChain->SetBranchAddress("isAccidental", &isAccidental, &b_isAccidental);
   fChain->SetBranchAddress("isMultiloop", &isMultiloop, &b_isMultiloop);
   fChain->SetBranchAddress("charge", &charge, &b_charge);
   fChain->SetBranchAddress("pid", &pid, &b_pid);
   fChain->SetBranchAddress("chisqr", &chisqr, &b_chisqr);
   fChain->SetBranchAddress("pval", &pval, &b_pval);
   fChain->SetBranchAddress("helix_cx", &helix_cx, &b_helix_cx);
   fChain->SetBranchAddress("helix_cy", &helix_cy, &b_helix_cy);
   fChain->SetBranchAddress("helix_z0", &helix_z0, &b_helix_z0);
   fChain->SetBranchAddress("helix_r", &helix_r, &b_helix_r);
   fChain->SetBranchAddress("helix_dz", &helix_dz, &b_helix_dz);
   fChain->SetBranchAddress("dE", &dE, &b_dE);
   fChain->SetBranchAddress("dEdx", &dEdx, &b_dEdx);
   fChain->SetBranchAddress("mom0", &mom0, &b_mom0);
   fChain->SetBranchAddress("path", &path, &b_path);
   fChain->SetBranchAddress("isElectron", &isElectron, &b_isElectron);
   fChain->SetBranchAddress("nsigma_triton", &nsigma_triton, &b_nsigma_triton);
   fChain->SetBranchAddress("nsigma_deutron", &nsigma_deutron, &b_nsigma_deutron);
   fChain->SetBranchAddress("nsigma_proton", &nsigma_proton, &b_nsigma_proton);
   fChain->SetBranchAddress("nsigma_kaon", &nsigma_kaon, &b_nsigma_kaon);
   fChain->SetBranchAddress("nsigma_pion", &nsigma_pion, &b_nsigma_pion);
   fChain->SetBranchAddress("nsigma_electron", &nsigma_electron, &b_nsigma_electron);
   fChain->SetBranchAddress("hitlayer", &hitlayer, &b_hitlayer);
   fChain->SetBranchAddress("hitpos_x", &hitpos_x, &b_hitpos_x);
   fChain->SetBranchAddress("hitpos_y", &hitpos_y, &b_hitpos_y);
   fChain->SetBranchAddress("hitpos_z", &hitpos_z, &b_hitpos_z);
   fChain->SetBranchAddress("calpos_x", &calpos_x, &b_calpos_x);
   fChain->SetBranchAddress("calpos_y", &calpos_y, &b_calpos_y);
   fChain->SetBranchAddress("calpos_z", &calpos_z, &b_calpos_z);
   fChain->SetBranchAddress("residual", &residual, &b_residual);
   fChain->SetBranchAddress("residual_x", &residual_x, &b_residual_x);
   fChain->SetBranchAddress("residual_y", &residual_y, &b_residual_y);
   fChain->SetBranchAddress("residual_z", &residual_z, &b_residual_z);
   fChain->SetBranchAddress("resolution_x", &resolution_x, &b_resolution_x);
   fChain->SetBranchAddress("resolution_y", &resolution_y, &b_resolution_y);
   fChain->SetBranchAddress("resolution_z", &resolution_z, &b_resolution_z);
   fChain->SetBranchAddress("pathhit", &pathhit, &b_pathhit);
   fChain->SetBranchAddress("alpha", &alpha, &b_alpha);
   fChain->SetBranchAddress("track_cluster_de", &track_cluster_de, &b_track_cluster_de);
   fChain->SetBranchAddress("track_cluster_size", &track_cluster_size, &b_track_cluster_size);
   fChain->SetBranchAddress("track_cluster_mrow", &track_cluster_mrow, &b_track_cluster_mrow);
   fChain->SetBranchAddress("track_cluster_de_center", &track_cluster_de_center, &b_track_cluster_de_center);
   fChain->SetBranchAddress("track_cluster_x_center", &track_cluster_x_center, &b_track_cluster_x_center);
   fChain->SetBranchAddress("track_cluster_y_center", &track_cluster_y_center, &b_track_cluster_y_center);
   fChain->SetBranchAddress("track_cluster_z_center", &track_cluster_z_center, &b_track_cluster_z_center);
   fChain->SetBranchAddress("track_cluster_row_center", &track_cluster_row_center, &b_track_cluster_row_center);
   fChain->SetBranchAddress("helix_t", &helix_t, &b_helix_t);
   fChain->SetBranchAddress("chargeIndistinguishable", &chargeIndistinguishable, &b_chargeIndistinguishable);
   fChain->SetBranchAddress("chisqr_inverted", &chisqr_inverted, &b_chisqr_inverted);
   fChain->SetBranchAddress("pval_inverted", &pval_inverted, &b_pval_inverted);
   fChain->SetBranchAddress("helix_cx_inverted", &helix_cx_inverted, &b_helix_cx_inverted);
   fChain->SetBranchAddress("helix_cy_inverted", &helix_cy_inverted, &b_helix_cy_inverted);
   fChain->SetBranchAddress("helix_z0_inverted", &helix_z0_inverted, &b_helix_z0_inverted);
   fChain->SetBranchAddress("helix_r_inverted", &helix_r_inverted, &b_helix_r_inverted);
   fChain->SetBranchAddress("helix_dz_inverted", &helix_dz_inverted, &b_helix_dz_inverted);
   fChain->SetBranchAddress("mom0_inverted", &mom0_inverted, &b_mom0_inverted);
   fChain->SetBranchAddress("pid_inverted", &pid_inverted, &b_pid_inverted);
   fChain->SetBranchAddress("nvtxTpc", &nvtxTpc, &b_nvtxTpc);
   fChain->SetBranchAddress("vtx_x", &vtx_x, &b_vtx_x);
   fChain->SetBranchAddress("vtx_y", &vtx_y, &b_vtx_y);
   fChain->SetBranchAddress("vtx_z", &vtx_z, &b_vtx_z);
   fChain->SetBranchAddress("vtx_dist", &vtx_dist, &b_vtx_dist);
   fChain->SetBranchAddress("vtx_angle", &vtx_angle, &b_vtx_angle);
   fChain->SetBranchAddress("vtxid", &vtxid, &b_vtxid);
   fChain->SetBranchAddress("vtxmom_theta", &vtxmom_theta, &b_vtxmom_theta);
   fChain->SetBranchAddress("vtxpos_x", &vtxpos_x, &b_vtxpos_x);
   fChain->SetBranchAddress("vtxpos_y", &vtxpos_y, &b_vtxpos_y);
   fChain->SetBranchAddress("vtxpos_z", &vtxpos_z, &b_vtxpos_z);
   fChain->SetBranchAddress("vtxmom_x", &vtxmom_x, &b_vtxmom_x);
   fChain->SetBranchAddress("vtxmom_y", &vtxmom_y, &b_vtxmom_y);
   fChain->SetBranchAddress("vtxmom_z", &vtxmom_z, &b_vtxmom_z);
   fChain->SetBranchAddress("GFcharge", &GFcharge, &b_GFcharge);
   fChain->SetBranchAddress("GFchisqr", &GFchisqr, &b_GFchisqr);
   fChain->SetBranchAddress("GFtof", &GFtof, &b_GFtof);
   fChain->SetBranchAddress("GFpval", &GFpval, &b_GFpval);
   fChain->SetBranchAddress("GFfitstatus", &GFfitstatus, &b_GFfitstatus);
   fChain->SetBranchAddress("GFpdgcode", &GFpdgcode, &b_GFpdgcode);
   fChain->SetBranchAddress("GFnhtrack", &GFnhtrack, &b_GFnhtrack);
   fChain->SetBranchAddress("GFlayer", &GFlayer, &b_GFlayer);
   fChain->SetBranchAddress("GFpos_x", &GFpos_x, &b_GFpos_x);
   fChain->SetBranchAddress("GFpos_y", &GFpos_y, &b_GFpos_y);
   fChain->SetBranchAddress("GFpos_z", &GFpos_z, &b_GFpos_z);
   fChain->SetBranchAddress("GFmom", &GFmom, &b_GFmom);
   fChain->SetBranchAddress("GFmom_x", &GFmom_x, &b_GFmom_x);
   fChain->SetBranchAddress("GFmom_y", &GFmom_y, &b_GFmom_y);
   fChain->SetBranchAddress("GFmom_z", &GFmom_z, &b_GFmom_z);
   fChain->SetBranchAddress("GFresidual_x", &GFresidual_x, &b_GFresidual_x);
   fChain->SetBranchAddress("GFresidual_y", &GFresidual_y, &b_GFresidual_y);
   fChain->SetBranchAddress("GFresidual_z", &GFresidual_z, &b_GFresidual_z);
   fChain->SetBranchAddress("GFresidual_p", &GFresidual_p, &b_GFresidual_p);
   fChain->SetBranchAddress("GFresidual_px", &GFresidual_px, &b_GFresidual_px);
   fChain->SetBranchAddress("GFresidual_py", &GFresidual_py, &b_GFresidual_py);
   fChain->SetBranchAddress("GFresidual_pz", &GFresidual_pz, &b_GFresidual_pz);
   fChain->SetBranchAddress("GFntTpc_target", &GFntTpc_target, &b_GFntTpc_target);
   fChain->SetBranchAddress("GFprodvtx_x", &GFprodvtx_x, &b_GFprodvtx_x);
   fChain->SetBranchAddress("GFprodvtx_y", &GFprodvtx_y, &b_GFprodvtx_y);
   fChain->SetBranchAddress("GFprodvtx_z", &GFprodvtx_z, &b_GFprodvtx_z);
   fChain->SetBranchAddress("GFinside", &GFinside, &b_GFinside);
   fChain->SetBranchAddress("GFtracklen", &GFtracklen, &b_GFtracklen);
   fChain->SetBranchAddress("GFtrack2vtxdist", &GFtrack2vtxdist, &b_GFtrack2vtxdist);
   fChain->SetBranchAddress("GFcalctof", &GFcalctof, &b_GFcalctof);
   fChain->SetBranchAddress("GFsegHtof", &GFsegHtof, &b_GFsegHtof);
   fChain->SetBranchAddress("GFtofHtof", &GFtofHtof, &b_GFtofHtof);
   fChain->SetBranchAddress("GFtdiffHtof", &GFtdiffHtof, &b_GFtdiffHtof);
   fChain->SetBranchAddress("GFposHtof", &GFposHtof, &b_GFposHtof);
   fChain->SetBranchAddress("GFposx", &GFposx, &b_GFposx);
   fChain->SetBranchAddress("GFposy", &GFposy, &b_GFposy);
   fChain->SetBranchAddress("GFposz", &GFposz, &b_GFposz);
   fChain->SetBranchAddress("GFinvbeta", &GFinvbeta, &b_GFinvbeta);
   fChain->SetBranchAddress("GFm2", &GFm2, &b_GFm2);
   fChain->SetBranchAddress("nsigma_tritonHtof", &nsigma_tritonHtof, &b_nsigma_tritonHtof);
   fChain->SetBranchAddress("nsigma_deutronHtof", &nsigma_deutronHtof, &b_nsigma_deutronHtof);
   fChain->SetBranchAddress("nsigma_protonHtof", &nsigma_protonHtof, &b_nsigma_protonHtof);
   fChain->SetBranchAddress("nsigma_kaonHtof", &nsigma_kaonHtof, &b_nsigma_kaonHtof);
   fChain->SetBranchAddress("nsigma_pionHtof", &nsigma_pionHtof, &b_nsigma_pionHtof);
   fChain->SetBranchAddress("nsigma_electronHtof", &nsigma_electronHtof, &b_nsigma_electronHtof);
   fChain->SetBranchAddress("isLambda", &isLambda, &b_isLambda);
   fChain->SetBranchAddress("ncombiLambda", &ncombiLambda, &b_ncombiLambda);
   fChain->SetBranchAddress("distLambda", &distLambda, &b_distLambda);
   fChain->SetBranchAddress("angleLambda", &angleLambda, &b_angleLambda);
   fChain->SetBranchAddress("bestmassLambda", &bestmassLambda, &b_bestmassLambda);
   fChain->SetBranchAddress("massLambda", &massLambda, &b_massLambda);
   fChain->SetBranchAddress("vtxLambda_x", &vtxLambda_x, &b_vtxLambda_x);
   fChain->SetBranchAddress("vtxLambda_y", &vtxLambda_y, &b_vtxLambda_y);
   fChain->SetBranchAddress("vtxLambda_z", &vtxLambda_z, &b_vtxLambda_z);
   fChain->SetBranchAddress("momLambda", &momLambda, &b_momLambda);
   fChain->SetBranchAddress("momLambda_x", &momLambda_x, &b_momLambda_x);
   fChain->SetBranchAddress("momLambda_y", &momLambda_y, &b_momLambda_y);
   fChain->SetBranchAddress("momLambda_z", &momLambda_z, &b_momLambda_z);
   fChain->SetBranchAddress("decaysidLambda", &decaysidLambda, &b_decaysidLambda);
   fChain->SetBranchAddress("decaysmomLambda", &decaysmomLambda, &b_decaysmomLambda);
   fChain->SetBranchAddress("decaysmomLambda_x", &decaysmomLambda_x, &b_decaysmomLambda_x);
   fChain->SetBranchAddress("decaysmomLambda_y", &decaysmomLambda_y, &b_decaysmomLambda_y);
   fChain->SetBranchAddress("decaysmomLambda_z", &decaysmomLambda_z, &b_decaysmomLambda_z);
   fChain->SetBranchAddress("nvtxTpcClustered", &nvtxTpcClustered, &b_nvtxTpcClustered);
   fChain->SetBranchAddress("clusteredVtx_x", &clusteredVtx_x, &b_clusteredVtx_x);
   fChain->SetBranchAddress("clusteredVtx_y", &clusteredVtx_y, &b_clusteredVtx_y);
   fChain->SetBranchAddress("clusteredVtx_z", &clusteredVtx_z, &b_clusteredVtx_z);
   fChain->SetBranchAddress("clusteredVtxid", &clusteredVtxid, &b_clusteredVtxid);
   fChain->SetBranchAddress("ncombiReconFailedLambda", &ncombiReconFailedLambda, &b_ncombiReconFailedLambda);
   fChain->SetBranchAddress("ReconFailedLambdaPId", &ReconFailedLambdaPId, &b_ReconFailedLambdaPId);
   fChain->SetBranchAddress("ReconFailedLambdaPiId", &ReconFailedLambdaPiId, &b_ReconFailedLambdaPiId);
   fChain->SetBranchAddress("ReconFailedLambdaMass", &ReconFailedLambdaMass, &b_ReconFailedLambdaMass);
   fChain->SetBranchAddress("ReconFailedLambdaDecayVtx_x", &ReconFailedLambdaDecayVtx_x, &b_ReconFailedLambdaDecayVtx_x);
   fChain->SetBranchAddress("ReconFailedLambdaDecayVtx_y", &ReconFailedLambdaDecayVtx_y, &b_ReconFailedLambdaDecayVtx_y);
   fChain->SetBranchAddress("ReconFailedLambdaDecayVtx_z", &ReconFailedLambdaDecayVtx_z, &b_ReconFailedLambdaDecayVtx_z);
   fChain->SetBranchAddress("ReconFailedLambdaMom", &ReconFailedLambdaMom, &b_ReconFailedLambdaMom);
   fChain->SetBranchAddress("ReconFailedLambdaMom_x", &ReconFailedLambdaMom_x, &b_ReconFailedLambdaMom_x);
   fChain->SetBranchAddress("ReconFailedLambdaMom_y", &ReconFailedLambdaMom_y, &b_ReconFailedLambdaMom_y);
   fChain->SetBranchAddress("ReconFailedLambdaMom_z", &ReconFailedLambdaMom_z, &b_ReconFailedLambdaMom_z);
   fChain->SetBranchAddress("ReconFailedLambdaPMom", &ReconFailedLambdaPMom, &b_ReconFailedLambdaPMom);
   fChain->SetBranchAddress("ReconFailedLambdaPMom_x", &ReconFailedLambdaPMom_x, &b_ReconFailedLambdaPMom_x);
   fChain->SetBranchAddress("ReconFailedLambdaPMom_y", &ReconFailedLambdaPMom_y, &b_ReconFailedLambdaPMom_y);
   fChain->SetBranchAddress("ReconFailedLambdaPMom_z", &ReconFailedLambdaPMom_z, &b_ReconFailedLambdaPMom_z);
   fChain->SetBranchAddress("ReconFailedLambdaPiMom", &ReconFailedLambdaPiMom, &b_ReconFailedLambdaPiMom);
   fChain->SetBranchAddress("ReconFailedLambdaPiMom_x", &ReconFailedLambdaPiMom_x, &b_ReconFailedLambdaPiMom_x);
   fChain->SetBranchAddress("ReconFailedLambdaPiMom_y", &ReconFailedLambdaPiMom_y, &b_ReconFailedLambdaPiMom_y);
   fChain->SetBranchAddress("ReconFailedLambdaPiMom_z", &ReconFailedLambdaPiMom_z, &b_ReconFailedLambdaPiMom_z);
   fChain->SetBranchAddress("ReconFailedLambdavtxCloseDist", &ReconFailedLambdavtxCloseDist, &b_ReconFailedLambdavtxCloseDist);
   fChain->SetBranchAddress("ncombiPiPair", &ncombiPiPair, &b_ncombiPiPair);
   fChain->SetBranchAddress("PiPairPipId", &PiPairPipId, &b_PiPairPipId);
   fChain->SetBranchAddress("PiPairPimId", &PiPairPimId, &b_PiPairPimId);
   fChain->SetBranchAddress("PiPairPipMom", &PiPairPipMom, &b_PiPairPipMom);
   fChain->SetBranchAddress("PiPairPipMom_x", &PiPairPipMom_x, &b_PiPairPipMom_x);
   fChain->SetBranchAddress("PiPairPipMom_y", &PiPairPipMom_y, &b_PiPairPipMom_y);
   fChain->SetBranchAddress("PiPairPipMom_z", &PiPairPipMom_z, &b_PiPairPipMom_z);
   fChain->SetBranchAddress("PiPairPimMom", &PiPairPimMom, &b_PiPairPimMom);
   fChain->SetBranchAddress("PiPairPimMom_x", &PiPairPimMom_x, &b_PiPairPimMom_x);
   fChain->SetBranchAddress("PiPairPimMom_y", &PiPairPimMom_y, &b_PiPairPimMom_y);
   fChain->SetBranchAddress("PiPairPimMom_z", &PiPairPimMom_z, &b_PiPairPimMom_z);
   fChain->SetBranchAddress("PiPairMom", &PiPairMom, &b_PiPairMom);
   fChain->SetBranchAddress("PiPairMom_x", &PiPairMom_x, &b_PiPairMom_x);
   fChain->SetBranchAddress("PiPairMom_y", &PiPairMom_y, &b_PiPairMom_y);
   fChain->SetBranchAddress("PiPairMom_z", &PiPairMom_z, &b_PiPairMom_z);
   fChain->SetBranchAddress("PiPairReconLambdaMass", &PiPairReconLambdaMass, &b_PiPairReconLambdaMass);
   fChain->SetBranchAddress("PiPairReconMass", &PiPairReconMass, &b_PiPairReconMass);
   fChain->SetBranchAddress("PiPairCloseDist", &PiPairCloseDist, &b_PiPairCloseDist);
   fChain->SetBranchAddress("ntK18", &ntK18, &b_ntK18);
   fChain->SetBranchAddress("pK18", &pK18, &b_pK18);
   fChain->SetBranchAddress("chisqrK18", &chisqrK18, &b_chisqrK18);
   fChain->SetBranchAddress("xtgtK18", &xtgtK18, &b_xtgtK18);
   fChain->SetBranchAddress("ytgtK18", &ytgtK18, &b_ytgtK18);
   fChain->SetBranchAddress("utgtK18", &utgtK18, &b_utgtK18);
   fChain->SetBranchAddress("vtgtK18", &vtgtK18, &b_vtgtK18);
   fChain->SetBranchAddress("isgoodTPCK18", &isgoodTPCK18, &b_isgoodTPCK18);
   fChain->SetBranchAddress("chisqrTPCK18", &chisqrTPCK18, &b_chisqrTPCK18);
   fChain->SetBranchAddress("qTPCK18", &qTPCK18, &b_qTPCK18);
   fChain->SetBranchAddress("pTPCK18", &pTPCK18, &b_pTPCK18);
   fChain->SetBranchAddress("xtgtTPCK18", &xtgtTPCK18, &b_xtgtTPCK18);
   fChain->SetBranchAddress("ytgtTPCK18", &ytgtTPCK18, &b_ytgtTPCK18);
   fChain->SetBranchAddress("utgtTPCK18", &utgtTPCK18, &b_utgtTPCK18);
   fChain->SetBranchAddress("vtgtTPCK18", &vtgtTPCK18, &b_vtgtTPCK18);
   fChain->SetBranchAddress("thetaTPCK18", &thetaTPCK18, &b_thetaTPCK18);
   fChain->SetBranchAddress("lhtofTPCK18", &lhtofTPCK18, &b_lhtofTPCK18);
   fChain->SetBranchAddress("xhtofTPCK18", &xhtofTPCK18, &b_xhtofTPCK18);
   fChain->SetBranchAddress("yhtofTPCK18", &yhtofTPCK18, &b_yhtofTPCK18);
   fChain->SetBranchAddress("lvpTPCK18", &lvpTPCK18, &b_lvpTPCK18);
   fChain->SetBranchAddress("xvpTPCK18", &xvpTPCK18, &b_xvpTPCK18);
   fChain->SetBranchAddress("yvpTPCK18", &yvpTPCK18, &b_yvpTPCK18);
   fChain->SetBranchAddress("ntKurama", &ntKurama, &b_ntKurama);
   fChain->SetBranchAddress("pKurama", &pKurama, &b_pKurama);
   fChain->SetBranchAddress("qKurama", &qKurama, &b_qKurama);
   fChain->SetBranchAddress("chisqrKurama", &chisqrKurama, &b_chisqrKurama);
   fChain->SetBranchAddress("m2Kurama", &m2Kurama, &b_m2Kurama);
   fChain->SetBranchAddress("xtgtKurama", &xtgtKurama, &b_xtgtKurama);
   fChain->SetBranchAddress("ytgtKurama", &ytgtKurama, &b_ytgtKurama);
   fChain->SetBranchAddress("utgtKurama", &utgtKurama, &b_utgtKurama);
   fChain->SetBranchAddress("vtgtKurama", &vtgtKurama, &b_vtgtKurama);
   fChain->SetBranchAddress("thetaKurama", &thetaKurama, &b_thetaKurama);
   fChain->SetBranchAddress("pathwcKurama", &pathwcKurama, &b_pathwcKurama);
   fChain->SetBranchAddress("xin", &xin, &b_xin);
   fChain->SetBranchAddress("yin", &yin, &b_yin);
   fChain->SetBranchAddress("zin", &zin, &b_zin);
   fChain->SetBranchAddress("pxin", &pxin, &b_pxin);
   fChain->SetBranchAddress("pyin", &pyin, &b_pyin);
   fChain->SetBranchAddress("pzin", &pzin, &b_pzin);
   fChain->SetBranchAddress("isgoodTPCKurama", &isgoodTPCKurama, &b_isgoodTPCKurama);
   fChain->SetBranchAddress("kflagTPCKurama", &kflagTPCKurama, &b_kflagTPCKurama);
   fChain->SetBranchAddress("chisqrTPCKurama", &chisqrTPCKurama, &b_chisqrTPCKurama);
   fChain->SetBranchAddress("pTPCKurama", &pTPCKurama, &b_pTPCKurama);
   fChain->SetBranchAddress("qTPCKurama", &qTPCKurama, &b_qTPCKurama);
   fChain->SetBranchAddress("m2TPCKurama", &m2TPCKurama, &b_m2TPCKurama);
   fChain->SetBranchAddress("xtgtTPCKurama", &xtgtTPCKurama, &b_xtgtTPCKurama);
   fChain->SetBranchAddress("ytgtTPCKurama", &ytgtTPCKurama, &b_ytgtTPCKurama);
   fChain->SetBranchAddress("utgtTPCKurama", &utgtTPCKurama, &b_utgtTPCKurama);
   fChain->SetBranchAddress("vtgtTPCKurama", &vtgtTPCKurama, &b_vtgtTPCKurama);
   fChain->SetBranchAddress("thetaTPCKurama", &thetaTPCKurama, &b_thetaTPCKurama);
   fChain->SetBranchAddress("pathTPCKurama", &pathTPCKurama, &b_pathTPCKurama);
   fChain->SetBranchAddress("lhtofTPCKurama", &lhtofTPCKurama, &b_lhtofTPCKurama);
   fChain->SetBranchAddress("xhtofTPCKurama", &xhtofTPCKurama, &b_xhtofTPCKurama);
   fChain->SetBranchAddress("yhtofTPCKurama", &yhtofTPCKurama, &b_yhtofTPCKurama);
   fChain->SetBranchAddress("lvpTPCKurama", &lvpTPCKurama, &b_lvpTPCKurama);
   fChain->SetBranchAddress("xvpTPCKurama", &xvpTPCKurama, &b_xvpTPCKurama);
   fChain->SetBranchAddress("yvpTPCKurama", &yvpTPCKurama, &b_yvpTPCKurama);
   fChain->SetBranchAddress("isgoodTPC", &isgoodTPC, &b_isgoodTPC);
   fChain->SetBranchAddress("insideTPC", &insideTPC, &b_insideTPC);
   fChain->SetBranchAddress("vtxTPC", &vtxTPC, &b_vtxTPC);
   fChain->SetBranchAddress("vtyTPC", &vtyTPC, &b_vtyTPC);
   fChain->SetBranchAddress("vtzTPC", &vtzTPC, &b_vtzTPC);
   fChain->SetBranchAddress("closeDistTPC", &closeDistTPC, &b_closeDistTPC);
   fChain->SetBranchAddress("MissMassTPC", &MissMassTPC, &b_MissMassTPC);
   fChain->SetBranchAddress("MissMassCorrTPC", &MissMassCorrTPC, &b_MissMassCorrTPC);
   fChain->SetBranchAddress("MissMassCorrDETPC", &MissMassCorrDETPC, &b_MissMassCorrDETPC);
   fChain->SetBranchAddress("pOrgTPC", &pOrgTPC, &b_pOrgTPC);
   fChain->SetBranchAddress("pCalcTPC", &pCalcTPC, &b_pCalcTPC);
   fChain->SetBranchAddress("pCorrTPC", &pCorrTPC, &b_pCorrTPC);
   fChain->SetBranchAddress("pCorrDETPC", &pCorrDETPC, &b_pCorrDETPC);
   fChain->SetBranchAddress("xbTPC", &xbTPC, &b_xbTPC);
   fChain->SetBranchAddress("ybTPC", &ybTPC, &b_ybTPC);
   fChain->SetBranchAddress("ubTPC", &ubTPC, &b_ubTPC);
   fChain->SetBranchAddress("vbTPC", &vbTPC, &b_vbTPC);
   fChain->SetBranchAddress("xsTPC", &xsTPC, &b_xsTPC);
   fChain->SetBranchAddress("ysTPC", &ysTPC, &b_ysTPC);
   fChain->SetBranchAddress("usTPC", &usTPC, &b_usTPC);
   fChain->SetBranchAddress("vsTPC", &vsTPC, &b_vsTPC);
   fChain->SetBranchAddress("nKK", &nKK, &b_nKK);
   fChain->SetBranchAddress("Kflag", &Kflag, &b_Kflag);
   fChain->SetBranchAddress("MissMass", &MissMass, &b_MissMass);
   fChain->SetBranchAddress("MissMassCorr", &MissMassCorr, &b_MissMassCorr);
   fChain->SetBranchAddress("MissMassCorrDE", &MissMassCorrDE, &b_MissMassCorrDE);
   fChain->SetBranchAddress("vtx", &vtx, &b_vtx);
   fChain->SetBranchAddress("vty", &vty, &b_vty);
   fChain->SetBranchAddress("vtz", &vtz, &b_vtz);
   fChain->SetBranchAddress("pOrg", &pOrg, &b_pOrg);
   fChain->SetBranchAddress("pCalc", &pCalc, &b_pCalc);
   fChain->SetBranchAddress("pCorr", &pCorr, &b_pCorr);
   fChain->SetBranchAddress("pCorrDE", &pCorrDE, &b_pCorrDE);
   fChain->SetBranchAddress("xb", &xb, &b_xb);
   fChain->SetBranchAddress("yb", &yb, &b_yb);
   fChain->SetBranchAddress("ub", &ub, &b_ub);
   fChain->SetBranchAddress("vb", &vb, &b_vb);
   fChain->SetBranchAddress("xs", &xs, &b_xs);
   fChain->SetBranchAddress("ys", &ys, &b_ys);
   fChain->SetBranchAddress("us", &us, &b_us);
   fChain->SetBranchAddress("vs", &vs, &b_vs);
   fChain->SetBranchAddress("KmMom_x", &KmMom_x, &b_KmMom_x);
   fChain->SetBranchAddress("KmMom_y", &KmMom_y, &b_KmMom_y);
   fChain->SetBranchAddress("KmMom_z", &KmMom_z, &b_KmMom_z);
   fChain->SetBranchAddress("KpMom_x", &KpMom_x, &b_KpMom_x);
   fChain->SetBranchAddress("KpMom_y", &KpMom_y, &b_KpMom_y);
   fChain->SetBranchAddress("KpMom_z", &KpMom_z, &b_KpMom_z);
   fChain->SetBranchAddress("BE", &BE, &b_BE);
   fChain->SetBranchAddress("BETPC", &BETPC, &b_BETPC);
   fChain->SetBranchAddress("BE_LL", &BE_LL, &b_BE_LL);
   fChain->SetBranchAddress("BETPC_LL", &BETPC_LL, &b_BETPC_LL);
   fChain->SetBranchAddress("GFKKXiProductionVtx_x", &GFKKXiProductionVtx_x, &b_GFKKXiProductionVtx_x);
   fChain->SetBranchAddress("GFKKXiProductionVtx_y", &GFKKXiProductionVtx_y, &b_GFKKXiProductionVtx_y);
   fChain->SetBranchAddress("GFKKXiProductionVtx_z", &GFKKXiProductionVtx_z, &b_GFKKXiProductionVtx_z);
   fChain->SetBranchAddress("GFKKLLProductionVtx_x", &GFKKLLProductionVtx_x, &b_GFKKLLProductionVtx_x);
   fChain->SetBranchAddress("GFKKLLProductionVtx_y", &GFKKLLProductionVtx_y, &b_GFKKLLProductionVtx_y);
   fChain->SetBranchAddress("GFKKLLProductionVtx_z", &GFKKLLProductionVtx_z, &b_GFKKLLProductionVtx_z);
   fChain->SetBranchAddress("GFKKLProductionVtx_x1", &GFKKLProductionVtx_x1, &b_GFKKLProductionVtx_x1);
   fChain->SetBranchAddress("GFKKLProductionVtx_y1", &GFKKLProductionVtx_y1, &b_GFKKLProductionVtx_y1);
   fChain->SetBranchAddress("GFKKLProductionVtx_z1", &GFKKLProductionVtx_z1, &b_GFKKLProductionVtx_z1);
   fChain->SetBranchAddress("GFKKLProductionVtx_x2", &GFKKLProductionVtx_x2, &b_GFKKLProductionVtx_x2);
   fChain->SetBranchAddress("GFKKLProductionVtx_y2", &GFKKLProductionVtx_y2, &b_GFKKLProductionVtx_y2);
   fChain->SetBranchAddress("GFKKLProductionVtx_z2", &GFKKLProductionVtx_z2, &b_GFKKLProductionVtx_z2);
   fChain->SetBranchAddress("GFKKLProductionVtx_x", &GFKKLProductionVtx_x, &b_GFKKLProductionVtx_x);
   fChain->SetBranchAddress("GFKKLProductionVtx_y", &GFKKLProductionVtx_y, &b_GFKKLProductionVtx_y);
   fChain->SetBranchAddress("GFKKLProductionVtx_z", &GFKKLProductionVtx_z, &b_GFKKLProductionVtx_z);
   fChain->SetBranchAddress("Xiflag", &Xiflag, &b_Xiflag);
   fChain->SetBranchAddress("XiPflag", &XiPflag, &b_XiPflag);
   fChain->SetBranchAddress("XiMass", &XiMass, &b_XiMass);
   fChain->SetBranchAddress("XiDecayVtx_x", &XiDecayVtx_x, &b_XiDecayVtx_x);
   fChain->SetBranchAddress("XiDecayVtx_y", &XiDecayVtx_y, &b_XiDecayVtx_y);
   fChain->SetBranchAddress("XiDecayVtx_z", &XiDecayVtx_z, &b_XiDecayVtx_z);
   fChain->SetBranchAddress("XiMom", &XiMom, &b_XiMom);
   fChain->SetBranchAddress("XiMom_x", &XiMom_x, &b_XiMom_x);
   fChain->SetBranchAddress("XiMom_y", &XiMom_y, &b_XiMom_y);
   fChain->SetBranchAddress("XiMom_z", &XiMom_z, &b_XiMom_z);
   fChain->SetBranchAddress("XiVtxCloseDist", &XiVtxCloseDist, &b_XiVtxCloseDist);
   fChain->SetBranchAddress("LambdaMass", &LambdaMass, &b_LambdaMass);
   fChain->SetBranchAddress("LambdaDecayVtx_x", &LambdaDecayVtx_x, &b_LambdaDecayVtx_x);
   fChain->SetBranchAddress("LambdaDecayVtx_y", &LambdaDecayVtx_y, &b_LambdaDecayVtx_y);
   fChain->SetBranchAddress("LambdaDecayVtx_z", &LambdaDecayVtx_z, &b_LambdaDecayVtx_z);
   fChain->SetBranchAddress("LambdaMom", &LambdaMom, &b_LambdaMom);
   fChain->SetBranchAddress("LambdaMom_x", &LambdaMom_x, &b_LambdaMom_x);
   fChain->SetBranchAddress("LambdaMom_y", &LambdaMom_y, &b_LambdaMom_y);
   fChain->SetBranchAddress("LambdaMom_z", &LambdaMom_z, &b_LambdaMom_z);
   fChain->SetBranchAddress("LambdaVtxCloseDist", &LambdaVtxCloseDist, &b_LambdaVtxCloseDist);
   fChain->SetBranchAddress("XiTarget_x", &XiTarget_x, &b_XiTarget_x);
   fChain->SetBranchAddress("XiTarget_y", &XiTarget_y, &b_XiTarget_y);
   fChain->SetBranchAddress("XiTarget_z", &XiTarget_z, &b_XiTarget_z);
   fChain->SetBranchAddress("XiTargetMom", &XiTargetMom, &b_XiTargetMom);
   fChain->SetBranchAddress("XiTargetMom_x", &XiTargetMom_x, &b_XiTargetMom_x);
   fChain->SetBranchAddress("XiTargetMom_y", &XiTargetMom_y, &b_XiTargetMom_y);
   fChain->SetBranchAddress("XiTargetMom_z", &XiTargetMom_z, &b_XiTargetMom_z);
   fChain->SetBranchAddress("XiTargetCloseDist", &XiTargetCloseDist, &b_XiTargetCloseDist);
   fChain->SetBranchAddress("GFXiMass", &GFXiMass, &b_GFXiMass);
   fChain->SetBranchAddress("GFXiDecayVtx_x", &GFXiDecayVtx_x, &b_GFXiDecayVtx_x);
   fChain->SetBranchAddress("GFXiDecayVtx_y", &GFXiDecayVtx_y, &b_GFXiDecayVtx_y);
   fChain->SetBranchAddress("GFXiDecayVtx_z", &GFXiDecayVtx_z, &b_GFXiDecayVtx_z);
   fChain->SetBranchAddress("GFXiMom", &GFXiMom, &b_GFXiMom);
   fChain->SetBranchAddress("GFXiMom_x", &GFXiMom_x, &b_GFXiMom_x);
   fChain->SetBranchAddress("GFXiMom_y", &GFXiMom_y, &b_GFXiMom_y);
   fChain->SetBranchAddress("GFXiMom_z", &GFXiMom_z, &b_GFXiMom_z);
   fChain->SetBranchAddress("GFXiVtxCloseDist", &GFXiVtxCloseDist, &b_GFXiVtxCloseDist);
   fChain->SetBranchAddress("GFXiKKVtx_x", &GFXiKKVtx_x, &b_GFXiKKVtx_x);
   fChain->SetBranchAddress("GFXiKKVtx_y", &GFXiKKVtx_y, &b_GFXiKKVtx_y);
   fChain->SetBranchAddress("GFXiKKVtx_z", &GFXiKKVtx_z, &b_GFXiKKVtx_z);
   fChain->SetBranchAddress("GFXiKKVtxMom", &GFXiKKVtxMom, &b_GFXiKKVtxMom);
   fChain->SetBranchAddress("GFXiKKVtxMom_x", &GFXiKKVtxMom_x, &b_GFXiKKVtxMom_x);
   fChain->SetBranchAddress("GFXiKKVtxMom_y", &GFXiKKVtxMom_y, &b_GFXiKKVtxMom_y);
   fChain->SetBranchAddress("GFXiKKVtxMom_z", &GFXiKKVtxMom_z, &b_GFXiKKVtxMom_z);
   fChain->SetBranchAddress("GFXiKKVtxCloseDist", &GFXiKKVtxCloseDist, &b_GFXiKKVtxCloseDist);
   fChain->SetBranchAddress("GFXiProductionVtx_x", &GFXiProductionVtx_x, &b_GFXiProductionVtx_x);
   fChain->SetBranchAddress("GFXiProductionVtx_y", &GFXiProductionVtx_y, &b_GFXiProductionVtx_y);
   fChain->SetBranchAddress("GFXiProductionVtx_z", &GFXiProductionVtx_z, &b_GFXiProductionVtx_z);
   fChain->SetBranchAddress("GFXiProductionVtxMom", &GFXiProductionVtxMom, &b_GFXiProductionVtxMom);
   fChain->SetBranchAddress("GFXiProductionVtxMom_x", &GFXiProductionVtxMom_x, &b_GFXiProductionVtxMom_x);
   fChain->SetBranchAddress("GFXiProductionVtxMom_y", &GFXiProductionVtxMom_y, &b_GFXiProductionVtxMom_y);
   fChain->SetBranchAddress("GFXiProductionVtxMom_z", &GFXiProductionVtxMom_z, &b_GFXiProductionVtxMom_z);
   fChain->SetBranchAddress("GFXiProductionVtxCloseDist", &GFXiProductionVtxCloseDist, &b_GFXiProductionVtxCloseDist);
   fChain->SetBranchAddress("GFXiTrackLen", &GFXiTrackLen, &b_GFXiTrackLen);
   fChain->SetBranchAddress("GFXiTof", &GFXiTof, &b_GFXiTof);
   fChain->SetBranchAddress("GFXiMomLoss", &GFXiMomLoss, &b_GFXiMomLoss);
   fChain->SetBranchAddress("GFXiExcitation", &GFXiExcitation, &b_GFXiExcitation);
   fChain->SetBranchAddress("GFXiTarget_x", &GFXiTarget_x, &b_GFXiTarget_x);
   fChain->SetBranchAddress("GFXiTarget_y", &GFXiTarget_y, &b_GFXiTarget_y);
   fChain->SetBranchAddress("GFXiTarget_z", &GFXiTarget_z, &b_GFXiTarget_z);
   fChain->SetBranchAddress("GFXiTargetMom", &GFXiTargetMom, &b_GFXiTargetMom);
   fChain->SetBranchAddress("GFXiTargetMom_x", &GFXiTargetMom_x, &b_GFXiTargetMom_x);
   fChain->SetBranchAddress("GFXiTargetMom_y", &GFXiTargetMom_y, &b_GFXiTargetMom_y);
   fChain->SetBranchAddress("GFXiTargetMom_z", &GFXiTargetMom_z, &b_GFXiTargetMom_z);
   fChain->SetBranchAddress("GFXiTargetCloseDist", &GFXiTargetCloseDist, &b_GFXiTargetCloseDist);
   fChain->SetBranchAddress("GFXiTargetCenter_x", &GFXiTargetCenter_x, &b_GFXiTargetCenter_x);
   fChain->SetBranchAddress("GFXiTargetCenter_y", &GFXiTargetCenter_y, &b_GFXiTargetCenter_y);
   fChain->SetBranchAddress("GFXiTargetCenter_z", &GFXiTargetCenter_z, &b_GFXiTargetCenter_z);
   fChain->SetBranchAddress("GFXiTargetCenterMom", &GFXiTargetCenterMom, &b_GFXiTargetCenterMom);
   fChain->SetBranchAddress("GFXiTargetCenterMom_x", &GFXiTargetCenterMom_x, &b_GFXiTargetCenterMom_x);
   fChain->SetBranchAddress("GFXiTargetCenterMom_y", &GFXiTargetCenterMom_y, &b_GFXiTargetCenterMom_y);
   fChain->SetBranchAddress("GFXiTargetCenterMom_z", &GFXiTargetCenterMom_z, &b_GFXiTargetCenterMom_z);
   fChain->SetBranchAddress("GFXiTargetCenterCloseDist", &GFXiTargetCenterCloseDist, &b_GFXiTargetCenterCloseDist);
   fChain->SetBranchAddress("LLflag", &LLflag, &b_LLflag);
   fChain->SetBranchAddress("LambdaLambdaVtx_x", &LambdaLambdaVtx_x, &b_LambdaLambdaVtx_x);
   fChain->SetBranchAddress("LambdaLambdaVtx_y", &LambdaLambdaVtx_y, &b_LambdaLambdaVtx_y);
   fChain->SetBranchAddress("LambdaLambdaVtx_z", &LambdaLambdaVtx_z, &b_LambdaLambdaVtx_z);
   fChain->SetBranchAddress("LambdaLambdaCloseDist", &LambdaLambdaCloseDist, &b_LambdaLambdaCloseDist);
   fChain->SetBranchAddress("LambdaMass1", &LambdaMass1, &b_LambdaMass1);
   fChain->SetBranchAddress("LambdaDecayVtx_x1", &LambdaDecayVtx_x1, &b_LambdaDecayVtx_x1);
   fChain->SetBranchAddress("LambdaDecayVtx_y1", &LambdaDecayVtx_y1, &b_LambdaDecayVtx_y1);
   fChain->SetBranchAddress("LambdaDecayVtx_z1", &LambdaDecayVtx_z1, &b_LambdaDecayVtx_z1);
   fChain->SetBranchAddress("LambdaMom1", &LambdaMom1, &b_LambdaMom1);
   fChain->SetBranchAddress("LambdaMom_x1", &LambdaMom_x1, &b_LambdaMom_x1);
   fChain->SetBranchAddress("LambdaMom_y1", &LambdaMom_y1, &b_LambdaMom_y1);
   fChain->SetBranchAddress("LambdaMom_z1", &LambdaMom_z1, &b_LambdaMom_z1);
   fChain->SetBranchAddress("LambdaVtxCloseDist1", &LambdaVtxCloseDist1, &b_LambdaVtxCloseDist1);
   fChain->SetBranchAddress("LambdaTargetCloseDist1", &LambdaTargetCloseDist1, &b_LambdaTargetCloseDist1);
   fChain->SetBranchAddress("LambdaTargetCloseVtx_x1", &LambdaTargetCloseVtx_x1, &b_LambdaTargetCloseVtx_x1);
   fChain->SetBranchAddress("LambdaTargetCloseVtx_y1", &LambdaTargetCloseVtx_y1, &b_LambdaTargetCloseVtx_y1);
   fChain->SetBranchAddress("LambdaTargetCloseVtx_z1", &LambdaTargetCloseVtx_z1, &b_LambdaTargetCloseVtx_z1);
   fChain->SetBranchAddress("LambdaMass2", &LambdaMass2, &b_LambdaMass2);
   fChain->SetBranchAddress("LambdaDecayVtx_x2", &LambdaDecayVtx_x2, &b_LambdaDecayVtx_x2);
   fChain->SetBranchAddress("LambdaDecayVtx_y2", &LambdaDecayVtx_y2, &b_LambdaDecayVtx_y2);
   fChain->SetBranchAddress("LambdaDecayVtx_z2", &LambdaDecayVtx_z2, &b_LambdaDecayVtx_z2);
   fChain->SetBranchAddress("LambdaMom2", &LambdaMom2, &b_LambdaMom2);
   fChain->SetBranchAddress("LambdaMom_x2", &LambdaMom_x2, &b_LambdaMom_x2);
   fChain->SetBranchAddress("LambdaMom_y2", &LambdaMom_y2, &b_LambdaMom_y2);
   fChain->SetBranchAddress("LambdaMom_z2", &LambdaMom_z2, &b_LambdaMom_z2);
   fChain->SetBranchAddress("LambdaVtxCloseDist2", &LambdaVtxCloseDist2, &b_LambdaVtxCloseDist2);
   fChain->SetBranchAddress("LambdaTargetCloseDist2", &LambdaTargetCloseDist2, &b_LambdaTargetCloseDist2);
   fChain->SetBranchAddress("LambdaTargetCloseVtx_x2", &LambdaTargetCloseVtx_x2, &b_LambdaTargetCloseVtx_x2);
   fChain->SetBranchAddress("LambdaTargetCloseVtx_y2", &LambdaTargetCloseVtx_y2, &b_LambdaTargetCloseVtx_y2);
   fChain->SetBranchAddress("LambdaTargetCloseVtx_z2", &LambdaTargetCloseVtx_z2, &b_LambdaTargetCloseVtx_z2);
   fChain->SetBranchAddress("GFLLExcitation", &GFLLExcitation, &b_GFLLExcitation);
   fChain->SetBranchAddress("GFLambdaMass1", &GFLambdaMass1, &b_GFLambdaMass1);
   fChain->SetBranchAddress("GFLambdaDecayVtx_x1", &GFLambdaDecayVtx_x1, &b_GFLambdaDecayVtx_x1);
   fChain->SetBranchAddress("GFLambdaDecayVtx_y1", &GFLambdaDecayVtx_y1, &b_GFLambdaDecayVtx_y1);
   fChain->SetBranchAddress("GFLambdaDecayVtx_z1", &GFLambdaDecayVtx_z1, &b_GFLambdaDecayVtx_z1);
   fChain->SetBranchAddress("GFLambdaMom1", &GFLambdaMom1, &b_GFLambdaMom1);
   fChain->SetBranchAddress("GFLambdaMom_x1", &GFLambdaMom_x1, &b_GFLambdaMom_x1);
   fChain->SetBranchAddress("GFLambdaMom_y1", &GFLambdaMom_y1, &b_GFLambdaMom_y1);
   fChain->SetBranchAddress("GFLambdaMom_z1", &GFLambdaMom_z1, &b_GFLambdaMom_z1);
   fChain->SetBranchAddress("GFLambdaVtxCloseDist1", &GFLambdaVtxCloseDist1, &b_GFLambdaVtxCloseDist1);
   fChain->SetBranchAddress("GFLambdaMass2", &GFLambdaMass2, &b_GFLambdaMass2);
   fChain->SetBranchAddress("GFLambdaDecayVtx_x2", &GFLambdaDecayVtx_x2, &b_GFLambdaDecayVtx_x2);
   fChain->SetBranchAddress("GFLambdaDecayVtx_y2", &GFLambdaDecayVtx_y2, &b_GFLambdaDecayVtx_y2);
   fChain->SetBranchAddress("GFLambdaDecayVtx_z2", &GFLambdaDecayVtx_z2, &b_GFLambdaDecayVtx_z2);
   fChain->SetBranchAddress("GFLambdaMom2", &GFLambdaMom2, &b_GFLambdaMom2);
   fChain->SetBranchAddress("GFLambdaMom_x2", &GFLambdaMom_x2, &b_GFLambdaMom_x2);
   fChain->SetBranchAddress("GFLambdaMom_y2", &GFLambdaMom_y2, &b_GFLambdaMom_y2);
   fChain->SetBranchAddress("GFLambdaMom_z2", &GFLambdaMom_z2, &b_GFLambdaMom_z2);
   fChain->SetBranchAddress("GFLambdaVtxCloseDist2", &GFLambdaVtxCloseDist2, &b_GFLambdaVtxCloseDist2);
   fChain->SetBranchAddress("GFLambdaMass1_Alt", &GFLambdaMass1_Alt, &b_GFLambdaMass1_Alt);
   fChain->SetBranchAddress("GFLambdaDecayVtx_x1_Alt", &GFLambdaDecayVtx_x1_Alt, &b_GFLambdaDecayVtx_x1_Alt);
   fChain->SetBranchAddress("GFLambdaDecayVtx_y1_Alt", &GFLambdaDecayVtx_y1_Alt, &b_GFLambdaDecayVtx_y1_Alt);
   fChain->SetBranchAddress("GFLambdaDecayVtx_z1_Alt", &GFLambdaDecayVtx_z1_Alt, &b_GFLambdaDecayVtx_z1_Alt);
   fChain->SetBranchAddress("GFLambdaMom1_Alt", &GFLambdaMom1_Alt, &b_GFLambdaMom1_Alt);
   fChain->SetBranchAddress("GFLambdaMom_x1_Alt", &GFLambdaMom_x1_Alt, &b_GFLambdaMom_x1_Alt);
   fChain->SetBranchAddress("GFLambdaMom_y1_Alt", &GFLambdaMom_y1_Alt, &b_GFLambdaMom_y1_Alt);
   fChain->SetBranchAddress("GFLambdaMom_z1_Alt", &GFLambdaMom_z1_Alt, &b_GFLambdaMom_z1_Alt);
   fChain->SetBranchAddress("GFLambdaVtxCloseDist1_Alt", &GFLambdaVtxCloseDist1_Alt, &b_GFLambdaVtxCloseDist1_Alt);
   fChain->SetBranchAddress("GFLambdaTargetCloseDist1_Alt", &GFLambdaTargetCloseDist1_Alt, &b_GFLambdaTargetCloseDist1_Alt);
   fChain->SetBranchAddress("GFLambdaTarget_x1_Alt", &GFLambdaTarget_x1_Alt, &b_GFLambdaTarget_x1_Alt);
   fChain->SetBranchAddress("GFLambdaTarget_y1_Alt", &GFLambdaTarget_y1_Alt, &b_GFLambdaTarget_y1_Alt);
   fChain->SetBranchAddress("GFLambdaTarget_z1_Alt", &GFLambdaTarget_z1_Alt, &b_GFLambdaTarget_z1_Alt);
   fChain->SetBranchAddress("GFLambdaTargetCenterCloseDist1_Alt", &GFLambdaTargetCenterCloseDist1_Alt, &b_GFLambdaTargetCenterCloseDist1_Alt);
   fChain->SetBranchAddress("GFLambdaTargetCenter_x1_Alt", &GFLambdaTargetCenter_x1_Alt, &b_GFLambdaTargetCenter_x1_Alt);
   fChain->SetBranchAddress("GFLambdaTargetCenter_y1_Alt", &GFLambdaTargetCenter_y1_Alt, &b_GFLambdaTargetCenter_y1_Alt);
   fChain->SetBranchAddress("GFLambdaTargetCenter_z1_Alt", &GFLambdaTargetCenter_z1_Alt, &b_GFLambdaTargetCenter_z1_Alt);
   fChain->SetBranchAddress("GFLambdaMass2_Alt", &GFLambdaMass2_Alt, &b_GFLambdaMass2_Alt);
   fChain->SetBranchAddress("GFLambdaDecayVtx_x2_Alt", &GFLambdaDecayVtx_x2_Alt, &b_GFLambdaDecayVtx_x2_Alt);
   fChain->SetBranchAddress("GFLambdaDecayVtx_y2_Alt", &GFLambdaDecayVtx_y2_Alt, &b_GFLambdaDecayVtx_y2_Alt);
   fChain->SetBranchAddress("GFLambdaDecayVtx_z2_Alt", &GFLambdaDecayVtx_z2_Alt, &b_GFLambdaDecayVtx_z2_Alt);
   fChain->SetBranchAddress("GFLambdaMom2_Alt", &GFLambdaMom2_Alt, &b_GFLambdaMom2_Alt);
   fChain->SetBranchAddress("GFLambdaMom_x2_Alt", &GFLambdaMom_x2_Alt, &b_GFLambdaMom_x2_Alt);
   fChain->SetBranchAddress("GFLambdaMom_y2_Alt", &GFLambdaMom_y2_Alt, &b_GFLambdaMom_y2_Alt);
   fChain->SetBranchAddress("GFLambdaMom_z2_Alt", &GFLambdaMom_z2_Alt, &b_GFLambdaMom_z2_Alt);
   fChain->SetBranchAddress("GFLambdaVtxCloseDist2_Alt", &GFLambdaVtxCloseDist2_Alt, &b_GFLambdaVtxCloseDist2_Alt);
   fChain->SetBranchAddress("GFLambdaTargetCloseDist2_Alt", &GFLambdaTargetCloseDist2_Alt, &b_GFLambdaTargetCloseDist2_Alt);
   fChain->SetBranchAddress("GFLambdaTarget_x2_Alt", &GFLambdaTarget_x2_Alt, &b_GFLambdaTarget_x2_Alt);
   fChain->SetBranchAddress("GFLambdaTarget_y2_Alt", &GFLambdaTarget_y2_Alt, &b_GFLambdaTarget_y2_Alt);
   fChain->SetBranchAddress("GFLambdaTarget_z2_Alt", &GFLambdaTarget_z2_Alt, &b_GFLambdaTarget_z2_Alt);
   fChain->SetBranchAddress("GFLambdaTargetCenterCloseDist2_Alt", &GFLambdaTargetCenterCloseDist2_Alt, &b_GFLambdaTargetCenterCloseDist2_Alt);
   fChain->SetBranchAddress("GFLambdaTargetCenter_x2_Alt", &GFLambdaTargetCenter_x2_Alt, &b_GFLambdaTargetCenter_x2_Alt);
   fChain->SetBranchAddress("GFLambdaTargetCenter_y2_Alt", &GFLambdaTargetCenter_y2_Alt, &b_GFLambdaTargetCenter_y2_Alt);
   fChain->SetBranchAddress("GFLambdaTargetCenter_z2_Alt", &GFLambdaTargetCenter_z2_Alt, &b_GFLambdaTargetCenter_z2_Alt);
   fChain->SetBranchAddress("GFLambdaLambdaVtx_x", &GFLambdaLambdaVtx_x, &b_GFLambdaLambdaVtx_x);
   fChain->SetBranchAddress("GFLambdaLambdaVtx_y", &GFLambdaLambdaVtx_y, &b_GFLambdaLambdaVtx_y);
   fChain->SetBranchAddress("GFLambdaLambdaVtx_z", &GFLambdaLambdaVtx_z, &b_GFLambdaLambdaVtx_z);
   fChain->SetBranchAddress("GFLambdaLambdaCloseDist", &GFLambdaLambdaCloseDist, &b_GFLambdaLambdaCloseDist);
   fChain->SetBranchAddress("GFLambdaTarget_x1", &GFLambdaTarget_x1, &b_GFLambdaTarget_x1);
   fChain->SetBranchAddress("GFLambdaTarget_y1", &GFLambdaTarget_y1, &b_GFLambdaTarget_y1);
   fChain->SetBranchAddress("GFLambdaTarget_z1", &GFLambdaTarget_z1, &b_GFLambdaTarget_z1);
   fChain->SetBranchAddress("GFLambdaTargetCloseDist1", &GFLambdaTargetCloseDist1, &b_GFLambdaTargetCloseDist1);
   fChain->SetBranchAddress("GFLambdaTarget_x2", &GFLambdaTarget_x2, &b_GFLambdaTarget_x2);
   fChain->SetBranchAddress("GFLambdaTarget_y2", &GFLambdaTarget_y2, &b_GFLambdaTarget_y2);
   fChain->SetBranchAddress("GFLambdaTarget_z2", &GFLambdaTarget_z2, &b_GFLambdaTarget_z2);
   fChain->SetBranchAddress("GFLambdaTargetCloseDist2", &GFLambdaTargetCloseDist2, &b_GFLambdaTargetCloseDist2);
   fChain->SetBranchAddress("GFLambdaTargetCenter_x1", &GFLambdaTargetCenter_x1, &b_GFLambdaTargetCenter_x1);
   fChain->SetBranchAddress("GFLambdaTargetCenter_y1", &GFLambdaTargetCenter_y1, &b_GFLambdaTargetCenter_y1);
   fChain->SetBranchAddress("GFLambdaTargetCenter_z1", &GFLambdaTargetCenter_z1, &b_GFLambdaTargetCenter_z1);
   fChain->SetBranchAddress("GFLambdaTargetCenterCloseDist1", &GFLambdaTargetCenterCloseDist1, &b_GFLambdaTargetCenterCloseDist1);
   fChain->SetBranchAddress("GFLambdaTargetCenter_x2", &GFLambdaTargetCenter_x2, &b_GFLambdaTargetCenter_x2);
   fChain->SetBranchAddress("GFLambdaTargetCenter_y2", &GFLambdaTargetCenter_y2, &b_GFLambdaTargetCenter_y2);
   fChain->SetBranchAddress("GFLambdaTargetCenter_z2", &GFLambdaTargetCenter_z2, &b_GFLambdaTargetCenter_z2);
   fChain->SetBranchAddress("GFLambdaTargetCenterCloseDist2", &GFLambdaTargetCenterCloseDist2, &b_GFLambdaTargetCenterCloseDist2);
   fChain->SetBranchAddress("GFLambdaProductionVtx_x1", &GFLambdaProductionVtx_x1, &b_GFLambdaProductionVtx_x1);
   fChain->SetBranchAddress("GFLambdaProductionVtx_y1", &GFLambdaProductionVtx_y1, &b_GFLambdaProductionVtx_y1);
   fChain->SetBranchAddress("GFLambdaProductionVtx_z1", &GFLambdaProductionVtx_z1, &b_GFLambdaProductionVtx_z1);
   fChain->SetBranchAddress("GFLambdaProductionVtxCloseDist1", &GFLambdaProductionVtxCloseDist1, &b_GFLambdaProductionVtxCloseDist1);
   fChain->SetBranchAddress("GFLambdaTrackLen1", &GFLambdaTrackLen1, &b_GFLambdaTrackLen1);
   fChain->SetBranchAddress("GFLambdaTof1", &GFLambdaTof1, &b_GFLambdaTof1);
   fChain->SetBranchAddress("GFLambdaProductionVtx_x2", &GFLambdaProductionVtx_x2, &b_GFLambdaProductionVtx_x2);
   fChain->SetBranchAddress("GFLambdaProductionVtx_y2", &GFLambdaProductionVtx_y2, &b_GFLambdaProductionVtx_y2);
   fChain->SetBranchAddress("GFLambdaProductionVtx_z2", &GFLambdaProductionVtx_z2, &b_GFLambdaProductionVtx_z2);
   fChain->SetBranchAddress("GFLambdaProductionVtxCloseDist2", &GFLambdaProductionVtxCloseDist2, &b_GFLambdaProductionVtxCloseDist2);
   fChain->SetBranchAddress("GFLambdaTrackLen2", &GFLambdaTrackLen2, &b_GFLambdaTrackLen2);
   fChain->SetBranchAddress("GFLambdaTof2", &GFLambdaTof2, &b_GFLambdaTof2);
   fChain->SetBranchAddress("KFLLExcitation", &KFLLExcitation, &b_KFLLExcitation);
   fChain->SetBranchAddress("KFKKLLProductionVtxChisqr", &KFKKLLProductionVtxChisqr, &b_KFKKLLProductionVtxChisqr);
   fChain->SetBranchAddress("KFKKLLProductionVtx_x", &KFKKLLProductionVtx_x, &b_KFKKLLProductionVtx_x);
   fChain->SetBranchAddress("KFKKLLProductionVtx_y", &KFKKLLProductionVtx_y, &b_KFKKLLProductionVtx_y);
   fChain->SetBranchAddress("KFKKLLProductionVtx_z", &KFKKLLProductionVtx_z, &b_KFKKLLProductionVtx_z);
   fChain->SetBranchAddress("KFKKLProductionVtx_x1", &KFKKLProductionVtx_x1, &b_KFKKLProductionVtx_x1);
   fChain->SetBranchAddress("KFKKLProductionVtx_y1", &KFKKLProductionVtx_y1, &b_KFKKLProductionVtx_y1);
   fChain->SetBranchAddress("KFKKLProductionVtx_z1", &KFKKLProductionVtx_z1, &b_KFKKLProductionVtx_z1);
   fChain->SetBranchAddress("KFKKLProductionVtx_x2", &KFKKLProductionVtx_x2, &b_KFKKLProductionVtx_x2);
   fChain->SetBranchAddress("KFKKLProductionVtx_y2", &KFKKLProductionVtx_y2, &b_KFKKLProductionVtx_y2);
   fChain->SetBranchAddress("KFKKLProductionVtx_z2", &KFKKLProductionVtx_z2, &b_KFKKLProductionVtx_z2);
   fChain->SetBranchAddress("KFLambdaLambdaVtx_x", &KFLambdaLambdaVtx_x, &b_KFLambdaLambdaVtx_x);
   fChain->SetBranchAddress("KFLambdaLambdaVtx_y", &KFLambdaLambdaVtx_y, &b_KFLambdaLambdaVtx_y);
   fChain->SetBranchAddress("KFLambdaLambdaVtx_z", &KFLambdaLambdaVtx_z, &b_KFLambdaLambdaVtx_z);
   fChain->SetBranchAddress("KFLambdaLambdaCloseDist", &KFLambdaLambdaCloseDist, &b_KFLambdaLambdaCloseDist);
   fChain->SetBranchAddress("KFLambdaProductionVtx_x1", &KFLambdaProductionVtx_x1, &b_KFLambdaProductionVtx_x1);
   fChain->SetBranchAddress("KFLambdaProductionVtx_y1", &KFLambdaProductionVtx_y1, &b_KFLambdaProductionVtx_y1);
   fChain->SetBranchAddress("KFLambdaProductionVtx_z1", &KFLambdaProductionVtx_z1, &b_KFLambdaProductionVtx_z1);
   fChain->SetBranchAddress("KFLambdaProductionVtxCloseDist1", &KFLambdaProductionVtxCloseDist1, &b_KFLambdaProductionVtxCloseDist1);
   fChain->SetBranchAddress("KFLambdaTrackLen1", &KFLambdaTrackLen1, &b_KFLambdaTrackLen1);
   fChain->SetBranchAddress("KFLambdaTof1", &KFLambdaTof1, &b_KFLambdaTof1);
   fChain->SetBranchAddress("KFLambdaProductionVtx_x2", &KFLambdaProductionVtx_x2, &b_KFLambdaProductionVtx_x2);
   fChain->SetBranchAddress("KFLambdaProductionVtx_y2", &KFLambdaProductionVtx_y2, &b_KFLambdaProductionVtx_y2);
   fChain->SetBranchAddress("KFLambdaProductionVtx_z2", &KFLambdaProductionVtx_z2, &b_KFLambdaProductionVtx_z2);
   fChain->SetBranchAddress("KFLambdaProductionVtxCloseDist2", &KFLambdaProductionVtxCloseDist2, &b_KFLambdaProductionVtxCloseDist2);
   fChain->SetBranchAddress("KFLambdaTrackLen2", &KFLambdaTrackLen2, &b_KFLambdaTrackLen2);
   fChain->SetBranchAddress("KFLambdaTof2", &KFLambdaTof2, &b_KFLambdaTof2);
   fChain->SetBranchAddress("KFLambdaPval1", &KFLambdaPval1, &b_KFLambdaPval1);
   fChain->SetBranchAddress("KFLambdaPull1", &KFLambdaPull1, &b_KFLambdaPull1);
   fChain->SetBranchAddress("KFLambdaChisqr1", &KFLambdaChisqr1, &b_KFLambdaChisqr1);
   fChain->SetBranchAddress("KFLambdaMom1", &KFLambdaMom1, &b_KFLambdaMom1);
   fChain->SetBranchAddress("KFLambdaMom_x1", &KFLambdaMom_x1, &b_KFLambdaMom_x1);
   fChain->SetBranchAddress("KFLambdaMom_y1", &KFLambdaMom_y1, &b_KFLambdaMom_y1);
   fChain->SetBranchAddress("KFLambdaMom_z1", &KFLambdaMom_z1, &b_KFLambdaMom_z1);
   fChain->SetBranchAddress("KFLambdaPval2", &KFLambdaPval2, &b_KFLambdaPval2);
   fChain->SetBranchAddress("KFLambdaPull2", &KFLambdaPull2, &b_KFLambdaPull2);
   fChain->SetBranchAddress("KFLambdaChisqr2", &KFLambdaChisqr2, &b_KFLambdaChisqr2);
   fChain->SetBranchAddress("KFLambdaMom2", &KFLambdaMom2, &b_KFLambdaMom2);
   fChain->SetBranchAddress("KFLambdaMom_x2", &KFLambdaMom_x2, &b_KFLambdaMom_x2);
   fChain->SetBranchAddress("KFLambdaMom_y2", &KFLambdaMom_y2, &b_KFLambdaMom_y2);
   fChain->SetBranchAddress("KFLambdaMom_z2", &KFLambdaMom_z2, &b_KFLambdaMom_z2);
   fChain->SetBranchAddress("Lflag", &Lflag, &b_Lflag);
   fChain->SetBranchAddress("LambdaTargetCloseDist", &LambdaTargetCloseDist, &b_LambdaTargetCloseDist);
   fChain->SetBranchAddress("LambdaTargetCloseVtx_x", &LambdaTargetCloseVtx_x, &b_LambdaTargetCloseVtx_x);
   fChain->SetBranchAddress("LambdaTargetCloseVtx_y", &LambdaTargetCloseVtx_y, &b_LambdaTargetCloseVtx_y);
   fChain->SetBranchAddress("LambdaTargetCloseVtx_z", &LambdaTargetCloseVtx_z, &b_LambdaTargetCloseVtx_z);
   fChain->SetBranchAddress("GFLambdaMass", &GFLambdaMass, &b_GFLambdaMass);
   fChain->SetBranchAddress("GFLambdaDecayVtx_x", &GFLambdaDecayVtx_x, &b_GFLambdaDecayVtx_x);
   fChain->SetBranchAddress("GFLambdaDecayVtx_y", &GFLambdaDecayVtx_y, &b_GFLambdaDecayVtx_y);
   fChain->SetBranchAddress("GFLambdaDecayVtx_z", &GFLambdaDecayVtx_z, &b_GFLambdaDecayVtx_z);
   fChain->SetBranchAddress("GFLambdaMom", &GFLambdaMom, &b_GFLambdaMom);
   fChain->SetBranchAddress("GFLambdaMom_x", &GFLambdaMom_x, &b_GFLambdaMom_x);
   fChain->SetBranchAddress("GFLambdaMom_y", &GFLambdaMom_y, &b_GFLambdaMom_y);
   fChain->SetBranchAddress("GFLambdaMom_z", &GFLambdaMom_z, &b_GFLambdaMom_z);
   fChain->SetBranchAddress("GFLambdaVtxCloseDist", &GFLambdaVtxCloseDist, &b_GFLambdaVtxCloseDist);
   fChain->SetBranchAddress("GFLambdaTargetCloseDist", &GFLambdaTargetCloseDist, &b_GFLambdaTargetCloseDist);
   fChain->SetBranchAddress("GFLambdaTarget_x", &GFLambdaTarget_x, &b_GFLambdaTarget_x);
   fChain->SetBranchAddress("GFLambdaTarget_y", &GFLambdaTarget_y, &b_GFLambdaTarget_y);
   fChain->SetBranchAddress("GFLambdaTarget_z", &GFLambdaTarget_z, &b_GFLambdaTarget_z);
   fChain->SetBranchAddress("GFLambdaTargetCenter_x", &GFLambdaTargetCenter_x, &b_GFLambdaTargetCenter_x);
   fChain->SetBranchAddress("GFLambdaTargetCenter_y", &GFLambdaTargetCenter_y, &b_GFLambdaTargetCenter_y);
   fChain->SetBranchAddress("GFLambdaTargetCenter_z", &GFLambdaTargetCenter_z, &b_GFLambdaTargetCenter_z);
   fChain->SetBranchAddress("GFLambdaTargetCenterCloseDist", &GFLambdaTargetCenterCloseDist, &b_GFLambdaTargetCenterCloseDist);
   fChain->SetBranchAddress("GFLambdaProductionVtx_x", &GFLambdaProductionVtx_x, &b_GFLambdaProductionVtx_x);
   fChain->SetBranchAddress("GFLambdaProductionVtx_y", &GFLambdaProductionVtx_y, &b_GFLambdaProductionVtx_y);
   fChain->SetBranchAddress("GFLambdaProductionVtx_z", &GFLambdaProductionVtx_z, &b_GFLambdaProductionVtx_z);
   fChain->SetBranchAddress("GFLambdaProductionVtxCloseDist", &GFLambdaProductionVtxCloseDist, &b_GFLambdaProductionVtxCloseDist);
   fChain->SetBranchAddress("GFLambdaTrackLen", &GFLambdaTrackLen, &b_GFLambdaTrackLen);
   fChain->SetBranchAddress("GFLambdaTof", &GFLambdaTof, &b_GFLambdaTof);
   fChain->SetBranchAddress("KFKKLProductionVtx_x", &KFKKLProductionVtx_x, &b_KFKKLProductionVtx_x);
   fChain->SetBranchAddress("KFKKLProductionVtx_y", &KFKKLProductionVtx_y, &b_KFKKLProductionVtx_y);
   fChain->SetBranchAddress("KFKKLProductionVtx_z", &KFKKLProductionVtx_z, &b_KFKKLProductionVtx_z);
   fChain->SetBranchAddress("KFLambdaProductionVtx_x", &KFLambdaProductionVtx_x, &b_KFLambdaProductionVtx_x);
   fChain->SetBranchAddress("KFLambdaProductionVtx_y", &KFLambdaProductionVtx_y, &b_KFLambdaProductionVtx_y);
   fChain->SetBranchAddress("KFLambdaProductionVtx_z", &KFLambdaProductionVtx_z, &b_KFLambdaProductionVtx_z);
   fChain->SetBranchAddress("KFLambdaProductionVtxCloseDist", &KFLambdaProductionVtxCloseDist, &b_KFLambdaProductionVtxCloseDist);
   fChain->SetBranchAddress("KFLambdaTrackLen", &KFLambdaTrackLen, &b_KFLambdaTrackLen);
   fChain->SetBranchAddress("KFLambdaTof", &KFLambdaTof, &b_KFLambdaTof);
   fChain->SetBranchAddress("LPhiflag", &LPhiflag, &b_LPhiflag);
   fChain->SetBranchAddress("PhiMass", &PhiMass, &b_PhiMass);
   fChain->SetBranchAddress("PhiDecayVtx_x", &PhiDecayVtx_x, &b_PhiDecayVtx_x);
   fChain->SetBranchAddress("PhiDecayVtx_y", &PhiDecayVtx_y, &b_PhiDecayVtx_y);
   fChain->SetBranchAddress("PhiDecayVtx_z", &PhiDecayVtx_z, &b_PhiDecayVtx_z);
   fChain->SetBranchAddress("PhiCosKK", &PhiCosKK, &b_PhiCosKK);
   fChain->SetBranchAddress("PhiMom", &PhiMom, &b_PhiMom);
   fChain->SetBranchAddress("PhiMom_x", &PhiMom_x, &b_PhiMom_x);
   fChain->SetBranchAddress("PhiMom_y", &PhiMom_y, &b_PhiMom_y);
   fChain->SetBranchAddress("PhiMom_z", &PhiMom_z, &b_PhiMom_z);
   fChain->SetBranchAddress("PhiVtxCloseDist", &PhiVtxCloseDist, &b_PhiVtxCloseDist);
   fChain->SetBranchAddress("PhiDecaysTrackId", &PhiDecaysTrackId, &b_PhiDecaysTrackId);
   fChain->SetBranchAddress("PhiDecaysMom", &PhiDecaysMom, &b_PhiDecaysMom);
   fChain->SetBranchAddress("PhiDecaysMom_x", &PhiDecaysMom_x, &b_PhiDecaysMom_x);
   fChain->SetBranchAddress("PhiDecaysMom_y", &PhiDecaysMom_y, &b_PhiDecaysMom_y);
   fChain->SetBranchAddress("PhiDecaysMom_z", &PhiDecaysMom_z, &b_PhiDecaysMom_z);
   fChain->SetBranchAddress("GFPhiMass", &GFPhiMass, &b_GFPhiMass);
   fChain->SetBranchAddress("GFPhiDecayVtx_x", &GFPhiDecayVtx_x, &b_GFPhiDecayVtx_x);
   fChain->SetBranchAddress("GFPhiDecayVtx_y", &GFPhiDecayVtx_y, &b_GFPhiDecayVtx_y);
   fChain->SetBranchAddress("GFPhiDecayVtx_z", &GFPhiDecayVtx_z, &b_GFPhiDecayVtx_z);
   fChain->SetBranchAddress("GFPhiCosKK", &GFPhiCosKK, &b_GFPhiCosKK);
   fChain->SetBranchAddress("GFPhiMom", &GFPhiMom, &b_GFPhiMom);
   fChain->SetBranchAddress("GFPhiMom_x", &GFPhiMom_x, &b_GFPhiMom_x);
   fChain->SetBranchAddress("GFPhiMom_y", &GFPhiMom_y, &b_GFPhiMom_y);
   fChain->SetBranchAddress("GFPhiMom_z", &GFPhiMom_z, &b_GFPhiMom_z);
   fChain->SetBranchAddress("GFPhiVtxCloseDist", &GFPhiVtxCloseDist, &b_GFPhiVtxCloseDist);
   fChain->SetBranchAddress("GFPhiProductionVtxCloseDist", &GFPhiProductionVtxCloseDist, &b_GFPhiProductionVtxCloseDist);
   fChain->SetBranchAddress("GFPhiKmMassSquare", &GFPhiKmMassSquare, &b_GFPhiKmMassSquare);
   fChain->SetBranchAddress("GFPhiKmInvbeta", &GFPhiKmInvbeta, &b_GFPhiKmInvbeta);
   fChain->SetBranchAddress("GFPhiKpMassSquare", &GFPhiKpMassSquare, &b_GFPhiKpMassSquare);
   fChain->SetBranchAddress("GFPhiKpInvbeta", &GFPhiKpInvbeta, &b_GFPhiKpInvbeta);
   fChain->SetBranchAddress("GFPhiDecaysMom", &GFPhiDecaysMom, &b_GFPhiDecaysMom);
   fChain->SetBranchAddress("GFPhiDecaysMom_x", &GFPhiDecaysMom_x, &b_GFPhiDecaysMom_x);
   fChain->SetBranchAddress("GFPhiDecaysMom_y", &GFPhiDecaysMom_y, &b_GFPhiDecaysMom_y);
   fChain->SetBranchAddress("GFPhiDecaysMom_z", &GFPhiDecaysMom_z, &b_GFPhiDecaysMom_z);
   fChain->SetBranchAddress("GFPhiMassKurama", &GFPhiMassKurama, &b_GFPhiMassKurama);
   fChain->SetBranchAddress("GFPhiMomKurama", &GFPhiMomKurama, &b_GFPhiMomKurama);
   fChain->SetBranchAddress("GFPhiMomKurama_x", &GFPhiMomKurama_x, &b_GFPhiMomKurama_x);
   fChain->SetBranchAddress("GFPhiMomKurama_y", &GFPhiMomKurama_y, &b_GFPhiMomKurama_y);
   fChain->SetBranchAddress("GFPhiMomKurama_z", &GFPhiMomKurama_z, &b_GFPhiMomKurama_z);
   fChain->SetBranchAddress("LPiflag", &LPiflag, &b_LPiflag);
   fChain->SetBranchAddress("PiPiflag", &PiPiflag, &b_PiPiflag);
   fChain->SetBranchAddress("Pimflag", &Pimflag, &b_Pimflag);
   fChain->SetBranchAddress("Emptyflag", &Emptyflag, &b_Emptyflag);
   fChain->SetBranchAddress("GFXiDecaysNhit", &GFXiDecaysNhit, &b_GFXiDecaysNhit);
   fChain->SetBranchAddress("GFXiDecaysChisqr", &GFXiDecaysChisqr, &b_GFXiDecaysChisqr);
   fChain->SetBranchAddress("GFXiDecaysCharge", &GFXiDecaysCharge, &b_GFXiDecaysCharge);
   fChain->SetBranchAddress("GFXiDecaysTof", &GFXiDecaysTof, &b_GFXiDecaysTof);
   fChain->SetBranchAddress("GFXiDecaysPval", &GFXiDecaysPval, &b_GFXiDecaysPval);
   fChain->SetBranchAddress("GFXiDecaysPdgcode", &GFXiDecaysPdgcode, &b_GFXiDecaysPdgcode);
   fChain->SetBranchAddress("GFXiDecaysHtofId", &GFXiDecaysHtofId, &b_GFXiDecaysHtofId);
   fChain->SetBranchAddress("GFXiDecaysTrackLen", &GFXiDecaysTrackLen, &b_GFXiDecaysTrackLen);
   fChain->SetBranchAddress("GFXiDecaysTrackTof", &GFXiDecaysTrackTof, &b_GFXiDecaysTrackTof);
   fChain->SetBranchAddress("GFXiDecaysMassSquare", &GFXiDecaysMassSquare, &b_GFXiDecaysMassSquare);
   fChain->SetBranchAddress("GFXiDecaysInvbeta", &GFXiDecaysInvbeta, &b_GFXiDecaysInvbeta);
   fChain->SetBranchAddress("GFXiDecaysMom", &GFXiDecaysMom, &b_GFXiDecaysMom);
   fChain->SetBranchAddress("GFXiDecaysMom_x", &GFXiDecaysMom_x, &b_GFXiDecaysMom_x);
   fChain->SetBranchAddress("GFXiDecaysMom_y", &GFXiDecaysMom_y, &b_GFXiDecaysMom_y);
   fChain->SetBranchAddress("GFXiDecaysMom_z", &GFXiDecaysMom_z, &b_GFXiDecaysMom_z);
   fChain->SetBranchAddress("GFXiDecaysMomCM", &GFXiDecaysMomCM, &b_GFXiDecaysMomCM);
   fChain->SetBranchAddress("GFXiDecaysMomCM_x", &GFXiDecaysMomCM_x, &b_GFXiDecaysMomCM_x);
   fChain->SetBranchAddress("GFXiDecaysMomCM_y", &GFXiDecaysMomCM_y, &b_GFXiDecaysMomCM_y);
   fChain->SetBranchAddress("GFXiDecaysMomCM_z", &GFXiDecaysMomCM_z, &b_GFXiDecaysMomCM_z);
   fChain->SetBranchAddress("GFXiDecaysMomLoss", &GFXiDecaysMomLoss, &b_GFXiDecaysMomLoss);
   fChain->SetBranchAddress("GFXiDecaysELoss", &GFXiDecaysELoss, &b_GFXiDecaysELoss);
   fChain->SetBranchAddress("XiDecaysTrackId", &XiDecaysTrackId, &b_XiDecaysTrackId);
   fChain->SetBranchAddress("XiDecaysMom", &XiDecaysMom, &b_XiDecaysMom);
   fChain->SetBranchAddress("XiDecaysMom_x", &XiDecaysMom_x, &b_XiDecaysMom_x);
   fChain->SetBranchAddress("XiDecaysMom_y", &XiDecaysMom_y, &b_XiDecaysMom_y);
   fChain->SetBranchAddress("XiDecaysMom_z", &XiDecaysMom_z, &b_XiDecaysMom_z);
   fChain->SetBranchAddress("XiDecaysMomCM", &XiDecaysMomCM, &b_XiDecaysMomCM);
   fChain->SetBranchAddress("XiDecaysMomCM_x", &XiDecaysMomCM_x, &b_XiDecaysMomCM_x);
   fChain->SetBranchAddress("XiDecaysMomCM_y", &XiDecaysMomCM_y, &b_XiDecaysMomCM_y);
   fChain->SetBranchAddress("XiDecaysMomCM_z", &XiDecaysMomCM_z, &b_XiDecaysMomCM_z);
   fChain->SetBranchAddress("GFLLDecaysNhit", &GFLLDecaysNhit, &b_GFLLDecaysNhit);
   fChain->SetBranchAddress("GFLLDecaysChisqr", &GFLLDecaysChisqr, &b_GFLLDecaysChisqr);
   fChain->SetBranchAddress("GFLLDecaysCharge", &GFLLDecaysCharge, &b_GFLLDecaysCharge);
   fChain->SetBranchAddress("GFLLDecaysTof", &GFLLDecaysTof, &b_GFLLDecaysTof);
   fChain->SetBranchAddress("GFLLDecaysPval", &GFLLDecaysPval, &b_GFLLDecaysPval);
   fChain->SetBranchAddress("GFLLDecaysPdgcode", &GFLLDecaysPdgcode, &b_GFLLDecaysPdgcode);
   fChain->SetBranchAddress("GFLLDecaysHtofId", &GFLLDecaysHtofId, &b_GFLLDecaysHtofId);
   fChain->SetBranchAddress("GFLLDecaysTrackLen", &GFLLDecaysTrackLen, &b_GFLLDecaysTrackLen);
   fChain->SetBranchAddress("GFLLDecaysTrackTof", &GFLLDecaysTrackTof, &b_GFLLDecaysTrackTof);
   fChain->SetBranchAddress("GFLLDecaysMassSquare", &GFLLDecaysMassSquare, &b_GFLLDecaysMassSquare);
   fChain->SetBranchAddress("GFLLDecaysInvbeta", &GFLLDecaysInvbeta, &b_GFLLDecaysInvbeta);
   fChain->SetBranchAddress("GFLLDecaysMom", &GFLLDecaysMom, &b_GFLLDecaysMom);
   fChain->SetBranchAddress("GFLLDecaysMom_x", &GFLLDecaysMom_x, &b_GFLLDecaysMom_x);
   fChain->SetBranchAddress("GFLLDecaysMom_y", &GFLLDecaysMom_y, &b_GFLLDecaysMom_y);
   fChain->SetBranchAddress("GFLLDecaysMom_z", &GFLLDecaysMom_z, &b_GFLLDecaysMom_z);
   fChain->SetBranchAddress("GFLLDecaysMomCM", &GFLLDecaysMomCM, &b_GFLLDecaysMomCM);
   fChain->SetBranchAddress("GFLLDecaysMomCM_x", &GFLLDecaysMomCM_x, &b_GFLLDecaysMomCM_x);
   fChain->SetBranchAddress("GFLLDecaysMomCM_y", &GFLLDecaysMomCM_y, &b_GFLLDecaysMomCM_y);
   fChain->SetBranchAddress("GFLLDecaysMomCM_z", &GFLLDecaysMomCM_z, &b_GFLLDecaysMomCM_z);
   fChain->SetBranchAddress("GFLLDecaysMomLoss", &GFLLDecaysMomLoss, &b_GFLLDecaysMomLoss);
   fChain->SetBranchAddress("GFLLDecaysELoss", &GFLLDecaysELoss, &b_GFLLDecaysELoss);
   fChain->SetBranchAddress("LLDecaysTrackId", &LLDecaysTrackId, &b_LLDecaysTrackId);
   fChain->SetBranchAddress("LLDecaysMom", &LLDecaysMom, &b_LLDecaysMom);
   fChain->SetBranchAddress("LLDecaysMom_x", &LLDecaysMom_x, &b_LLDecaysMom_x);
   fChain->SetBranchAddress("LLDecaysMom_y", &LLDecaysMom_y, &b_LLDecaysMom_y);
   fChain->SetBranchAddress("LLDecaysMom_z", &LLDecaysMom_z, &b_LLDecaysMom_z);
   fChain->SetBranchAddress("LLDecaysMomCM", &LLDecaysMomCM, &b_LLDecaysMomCM);
   fChain->SetBranchAddress("LLDecaysMomCM_x", &LLDecaysMomCM_x, &b_LLDecaysMomCM_x);
   fChain->SetBranchAddress("LLDecaysMomCM_y", &LLDecaysMomCM_y, &b_LLDecaysMomCM_y);
   fChain->SetBranchAddress("LLDecaysMomCM_z", &LLDecaysMomCM_z, &b_LLDecaysMomCM_z);
   fChain->SetBranchAddress("GFDecaysNhit", &GFDecaysNhit, &b_GFDecaysNhit);
   fChain->SetBranchAddress("GFDecaysChisqr", &GFDecaysChisqr, &b_GFDecaysChisqr);
   fChain->SetBranchAddress("GFDecaysCharge", &GFDecaysCharge, &b_GFDecaysCharge);
   fChain->SetBranchAddress("GFDecaysTof", &GFDecaysTof, &b_GFDecaysTof);
   fChain->SetBranchAddress("GFDecaysPval", &GFDecaysPval, &b_GFDecaysPval);
   fChain->SetBranchAddress("GFDecaysPdgcode", &GFDecaysPdgcode, &b_GFDecaysPdgcode);
   fChain->SetBranchAddress("GFDecaysHtofId", &GFDecaysHtofId, &b_GFDecaysHtofId);
   fChain->SetBranchAddress("GFDecaysTrackLen", &GFDecaysTrackLen, &b_GFDecaysTrackLen);
   fChain->SetBranchAddress("GFDecaysTrackTof", &GFDecaysTrackTof, &b_GFDecaysTrackTof);
   fChain->SetBranchAddress("GFDecaysMassSquare", &GFDecaysMassSquare, &b_GFDecaysMassSquare);
   fChain->SetBranchAddress("GFDecaysInvbeta", &GFDecaysInvbeta, &b_GFDecaysInvbeta);
   fChain->SetBranchAddress("GFDecaysMom", &GFDecaysMom, &b_GFDecaysMom);
   fChain->SetBranchAddress("GFDecaysMom_x", &GFDecaysMom_x, &b_GFDecaysMom_x);
   fChain->SetBranchAddress("GFDecaysMom_y", &GFDecaysMom_y, &b_GFDecaysMom_y);
   fChain->SetBranchAddress("GFDecaysMom_z", &GFDecaysMom_z, &b_GFDecaysMom_z);
   fChain->SetBranchAddress("GFDecaysMomCM", &GFDecaysMomCM, &b_GFDecaysMomCM);
   fChain->SetBranchAddress("GFDecaysMomCM_x", &GFDecaysMomCM_x, &b_GFDecaysMomCM_x);
   fChain->SetBranchAddress("GFDecaysMomCM_y", &GFDecaysMomCM_y, &b_GFDecaysMomCM_y);
   fChain->SetBranchAddress("GFDecaysMomCM_z", &GFDecaysMomCM_z, &b_GFDecaysMomCM_z);
   fChain->SetBranchAddress("GFDecaysMomLoss", &GFDecaysMomLoss, &b_GFDecaysMomLoss);
   fChain->SetBranchAddress("GFDecaysELoss", &GFDecaysELoss, &b_GFDecaysELoss);
   fChain->SetBranchAddress("DecaysTrackId", &DecaysTrackId, &b_DecaysTrackId);
   fChain->SetBranchAddress("DecaysMom", &DecaysMom, &b_DecaysMom);
   fChain->SetBranchAddress("DecaysMom_x", &DecaysMom_x, &b_DecaysMom_x);
   fChain->SetBranchAddress("DecaysMom_y", &DecaysMom_y, &b_DecaysMom_y);
   fChain->SetBranchAddress("DecaysMom_z", &DecaysMom_z, &b_DecaysMom_z);
   fChain->SetBranchAddress("DecaysMomCM", &DecaysMomCM, &b_DecaysMomCM);
   fChain->SetBranchAddress("DecaysMomCM_x", &DecaysMomCM_x, &b_DecaysMomCM_x);
   fChain->SetBranchAddress("DecaysMomCM_y", &DecaysMomCM_y, &b_DecaysMomCM_y);
   fChain->SetBranchAddress("DecaysMomCM_z", &DecaysMomCM_z, &b_DecaysMomCM_z);
   fChain->SetBranchAddress("AccidentalMultiplicity", &AccidentalMultiplicity, &b_AccidentalMultiplicity);
   fChain->SetBranchAddress("AccidentalTrackId", &AccidentalTrackId, &b_AccidentalTrackId);
   fChain->SetBranchAddress("LLResidualsMultiplicity", &LLResidualsMultiplicity, &b_LLResidualsMultiplicity);
   fChain->SetBranchAddress("LLpMultiplicity", &LLpMultiplicity, &b_LLpMultiplicity);
   fChain->SetBranchAddress("LLpipMultiplicity", &LLpipMultiplicity, &b_LLpipMultiplicity);
   fChain->SetBranchAddress("LLpimMultiplicity", &LLpimMultiplicity, &b_LLpimMultiplicity);
   fChain->SetBranchAddress("LLepMultiplicity", &LLepMultiplicity, &b_LLepMultiplicity);
   fChain->SetBranchAddress("LLemMultiplicity", &LLemMultiplicity, &b_LLemMultiplicity);
   fChain->SetBranchAddress("LLppipMultiplicity", &LLppipMultiplicity, &b_LLppipMultiplicity);
   fChain->SetBranchAddress("LLResidualsTrackId", &LLResidualsTrackId, &b_LLResidualsTrackId);
   fChain->SetBranchAddress("LLResidualsMassSquare", &LLResidualsMassSquare, &b_LLResidualsMassSquare);
   fChain->SetBranchAddress("LLResidualsInvbeta", &LLResidualsInvbeta, &b_LLResidualsInvbeta);
   fChain->SetBranchAddress("LLResidualsCloseDistTgt", &LLResidualsCloseDistTgt, &b_LLResidualsCloseDistTgt);
   fChain->SetBranchAddress("LLResidualsGFProductionVtxCloseDist", &LLResidualsGFProductionVtxCloseDist, &b_LLResidualsGFProductionVtxCloseDist);
   fChain->SetBranchAddress("LLResidualsKFProductionVtxCloseDist", &LLResidualsKFProductionVtxCloseDist, &b_LLResidualsKFProductionVtxCloseDist);
   fChain->SetBranchAddress("LLResidualsMom", &LLResidualsMom, &b_LLResidualsMom);
   fChain->SetBranchAddress("LLResidualsMom_x", &LLResidualsMom_x, &b_LLResidualsMom_x);
   fChain->SetBranchAddress("LLResidualsMom_y", &LLResidualsMom_y, &b_LLResidualsMom_y);
   fChain->SetBranchAddress("LLResidualsMom_z", &LLResidualsMom_z, &b_LLResidualsMom_z);
   fChain->SetBranchAddress("LLResidualsCharge", &LLResidualsCharge, &b_LLResidualsCharge);
   fChain->SetBranchAddress("XiResidualsMultiplicity", &XiResidualsMultiplicity, &b_XiResidualsMultiplicity);
   fChain->SetBranchAddress("XipMultiplicity", &XipMultiplicity, &b_XipMultiplicity);
   fChain->SetBranchAddress("XipipMultiplicity", &XipipMultiplicity, &b_XipipMultiplicity);
   fChain->SetBranchAddress("XipimMultiplicity", &XipimMultiplicity, &b_XipimMultiplicity);
   fChain->SetBranchAddress("XiepMultiplicity", &XiepMultiplicity, &b_XiepMultiplicity);
   fChain->SetBranchAddress("XiemMultiplicity", &XiemMultiplicity, &b_XiemMultiplicity);
   fChain->SetBranchAddress("XippipMultiplicity", &XippipMultiplicity, &b_XippipMultiplicity);
   fChain->SetBranchAddress("XiResidualsTrackId", &XiResidualsTrackId, &b_XiResidualsTrackId);
   fChain->SetBranchAddress("XiResidualsMassSquare", &XiResidualsMassSquare, &b_XiResidualsMassSquare);
   fChain->SetBranchAddress("XiResidualsInvbeta", &XiResidualsInvbeta, &b_XiResidualsInvbeta);
   fChain->SetBranchAddress("XiResidualsCloseDistTgt", &XiResidualsCloseDistTgt, &b_XiResidualsCloseDistTgt);
   fChain->SetBranchAddress("XiResidualsGFProductionVtxCloseDist", &XiResidualsGFProductionVtxCloseDist, &b_XiResidualsGFProductionVtxCloseDist);
   fChain->SetBranchAddress("XiResidualsKFProductionVtxCloseDist", &XiResidualsKFProductionVtxCloseDist, &b_XiResidualsKFProductionVtxCloseDist);
   fChain->SetBranchAddress("XiResidualsMom", &XiResidualsMom, &b_XiResidualsMom);
   fChain->SetBranchAddress("XiResidualsMom_x", &XiResidualsMom_x, &b_XiResidualsMom_x);
   fChain->SetBranchAddress("XiResidualsMom_y", &XiResidualsMom_y, &b_XiResidualsMom_y);
   fChain->SetBranchAddress("XiResidualsMom_z", &XiResidualsMom_z, &b_XiResidualsMom_z);
   fChain->SetBranchAddress("XiResidualsCharge", &XiResidualsCharge, &b_XiResidualsCharge);
   fChain->SetBranchAddress("ResidualsMultiplicity", &ResidualsMultiplicity, &b_ResidualsMultiplicity);
   fChain->SetBranchAddress("pMultiplicity", &pMultiplicity, &b_pMultiplicity);
   fChain->SetBranchAddress("pipMultiplicity", &pipMultiplicity, &b_pipMultiplicity);
   fChain->SetBranchAddress("pimMultiplicity", &pimMultiplicity, &b_pimMultiplicity);
   fChain->SetBranchAddress("epMultiplicity", &epMultiplicity, &b_epMultiplicity);
   fChain->SetBranchAddress("emMultiplicity", &emMultiplicity, &b_emMultiplicity);
   fChain->SetBranchAddress("ppipMultiplicity", &ppipMultiplicity, &b_ppipMultiplicity);
   fChain->SetBranchAddress("ResidualsTrackId", &ResidualsTrackId, &b_ResidualsTrackId);
   fChain->SetBranchAddress("ResidualsMassSquare", &ResidualsMassSquare, &b_ResidualsMassSquare);
   fChain->SetBranchAddress("ResidualsInvbeta", &ResidualsInvbeta, &b_ResidualsInvbeta);
   fChain->SetBranchAddress("ResidualsCloseDistTgt", &ResidualsCloseDistTgt, &b_ResidualsCloseDistTgt);
   fChain->SetBranchAddress("ResidualsGFProductionVtxCloseDist", &ResidualsGFProductionVtxCloseDist, &b_ResidualsGFProductionVtxCloseDist);
   fChain->SetBranchAddress("ResidualsKFProductionVtxCloseDist", &ResidualsKFProductionVtxCloseDist, &b_ResidualsKFProductionVtxCloseDist);
   fChain->SetBranchAddress("ResidualsMom", &ResidualsMom, &b_ResidualsMom);
   fChain->SetBranchAddress("ResidualsMom_x", &ResidualsMom_x, &b_ResidualsMom_x);
   fChain->SetBranchAddress("ResidualsMom_y", &ResidualsMom_y, &b_ResidualsMom_y);
   fChain->SetBranchAddress("ResidualsMom_z", &ResidualsMom_z, &b_ResidualsMom_z);
   fChain->SetBranchAddress("ResidualsCharge", &ResidualsCharge, &b_ResidualsCharge);
   fChain->SetBranchAddress("KFLambdaMomPpi", &KFLambdaMomPpi, &b_KFLambdaMomPpi);
   fChain->SetBranchAddress("KFLambdaMomPpi_x", &KFLambdaMomPpi_x, &b_KFLambdaMomPpi_x);
   fChain->SetBranchAddress("KFLambdaMomPpi_y0", &KFLambdaMomPpi_y0, &b_KFLambdaMomPpi_y0);
   fChain->SetBranchAddress("KFLambdaMomPpi_z0", &KFLambdaMomPpi_z0, &b_KFLambdaMomPpi_z0);
   fChain->SetBranchAddress("KFLambdaMom", &KFLambdaMom, &b_KFLambdaMom);
   fChain->SetBranchAddress("KFLambdaMom_x", &KFLambdaMom_x, &b_KFLambdaMom_x);
   fChain->SetBranchAddress("KFLambdaMom_y", &KFLambdaMom_y, &b_KFLambdaMom_y);
   fChain->SetBranchAddress("KFLambdaMom_z", &KFLambdaMom_z, &b_KFLambdaMom_z);
   fChain->SetBranchAddress("KFLambdaChisqr", &KFLambdaChisqr, &b_KFLambdaChisqr);
   fChain->SetBranchAddress("KFLambdaPval", &KFLambdaPval, &b_KFLambdaPval);
   fChain->SetBranchAddress("KFLambdaPull", &KFLambdaPull, &b_KFLambdaPull);
   fChain->SetBranchAddress("KFLambdaCovMatrix", &KFLambdaCovMatrix, &b_KFLambdaCovMatrix);
   fChain->SetBranchAddress("KFXiVtxCloseDist", &KFXiVtxCloseDist, &b_KFXiVtxCloseDist);
   fChain->SetBranchAddress("KFXiMom", &KFXiMom, &b_KFXiMom);
   fChain->SetBranchAddress("KFXiMom_x", &KFXiMom_x, &b_KFXiMom_x);
   fChain->SetBranchAddress("KFXiMom_y", &KFXiMom_y, &b_KFXiMom_y);
   fChain->SetBranchAddress("KFXiMom_z", &KFXiMom_z, &b_KFXiMom_z);
   fChain->SetBranchAddress("KFXiChisqr", &KFXiChisqr, &b_KFXiChisqr);
   fChain->SetBranchAddress("KFXiPval", &KFXiPval, &b_KFXiPval);
   fChain->SetBranchAddress("KFXiMass", &KFXiMass, &b_KFXiMass);
   fChain->SetBranchAddress("KFXiDecayVtx_x", &KFXiDecayVtx_x, &b_KFXiDecayVtx_x);
   fChain->SetBranchAddress("KFXiDecayVtx_y", &KFXiDecayVtx_y, &b_KFXiDecayVtx_y);
   fChain->SetBranchAddress("KFXiDecayVtx_z", &KFXiDecayVtx_z, &b_KFXiDecayVtx_z);
   fChain->SetBranchAddress("KFXiPull", &KFXiPull, &b_KFXiPull);
   fChain->SetBranchAddress("KFXiCovMatrix", &KFXiCovMatrix, &b_KFXiCovMatrix);
   fChain->SetBranchAddress("KFKKXiProductionVtxChisqr", &KFKKXiProductionVtxChisqr, &b_KFKKXiProductionVtxChisqr);
   fChain->SetBranchAddress("KFKKXiProductionVtx_x", &KFKKXiProductionVtx_x, &b_KFKKXiProductionVtx_x);
   fChain->SetBranchAddress("KFKKXiProductionVtx_y", &KFKKXiProductionVtx_y, &b_KFKKXiProductionVtx_y);
   fChain->SetBranchAddress("KFKKXiProductionVtx_z", &KFKKXiProductionVtx_z, &b_KFKKXiProductionVtx_z);
   fChain->SetBranchAddress("KFXiProductionVtx_x", &KFXiProductionVtx_x, &b_KFXiProductionVtx_x);
   fChain->SetBranchAddress("KFXiProductionVtx_y", &KFXiProductionVtx_y, &b_KFXiProductionVtx_y);
   fChain->SetBranchAddress("KFXiProductionVtx_z", &KFXiProductionVtx_z, &b_KFXiProductionVtx_z);
   fChain->SetBranchAddress("KFXiProductionVtxMom", &KFXiProductionVtxMom, &b_KFXiProductionVtxMom);
   fChain->SetBranchAddress("KFXiProductionVtxMom_x", &KFXiProductionVtxMom_x, &b_KFXiProductionVtxMom_x);
   fChain->SetBranchAddress("KFXiProductionVtxMom_y", &KFXiProductionVtxMom_y, &b_KFXiProductionVtxMom_y);
   fChain->SetBranchAddress("KFXiProductionVtxMom_z", &KFXiProductionVtxMom_z, &b_KFXiProductionVtxMom_z);
   fChain->SetBranchAddress("KFXiProductionVtxCloseDist", &KFXiProductionVtxCloseDist, &b_KFXiProductionVtxCloseDist);
   fChain->SetBranchAddress("KFXiTrackLen", &KFXiTrackLen, &b_KFXiTrackLen);
   fChain->SetBranchAddress("KFXiTof", &KFXiTof, &b_KFXiTof);
   fChain->SetBranchAddress("KFXiMomLoss", &KFXiMomLoss, &b_KFXiMomLoss);
   fChain->SetBranchAddress("KFXiExcitation", &KFXiExcitation, &b_KFXiExcitation);
   fChain->SetBranchAddress("KFKpXiProductionVtx_x", &KFKpXiProductionVtx_x, &b_KFKpXiProductionVtx_x);
   fChain->SetBranchAddress("KFKpXiProductionVtx_y", &KFKpXiProductionVtx_y, &b_KFKpXiProductionVtx_y);
   fChain->SetBranchAddress("KFKpXiProductionVtx_z", &KFKpXiProductionVtx_z, &b_KFKpXiProductionVtx_z);
   fChain->SetBranchAddress("KFXiProductionVtx_x_KpXi", &KFXiProductionVtx_x_KpXi, &b_KFXiProductionVtx_x_KpXi);
   fChain->SetBranchAddress("KFXiProductionVtx_y_KpXi", &KFXiProductionVtx_y_KpXi, &b_KFXiProductionVtx_y_KpXi);
   fChain->SetBranchAddress("KFXiProductionVtx_z_KpXi", &KFXiProductionVtx_z_KpXi, &b_KFXiProductionVtx_z_KpXi);
   fChain->SetBranchAddress("KFXiProductionVtxMom_KpXi", &KFXiProductionVtxMom_KpXi, &b_KFXiProductionVtxMom_KpXi);
   fChain->SetBranchAddress("KFXiProductionVtxMom_x_KpXi", &KFXiProductionVtxMom_x_KpXi, &b_KFXiProductionVtxMom_x_KpXi);
   fChain->SetBranchAddress("KFXiProductionVtxMom_y_KpXi", &KFXiProductionVtxMom_y_KpXi, &b_KFXiProductionVtxMom_y_KpXi);
   fChain->SetBranchAddress("KFXiProductionVtxMom_z_KpXi", &KFXiProductionVtxMom_z_KpXi, &b_KFXiProductionVtxMom_z_KpXi);
   fChain->SetBranchAddress("KFXiProductionVtxCloseDist_KpXi", &KFXiProductionVtxCloseDist_KpXi, &b_KFXiProductionVtxCloseDist_KpXi);
   fChain->SetBranchAddress("KFXiProductionVtx_x_KK", &KFXiProductionVtx_x_KK, &b_KFXiProductionVtx_x_KK);
   fChain->SetBranchAddress("KFXiProductionVtx_y_KK", &KFXiProductionVtx_y_KK, &b_KFXiProductionVtx_y_KK);
   fChain->SetBranchAddress("KFXiProductionVtx_z_KK", &KFXiProductionVtx_z_KK, &b_KFXiProductionVtx_z_KK);
   fChain->SetBranchAddress("KFXiProductionVtxMom_KK", &KFXiProductionVtxMom_KK, &b_KFXiProductionVtxMom_KK);
   fChain->SetBranchAddress("KFXiProductionVtxMom_x_KK", &KFXiProductionVtxMom_x_KK, &b_KFXiProductionVtxMom_x_KK);
   fChain->SetBranchAddress("KFXiProductionVtxMom_y_KK", &KFXiProductionVtxMom_y_KK, &b_KFXiProductionVtxMom_y_KK);
   fChain->SetBranchAddress("KFXiProductionVtxMom_z_KK", &KFXiProductionVtxMom_z_KK, &b_KFXiProductionVtxMom_z_KK);
   fChain->SetBranchAddress("KFXiProductionVtxCloseDist_KK", &KFXiProductionVtxCloseDist_KK, &b_KFXiProductionVtxCloseDist_KK);
   fChain->SetBranchAddress("KFLLDecaysMom", &KFLLDecaysMom, &b_KFLLDecaysMom);
   fChain->SetBranchAddress("KFLLDecaysMom_x", &KFLLDecaysMom_x, &b_KFLLDecaysMom_x);
   fChain->SetBranchAddress("KFLLDecaysMom_y", &KFLLDecaysMom_y, &b_KFLLDecaysMom_y);
   fChain->SetBranchAddress("KFLLDecaysMom_z", &KFLLDecaysMom_z, &b_KFLLDecaysMom_z);
   fChain->SetBranchAddress("KFLLDecaysMomCM", &KFLLDecaysMomCM, &b_KFLLDecaysMomCM);
   fChain->SetBranchAddress("KFLLDecaysMomCM_x", &KFLLDecaysMomCM_x, &b_KFLLDecaysMomCM_x);
   fChain->SetBranchAddress("KFLLDecaysMomCM_y", &KFLLDecaysMomCM_y, &b_KFLLDecaysMomCM_y);
   fChain->SetBranchAddress("KFLLDecaysMomCM_z", &KFLLDecaysMomCM_z, &b_KFLLDecaysMomCM_z);
   fChain->SetBranchAddress("KFXiDecaysMom", &KFXiDecaysMom, &b_KFXiDecaysMom);
   fChain->SetBranchAddress("KFXiDecaysMom_x", &KFXiDecaysMom_x, &b_KFXiDecaysMom_x);
   fChain->SetBranchAddress("KFXiDecaysMom_y", &KFXiDecaysMom_y, &b_KFXiDecaysMom_y);
   fChain->SetBranchAddress("KFXiDecaysMom_z", &KFXiDecaysMom_z, &b_KFXiDecaysMom_z);
   fChain->SetBranchAddress("KFXiDecaysMomCM", &KFXiDecaysMomCM, &b_KFXiDecaysMomCM);
   fChain->SetBranchAddress("KFXiDecaysMomCM_x", &KFXiDecaysMomCM_x, &b_KFXiDecaysMomCM_x);
   fChain->SetBranchAddress("KFXiDecaysMomCM_y", &KFXiDecaysMomCM_y, &b_KFXiDecaysMomCM_y);
   fChain->SetBranchAddress("KFXiDecaysMomCM_z", &KFXiDecaysMomCM_z, &b_KFXiDecaysMomCM_z);
   fChain->SetBranchAddress("KFDecaysMom", &KFDecaysMom, &b_KFDecaysMom);
   fChain->SetBranchAddress("KFDecaysMom_x", &KFDecaysMom_x, &b_KFDecaysMom_x);
   fChain->SetBranchAddress("KFDecaysMom_y", &KFDecaysMom_y, &b_KFDecaysMom_y);
   fChain->SetBranchAddress("KFDecaysMom_z", &KFDecaysMom_z, &b_KFDecaysMom_z);
   fChain->SetBranchAddress("KFDecaysMomCM", &KFDecaysMomCM, &b_KFDecaysMomCM);
   fChain->SetBranchAddress("KFDecaysMomCM_x", &KFDecaysMomCM_x, &b_KFDecaysMomCM_x);
   fChain->SetBranchAddress("KFDecaysMomCM_y", &KFDecaysMomCM_y, &b_KFDecaysMomCM_y);
   fChain->SetBranchAddress("KFDecaysMomCM_z", &KFDecaysMomCM_z, &b_KFDecaysMomCM_z);
   fChain->SetBranchAddress("LLnGamma", &LLnGamma, &b_LLnGamma);
   fChain->SetBranchAddress("LLGammaEpTrackId", &LLGammaEpTrackId, &b_LLGammaEpTrackId);
   fChain->SetBranchAddress("LLGammaEmTrackId", &LLGammaEmTrackId, &b_LLGammaEmTrackId);
   fChain->SetBranchAddress("LLGammaMomId", &LLGammaMomId, &b_LLGammaMomId);
   fChain->SetBranchAddress("LLGammaEpMom_x", &LLGammaEpMom_x, &b_LLGammaEpMom_x);
   fChain->SetBranchAddress("LLGammaEpMom_y", &LLGammaEpMom_y, &b_LLGammaEpMom_y);
   fChain->SetBranchAddress("LLGammaEpMom_z", &LLGammaEpMom_z, &b_LLGammaEpMom_z);
   fChain->SetBranchAddress("LLGammaEmMom", &LLGammaEmMom, &b_LLGammaEmMom);
   fChain->SetBranchAddress("LLGammaEmMom_x", &LLGammaEmMom_x, &b_LLGammaEmMom_x);
   fChain->SetBranchAddress("LLGammaEmMom_y", &LLGammaEmMom_y, &b_LLGammaEmMom_y);
   fChain->SetBranchAddress("LLGammaEmMom_z", &LLGammaEmMom_z, &b_LLGammaEmMom_z);
   fChain->SetBranchAddress("LLGammaMom", &LLGammaMom, &b_LLGammaMom);
   fChain->SetBranchAddress("LLGammaMom_x", &LLGammaMom_x, &b_LLGammaMom_x);
   fChain->SetBranchAddress("LLGammaMom_y", &LLGammaMom_y, &b_LLGammaMom_y);
   fChain->SetBranchAddress("LLGammaMom_z", &LLGammaMom_z, &b_LLGammaMom_z);
   fChain->SetBranchAddress("LLGammaVtxCloseDist", &LLGammaVtxCloseDist, &b_LLGammaVtxCloseDist);
   fChain->SetBranchAddress("LLGammaDecayVtx_x", &LLGammaDecayVtx_x, &b_LLGammaDecayVtx_x);
   fChain->SetBranchAddress("LLGammaDecayVtx_y", &LLGammaDecayVtx_y, &b_LLGammaDecayVtx_y);
   fChain->SetBranchAddress("LLGammaDecayVtx_z", &LLGammaDecayVtx_z, &b_LLGammaDecayVtx_z);
   fChain->SetBranchAddress("XinGamma", &XinGamma, &b_XinGamma);
   fChain->SetBranchAddress("XiGammaEpTrackId", &XiGammaEpTrackId, &b_XiGammaEpTrackId);
   fChain->SetBranchAddress("XiGammaEmTrackId", &XiGammaEmTrackId, &b_XiGammaEmTrackId);
   fChain->SetBranchAddress("XiGammaMomId", &XiGammaMomId, &b_XiGammaMomId);
   fChain->SetBranchAddress("XiGammaEpMom_x", &XiGammaEpMom_x, &b_XiGammaEpMom_x);
   fChain->SetBranchAddress("XiGammaEpMom_y", &XiGammaEpMom_y, &b_XiGammaEpMom_y);
   fChain->SetBranchAddress("XiGammaEpMom_z", &XiGammaEpMom_z, &b_XiGammaEpMom_z);
   fChain->SetBranchAddress("XiGammaEmMom", &XiGammaEmMom, &b_XiGammaEmMom);
   fChain->SetBranchAddress("XiGammaEmMom_x", &XiGammaEmMom_x, &b_XiGammaEmMom_x);
   fChain->SetBranchAddress("XiGammaEmMom_y", &XiGammaEmMom_y, &b_XiGammaEmMom_y);
   fChain->SetBranchAddress("XiGammaEmMom_z", &XiGammaEmMom_z, &b_XiGammaEmMom_z);
   fChain->SetBranchAddress("XiGammaMom", &XiGammaMom, &b_XiGammaMom);
   fChain->SetBranchAddress("XiGammaMom_x", &XiGammaMom_x, &b_XiGammaMom_x);
   fChain->SetBranchAddress("XiGammaMom_y", &XiGammaMom_y, &b_XiGammaMom_y);
   fChain->SetBranchAddress("XiGammaMom_z", &XiGammaMom_z, &b_XiGammaMom_z);
   fChain->SetBranchAddress("XiGammaVtxCloseDist", &XiGammaVtxCloseDist, &b_XiGammaVtxCloseDist);
   fChain->SetBranchAddress("XiGammaDecayVtx_x", &XiGammaDecayVtx_x, &b_XiGammaDecayVtx_x);
   fChain->SetBranchAddress("XiGammaDecayVtx_y", &XiGammaDecayVtx_y, &b_XiGammaDecayVtx_y);
   fChain->SetBranchAddress("XiGammaDecayVtx_z", &XiGammaDecayVtx_z, &b_XiGammaDecayVtx_z);
   fChain->SetBranchAddress("nGamma", &nGamma, &b_nGamma);
   fChain->SetBranchAddress("GammaEpTrackId", &GammaEpTrackId, &b_GammaEpTrackId);
   fChain->SetBranchAddress("GammaEmTrackId", &GammaEmTrackId, &b_GammaEmTrackId);
   fChain->SetBranchAddress("GammaMomId", &GammaMomId, &b_GammaMomId);
   fChain->SetBranchAddress("GammaEpMom_x", &GammaEpMom_x, &b_GammaEpMom_x);
   fChain->SetBranchAddress("GammaEpMom_y", &GammaEpMom_y, &b_GammaEpMom_y);
   fChain->SetBranchAddress("GammaEpMom_z", &GammaEpMom_z, &b_GammaEpMom_z);
   fChain->SetBranchAddress("GammaEmMom", &GammaEmMom, &b_GammaEmMom);
   fChain->SetBranchAddress("GammaEmMom_x", &GammaEmMom_x, &b_GammaEmMom_x);
   fChain->SetBranchAddress("GammaEmMom_y", &GammaEmMom_y, &b_GammaEmMom_y);
   fChain->SetBranchAddress("GammaEmMom_z", &GammaEmMom_z, &b_GammaEmMom_z);
   fChain->SetBranchAddress("GammaMom", &GammaMom, &b_GammaMom);
   fChain->SetBranchAddress("GammaMom_x", &GammaMom_x, &b_GammaMom_x);
   fChain->SetBranchAddress("GammaMom_y", &GammaMom_y, &b_GammaMom_y);
   fChain->SetBranchAddress("GammaMom_z", &GammaMom_z, &b_GammaMom_z);
   fChain->SetBranchAddress("GammaVtxCloseDist", &GammaVtxCloseDist, &b_GammaVtxCloseDist);
   fChain->SetBranchAddress("GammaDecayVtx_x", &GammaDecayVtx_x, &b_GammaDecayVtx_x);
   fChain->SetBranchAddress("GammaDecayVtx_y", &GammaDecayVtx_y, &b_GammaDecayVtx_y);
   fChain->SetBranchAddress("GammaDecayVtx_z", &GammaDecayVtx_z, &b_GammaDecayVtx_z);
   Notify();
}

bool tpc::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void tpc::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tpc::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef tpc_cxx
