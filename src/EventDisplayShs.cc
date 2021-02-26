// -*- C++ -*-

#include "EventDisplayShs.hh"

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include <TApplication.h>
#include <TBRIK.h>
#include <TArc.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TEnv.h>
#include <TFile.h>
#include <TGeometry.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH2Poly.h>
#include <TF1.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TMarker3DBox.h>
#include <TMixture.h>
#include <TNode.h>
#include <TPad.h>
#include <TSystem.h>
#include <TPave.h>
#include <TPaveLabel.h>
#include <TPaveText.h>
#include <TPolyLine.h>
#include <TPolyLine3D.h>
#include <TPolyMarker.h>
#include <TPolyMarker3D.h>
#include <TROOT.h>
#include <TRint.h>
#include <TRotMatrix.h>
#include <TStyle.h>
#include <TTRD1.h>
#include <TTRD2.h>
#include <TTUBS.h>
#include <TTUBS.h>
#include <TView.h>
#include <TLine.h>

#include <std_ostream.hh>

#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "DeleteUtility.hh"
#include "HodoParamMan.hh"
#include "DCTdcCalibMan.hh"
#include "TPCPadHelper.hh"

#define BH2        1
#define BcOut      1
#define KURAMA     1
#define SdcIn      1
#define SdcOut     1
#define SCH        1
#define FBT        1
#define TOF        1
#define Vertex     1
#define Hist       1
#define Hist_Timing 1
#define Hist_SdcOut 1
#define Hist_BcIn   1

namespace
{
  const DCGeomMan& gGeom = DCGeomMan::GetInstance();
  const int& IdBH2 = gGeom.DetectorId("BH2");
  const int& IdBH1 = gGeom.DetectorId("BH1");
  const int& IdSFT_U = gGeom.DetectorId("SFT-U");
  const int& IdSFT_V = gGeom.DetectorId("SFT-V");
  const int& IdSFT_X = gGeom.DetectorId("SFT-X");
  const int& IdSCH = gGeom.DetectorId("SCH");
  const int& IdFBT1_D1 = gGeom.DetectorId("FBT1-DX1");
  const int& IdFBT1_U1 = gGeom.DetectorId("FBT1-UX1");
  const int& IdFBT1_D2 = gGeom.DetectorId("FBT1-DX2");
  const int& IdFBT1_U2 = gGeom.DetectorId("FBT1-UX2");
  const int& IdFBT2_D1 = gGeom.DetectorId("FBT2-DX1");
  const int& IdFBT2_U1 = gGeom.DetectorId("FBT2-UX1");
  const int& IdFBT2_D2 = gGeom.DetectorId("FBT2-DX2");
  const int& IdFBT2_U2 = gGeom.DetectorId("FBT2-UX2");
  const int& IdTOF = gGeom.DetectorId("TOF");
  const int& IdTarget    = gGeom.DetectorId("Target");
  const double& zTarget    = gGeom.LocalZ("Target");
  const double& zK18Target = gGeom.LocalZ("K18Target");
  const double& gzK18Target = gGeom.GlobalZ("K18Target");
  //const double& gxK18Target = gGeom.GetGlobalPosition("K18Target").x();
  const double& gxK18Target = -240.;
  const double& zBFT    = gGeom.LocalZ("BFT");

  //const double BeamAxis = -150.; //E07
  const double BeamAxis = -240.; //E40
#if Vertex
  const double MinX = -50.;
  const double MaxX =  50.;
  const double MinY = -50.;
  const double MaxY =  50.;
  const double MinZ = -25.;
#endif
  const double MaxZ =  50.;

  const HodoParamMan& gHodo = HodoParamMan::GetInstance();
  const DCTdcCalibMan& gTdc = DCTdcCalibMan::GetInstance();

//______________________________________________________________________________
// local function
[[maybe_unused]] void
ConstructionDone( const char* name, std::ostream& ost=hddaq::cout )
{
  const std::size_t n = 20;
  std::size_t s = std::string(name).size();
  ost << " " << name << " ";
  for( std::size_t i=0; i<n-s; ++i )
    ost << ".";
  ost << " done" << std::endl;
}

}

//______________________________________________________________________________
EventDisplayShs::EventDisplayShs( void )
  : m_is_ready(false),
    m_theApp(nullptr),
    m_geometry(nullptr),
    m_node(nullptr),
    m_canvas(nullptr),
    m_tpc_adc(nullptr),
    m_tpc_tdc(nullptr),
    m_tpc_adc2d(nullptr),
    m_tpc_tdc2d(nullptr)
{
}

//______________________________________________________________________________
EventDisplayShs::~EventDisplayShs( void )
{
}

//______________________________________________________________________________
void
EventDisplayShs::FillTPCADC(Int_t layer, Int_t row, Double_t adc)
{
  Int_t pad = tpc::GetPadId(layer, row);
  m_tpc_adc->Fill(adc);
  m_tpc_adc2d->SetBinContent(pad + 1, adc);
}

//______________________________________________________________________________
void
EventDisplayShs::FillTPCTDC(Int_t layer, Int_t row, Double_t tdc)
{
  Int_t pad = tpc::GetPadId(layer, row);
  m_tpc_tdc->Fill(tdc);
  m_tpc_tdc2d->SetBinContent(pad + 1, tdc);
}

//______________________________________________________________________________
Bool_t
EventDisplayShs::Initialize( void )
{
  if (m_is_ready) {
    hddaq::cerr << FUNC_NAME << " already initialied" << std::endl;
    return false;
  }

  m_theApp = new TApplication("App", 0, 0);

#if ROOT_VERSION_CODE > ROOT_VERSION(6,4,0)
  gStyle->SetPalette(kBird);
#endif
  gStyle->SetNumberContours(255);

  m_canvas = new TCanvas("Event Display", "Event Display", 800, 800);

  m_tpc_adc = new TH1D("h_tpc_adc", "TPC ADC", 4096, 0, 4096);
  m_tpc_tdc = new TH1D("h_tpc_adc", "TPC ADC",
                       NumOfTimeBucket, 0, NumOfTimeBucket);
  const Double_t MinX = -300.;
  const Double_t MaxX =  300.;
  const Double_t MinZ = -300.;
  const Double_t MaxZ =  300.;
  m_tpc_adc2d = new TH2Poly("h_tpc_adc2d", "TPC ADC;Z;X", MinZ, MaxZ, MinX, MaxX);
  m_tpc_tdc2d = new TH2Poly("h_tpc_tdc2d", "TPC TDC;Z;X", MinZ, MaxZ, MinX, MaxX);

  Double_t X[5];
  Double_t Y[5];
  for (Int_t l=0; l<NumOfLayersTPC; ++l) {
    Double_t pLength = tpc::padParameter[l][5];
    Double_t st      = (180.-(360./tpc::padParameter[l][3]) *
                        tpc::padParameter[l][1]/2.);
    Double_t sTheta  = (-1+st/180.)*TMath::Pi();
    Double_t dTheta  = (360./tpc::padParameter[l][3])/180.*TMath::Pi();
    Double_t cRad    = tpc::padParameter[l][2];
    Int_t    nPad    = tpc::padParameter[l][1];
    for (Int_t j=0; j<nPad; ++j) {
      X[1] = (cRad+(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
      X[2] = (cRad+(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
      X[3] = (cRad-(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
      X[4] = (cRad-(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
      X[0] = X[4];
      Y[1] = (cRad+(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
      Y[2] = (cRad+(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
      Y[3] = (cRad-(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
      Y[4] = (cRad-(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
      Y[0] = Y[4];
      for (Int_t k=0; k<5; ++k) X[k] -=143;
      m_tpc_adc2d->AddBin(5, X, Y);
      m_tpc_tdc2d->AddBin(5, X, Y);
    }
  }

  m_tpc_adc2d->SetStats(0);
  m_tpc_tdc2d->SetStats(0);
  m_canvas->Divide(2, 2);
  m_canvas->cd(1);
  m_tpc_adc2d->Draw("colz");
  m_canvas->cd(2);
  m_tpc_tdc2d->Draw("colz");
  m_canvas->cd(3);
  m_tpc_adc->Draw("colz");
  m_canvas->cd(4);
  m_tpc_tdc->Draw("colz");

  m_canvas->Update();

  m_is_ready = true;
  return m_is_ready;
}

//______________________________________________________________________________
void
EventDisplayShs::Reset( void )
{
  m_tpc_adc2d->Reset("");
  m_tpc_tdc2d->Reset("");
}

//______________________________________________________________________________
Int_t
EventDisplayShs::GetCommand( void ) const
{
  char ch;
  char data[100];
  static int stat   = 0;
  static int Nevent = 0;
  static int ev     = 0;

  const int kSkip = 1;
  const int kNormal = 0;

  if (stat == 1 && Nevent > 0 && ev<Nevent) {
    ev++;
    return kSkip;
  }
  if (ev==Nevent) {
    stat=0;
    ev=0;
  }

  if (stat == 0) {
    hddaq::cout << std::endl
                << "q: Quit." << std::endl
                << "n: Process Nevent." << std::endl
                << "p: Run TApplication." << std::endl
                << std::endl
                << "q|n|p> " << std::flush;
    std::scanf("%c",&ch);
    if (ch!='\n')
      while(getchar() != '\n');
    switch (ch) {
    case 'q':
      std::exit(EXIT_SUCCESS);
    case 'n':
      stat = 1;
      do {
	printf("event#>");
	scanf("%s",data);
      } while ( (Nevent=atoi(data))<=0 );
      //hddaq::cout << "Continue " << Nevent << "event" << std::endl;
      hddaq::cout << "Skip " << Nevent << "event" << std::endl;
      break;
    case 'p':
      m_theApp->Run(kTRUE);
      break;
    }
  }

  if (stat == 1)
    return kSkip;

  return kNormal;
}

//______________________________________________________________________________
void
EventDisplayShs::Run(Bool_t flag)
{
  hddaq::cout << FUNC_NAME << " TApplication is running" << std::endl;
  m_theApp->Run(flag);
}

//______________________________________________________________________________
void
EventDisplayShs::Update(void)
{
  TIter canvas_iterator(gROOT->GetListOfCanvases());
  while (true) {
    auto canvas = dynamic_cast<TCanvas*>(canvas_iterator.Next());
    if (!canvas) { break; }
    canvas->UseCurrentStyle();
    canvas->Modified();
    canvas->Update();
  }
}
