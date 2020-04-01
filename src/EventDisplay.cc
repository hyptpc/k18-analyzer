/**
 *  file: EventDisplay.cc
 *  date: 2017.04.10
 *  note: Unit is mm.
 *
 */

#include "EventDisplay.hh"

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
#include "MathTools.hh"
#include "DeleteUtility.hh"
#include "HodoParamMan.hh"
#include "DCTdcCalibMan.hh"

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
#define CATCH      1

namespace
{
  const std::string& class_name("EventDisplay");
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
// FBT
//const double p0_FBT = gHodo.GetOffset(DetIdFBT1, 0, 0, 0);
//const double p1_FBT = gHodo.GetGain(DetIdFBT1, 0, 0, 0);
//const double p0_FBT = 526.;
//const double p1_FBT = -1;
  //const double offsetCATCH = 180;
  const double offsetCATCH = 155;
  const double offsetBGO   = 60; // offset from CFT
}

//______________________________________________________________________________
EventDisplay::EventDisplay( void )
  : m_is_ready(false),
    m_theApp(0),
    m_geometry(0),
    m_node(0),
    m_canvas(0),
    m_canvas_vertex(0),
    m_canvas_hist(0),
    m_canvas_hist2(0),
    m_canvas_hist3(0),
    m_canvas_hist4(0),
    m_canvas_hist5(0),
    m_canvas_hist6(0),
    m_canvas_hist7(0),
    m_canvas_hist8(0),
    m_canvas_hist9(0),
    m_canvas_catch(0),
    m_canvas_catch3d(0),
    m_hist_vertex_x(0),
    m_hist_vertex_y(0),
    m_hist_p(0),
    m_hist_m2(0),
    m_hist_missmass(0),
    m_hist_bh1(0),
    m_hist_bft(0),
    m_hist_bft_p(0),
    m_BH1box_cont(0),
    m_hist_bcIn(0),
    m_BH2box_cont(0),
    m_hist_bcOut(0),
    m_hist_bh2(0),
    m_hist_bcOut_sdcIn(0),
    m_hist_sdcIn_predict(0),
    m_hist_sdcIn_predict2(0),
    m_TargetXZ_box2(0),
    m_TargetYZ_box2(0),
    m_hist_sft_x(0),
    m_hist_sft_u(0),
    m_hist_sft_v(0),
    m_hist_sch(0),
    m_hist_tof(0),
    m_hist_sdc1(0),
    m_hist_sdc1p(0),
    m_hist_sdc2_l(0),
    m_hist_sdc2_t(0),
    m_hist_sdc2p_l(0),
    m_hist_sdc2p_t(0),
    m_hist_sdc2y_l(0),
    m_hist_sdc2y_t(0),
    m_hist_sdc2yp_l(0),
    m_hist_sdc2yp_t(0),
    m_hist_sdc3_l(0),
    m_hist_sdc3_t(0),
    m_hist_sdc3p_l(0),
    m_hist_sdc3p_t(0),
    m_hist_sdc3y_l(0),
    m_hist_sdc3y_t(0),
    m_hist_sdc3yp_l(0),
    m_hist_sdc3yp_t(0),
    m_hist_bc3(0),
    m_hist_bc3p(0),
    m_hist_bc3u(0),
    m_hist_bc3up(0),
    m_hist_bc3v(0),
    m_hist_bc3vp(0),
    m_hist_bc4(0),
    m_hist_bc4p(0),
    m_hist_bc4u(0),
    m_hist_bc4up(0),
    m_hist_bc4v(0),
    m_hist_bc4vp(0),
    m_hist_bc3_time(0),
    m_hist_bc3p_time(0),
    m_hist_bc4_time(0),
    m_hist_bc4p_time(0),
    m_hist_fbt1u(0),
    m_hist_fbt1up(0),
    m_hist_fbt2u(0),
    m_hist_fbt2up(0),
    m_hist_fbt1d(0),
    m_hist_fbt1dp(0),
    m_hist_fbt2d(0),
    m_hist_fbt2dp(0),
    m_hist_cft1_l(0),
    m_hist_cft1_t(0),
    m_hist_cft1_hi(0),
    m_hist_cft1_lo(0),
    m_hist_cft2_l(0),
    m_hist_cft2_t(0),
    m_hist_cft2_hi(0),
    m_hist_cft2_lo(0),
    m_hist_cft3_l(0),
    m_hist_cft3_t(0),
    m_hist_cft3_hi(0),
    m_hist_cft3_lo(0),
    m_hist_cft4_l(0),
    m_hist_cft4_t(0),
    m_hist_cft4_hi(0),
    m_hist_cft4_lo(0),
    m_hist_cft5_l(0),
    m_hist_cft5_t(0),
    m_hist_cft5_hi(0),
    m_hist_cft5_lo(0),
    m_hist_cft6_l(0),
    m_hist_cft6_t(0),
    m_hist_cft6_hi(0),
    m_hist_cft6_lo(0),
    m_hist_cft7_l(0),
    m_hist_cft7_t(0),
    m_hist_cft7_hi(0),
    m_hist_cft7_lo(0),
    m_hist_cft8_l(0),
    m_hist_cft8_t(0),
    m_hist_cft8_hi(0),
    m_hist_cft8_lo(0),
    m_hist_bgo(0),
    m_hist_piid_l(0),
    m_hist_piid_t(0),
    m_target_node(0),
    m_kurama_inner_node(0),
    m_kurama_outer_node(0),
    m_BH2wall_node(0),
    m_SCHwall_node(0),
    m_FBT1d1wall_node(0),
    m_FBT1u1wall_node(0),
    m_FBT1d2wall_node(0),
    m_FBT1u2wall_node(0),
    m_FBT2d1wall_node(0),
    m_FBT2u1wall_node(0),
    m_FBT2d2wall_node(0),
    m_FBT2u2wall_node(0),
    m_TOFwall_node(0),
    m_init_step_mark(0),
    m_kurama_step_mark(0),
    m_TargetXZ_box(0),
    m_TargetYZ_box(0),
    m_VertexPointXZ(0),
    m_VertexPointYZ(0),
    m_KuramaMarkVertexX(0),
    m_KuramaMarkVertexY(0),
    m_MissMomXZ_line(0),
    m_MissMomYZ_line(0),
    m_geometry_catch(0),
    m_node_catch(0),
    m_Tgt_Arc(0),
    m_CFRP_Arc(0),
    m_hbase_catch(0)
{
}

//______________________________________________________________________________
EventDisplay::~EventDisplay( void )
{
  del::DeleteObject( m_theApp );
}

//______________________________________________________________________________
// local function
void
ConstructionDone( const char* name, std::ostream& ost=hddaq::cout )
{
  const std::size_t n = 20;
  std::size_t s = std::string(name).size();
  ost << " " << name << " ";
  for( std::size_t i=0; i<n-s; ++i )
    ost << ".";
  ost << " done" << std::endl;
}

//______________________________________________________________________________
bool
EventDisplay::Initialize( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_ready ){
    hddaq::cerr << "#W " << func_name
		<< " already initialied" << std::endl;
    return false;
  }

  m_theApp = new TApplication( "App", 0, 0 );

#if ROOT_VERSION_CODE > ROOT_VERSION(6,4,0)
  gStyle->SetPalette(kCool);
#endif
  gStyle->SetNumberContours(255);

  m_geometry = new TGeometry( "evdisp","K1.8 Event Display" );

  ThreeVector worldSize( 1000., 1000., 1000. ); /*mm*/
  new TBRIK( "world", "world", "void",
	     worldSize.x(), worldSize.y(), worldSize.z() );

  m_node = new TNode( "node", "node", "world", 0., 0., 0. );
  m_geometry->GetNode("node")->SetVisibility(0);

#if BH2
  ConstructBH2();
#endif

  ConstructTarget();

#if KURAMA
  ConstructKURAMA();
#endif

#if BcOut
  ConstructBcOut();
#endif

#if SdcIn
  ConstructSdcIn();
#endif

#if SdcOut
  ConstructSdcOut();
#endif

#if SCH
  ConstructSCH();
#endif

#if FBT
  ConstructFBT();
#endif

#if TOF
  ConstructTOF();
#endif

  m_canvas = new TCanvas( "canvas", "K1.8 Event Display",
			  1000, 1000 );
  m_canvas->Divide(1,2);
  //m_canvas->cd(1)->SetPad( 0.00, 0.92, 1.00, 1.00 );
  //m_canvas->cd(2)->SetPad( 0.00, 0.00, 1.00, 0.92 );
  m_canvas->cd(1)->SetPad( 0.00, 0.72, 1.00, 1.00 );
  m_canvas->cd(2)->SetPad( 0.00, 0.00, 1.00, 0.72 );
  m_canvas->cd(2);

  m_geometry->Draw();

  gPad->SetPhi( 180. );
  gPad->SetTheta( -20. );
  gPad->GetView()->ZoomIn();

  m_canvas->Update();

#if Vertex

  m_canvas_vertex = new TCanvas( "canvas_vertex", "K1.8 Event Display (Vertex)",
				 1000, 800 );
  m_canvas_vertex->Divide(1,2);
  m_canvas_vertex->cd(1);
  gPad->DrawFrame( MinZ, MinX, MaxZ, MaxX, "Vertex XZ projection" );

  m_TargetXZ_box->Draw("L");

  gPad->Update();

  m_canvas_vertex->cd(2);
  gPad->DrawFrame( MinZ, MinY, MaxZ, MaxY, "Vertex YZ projection" );
  m_TargetYZ_box->Draw("L");
  gPad->Update();
#endif

#if Hist
  m_canvas_hist = new TCanvas( "canvas_hist", "EventDisplay Hist",
			       400, 800 );
  m_canvas_hist->Divide(1,3);
  m_hist_p  = new TH1F( "hist1", "Momentum", 100, 0., 3. );
  m_hist_m2 = new TH1F( "hist2", "Mass Square", 200, -0.5, 1.5 );
  m_hist_missmass = new TH1F( "hist3", "Missing Mass", 100, 1, 1.4 );
  m_canvas_hist->cd(1)->SetGrid();
  m_hist_p->Draw();
  m_canvas_hist->cd(2)->SetGrid();
  m_hist_m2->Draw();
  m_canvas_hist->cd(3)->SetGrid();
  m_hist_missmass->Draw();
  // gStyle->SetOptTitle(0);
  gStyle->SetOptStat(1111110);
#endif

#if Hist_Timing
  m_canvas_hist2 = new TCanvas( "canvas_hist2", "EventDisplay Detector Timing",
			       800, 1000 );
  m_canvas_hist2->Divide(3,3);
  m_hist_bh2  = new TH2F( "hist_bh2", "BH2", NumOfSegBH2, 0., NumOfSegBH2, 500, -500, 500 );
  m_hist_bh2->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_sft_x  = new TH2F( "hist_sft_x", "SFT_X", NumOfSegSFT_X, 0., NumOfSegSFT_X, 500, -500, 500 );
  m_hist_sft_x->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_sft_u  = new TH2F( "hist_sft_u", "SFT_U", NumOfSegSFT_UV, 0., NumOfSegSFT_UV, 500, -500, 500 );
  m_hist_sft_u->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_sft_v  = new TH2F( "hist_sft_v", "SFT_V", NumOfSegSFT_UV, 0., NumOfSegSFT_UV, 500, -500, 500 );
  m_hist_sft_v->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_sch = new TH2F( "hist_sch", "SCH", NumOfSegSCH, 0, NumOfSegSCH, 500, -500, 500 );
  m_hist_sch->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_tof = new TH2F( "hist_tof", "TOF", NumOfSegTOF, 0, NumOfSegTOF, 500, -500, 500 );
  m_hist_tof->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_sdc1 = new TH2F( "hist_sdc1", "SDC1", MaxWireSDC1, 0, MaxWireSDC1, 500, -500, 500 );
  m_hist_sdc1->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_sdc1p = new TH2F( "hist_sdc1p", "SDC1 Xp", MaxWireSDC1, 0, MaxWireSDC1, 500, -500, 500 );

  m_hist_sdc1p->SetFillColor(kBlack);
  m_canvas_hist2->cd(1)->SetGrid();
  m_hist_bh2->Draw("box");
  m_canvas_hist2->cd(2)->SetGrid();
  m_hist_sft_x->Draw("box");
  m_canvas_hist2->cd(4)->SetGrid();
  m_hist_sft_u->Draw("box");
  m_canvas_hist2->cd(5)->SetGrid();
  m_hist_sft_v->Draw("box");
  m_canvas_hist2->cd(7)->SetGrid();
  m_hist_sch->Draw("box");
  m_canvas_hist2->cd(8)->SetGrid();
  m_hist_tof->Draw("box");
  m_canvas_hist2->cd(3)->SetGrid();
  m_hist_sdc1->Draw("box");
  m_hist_sdc1p->Draw("samebox");

  m_canvas_hist3 = new TCanvas( "canvas_hist3", "EventDisplay Detector Timing",
			       800, 1000 );
  m_canvas_hist3->Divide(3,2);

  m_hist_bc3 = new TH2F( "hist_bc3", "BC3 X", MaxWireBC3, 0, MaxWireBC3, 500, -500, 500 );
  m_hist_bc3p = new TH2F( "hist_bc3p", "BC3 Xp", MaxWireBC3, 0, MaxWireBC3, 500, -500, 500 );
  m_hist_bc3p->SetFillColor(kBlack);
  m_hist_bc3->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bc3p->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_bc3_time = new TH1F( "hist_bc3_time", "BC3 X", 500, -500, 500 );
  m_hist_bc3p_time = new TH1F( "hist_bc3p_time", "BC3 X", 500, -500, 500 );
  m_hist_bc3p_time->SetFillStyle(3001);
  m_hist_bc3p_time->SetFillColor(kBlack);

  m_hist_bc3u = new TH2F( "hist_bc3u", "BC3 U", MaxWireBC3, 0, MaxWireBC3, 500, -500, 500 );
  m_hist_bc3up = new TH2F( "hist_bc3up", "BC3 Up", MaxWireBC3, 0, MaxWireBC3, 500, -500, 500 );
  m_hist_bc3up->SetFillColor(kBlack);
  m_hist_bc3u->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bc3up->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_bc3v = new TH2F( "hist_bc3v", "BC3 V", MaxWireBC3, 0, MaxWireBC3, 500, -500, 500 );
  m_hist_bc3vp = new TH2F( "hist_bc3vp", "BC3 Vp", MaxWireBC3, 0, MaxWireBC3, 500, -500, 500 );
  m_hist_bc3vp->SetFillColor(kBlack);
  m_hist_bc3v->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bc3vp->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_bc4 = new TH2F( "hist_bc4", "BC4 X", MaxWireBC4, 0, MaxWireBC4, 500, -500, 500 );
  m_hist_bc4p = new TH2F( "hist_bc4p", "BC4 Xp", MaxWireBC4, 0, MaxWireBC4, 500, -500, 500 );
  m_hist_bc4p->SetFillColor(kBlack);
  m_hist_bc4->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bc4p->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_bc4_time = new TH1F( "hist_bc4_time", "BC4 X", 500, -500, 500 );
  m_hist_bc4p_time = new TH1F( "hist_bc4p_time", "BC4 X", 500, -500, 500 );
  m_hist_bc4p_time->SetFillStyle(3001);
  m_hist_bc4p_time->SetFillColor(kBlack);

  m_hist_bc4u = new TH2F( "hist_bc4u", "BC4 U", MaxWireBC4, 0, MaxWireBC4, 500, -500, 500 );
  m_hist_bc4up = new TH2F( "hist_bc4up", "BC4 Up", MaxWireBC4, 0, MaxWireBC4, 500, -500, 500 );
  m_hist_bc4up->SetFillColor(kBlack);
  m_hist_bc4u->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bc4up->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_bc4v = new TH2F( "hist_bc4v", "BC4 V", MaxWireBC4, 0, MaxWireBC4, 500, -500, 500 );
  m_hist_bc4vp = new TH2F( "hist_bc4vp", "BC4 Vp", MaxWireBC4, 0, MaxWireBC4, 500, -500, 500 );
  m_hist_bc4vp->SetFillColor(kBlack);
  m_hist_bc4v->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bc4vp->GetYaxis()->SetRangeUser(-100, 100);

  m_canvas_hist3->cd(1)->SetGrid();
  m_hist_bc3->Draw("box");
  m_hist_bc3p->Draw("samebox");
  m_canvas_hist3->cd(2)->SetGrid();
  m_hist_bc3v->Draw("box");
  m_hist_bc3vp->Draw("samebox");
  m_canvas_hist3->cd(3)->SetGrid();
  m_hist_bc3u->Draw("box");
  m_hist_bc3up->Draw("samebox");
  m_canvas_hist3->cd(4)->SetGrid();
  m_hist_bc4->Draw("box");
  m_hist_bc4p->Draw("samebox");
  m_canvas_hist3->cd(5)->SetGrid();
  m_hist_bc4v->Draw("box");
  m_hist_bc4vp->Draw("samebox");
  m_canvas_hist3->cd(6)->SetGrid();
  m_hist_bc4u->Draw("box");
  m_hist_bc4up->Draw("samebox");
  /*
  m_canvas_hist3->cd(3)->SetGrid();
  m_hist_bc3_time->Draw("");
  m_hist_bc3p_time->Draw("same");
  m_canvas_hist3->cd(4)->SetGrid();
  m_hist_bc4_time->Draw();
  m_hist_bc4p_time->Draw("same");
  */
  // gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

#endif

#if Hist_SdcOut
  m_canvas_hist4 = new TCanvas( "canvas_hist4", "EventDisplay Detector Timing (SdcOut)", 800, 800);
  m_canvas_hist4->Divide(2,4);  

  m_hist_fbt1u = new TH2F( "hist_fbt1u", "FBT1U", MaxSegFBT1, 0, MaxSegFBT1, 500, -500, 500 );
  m_hist_fbt1up = new TH2F( "hist_fbt1u", "FBT1U p", MaxSegFBT1, 0, MaxSegFBT1, 500, -500, 500 );
  m_hist_fbt1up->SetFillColor(kBlack);

  m_hist_fbt1u->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_fbt1up->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_fbt2u = new TH2F( "hist_fbt2u", "FBT2U", MaxSegFBT2, 0, MaxSegFBT2, 500, -500, 500 );
  m_hist_fbt2up = new TH2F( "hist_fbt2u", "FBT2U p", MaxSegFBT2, 0, MaxSegFBT2, 500, -500, 500 );
  m_hist_fbt2up->SetFillColor(kBlack);

  m_hist_fbt2u->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_fbt2up->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_fbt1d = new TH2F( "hist_fbt1d", "FBT1D", MaxSegFBT1, 0, MaxSegFBT1, 500, -500, 500 );
  m_hist_fbt1dp = new TH2F( "hist_fbt1d", "FBT1D p", MaxSegFBT1, 0, MaxSegFBT1, 500, -500, 500 );
  m_hist_fbt1dp->SetFillColor(kBlack);

  m_hist_fbt1d->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_fbt1dp->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_fbt2d = new TH2F( "hist_fbt2d", "FBT2D", MaxSegFBT2, 0, MaxSegFBT2, 500, -500, 500 );
  m_hist_fbt2dp = new TH2F( "hist_fbt2d", "FBT2D p", MaxSegFBT2, 0, MaxSegFBT2, 500, -500, 500 );
  m_hist_fbt2dp->SetFillColor(kBlack);

  m_hist_fbt2d->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_fbt2dp->GetYaxis()->SetRangeUser(-100, 100);

  m_canvas_hist4->cd(1)->SetGrid();
  m_hist_fbt1u->Draw("box");
  m_canvas_hist4->cd(1)->SetGrid();
  m_hist_fbt1up->Draw("samebox");

  m_canvas_hist4->cd(2)->SetGrid();
  m_hist_fbt2u->Draw("box");
  m_canvas_hist4->cd(2)->SetGrid();
  m_hist_fbt2up->Draw("samebox");

  m_canvas_hist4->cd(3)->SetGrid();
  m_hist_fbt1d->Draw("box");
  m_canvas_hist4->cd(3)->SetGrid();
  m_hist_fbt1dp->Draw("samebox");

  m_canvas_hist4->cd(4)->SetGrid();
  m_hist_fbt2d->Draw("box");
  m_canvas_hist4->cd(4)->SetGrid();
  m_hist_fbt2dp->Draw("samebox");

  m_hist_sdc2_l = new TH2F( "hist_sdc2_l", "SDC2 (leading)", MaxWireSDC2, 0, MaxWireSDC2, 500, -500, 500 );
  m_hist_sdc2_t = new TH2F( "hist_sdc2_t", "SDC2 (trailing)", MaxWireSDC2, 0, MaxWireSDC2, 500, -500, 500 );
  m_hist_sdc2_t->SetFillColor(kRed);

  m_hist_sdc2p_l = new TH2F( "hist_sdc2p_l", "SDC2 Xp (leading)", MaxWireSDC2, 0, MaxWireSDC2, 500, -500, 500 );
  m_hist_sdc2p_t = new TH2F( "hist_sdc2p_t", "SDC2 Xp(trailing)", MaxWireSDC2, 0, MaxWireSDC2, 500, -500, 500 );
  m_hist_sdc2p_l->SetFillColor(kBlack);
  m_hist_sdc2p_t->SetFillColor(kGreen);

  m_hist_sdc2y_l = new TH2F( "hist_sdc2y_l", "SDC2 Y (leading)", MaxWireSDC2, 0, MaxWireSDC2, 500, -500, 500 );
  m_hist_sdc2y_t = new TH2F( "hist_sdc2y_t", "SDC2 Y (trailing)", MaxWireSDC2, 0, MaxWireSDC2, 500, -500, 500 );
  m_hist_sdc2y_t->SetFillColor(kRed);

  m_hist_sdc2yp_l = new TH2F( "hist_sdc2yp_l", "SDC2 Yp (leading)", MaxWireSDC2, 0, MaxWireSDC2, 500, -500, 500 );
  m_hist_sdc2yp_t = new TH2F( "hist_sdc2yp_t", "SDC2 Yp(trailing)", MaxWireSDC2, 0, MaxWireSDC2, 500, -500, 500 );
  m_hist_sdc2yp_l->SetFillColor(kBlack);
  m_hist_sdc2yp_t->SetFillColor(kGreen);


  m_hist_sdc3_l = new TH2F( "hist_sdc3_l", "SDC3 (leading)", MaxWireSDC3X, 0, MaxWireSDC3X, 500, -500, 500 );
  m_hist_sdc3_t = new TH2F( "hist_sdc3_t", "SDC3 (trailing)", MaxWireSDC3X, 0, MaxWireSDC3X, 500, -500, 500 );
  m_hist_sdc3_t->SetFillColor(kRed);

  m_hist_sdc3p_l = new TH2F( "hist_sdc3p_l", "SDC3 Xp(leading)", MaxWireSDC3X, 0, MaxWireSDC3X, 500, -500, 500 );
  m_hist_sdc3p_t = new TH2F( "hist_sdc3p_t", "SDC3 Xp(trailing)", MaxWireSDC3X, 0, MaxWireSDC3X, 500, -500, 500 );
  m_hist_sdc3p_l->SetFillColor(kBlack);
  m_hist_sdc3p_t->SetFillColor(kGreen);

  m_hist_sdc3y_l = new TH2F( "hist_sdc3y_l", "SDC3 Y (leading)", MaxWireSDC3Y, 0, MaxWireSDC3Y, 500, -500, 500 );
  m_hist_sdc3y_t = new TH2F( "hist_sdc3y_t", "SDC3 Y (trailing)", MaxWireSDC3Y, 0, MaxWireSDC3Y, 500, -500, 500 );
  m_hist_sdc3y_t->SetFillColor(kRed);

  m_hist_sdc3yp_l = new TH2F( "hist_sdc3yp_l", "SDC3 Yp(leading)", MaxWireSDC3Y, 0, MaxWireSDC3Y, 500, -500, 500 );
  m_hist_sdc3yp_t = new TH2F( "hist_sdc3yp_t", "SDC3 Yp(trailing)", MaxWireSDC3Y, 0, MaxWireSDC3Y, 500, -500, 500 );
  m_hist_sdc3yp_l->SetFillColor(kBlack);
  m_hist_sdc3yp_t->SetFillColor(kGreen);

  m_canvas_hist4->cd(5)->SetGrid();
  m_hist_sdc2_l->Draw("box");
  m_hist_sdc2_t->Draw("samebox");
  m_hist_sdc2p_l->Draw("samebox");
  m_hist_sdc2p_t->Draw("samebox");

  m_canvas_hist4->cd(6)->SetGrid();
  m_hist_sdc2y_l->Draw("box");
  m_hist_sdc2y_t->Draw("samebox");
  m_hist_sdc2yp_l->Draw("samebox");
  m_hist_sdc2yp_t->Draw("samebox");

  m_canvas_hist4->cd(7)->SetGrid();
  m_hist_sdc3_l->Draw("box");
  m_hist_sdc3_t->Draw("samebox");
  m_hist_sdc3p_l->Draw("samebox");
  m_hist_sdc3p_t->Draw("samebox");

  m_canvas_hist4->cd(8)->SetGrid();
  m_hist_sdc3y_l->Draw("box");
  m_hist_sdc3y_t->Draw("samebox");
  m_hist_sdc3yp_l->Draw("samebox");
  m_hist_sdc3yp_t->Draw("samebox");


#endif

#if Hist_BcIn
  m_canvas_hist5 = new TCanvas( "canvas_hist5", "EventDisplay Detector Timing (BcIn)",
			       800, 1000 );
  m_canvas_hist5->Divide(2,2);
  m_hist_bh1  = new TH2F( "hist_bh1", "BH1", NumOfSegBH1, 0., NumOfSegBH1, 500, -500, 500 );
  m_hist_bh1->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bft  = new TH2F( "hist_bft", "BFT", NumOfSegBFT, 0., NumOfSegBFT, 500, -500, 500 );
  m_hist_bft_p  = new TH2F( "hist_bft_p", "BFT prime", NumOfSegBFT, 0., NumOfSegBFT, 500, -500, 500 );
  m_hist_bft->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bft_p->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bft_p->SetFillColor(kBlack);

  m_hist_bcIn  = new TH2F( "hist_bcIn", "BcIn Tracking", 200, -100, 100, 200, -150, 50 );

  m_hist_bcOut  = new TH2F( "hist_bcOut", "BcOut Tracking", 200, -100, 100, 200, 0, 600 );


  m_canvas_hist5->cd(1)->SetGrid();
  m_hist_bh1->Draw("box");
  m_canvas_hist5->cd(2)->SetGrid();
  m_hist_bft->Draw("box");
  m_hist_bft_p->Draw("samebox");
  m_canvas_hist5->cd(3)->SetGrid();
  m_hist_bcIn->Draw("box");

  double Bh1SegX[NumOfSegBH1] = {30./2., 20./2., 16./2., 12./2., 8./2., 8./2., 8./2., 12./2., 16./2., 20./2., 30./2.};
  double Bh1SegY[NumOfSegBH1] = {5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2.};

  //double localPosBh1Z = -515.;
  double localPosBh1Z = gGeom.GetLocalZ( IdBH1 );
  double localPosBh1X_dX = 0.;
  double localPosBh1X[NumOfSegBH1] = {-70. + localPosBh1X_dX,
                                      -46. + localPosBh1X_dX,
                                      -29. + localPosBh1X_dX,
                                      -16. + localPosBh1X_dX,
                                      -7. + localPosBh1X_dX,
                                      0. + localPosBh1X_dX,
                                      7. + localPosBh1X_dX,
                                      16. + localPosBh1X_dX,
                                      29. + localPosBh1X_dX,
                                      46. + localPosBh1X_dX,
                                      70. + localPosBh1X_dX};
  double localPosBh1_dZ[NumOfSegBH1] = {4.5, -4.5, 4.5, -4.5, 4.5, -4.5, 4.5, -4.5, 4.5, -4.5, 4.5};

  for (int i=0; i<NumOfSegBH1; i++) {
    m_BH1box_cont.push_back( new TBox( localPosBh1X[i]-Bh1SegX[i], 
				       localPosBh1Z+localPosBh1_dZ[i]-Bh1SegY[i], 
				       localPosBh1X[i]+Bh1SegX[i],
				       localPosBh1Z+localPosBh1_dZ[i]+Bh1SegY[i])); 
  }
  for (int i=0; i<NumOfSegBH1; i++) {
    m_BH1box_cont[i]->SetFillColor(kWhite);
    m_BH1box_cont[i]->SetLineColor(kBlack);
    m_BH1box_cont[i]->Draw("L");
  }

  m_canvas_hist5->cd(4)->SetGrid();
  m_hist_bcOut->Draw("box");

  double Bh2SegX[NumOfSegBH2] = {35./2., 10./2., 7./2., 7./2., 7./2., 7./2., 10./2., 35./2.};
  double Bh2SegY[NumOfSegBH2] = {5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2.};

  double localPosBh2X[NumOfSegBH2] = {-41.5, -19.0, -10.5, -3.5, 3.5, 10.5, 19.0, 41.5};

  double localPosBh2Z = gGeom.GetLocalZ( IdBH2 );

  for (int i=0; i<NumOfSegBH2; i++) {
    m_BH2box_cont.push_back( new TBox( localPosBh2X[i]-Bh2SegX[i], 
				       localPosBh2Z-Bh2SegY[i], 
				       localPosBh2X[i]+Bh2SegX[i],
				       localPosBh2Z+Bh2SegY[i])); 
  }
  for (int i=0; i<NumOfSegBH2; i++) {
    m_BH2box_cont[i]->SetFillColor(kWhite);
    m_BH2box_cont[i]->SetLineColor(kBlack);
    m_BH2box_cont[i]->Draw("L");
  }

  gStyle->SetOptStat(0);


  m_canvas_hist6 = new TCanvas( "canvas_hist6", "EventDisplay Detector Timing (BcOut SdcIn)",
				800, 1000 );

  m_hist_sdcIn_predict  = new TH2F( "hist_sdcIn_predict", "BcOut-SdcIn Tracking", 200, -400, 400, 750, -3000, 0 );
  m_hist_sdcIn_predict->SetFillColor(kRed);
  m_hist_sdcIn_predict->SetLineColor(kRed);
  m_hist_sdcIn_predict->Draw("box");

  m_hist_sdcIn_predict2  = new TH2F( "hist_sdcIn_predict2", "BcOut-SdcIn Tracking", 200, -400, 400, 750, -3000, 0 );
  m_hist_sdcIn_predict2->SetFillColor(kMagenta);
  m_hist_sdcIn_predict2->SetLineColor(kMagenta);
  m_hist_sdcIn_predict2->Draw("samebox");

  m_hist_bcOut_sdcIn  = new TH2F( "hist_bcOut_sdcIn", "BcOut-SdcIn Tracking", 400, -400, 400, 1500, -3000, 0 );
  m_canvas_hist6->SetGrid();
  m_hist_bcOut_sdcIn->Draw("samebox");

  double globalPosBh2Z = gGeom.GetGlobalPosition( IdBH2 ).z();
  double globalPosBh2X = gGeom.GetGlobalPosition( IdBH2 ).x();

  for (int i=0; i<NumOfSegBH2; i++) {
    m_BH2box_cont2.push_back( new TBox( globalPosBh2X+localPosBh2X[i]-Bh2SegX[i], 
					globalPosBh2Z-Bh2SegY[i]-20, 
					globalPosBh2X+localPosBh2X[i]+Bh2SegX[i],
					globalPosBh2Z+Bh2SegY[i]+20)); 
  }
  for (int i=0; i<NumOfSegBH2; i++) {
    m_BH2box_cont2[i]->SetFillColor(kWhite);
    m_BH2box_cont2[i]->SetLineColor(kBlack);
    m_BH2box_cont2[i]->Draw("L");
  }

  double globalPosSchZ = gGeom.GetGlobalPosition( IdSCH ).z();

  double SchSegX = 11.5/2.;
  double SchSegY = 2.0/2.;
  for (int i=1; i<=NumOfSegSCH; i++) {
    double globalPosSchX = gGeom.CalcWirePosition( IdSCH, (double)i );
    m_SCHbox_cont.push_back( new TBox( globalPosSchX-SchSegX, 
				       globalPosSchZ-SchSegY-20, 
				       globalPosSchX+SchSegX,
				       globalPosSchZ+SchSegY+20)); 
  }
  for (int i=0; i<NumOfSegSCH; i++) {
    m_SCHbox_cont[i]->SetFillColor(kWhite);
    m_SCHbox_cont[i]->SetLineColor(kBlack);
    m_SCHbox_cont[i]->Draw("L");
  }


  double globalPosTarget_x = gGeom.GetGlobalPosition( IdTarget ).x();
  double globalPosTarget_y = gGeom.GetGlobalPosition( IdTarget ).y();
  double globalPosTarget_z = gGeom.GetGlobalPosition( IdTarget ).z();
  double target_r = 20.;
  double target_z = 300./2;
  m_TargetXZ_box2 =  new TBox( globalPosTarget_x-target_r, 
			       globalPosTarget_z-target_z, 
			       globalPosTarget_x+target_r, 
			       globalPosTarget_z+target_z); 

  m_TargetXZ_box2->SetFillColor(kWhite);
  m_TargetXZ_box2->SetLineColor(kBlack);
  m_TargetXZ_box2->Draw("L");

  m_TargetYZ_box2 =  new TBox( globalPosTarget_y-target_r, 
			       globalPosTarget_z-target_z, 
			       globalPosTarget_y+target_r, 
			       globalPosTarget_z+target_z); 

  m_TargetYZ_box2->SetFillColor(kWhite);
  m_TargetYZ_box2->SetLineColor(kBlack);
  m_TargetYZ_box2->Draw("L");

#endif

#if CATCH
  m_canvas_hist7 = new TCanvas( "canvas_hist7", "EventDisplay Detector Timing (CATCH)", 800, 1000 );
  m_canvas_hist7->Divide(4,3);  


  for (int layer=0; layer<8; layer++) {
    TH2 *hp_l=0, *hp_t=0;
    hp_l = new TH2F( Form( "hist_cft%d_l", layer+1 ),
		      Form( "CFT layer%d (Leading)", layer+1 ),
		      NumOfSegCFT[layer], 0, NumOfSegCFT[layer], 500, -500, 500 );
    hp_l->GetYaxis()->SetRangeUser(-100, 100);
    hp_t = new TH2F( Form( "hist_cft%d_t", layer+1 ),
		      Form( "CFT layer%d (Trailing)", layer+1 ),
		      NumOfSegCFT[layer], 0, NumOfSegCFT[layer], 500, -500, 500 );
    hp_t->GetYaxis()->SetRangeUser(-100, 100);
    hp_t->SetFillColor(kRed);

    m_canvas_hist7->cd(layer+1)->SetGrid();
    hp_l->Draw("box");
    hp_t->Draw("samebox");

    if (layer == 0) {
      m_hist_cft1_l = hp_l;
      m_hist_cft1_t = hp_t;
    } else if (layer == 1) {
      m_hist_cft2_l = hp_l;
      m_hist_cft2_t = hp_t;
    } else if (layer == 2) {
      m_hist_cft3_l = hp_l;
      m_hist_cft3_t = hp_t;
    } else if (layer == 3) {
      m_hist_cft4_l = hp_l;
      m_hist_cft4_t = hp_t;
    } else if (layer == 4) {
      m_hist_cft5_l = hp_l;
      m_hist_cft5_t = hp_t;
    } else if (layer == 5) {
      m_hist_cft6_l = hp_l;
      m_hist_cft6_t = hp_t;
    } else if (layer == 6) {
      m_hist_cft7_l = hp_l;
      m_hist_cft7_t = hp_t;
    } else if (layer == 7) {
      m_hist_cft8_l = hp_l;
      m_hist_cft8_t = hp_t;
    }

  }

  m_hist_bgo = new TH2F( "hist_bgo","BGO (Leading)",
		   NumOfSegBGO, 0, NumOfSegBGO, 500, -500, 500 );
  m_canvas_hist7->cd(9)->SetGrid();
  m_hist_bgo->Draw("box");

  m_hist_piid_l = new TH2F( "hist_piid_l","PiId (Leading)",
			   NumOfSegPiID, 0, NumOfSegPiID, 500, -500, 500 );
  m_hist_piid_t = new TH2F( "hist_piid_t","PiId (Leading)",
			   NumOfSegPiID, 0, NumOfSegPiID, 500, -500, 500 );
  m_canvas_hist7->cd(10)->SetGrid();
  m_hist_piid_l->Draw("box");
  m_hist_piid_t->Draw("samebox");


  m_canvas_hist8 = new TCanvas( "canvas_hist8", "EventDisplay Detector ADC (CATCH)", 800, 1000 );
  m_canvas_hist8->Divide(4,3);  

  for (int layer=0; layer<8; layer++) {
    TH2 *hp_hi=0, *hp_lo=0;
    hp_hi = new TH2F( Form( "hist_cft%d_hi", layer+1 ),
		      Form( "CFT layer%d (High Gain)", layer+1 ),
		      NumOfSegCFT[layer], 0, NumOfSegCFT[layer], 500, -50, 3450 );
    hp_lo = new TH2F( Form( "hist_cft%d_low", layer+1 ),
		      Form( "CFT layer%d (Low Gain)", layer+1 ),
		      NumOfSegCFT[layer], 0, NumOfSegCFT[layer], 500, -50, 3450 );
    hp_lo->SetFillColor(kRed);

    m_canvas_hist8->cd(layer+1)->SetGrid();
    hp_hi->Draw("box");
    hp_lo->Draw("samebox");

    if (layer == 0) {
      m_hist_cft1_hi = hp_hi;
      m_hist_cft1_lo = hp_lo;
    } else if (layer == 1) {
      m_hist_cft2_hi = hp_hi;
      m_hist_cft2_lo = hp_lo;
    } else if (layer == 2) {
      m_hist_cft3_hi = hp_hi;
      m_hist_cft3_lo = hp_lo;
    } else if (layer == 3) {
      m_hist_cft4_hi = hp_hi;
      m_hist_cft4_lo = hp_lo;
    } else if (layer == 4) {
      m_hist_cft5_hi = hp_hi;
      m_hist_cft5_lo = hp_lo;
    } else if (layer == 5) {
      m_hist_cft6_hi = hp_hi;
      m_hist_cft6_lo = hp_lo;
    } else if (layer == 6) {
      m_hist_cft7_hi = hp_hi;
      m_hist_cft7_lo = hp_lo;
    } else if (layer == 7) {
      m_hist_cft8_hi = hp_hi;
      m_hist_cft8_lo = hp_lo;
    }

  }

  m_canvas_hist9 = new TCanvas( "canvas_hist9", "EventDisplay Detector BGO waveform (CATCH)", 1000, 1000 );
  m_canvas_hist9->Divide(4,6);  
  for (int npad = 0; npad<24; npad++)
    m_canvas_hist9->cd(npad+1)->SetGrid();

  ConstructCATCH();
  ConstructCATCH3d();

#endif


  ResetVisibility();

  m_is_ready = true;
  return m_is_ready;
}

//______________________________________________________________________________
bool
EventDisplay::ConstructBH2( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  const int lid = gGeom.GetDetectorId("BH2");

  double rotMatBH2[9] = {};
  double BH2wallX = 120.0/2.; // X
  double BH2wallY =   6.0/2.; // Z
  double BH2wallZ =  40.0/2.; // Y
  double BH2SizeX[NumOfSegBH2] = { 120./2. }; // X
  double BH2SizeY[NumOfSegBH2] = {   6./2. }; // Z
  double BH2SizeZ[NumOfSegBH2] = {  40./2. }; // Y
  double BH2PosX[NumOfSegBH2]  = { 0./2. };
  double BH2PosY[NumOfSegBH2]  = { 0./2. };
  double BH2PosZ[NumOfSegBH2]  = { 0./2. };

  CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		 gGeom.GetRotAngle1( lid ),
		 gGeom.GetRotAngle2( lid ),
		 rotMatBH2 );

  new TRotMatrix( "rotBH2", "rotBH2", rotMatBH2 );
  const ThreeVector& BH2wallPos = gGeom.GetGlobalPosition( lid );
  new TBRIK( "BH2wall_brik", "BH2wall_brik", "void",
	     BH2wallX, BH2wallY, BH2wallZ );
  m_BH2wall_node = new TNode( "BH2wall_node", "BH2wall_node", "BH2wall_brik",
			      BH2wallPos.x(),
			      BH2wallPos.y(),
			      BH2wallPos.z(), "rotBH2", "void" );
  m_BH2wall_node->SetVisibility(0);
  m_BH2wall_node->cd();

  for( int i=0; i<NumOfSegBH2; ++i ){
    new TBRIK( Form( "BH2seg_brik_%d", i ),
	       Form( "BH2seg_brik_%d", i ),
	       "void", BH2SizeX[i], BH2SizeY[i], BH2SizeZ[i]);
    m_BH2seg_node.push_back( new TNode( Form( "BH2seg_node_%d", i ),
					Form( "BH2seg_node_%d", i ),
					Form( "BH2seg_brik_%d", i ),
					BH2PosX[i], BH2PosY[i], BH2PosZ[i] ) );
  }
  m_node->cd();
  ConstructionDone(__func__);
  return true;
}

//______________________________________________________________________________
bool
EventDisplay::ConstructKURAMA( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  double Matrix[9] = {};

  double inner_x = 1400.0/2.; // X
  double inner_y =  800.0/2.; // Z
  double inner_z =  800.0/2.; // Y

  double outer_x = 2200.0/2.; // X
  double outer_y =  800.0/2.; // Z
  double outer_z = 1540.0/2.; // Y

  double uguard_inner_x = 1600.0/2.; // X
  double uguard_inner_y =  100.0/2.; // Z
  double uguard_inner_z = 1940.0/2.; // Y
  double uguard_outer_x =  600.0/2.; // X
  double uguard_outer_y =  100.0/2.; // Z
  double uguard_outer_z =  300.0/2.; // Y

  double dguard_inner_x = 1600.0/2.; // X
  double dguard_inner_y =  100.0/2.; // Z
  double dguard_inner_z = 1940.0/2.; // Y
  double dguard_outer_x = 1100.0/2.; // X
  double dguard_outer_y =  100.0/2.; // Z
  double dguard_outer_z = 1100.0/2.; // Y

  CalcRotMatrix( 0., 0., 0., Matrix );

  new TRotMatrix( "rotKURAMA", "rotKURAMA", Matrix );

  new TBRIK( "kurama_inner_brik", "kurama_inner_brik",
	     "void", inner_x, inner_y, inner_z );

  new TBRIK( "kurama_outer_brik", "kurama_outer_brik",
	     "void", outer_x, outer_y, outer_z );

  new TBRIK( "uguard_inner_brik", "uguard_inner_brik",
	     "void", uguard_inner_x, uguard_inner_y, uguard_inner_z );

  new TBRIK( "uguard_outer_brik", "uguard_outer_brik",
	     "void", uguard_outer_x, uguard_outer_y, uguard_outer_z );

  new TBRIK( "dguard_inner_brik", "dguard_inner_brik",
	     "void", dguard_inner_x, dguard_inner_y, dguard_inner_z );

  new TBRIK( "dguard_outer_brik", "dguard_outer_brik",
	     "void", dguard_outer_x, dguard_outer_y, dguard_outer_z );


  m_kurama_inner_node = new TNode( "kurama_inner_node",
				   "kurama_inner_node",
				   "kurama_inner_brik",
				   0., 0., 0., "rotKURAMA", "void" );
  m_kurama_outer_node = new TNode( "kurama_outer_node",
				   "kurama_outer_node",
				   "kurama_outer_brik",
				   0., 0., 0., "rotKURAMA", "void" );

  TNode *uguard_inner = new TNode( "uguard_inner_node",
				   "uguard_inner_node",
				   "uguard_inner_brik",
				   0., 0., -820.,
				   "rotKURAMA", "void" );
  TNode *uguard_outer = new TNode( "uguard_outer_node",
				   "uguard_outer_node",
				   "uguard_outer_brik",
				   0., 0., -820.,
				   "rotKURAMA", "void" );

  TNode *dguard_inner = new TNode( "dguard_inner_node",
				   "dguard_inner_node",
				   "dguard_inner_brik",
				   0., 0., 820.,
				   "rotKURAMA", "void" );
  TNode *dguard_outer = new TNode( "dguard_outer_node",
				   "dguard_outer_node",
				   "dguard_outer_brik",
				   0., 0., 820.,
				   "rotKURAMA", "void" );

  const Color_t color = kBlack;

  m_kurama_inner_node->SetLineColor( color );
  m_kurama_outer_node->SetLineColor( color );
  uguard_inner->SetLineColor( color );
  uguard_outer->SetLineColor( color );
  dguard_inner->SetLineColor( color );
  dguard_outer->SetLineColor( color );

  m_node->cd();
  ConstructionDone(__func__);
  return true;
}

//______________________________________________________________________________
bool
EventDisplay::ConstructBcOut( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  const double wireL = 200.0;

  const double offsetZ = -3068.4;

  // BC3 X1
  {
    const int lid = gGeom.GetDetectorId("BC3-X1");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotBC3X1", "rotBC3X1", Matrix );
    new TTUBE( "BC3X1Tube", "BC3X1Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<=MaxWireBC3; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_BC3x1_node.push_back( new TNode( Form( "BC3x1_Node_%d", wire ),
					 Form( "BC3x1_Node_%d", wire ),
					 "BC3X1Tube",
					 localPos+BeamAxis,
					 wireGlobalPos.y(),
					 wireGlobalPos.z()+offsetZ,
					 "rotBC3X1", "void" ) );
    }
  }

  // BC3 X2
  {
    const int lid = gGeom.GetDetectorId("BC3-X2");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotBC3X2", "rotBC3X2", Matrix );
    new TTUBE( "BC3X2Tube", "BC3X2Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<=MaxWireBC3; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_BC3x2_node.push_back( new TNode( Form( "BC3x2_Node_%d", wire ),
					 Form( "BC3x2_Node_%d", wire ),
					 "BC3X2Tube",
					 localPos+BeamAxis,
					 wireGlobalPos.y(),
					 wireGlobalPos.z()+offsetZ,
					 "rotBC3X2", "void" ) );
    }
  }

  // BC3 V1
  {
    const int lid = gGeom.GetDetectorId("BC3-V1");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotBC3V1", "rotBC3V1", Matrix );
    new TTUBE( "BC3V1Tube", "BC3V1Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<=MaxWireBC3; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_BC3v1_node.push_back( new TNode( Form( "BC3v1_Node_%d", wire ),
					 Form( "BC3v1_Node_%d", wire ),
					 "BC3V1Tube",
					 localPos+BeamAxis,
					 wireGlobalPos.y(),
					 wireGlobalPos.z()+offsetZ,
					 "rotBC3V1", "void" ) );
    }
  }

  // BC3 V2
  {
    const int lid = gGeom.GetDetectorId("BC3-V2");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotBC3V2", "rotBC3V2", Matrix );
    new TTUBE( "BC3V2Tube", "BC3V2Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<=MaxWireBC3; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_BC3v2_node.push_back( new TNode( Form( "BC3v2_Node_%d", wire ),
					 Form( "BC3v2_Node_%d", wire ),
					 "BC3V2Tube",
					 localPos+BeamAxis,
					 wireGlobalPos.y(),
					 wireGlobalPos.z()+offsetZ,
					 "rotBC3V2", "void" ) );
    }
  }

  // BC3 U1
  {
    const int lid = gGeom.GetDetectorId("BC3-U1");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotBC3U1", "rotBC3U1", Matrix );
    new TTUBE( "BC3U1Tube", "BC3U1Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<=MaxWireBC3; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_BC3u1_node.push_back( new TNode( Form( "BC3u1_Node_%d", wire ),
					 Form( "BC3u1_Node_%d", wire ),
					 "BC3U1Tube",
					 localPos+BeamAxis,
					 wireGlobalPos.y(),
					 wireGlobalPos.z()+offsetZ,
					 "rotBC3U1", "void" ) );
    }
  }

  // BC3 U2
  {
    const int lid = gGeom.GetDetectorId("BC3-U2");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotBC3U2", "rotBC3U2", Matrix );
    new TTUBE( "BC3U2Tube", "BC3U2Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<=MaxWireBC3; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_BC3u2_node.push_back( new TNode( Form( "BC3u2_Node_%d", wire ),
					 Form( "BC3u2_Node_%d", wire ),
					 "BC3U2Tube",
					 localPos+BeamAxis,
					 wireGlobalPos.y(),
					 wireGlobalPos.z()+offsetZ,
					 "rotBC3U2", "void" ) );
    }
  }

  // BC4 U1
  {
    const int lid = gGeom.GetDetectorId("BC4-U1");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotBC4U1", "rotBC4U1", Matrix );
    new TTUBE( "BC4U1Tube", "BC4U1Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<=MaxWireBC4; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_BC4u1_node.push_back( new TNode( Form( "BC4u1_Node_%d", wire ),
					 Form( "BC4u1_Node_%d", wire ),
					 "BC4U1Tube",
					 localPos+BeamAxis,
					 wireGlobalPos.y(),
					 wireGlobalPos.z()+offsetZ,
					 "rotBC4U1", "void" ) );
    }
  }

  // BC4 U2
  {
    const int lid = gGeom.GetDetectorId("BC4-U2");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotBC4U2", "rotBC4U2", Matrix );
    new TTUBE( "BC4U2Tube", "BC4U2Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<=MaxWireBC4; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_BC4u2_node.push_back( new TNode( Form( "BC4u2_Node_%d", wire ),
					 Form( "BC4u2_Node_%d", wire ),
					 "BC4U2Tube",
					 localPos+BeamAxis,
					 wireGlobalPos.y(),
					 wireGlobalPos.z()+offsetZ,
					 "rotBC4U2", "void" ) );
    }
  }

  // BC4 V1
  {
    const int lid = gGeom.GetDetectorId("BC4-V1");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotBC4V1", "rotBC4V1", Matrix );
    new TTUBE( "BC4V1Tube", "BC4V1Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<=MaxWireBC4; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_BC4v1_node.push_back( new TNode( Form( "BC4v1_Node_%d", wire ),
					 Form( "BC4v1_Node_%d", wire ),
					 "BC4V1Tube",
					 localPos+BeamAxis,
					 wireGlobalPos.y(),
					 wireGlobalPos.z()+offsetZ,
					 "rotBC4V1", "void" ) );
    }
  }

  // BC4 V2
  {
    const int lid = gGeom.GetDetectorId("BC4-V2");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotBC4V2", "rotBC4V2", Matrix );
    new TTUBE( "BC4V2Tube", "BC4V2Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<=MaxWireBC4; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_BC4v2_node.push_back( new TNode( Form( "BC4v2_Node_%d", wire ),
					 Form( "BC4v2_Node_%d", wire ),
					 "BC4V2Tube",
					 localPos+BeamAxis,
					 wireGlobalPos.y(),
					 wireGlobalPos.z()+offsetZ,
					 "rotBC4V2", "void" ) );
    }
  }

  // BC4 X1
  {
    const int lid = gGeom.GetDetectorId("BC4-X1");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotBC4X1", "rotBC4X1", Matrix );
    new TTUBE( "BC4X1Tube", "BC4X1Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<=MaxWireBC4; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_BC4x1_node.push_back( new TNode( Form( "BC4x1_Node_%d", wire ),
					 Form( "BC4x1_Node_%d", wire ),
					 "BC4X1Tube",
					 localPos+BeamAxis,
					 wireGlobalPos.y(),
					 wireGlobalPos.z()+offsetZ,
					 "rotBC4X1", "void" ) );
    }
  }

  // BC4 X2
  {
    const int lid = gGeom.GetDetectorId("BC4-X2");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotBC4X2", "rotBC4X2", Matrix );
    new TTUBE( "BC4X2Tube", "BC4X2Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<=MaxWireBC4; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_BC4x2_node.push_back( new TNode( Form( "BC4x2_Node_%d", wire ),
					 Form( "BC4x2_Node_%d", wire ),
					 "BC4X2Tube",
					 localPos+BeamAxis,
					 wireGlobalPos.y(),
					 wireGlobalPos.z()+offsetZ,
					 "rotBC4X2", "void" ) );
    }
  }

  ConstructionDone(__func__);
  return true;
}

//______________________________________________________________________________
bool
EventDisplay::ConstructSdcIn( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  const double wireL = 200.0;

  // SDC1 V1
  {
    const int lid = gGeom.GetDetectorId("SDC1-V1");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotSDC1V1", "rotSDC1V1", Matrix );
    new TTUBE( "SDC1V1Tube", "SDC1V1Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<=MaxWireSDC1; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_SDC1v1_node.push_back( new TNode( Form( "SDC1v1_Node_%d", wire ),
					  Form( "SDC1v1_Node_%d", wire ),
					  "SDC1V1Tube",
					  wireGlobalPos.x()+localPos,
					  wireGlobalPos.y(),
					  wireGlobalPos.z(),
					  "rotSDC1V1", "void" ) );
    }
  }

  // SDC1 V2
  {
    const int lid = gGeom.GetDetectorId("SDC1-V2");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double Z    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotV2", "rotV2", Matrix );
    new TTUBE( "SDC1V2Tube", "SDC1V2Tube", "void", Rmin, Rmax, Z );
    for( int wire=1; wire<=MaxWireSDC1; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_SDC1v2_node.push_back( new TNode( Form( "SDC1v2_Node_%d", wire ),
					  Form( "SDC1v2_Node_%d", wire ),
					  "SDC1V2Tube",
					  wireGlobalPos.x()+localPos,
					  wireGlobalPos.y(),
					  wireGlobalPos.z(),
					  "rotV2", "void" ) );
    }
  }

  // SDC1 X1
  {
    const int lid = gGeom.GetDetectorId("SDC1-X1");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double Z    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotX1", "rotX1", Matrix );
    new TTUBE( "SDC1X1Tube", "SDC1X1Tube", "void", Rmin, Rmax, Z );
    for( int wire=1; wire<=MaxWireSDC1; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_SDC1x1_node.push_back( new TNode( Form( "SDC1x1_Node_%d", wire ),
					  Form( "SDC1x1_Node_%d", wire ),
					  "SDC1X1Tube",
					  wireGlobalPos.x()+localPos,
					  wireGlobalPos.y(),
					  wireGlobalPos.z(),
					  "rotX1", "void" ) );
    }
  }

  // SDC1 X2
  {
    const int lid = gGeom.GetDetectorId("SDC1-X2");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double Z    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotX2", "rotX2", Matrix );
    new TTUBE( "SDC1X2Tube", "SDC1X2Tube", "void", Rmin, Rmax, Z );
    for( int wire=1; wire<=MaxWireSDC1; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_SDC1x2_node.push_back( new TNode( Form( "SDC1x2_Node_%d", wire ),
					  Form( "SDC1x2_Node_%d", wire ),
					  "SDC1X2Tube",
					  wireGlobalPos.x()+localPos,
					  wireGlobalPos.y(),
					  wireGlobalPos.z(),
					  "rotX2", "void" ) );
    }
  }

  // SDC1 U1
  {
    const int lid = gGeom.GetDetectorId("SDC1-U1");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double Z    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotU1", "rotU1", Matrix );
    new TTUBE( "SDC1U1Tube", "SDC1U1Tube", "void", Rmin, Rmax, Z );
    for( int wire=1; wire<=MaxWireSDC1; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_SDC1u1_node.push_back( new TNode( Form( "SDC1u1_Node_%d", wire ),
					  Form( "SDC1u1_Node_%d", wire ),
					  "SDC1U1Tube",
					  wireGlobalPos.x()+localPos,
					  wireGlobalPos.y(),
					  wireGlobalPos.z(),
					  "rotU1", "void" ) );
    }
  }

  // SDC1 U2
  {
    const int lid = gGeom.GetDetectorId("SDC1-U2");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double Z    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotU2", "rotU2", Matrix );
    new TTUBE( "SDC1U2Tube", "SDC1U2Tube", "void", Rmin, Rmax, Z );
    for( int wire=1; wire<=MaxWireSDC1; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_SDC1u2_node.push_back( new TNode( Form( "SDC1u2_Node_%d", wire ),
					  Form( "SDC1u2_Node_%d", wire ),
					  "SDC1U2Tube",
					  wireGlobalPos.x()+localPos,
					  wireGlobalPos.y(),
					  wireGlobalPos.z(),
					  "rotU2", "void" ) );
    }
      
  }


  // SFT-U
  {
    const int lid = gGeom.GetDetectorId("SFT-U");      
    double Rmin = 0.0;
    double Rmax = 0.375;
    double Z    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotSFTU", "rotSFTU", Matrix );
    new TTUBE( "SFTUTube", "SFTUTube", "void", Rmin, Rmax, Z );

    for( int wire=1; wire<=NumOfSegSFT_UV; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_SFTu_node.push_back( new TNode( Form( "SFTu_Node_%d", wire ),
					Form( "SFTu_Node_%d", wire ),
					"SFTUTube",
					wireGlobalPos.x()+localPos,
					wireGlobalPos.y(),
					wireGlobalPos.z(),
					"rotSFTU", "void" ) );
    }

  }

  // SFT-V
  {
    const int lid = gGeom.GetDetectorId("SFT-V");      
    double Rmin = 0.0;
    double Rmax = 0.375;
    double Z    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotSFTV", "rotSFTV", Matrix );
    new TTUBE( "SFTVTube", "SFTVTube", "void", Rmin, Rmax, Z );

    for( int wire=1; wire<=NumOfSegSFT_UV; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_SFTv_node.push_back( new TNode( Form( "SFTv_Node_%d", wire ),
					Form( "SFTv_Node_%d", wire ),
					"SFTVTube",
					wireGlobalPos.x()+localPos,
					wireGlobalPos.y(),
					wireGlobalPos.z(),
					"rotSFTV", "void" ) );
    }

  }

  // SFT-X
  {
    const int lid = gGeom.GetDetectorId("SFT-X");      
    double Rmin = 0.0;
    double Rmax = 0.5;
    double Z    = wireL/cos(gGeom.GetTiltAngle(lid)*math::Deg2Rad())/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotSFTX", "rotSFTX", Matrix );
    new TTUBE( "SFTXTube", "SFTXTube", "void", Rmin, Rmax, Z );

    for( int wire=1; wire<=NumOfSegSFT_X; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition( lid );
      m_SFTx_node.push_back( new TNode( Form( "SFTx_Node_%d", wire ),
					Form( "SFTx_Node_%d", wire ),
					"SFTXTube",
					wireGlobalPos.x()+localPos,
					wireGlobalPos.y(),
					wireGlobalPos.z(),
					"rotSFTX", "void" ) );
    }

  }

  ConstructionDone(__func__);
  return true;
}

//______________________________________________________________________________
bool
EventDisplay::ConstructSdcOut( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  const double wireLSDC2  = 1152.0;
  const double wireLSDC3Y = 1920.0;
  const double wireLSDC3X = 1280.0;

  // SDC2 X1
  {
    const int lid = gGeom.GetDetectorId("SDC2-X1");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireLSDC2/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotX1", "rotX1", Matrix );
    new TTUBE( "SDC2X1Tube", "SDC2X1Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<= MaxWireSDC2; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireLocalPos( localPos, 0., 0. );
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos( lid, wireLocalPos );
      m_SDC2x1_node.push_back( new TNode( Form( "SDC2x1_Node_%d", wire ),
					  Form( "SDC2x1_Node_%d", wire ),
					  "SDC2X1Tube",
					  wireGlobalPos.x(),
					  wireGlobalPos.y(),
					  wireGlobalPos.z(),
					  "rotX1", "void" ) );
    }
  }

  // SDC2 X2
  {
    const int lid = gGeom.GetDetectorId("SDC2-X2");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireLSDC2/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotX2", "rotX2", Matrix );
    new TTUBE( "SDC2X2Tube", "SDC2X2Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<= MaxWireSDC2; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireLocalPos( localPos, 0., 0. );
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos( lid, wireLocalPos );
      m_SDC2x2_node.push_back( new TNode( Form( "SDC2x2_Node_%d", wire ),
					  Form( "SDC2x2_Node_%d", wire ),
					  "SDC2X2Tube",
					  wireGlobalPos.x(),
					  wireGlobalPos.y(),
					  wireGlobalPos.z(),
					  "rotX2", "void" ) );
    }
  }

  // SDC2 Y1
  {
    const int lid = gGeom.GetDetectorId("SDC2-Y1");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireLSDC2/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotY1", "rotY1", Matrix );
    new TTUBE( "SDC2Y1Tube", "SDC2Y1Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<= MaxWireSDC2; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireLocalPos( localPos, 0., 0. );
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos( lid, wireLocalPos );
      m_SDC2y1_node.push_back( new TNode( Form( "SDC2y1_Node_%d", wire ),
					  Form( "SDC2y1_Node_%d", wire ),
					  "SDC2Y1Tube",
					  wireGlobalPos.x(),
					  wireGlobalPos.y(),
					  wireGlobalPos.z(),
					  "rotY1", "void" ) );
    }
  }

  // SDC2 Y2
  {
    const int lid = gGeom.GetDetectorId("SDC2-Y2");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireLSDC2/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotY2", "rotY2", Matrix );
    new TTUBE( "SDC2Y2Tube", "SDC2Y2Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<= MaxWireSDC2; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireLocalPos( localPos, 0., 0. );
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos( lid, wireLocalPos );
      m_SDC2y2_node.push_back( new TNode( Form( "SDC2y2_Node_%d", wire ),
					  Form( "SDC2y2_Node_%d", wire ),
					  "SDC2Y2Tube",
					  wireGlobalPos.x(),
					  wireGlobalPos.y(),
					  wireGlobalPos.z(),
					  "rotY2", "void" ) );
    }
  }

  // SDC3 Y1
  {
    const int lid = gGeom.GetDetectorId("SDC3-Y1");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireLSDC3Y/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ), Matrix );
    new TRotMatrix( "rotY1", "rotY1", Matrix );
    new TTUBE( "SDC3Y1Tube", "SDC3Y1Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<= MaxWireSDC3Y; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireLocalPos( localPos, 0., 0. );
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos( lid, wireLocalPos );
      m_SDC3y1_node.push_back( new TNode( Form( "SDC3y1_Node_%d", wire ),
					  Form( "SDC3y1_Node_%d", wire ),
					  "SDC3Y1Tube",
					  wireGlobalPos.x(),
					  wireGlobalPos.y(),
					  wireGlobalPos.z(),
					  "rotY1", "void" ) );
    }
  }

  // SDC3 Y2
  {
    const int lid = gGeom.GetDetectorId("SDC3-Y2");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireLSDC3Y/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotY2", "rotY2", Matrix );
    new TTUBE( "SDC3Y2Tube", "SDC3Y2Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<= MaxWireSDC3Y; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireLocalPos( localPos, 0., 0. );
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos( lid, wireLocalPos );
      m_SDC3y2_node.push_back( new TNode( Form( "SDC3y2_Node_%d", wire ),
					  Form( "SDC3y2_Node_%d", wire ),
					  "SDC3Y2Tube",
					  wireGlobalPos.x(),
					  wireGlobalPos.y(),
					  wireGlobalPos.z(),
					  "rotY2", "void" ) );
    }
  }

  // SDC3 X1
  {
    const int lid = gGeom.GetDetectorId("SDC3-X1");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireLSDC3X/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotX1", "rotX1", Matrix );
    new TTUBE( "SDC3X1Tube", "SDC3X1Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<= MaxWireSDC3X; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireLocalPos( localPos, 0., 0. );
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos( lid, wireLocalPos );
      m_SDC3x1_node.push_back( new TNode( Form( "SDC3x1_Node_%d", wire ),
					  Form( "SDC3x1_Node_%d", wire ),
					  "SDC3X1Tube",
					  wireGlobalPos.x(),
					  wireGlobalPos.y(),
					  wireGlobalPos.z(),
					  "rotX1", "void" ) );
    }
  }

  // SDC3 X2
  {
    const int lid = gGeom.GetDetectorId("SDC3-X2");
    double Rmin = 0.0;
    double Rmax = 0.01;
    double L    = wireLSDC3X/2.;
    double Matrix[9] = {};
    CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		   gGeom.GetRotAngle1( lid ),
		   gGeom.GetRotAngle2( lid ),
		   Matrix );
    new TRotMatrix( "rotX2", "rotX2", Matrix );
    new TTUBE( "SDC3X2Tube", "SDC3X2Tube", "void", Rmin, Rmax, L );
    for( int wire=1; wire<= MaxWireSDC3X; ++wire ){
      double localPos = gGeom.CalcWirePosition( lid, wire );
      ThreeVector wireLocalPos( localPos, 0., 0. );
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos( lid, wireLocalPos );
      m_SDC3x2_node.push_back( new TNode( Form( "SDC3x2_Node_%d", wire ),
					  Form( "SDC3x2_Node_%d", wire ),
					  "SDC3X2Tube",
					  wireGlobalPos.x(),
					  wireGlobalPos.y(),
					  wireGlobalPos.z(),
					  "rotX2", "void" ) );
    }
  }

  ConstructionDone(__func__);
  return true;
}

//______________________________________________________________________________
bool
EventDisplay::ConstructTarget( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  double TargetX = 50.0/2.0;
  double TargetY = 30.0/2.0; // Z
  double TargetZ = 30.0/2.0; // Y
  double rotMatTarget[9];

  new TBRIK( "target_brik", "target_brik", "void",
	     TargetX, TargetY, TargetZ );

  CalcRotMatrix( 0., 0., 0., rotMatTarget );
  new TRotMatrix( "rotTarget", "rotTarget", rotMatTarget );

  const int lid = IdTarget;
  ThreeVector GlobalPos = gGeom.GetGlobalPosition( lid );
  m_target_node = new TNode( "target_node", "target_node", "target_brik",
			     GlobalPos.x(), GlobalPos.y(), GlobalPos.z(),
			     "rotTarget", "void");

#if Vertex
  m_TargetXZ_box = new TBox( -TargetY, -TargetX, TargetY, TargetX );
  m_TargetXZ_box->SetFillColor(kWhite);
  m_TargetXZ_box->SetLineColor(kBlack);
  m_TargetYZ_box = new TBox( -TargetY, -TargetZ, TargetY, TargetZ );
  m_TargetYZ_box->SetFillColor(kWhite);
  m_TargetYZ_box->SetLineColor(kBlack);
#endif

  ConstructionDone(__func__);
  return true;
}

//______________________________________________________________________________
bool
EventDisplay::ConstructSCH( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  const int lid = gGeom.GetDetectorId("SCH");

  double rotMatSCH[9] = {};
  double SCHwallX =   11.5/2.0*NumOfSegSCH; // X
  double SCHwallY =    2.0/2.0*2.0; // Z
  double SCHwallZ =  450.0/2.0; // Y

  double SCHSegX =   11.5/2.0; // X
  double SCHSegY =    2.0/2.0; // Z
  double SCHSegZ =  450.0/2.0; // Y

  double overlap = 1.0;

  CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		 gGeom.GetRotAngle1( lid ),
		 gGeom.GetRotAngle2( lid ),
		 rotMatSCH );

  new TRotMatrix( "rotSCH", "rotSCH", rotMatSCH );
  ThreeVector  SCHwallPos = gGeom.GetGlobalPosition( lid );
  double offset =  gGeom.CalcWirePosition(lid, (double)NumOfSegSCH/2.-0.5 );
  new TBRIK( "SCHwall_brik", "SCHwall_brik", "void",
	     SCHwallX, SCHwallY, SCHwallZ);
  m_SCHwall_node = new TNode( "SCHwall_node", "SCHwall_node", "SCHwall_brik",
			      SCHwallPos.x() + offset,
			      SCHwallPos.y(),
			      SCHwallPos.z(),
			      "rotSCH", "void");

  m_SCHwall_node->SetVisibility(0);
  m_SCHwall_node->cd();


  new TBRIK( "SCHseg_brik", "SCHseg_brik", "void",
	     SCHSegX, SCHSegY, SCHSegZ );
  for( int i=0; i<NumOfSegSCH; i++ ){
    ThreeVector schSegLocalPos( (-NumOfSegSCH/2.+i)*(SCHSegX*2.-overlap)+10.5/2.,
				(-(i%2)*2+1)*SCHSegY-(i%2)*2+1,
				0. );
    m_SCHseg_node.push_back( new TNode( Form( "SCHseg_node_%d", i ),
					Form( "SCHseg_node_%d", i ),
					"SCHseg_brik",
					schSegLocalPos.x(),
					schSegLocalPos.y(),
					schSegLocalPos.z() ) );
  }

  m_node->cd();
  ConstructionDone(__func__);
  return true;
}

//______________________________________________________________________________
bool
EventDisplay::ConstructFBT( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  double FBTSegX =   6.0/2.0; // X
  double FBTSegY =    2.0/2.0; // Z
  double FBTSegZ =  550.0/2.0; // Y

  new TBRIK( "FBTseg_brik", "FBTseg_brik", "void",
	     FBTSegX, FBTSegY, FBTSegZ );

  double overlap = 2.0;

  double FBT1wallX =   4.0/2.0*MaxSegFBT1 + 2.0/2.0; // X
  double FBT1wallY =    2.0/2.0*2.0; // Z
  double FBT1wallZ =  450.0/2.0; // Y

  const int lid1D1 = gGeom.GetDetectorId("FBT1-DX1");
  double rotMatFBT1D1[9] = {};
  CalcRotMatrix( gGeom.GetTiltAngle( lid1D1 ),
		 gGeom.GetRotAngle1( lid1D1 ),
		 gGeom.GetRotAngle2( lid1D1 ),
		 rotMatFBT1D1 );

  new TRotMatrix( "rotFBT1D1", "rotFBT1D1", rotMatFBT1D1 );
  ThreeVector  FBT1D1wallPos = gGeom.GetGlobalPosition( lid1D1 );
  double offset1d1 =  gGeom.CalcWirePosition(lid1D1, (double)MaxSegFBT1/2.-0.5 );
  std::cout << "ofset1d1 = " << offset1d1 << std::endl;
  new TBRIK( "FBT1D1wall_brik", "FBT1D1wall_brik", "void",
	     FBT1wallX, FBT1wallY, FBT1wallZ);
  m_FBT1d1wall_node = new TNode( "FBT1D1wall_node", "FBT1D1wall_node", 
				 "FBT1D1wall_brik",
				 FBT1D1wallPos.x()+offset1d1,
				 FBT1D1wallPos.y(),
				 FBT1D1wallPos.z(),
				 "rotFBT1D1", "void");

  m_FBT1d1wall_node->SetVisibility(0);
  m_FBT1d1wall_node->cd();


  for( int i=0; i<MaxSegFBT1; i++ ){
    ThreeVector FBT1D1SegLocalPos( (-MaxSegFBT1/2.+i)*((FBTSegX-overlap/2.)*2)+2.0/2.,
				    (-(i%2)*2+1)*FBTSegY-(i%2)*2+1,
				    0. );
    m_FBT1d1seg_node.push_back( new TNode( Form( "FBT1D1seg_node_%d", i ),
					   Form( "FBT1D1seg_node_%d", i ),
					   "FBTseg_brik",
					   FBT1D1SegLocalPos.x(),
					   FBT1D1SegLocalPos.y(),
					   FBT1D1SegLocalPos.z() ) );
  }

  m_node->cd();

  const int lid1U1 = gGeom.GetDetectorId("FBT1-UX1");
  double rotMatFBT1U1[9] = {};
  CalcRotMatrix( gGeom.GetTiltAngle( lid1U1 ),
		 gGeom.GetRotAngle1( lid1U1 ),
		 gGeom.GetRotAngle2( lid1U1 ),
		 rotMatFBT1U1 );

  new TRotMatrix( "rotFBT1U1", "rotFBT1U1", rotMatFBT1U1 );
  ThreeVector  FBT1U1wallPos = gGeom.GetGlobalPosition( lid1U1 );
  new TBRIK( "FBT1U1wall_brik", "FBT1U1wall_brik", "void",
	     FBT1wallX, FBT1wallY, FBT1wallZ);
  double offset1u1 =  gGeom.CalcWirePosition(lid1U1, (double)MaxSegFBT1/2.-0.5 );
  m_FBT1u1wall_node = new TNode( "FBT1U1wall_node", "FBT1U1wall_node", 
				 "FBT1U1wall_brik",
				 FBT1U1wallPos.x() + offset1u1,
				 FBT1U1wallPos.y(),
				 FBT1U1wallPos.z(),
				 "rotFBT1U1", "void");

  m_FBT1u1wall_node->SetVisibility(0);
  m_FBT1u1wall_node->cd();


  for( int i=0; i<MaxSegFBT1; i++ ){
    ThreeVector FBT1U1SegLocalPos( (-MaxSegFBT1/2.+i)*(FBTSegX-overlap/2.)*2+2./2.,
				   (-(i%2)*2+1)*FBTSegY-(i%2)*2+1,
				   0. );
    m_FBT1u1seg_node.push_back( new TNode( Form( "FBT1U1seg_node_%d", i ),
					   Form( "FBT1U1seg_node_%d", i ),
					   "FBTseg_brik",
					   FBT1U1SegLocalPos.x(),
					   FBT1U1SegLocalPos.y(),
					   FBT1U1SegLocalPos.z() ) );
  }

  m_node->cd();

  const int lid1D2 = gGeom.GetDetectorId("FBT1-DX2");
  double rotMatFBT1D2[9] = {};
  CalcRotMatrix( gGeom.GetTiltAngle( lid1D2 ),
		 gGeom.GetRotAngle1( lid1D2 ),
		 gGeom.GetRotAngle2( lid1D2 ),
		 rotMatFBT1D2 );

  new TRotMatrix( "rotFBT1D2", "rotFBT1D2", rotMatFBT1D2 );
  ThreeVector  FBT1D2wallPos = gGeom.GetGlobalPosition( lid1D2 );
  new TBRIK( "FBT1D2wall_brik", "FBT1D2wall_brik", "void",
	     FBT1wallX, FBT1wallY, FBT1wallZ);
  double offset1d2 =  gGeom.CalcWirePosition(lid1D2, (double)MaxSegFBT1/2.-0.5 );
  m_FBT1d2wall_node = new TNode( "FBT1D2wall_node", "FBT1D2wall_node", 
				 "FBT1D2wall_brik",
				 FBT1D2wallPos.x() + offset1d2,
				 FBT1D2wallPos.y(),
				 FBT1D2wallPos.z(),
				 "rotFBT1D2", "void");

  m_FBT1d2wall_node->SetVisibility(0);
  m_FBT1d2wall_node->cd();


  for( int i=0; i<MaxSegFBT1; i++ ){
    ThreeVector FBT1D2SegLocalPos( (-MaxSegFBT1/2.+i)*(FBTSegX-overlap/2.)*2+2.0/2.,
				    (-(i%2)*2+1)*FBTSegY-(i%2)*2+1,
				    0. );
    m_FBT1d2seg_node.push_back( new TNode( Form( "FBT1D2seg_node_%d", i ),
					   Form( "FBT1D2seg_node_%d", i ),
					   "FBTseg_brik",
					   FBT1D2SegLocalPos.x(),
					   FBT1D2SegLocalPos.y(),
					   FBT1D2SegLocalPos.z() ) );
  }

  m_node->cd();

  const int lid1U2 = gGeom.GetDetectorId("FBT1-UX2");
  double rotMatFBT1U2[9] = {};
  CalcRotMatrix( gGeom.GetTiltAngle( lid1U2 ),
		 gGeom.GetRotAngle1( lid1U2 ),
		 gGeom.GetRotAngle2( lid1U2 ),
		 rotMatFBT1U2 );

  new TRotMatrix( "rotFBT1U2", "rotFBT1U2", rotMatFBT1U2 );
  ThreeVector  FBT1U2wallPos = gGeom.GetGlobalPosition( lid1U2 );
  new TBRIK( "FBT1U2wall_brik", "FBT1U2wall_brik", "void",
	     FBT1wallX, FBT1wallY, FBT1wallZ);
  double offset1u2 =  gGeom.CalcWirePosition(lid1U2, (double)MaxSegFBT1/2.-0.5 );

  m_FBT1u2wall_node = new TNode( "FBT1U2wall_node", "FBT1U2wall_node", 
				 "FBT1U2wall_brik",
				 FBT1U2wallPos.x() + offset1u2,
				 FBT1U2wallPos.y(),
				 FBT1U2wallPos.z(),
				 "rotFBT1U2", "void");

  m_FBT1u2wall_node->SetVisibility(0);
  m_FBT1u2wall_node->cd();


  for( int i=0; i<MaxSegFBT1; i++ ){
    ThreeVector FBT1U2SegLocalPos( (-MaxSegFBT1/2.+i)*(FBTSegX-overlap/2)*2+2.0/2.,
				    (-(i%2)*2+1)*FBTSegY-(i%2)*2+1,
				    0. );
    m_FBT1u2seg_node.push_back( new TNode( Form( "FBT1U2seg_node_%d", i ),
					   Form( "FBT1U2seg_node_%d", i ),
					   "FBTseg_brik",
					   FBT1U2SegLocalPos.x(),
					   FBT1U2SegLocalPos.y(),
					   FBT1U2SegLocalPos.z() ) );
  }

  m_node->cd();

  double FBT2wallX =   11.5/2.0*MaxSegFBT2; // X
  double FBT2wallY =    2.0/2.0*2.0; // Z
  double FBT2wallZ =  450.0/2.0; // Y

  const int lid2D1 = gGeom.GetDetectorId("FBT2-DX1");
  double rotMatFBT2D1[9] = {};
  CalcRotMatrix( gGeom.GetTiltAngle( lid2D1 ),
		 gGeom.GetRotAngle1( lid2D1 ),
		 gGeom.GetRotAngle2( lid2D1 ),
		 rotMatFBT2D1 );

  new TRotMatrix( "rotFBT2D1", "rotFBT2D1", rotMatFBT2D1 );
  ThreeVector  FBT2D1wallPos = gGeom.GetGlobalPosition( lid2D1 );
  new TBRIK( "FBT2D1wall_brik", "FBT2D1wall_brik", "void",
	     FBT2wallX, FBT2wallY, FBT2wallZ);
  double offset2d1 =  gGeom.CalcWirePosition(lid2D1, (double)MaxSegFBT2/2.-0.5 );
  m_FBT2d1wall_node = new TNode( "FBT2D1wall_node", "FBT2D1wall_node", 
				 "FBT2D1wall_brik",
				 FBT2D1wallPos.x() + offset2d1,
				 FBT2D1wallPos.y(),
				 FBT2D1wallPos.z(),
				 "rotFBT2D1", "void");

  m_FBT2d1wall_node->SetVisibility(0);
  m_FBT2d1wall_node->cd();


  for( int i=0; i<MaxSegFBT2; i++ ){
    ThreeVector FBT2D1SegLocalPos( (-MaxSegFBT2/2.+i)*(FBTSegX-overlap/2)*2+2.0/2.,
				    (-(i%2)*2+1)*FBTSegY-(i%2)*2+1,
				    0. );
    m_FBT2d1seg_node.push_back( new TNode( Form( "FBT2D1seg_node_%d", i ),
					   Form( "FBT2D1seg_node_%d", i ),
					   "FBTseg_brik",
					   FBT2D1SegLocalPos.x(),
					   FBT2D1SegLocalPos.y(),
					   FBT2D1SegLocalPos.z() ) );
  }

  m_node->cd();

  const int lid2U1 = gGeom.GetDetectorId("FBT2-UX1");
  double rotMatFBT2U1[9] = {};
  CalcRotMatrix( gGeom.GetTiltAngle( lid2U1 ),
		 gGeom.GetRotAngle1( lid2U1 ),
		 gGeom.GetRotAngle2( lid2U1 ),
		 rotMatFBT2U1 );

  new TRotMatrix( "rotFBT2U1", "rotFBT2U1", rotMatFBT2U1 );
  ThreeVector  FBT2U1wallPos = gGeom.GetGlobalPosition( lid2U1 );
  new TBRIK( "FBT2U1wall_brik", "FBT2U1wall_brik", "void",
	     FBT2wallX, FBT2wallY, FBT2wallZ);
  double offset2u1 =  gGeom.CalcWirePosition(lid2U1, (double)MaxSegFBT2/2.-0.5 );
  m_FBT2u1wall_node = new TNode( "FBT2U1wall_node", "FBT2U1wall_node", 
				 "FBT2U1wall_brik",
				 FBT2U1wallPos.x() + offset2u1,
				 FBT2U1wallPos.y(),
				 FBT2U1wallPos.z(),
				 "rotFBT2U1", "void");

  m_FBT2u1wall_node->SetVisibility(0);
  m_FBT2u1wall_node->cd();


  for( int i=0; i<MaxSegFBT2; i++ ){
    ThreeVector FBT2U1SegLocalPos( (-MaxSegFBT2/2.+i)*(FBTSegX-overlap/2.)*2+2.0/2.,
				    (-(i%2)*2+1)*FBTSegY-(i%2)*2+1,
				    0. );
    m_FBT2u1seg_node.push_back( new TNode( Form( "FBT2U1seg_node_%d", i ),
					   Form( "FBT2U1seg_node_%d", i ),
					   "FBTseg_brik",
					   FBT2U1SegLocalPos.x(),
					   FBT2U1SegLocalPos.y(),
					   FBT2U1SegLocalPos.z() ) );
  }

  m_node->cd();

  const int lid2D2 = gGeom.GetDetectorId("FBT2-DX2");
  double rotMatFBT2D2[9] = {};
  CalcRotMatrix( gGeom.GetTiltAngle( lid2D2 ),
		 gGeom.GetRotAngle1( lid2D2 ),
		 gGeom.GetRotAngle2( lid2D2 ),
		 rotMatFBT2D2 );

  new TRotMatrix( "rotFBT2D2", "rotFBT2D2", rotMatFBT2D2 );
  ThreeVector  FBT2D2wallPos = gGeom.GetGlobalPosition( lid2D2 );
  new TBRIK( "FBT2D2wall_brik", "FBT2D2wall_brik", "void",
	     FBT2wallX, FBT2wallY, FBT2wallZ);
  double offset2d2 =  gGeom.CalcWirePosition(lid2D2, (double)MaxSegFBT2/2.-0.5 );
  m_FBT2d2wall_node = new TNode( "FBT2D2wall_node", "FBT2D2wall_node", 
				 "FBT2D2wall_brik",
				 FBT2D2wallPos.x() + offset2d2,
				 FBT2D2wallPos.y(),
				 FBT2D2wallPos.z(),
				 "rotFBT2D2", "void");

  m_FBT2d2wall_node->SetVisibility(0);
  m_FBT2d2wall_node->cd();


  for( int i=0; i<MaxSegFBT2; i++ ){
    ThreeVector FBT2D2SegLocalPos( (-MaxSegFBT2/2.+i)*(FBTSegX-overlap/2.)*2+2.0/2.,
				    (-(i%2)*2+1)*FBTSegY-(i%2)*2+1,
				    0. );
    m_FBT2d2seg_node.push_back( new TNode( Form( "FBT2D2seg_node_%d", i ),
					   Form( "FBT2D2seg_node_%d", i ),
					   "FBTseg_brik",
					   FBT2D2SegLocalPos.x(),
					   FBT2D2SegLocalPos.y(),
					   FBT2D2SegLocalPos.z() ) );
  }

  m_node->cd();

  const int lid2U2 = gGeom.GetDetectorId("FBT2-UX2");
  double rotMatFBT2U2[9] = {};
  CalcRotMatrix( gGeom.GetTiltAngle( lid2U2 ),
		 gGeom.GetRotAngle1( lid2U2 ),
		 gGeom.GetRotAngle2( lid2U2 ),
		 rotMatFBT2U2 );

  new TRotMatrix( "rotFBT2U2", "rotFBT2U2", rotMatFBT2U2 );
  ThreeVector  FBT2U2wallPos = gGeom.GetGlobalPosition( lid2U2 );
  new TBRIK( "FBT2U2wall_brik", "FBT2U2wall_brik", "void",
	     FBT2wallX, FBT2wallY, FBT2wallZ);
  double offset2u2 =  gGeom.CalcWirePosition(lid2U2, (double)MaxSegFBT2/2.-0.5 );
  m_FBT2u2wall_node = new TNode( "FBT2U2wall_node", "FBT2U2wall_node", 
				 "FBT2U2wall_brik",
				 FBT2U2wallPos.x() + offset2u2,
				 FBT2U2wallPos.y(),
				 FBT2U2wallPos.z(),
				 "rotFBT2U2", "void");

  m_FBT2u2wall_node->SetVisibility(0);
  m_FBT2u2wall_node->cd();


  for( int i=0; i<MaxSegFBT2; i++ ){
    ThreeVector FBT2U2SegLocalPos( (-MaxSegFBT2/2.+i)*(FBTSegX-overlap/2.)*2.+2.0/2.,
				    (-(i%2)*2+1)*FBTSegY-(i%2)*2+1,
				    0. );
    m_FBT2u2seg_node.push_back( new TNode( Form( "FBT2U2seg_node_%d", i ),
					   Form( "FBT2U2seg_node_%d", i ),
					   "FBTseg_brik",
					   FBT2U2SegLocalPos.x(),
					   FBT2U2SegLocalPos.y(),
					   FBT2U2SegLocalPos.z() ) );
  }

  m_node->cd();


  ConstructionDone(__func__);
  return true;
}

//______________________________________________________________________________
bool
EventDisplay::ConstructTOF( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  const int lid = gGeom.GetDetectorId("TOF");

  double rotMatTOF[9] = {};
  double TOFwallX =   80.0/2.0*NumOfSegTOF; // X
  double TOFwallY =   30.0/2.0*2.0; // Z
  double TOFwallZ = 1800.0/2.0; // Y

  double TOFSegX =   80.0/2.0; // X
  double TOFSegY =   30.0/2.0; // Z
  double TOFSegZ = 1800.0/2.0; // Y

  double overlap = 5.0;

  CalcRotMatrix( gGeom.GetTiltAngle( lid ),
		 gGeom.GetRotAngle1( lid ),
		 gGeom.GetRotAngle2( lid ),
		 rotMatTOF );

  new TRotMatrix( "rotTOF", "rotTOF", rotMatTOF );
  ThreeVector  TOFwallPos = gGeom.GetGlobalPosition( lid );
  double offset =  gGeom.CalcWirePosition(lid, (double)NumOfSegTOF/2.-0.5 );
  new TBRIK( "TOFwall_brik", "TOFwall_brik", "void",
	     TOFwallX, TOFwallY, TOFwallZ);
  m_TOFwall_node = new TNode( "TOFwall_node", "TOFwall_node", "TOFwall_brik",
			      TOFwallPos.x() + offset,
			      TOFwallPos.y(),
			      TOFwallPos.z(),
			      "rotTOF", "void");

  m_TOFwall_node->SetVisibility(0);
  m_TOFwall_node->cd();


  new TBRIK( "TOFseg_brik", "TOFseg_brik", "void",
	     TOFSegX, TOFSegY, TOFSegZ);
  for( int i=0; i<NumOfSegTOF; i++ ){
    ThreeVector tofSegLocalPos( (-NumOfSegTOF/2.+i)*(TOFSegX*2.-overlap)+75./2.,
				(-(i%2)*2+1)*TOFSegY-(i%2)*2+1,
				0. );
    m_TOFseg_node.push_back( new TNode( Form( "TOFseg_node_%d", i ),
					Form( "TOFseg_node_%d", i ),
					"TOFseg_brik",
					tofSegLocalPos.x(),
					tofSegLocalPos.y(),
					tofSegLocalPos.z() ) );
  }

  m_node->cd();
  ConstructionDone(__func__);
  return true;
}

//______________________________________________________________________________

bool EventDisplay::ConstructCATCH(void)
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  m_canvas_catch = new TCanvas( "canvas_catch", "CATCH Event Display",
				1600, 800 );
  m_canvas_catch->Divide(1,3);
  m_canvas_catch->cd(1)->SetPad( 0.001,0.001,0.499,0.999);
  m_canvas_catch->cd(2)->SetPad( 0.501,0.501,0.999,0.999);
  m_canvas_catch->cd(3)->SetPad( 0.501,0.001,0.999,0.499);

  m_canvas_catch->cd(1)->SetGrid(); 

  m_hbase_catch = new TH2F("hbase_catch","Event Display XY plane", 180, -180, 180, 180, -180, 180);
  m_hbase_catch->SetMaximum(200);
  m_hbase_catch->SetMinimum(-1);
  m_hbase_catch->Draw();

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  m_Tgt_Arc = new TArc(0, 0, 40./2);
  m_Tgt_Arc->SetLineColor(kCyan);  
  m_Tgt_Arc->SetFillStyle(0);
  m_Tgt_Arc->Draw("same");

  m_CFRP_Arc = new TArc(0, 0, 80./2);
  m_CFRP_Arc->SetLineColor(kBlack);  
  m_CFRP_Arc->SetFillStyle(0);
  m_CFRP_Arc->SetLineWidth(2);
  m_CFRP_Arc->Draw("same");

  ConstructCFT();
  ConstructBGO();
  ConstructPiID();

  m_canvas_catch->cd(2)->SetGrid(); 
  m_hbase_catch_zx = new TH2F("hbase_catch_zx","Event Display ZX plane", 600, -200, 400, 180, -180, 180);
  m_hbase_catch_zx->Draw();

  m_canvas_catch->cd(3)->SetGrid(); 
  m_hbase_catch_zy = new TH2F("hbase_catch_zy","Event Display ZY plane", 600, -200, 400, 180, -180, 180);
  m_hbase_catch_zy->Draw();


  return true;
}

//______________________________________________________________________________
bool EventDisplay::ConstructCFT(void)
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  m_canvas_catch->cd(1);

  for (int i=0; i<NumOfPlaneCFT; i++) {
    m_CFT_Arc_cont[i].reserve(NumOfSegCFT[i]);
    for (int seg=0; seg<NumOfSegCFT[i]; seg++) {
      double x, y;
      FiberPosPhi(i, seg, &x, &y);
      TArc *arc = new TArc(x, y, 0.75/2);	
      arc->SetLineColor(kBlack);
      arc->SetFillStyle(0);
      arc->Draw("same");
      m_CFT_Arc_cont[i].push_back(arc);
    }
  }

  return true;
}

//______________________________________________________________________________
bool EventDisplay::ConstructBGO(void)
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  m_canvas_catch->cd(1);

  int unit=0;

  for (int i=0; i<NumOfBGOUnit; i++) {
    double theta = (double)i*45.;

    for (int j=0; j<NumOfBGOInOneUnit; j++) {
      double x0 = RadiusOfBGOSurface+BGO_Y/2;
      double y0 = (double)(j-0.5)*BGO_X;

      double x1 = x0+BGO_Y/2;
      double y1 = y0+BGO_X/2;

      double x2 = x0-BGO_Y/2;
      double y2 = y0+BGO_X/2;

      double x3 = x0-BGO_Y/2;
      double y3 = y0-BGO_X/2;

      double x4 = x0+BGO_Y/2;
      double y4 = y0-BGO_X/2;

      ThreeVector pos1((x1*cos(theta*math::Deg2Rad()) - y1*sin(theta*math::Deg2Rad())),
		       (x1*sin(theta*math::Deg2Rad()) + y1*cos(theta*math::Deg2Rad())),
		       0);
      ThreeVector pos2((x2*cos(theta*math::Deg2Rad()) - y2*sin(theta*math::Deg2Rad())),
		       (x2*sin(theta*math::Deg2Rad()) + y2*cos(theta*math::Deg2Rad())),
		       0);
      ThreeVector pos3((x3*cos(theta*math::Deg2Rad()) - y3*sin(theta*math::Deg2Rad())),
		       (x3*sin(theta*math::Deg2Rad()) + y3*cos(theta*math::Deg2Rad())),
		       0);
      ThreeVector pos4((x4*cos(theta*math::Deg2Rad()) - y4*sin(theta*math::Deg2Rad())),
		       (x4*sin(theta*math::Deg2Rad()) + y4*cos(theta*math::Deg2Rad())),
		       0);

      TLine *l1 = new TLine(pos1.x(), pos1.y(), pos2.x(), pos2.y());
      TLine *l2 = new TLine(pos2.x(), pos2.y(), pos3.x(), pos3.y());
      TLine *l3 = new TLine(pos3.x(), pos3.y(), pos4.x(), pos4.y());
      TLine *l4 = new TLine(pos4.x(), pos4.y(), pos1.x(), pos1.y());
      l1->Draw("same");
      l2->Draw("same");
      l3->Draw("same");
      l4->Draw("same");

      unit = j+3*i;
      m_BGO_Line_cont[unit].push_back(l1);
      m_BGO_Line_cont[unit].push_back(l2);
      m_BGO_Line_cont[unit].push_back(l3);
      m_BGO_Line_cont[unit].push_back(l4);

      //std::cout << unit << std::endl;
    }
  }

  for (int i=0; i<NumOfBGOUnit; i++) {
    double theta = 22.5 + (double)i*45.;

    for (int j=0; j<NumOfBGOInOneUnit2; j++) {
      double x0 = RadiusOfBGOSurface2+BGO_Y/2;
      double y0 = (double)(j)*BGO_X;

      double x1 = x0+BGO_Y/2;
      double y1 = y0+BGO_X/2;

      double x2 = x0-BGO_Y/2;
      double y2 = y0+BGO_X/2;

      double x3 = x0-BGO_Y/2;
      double y3 = y0-BGO_X/2;

      double x4 = x0+BGO_Y/2;
      double y4 = y0-BGO_X/2;

      ThreeVector pos1((x1*cos(theta*math::Deg2Rad()) - y1*sin(theta*math::Deg2Rad())),
		       (x1*sin(theta*math::Deg2Rad()) + y1*cos(theta*math::Deg2Rad())),
		       0);
      ThreeVector pos2((x2*cos(theta*math::Deg2Rad()) - y2*sin(theta*math::Deg2Rad())),
		       (x2*sin(theta*math::Deg2Rad()) + y2*cos(theta*math::Deg2Rad())),
		       0);
      ThreeVector pos3((x3*cos(theta*math::Deg2Rad()) - y3*sin(theta*math::Deg2Rad())),
		       (x3*sin(theta*math::Deg2Rad()) + y3*cos(theta*math::Deg2Rad())),
		       0);
      ThreeVector pos4((x4*cos(theta*math::Deg2Rad()) - y4*sin(theta*math::Deg2Rad())),
		       (x4*sin(theta*math::Deg2Rad()) + y4*cos(theta*math::Deg2Rad())),
		       0);

      TLine *l1 = new TLine(pos1.x(), pos1.y(), pos2.x(), pos2.y());
      TLine *l2 = new TLine(pos2.x(), pos2.y(), pos3.x(), pos3.y());
      TLine *l3 = new TLine(pos3.x(), pos3.y(), pos4.x(), pos4.y());
      TLine *l4 = new TLine(pos4.x(), pos4.y(), pos1.x(), pos1.y());
      l1->Draw("same");
      l2->Draw("same");
      l3->Draw("same");
      l4->Draw("same");

      unit = j+NumOfBGOInOneUnit+3*i;
      m_BGO_Line_cont[unit].push_back(l1);
      m_BGO_Line_cont[unit].push_back(l2);
      m_BGO_Line_cont[unit].push_back(l3);
      m_BGO_Line_cont[unit].push_back(l4);

      //std::cout << unit << std::endl;
    }
  }

  return true;
}

//______________________________________________________________________________
bool EventDisplay::ConstructPiID(void)
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  m_canvas_catch->cd(1);

  int unit=0;

  for (int i=0; i<NumOfPiIDUnit; i++) {
    double theta = (double)i*45.;

    for (int j=0; j<NumOfPiIDInOneUnit; j++) {
      double x0 = RadiusOfPiIDSurface+PiID_Y/2;
      double y0 = (double)(j-1)*PiID_X;

      double x1 = x0+PiID_Y/2;
      double y1 = y0+PiID_X/2;

      double x2 = x0-PiID_Y/2;
      double y2 = y0+PiID_X/2;

      double x3 = x0-PiID_Y/2;
      double y3 = y0-PiID_X/2;

      double x4 = x0+PiID_Y/2;
      double y4 = y0-PiID_X/2;

      ThreeVector pos1((x1*cos(theta*math::Deg2Rad()) - y1*sin(theta*math::Deg2Rad())),
		       (x1*sin(theta*math::Deg2Rad()) + y1*cos(theta*math::Deg2Rad())),
		       0);
      ThreeVector pos2((x2*cos(theta*math::Deg2Rad()) - y2*sin(theta*math::Deg2Rad())),
		       (x2*sin(theta*math::Deg2Rad()) + y2*cos(theta*math::Deg2Rad())),
		       0);
      ThreeVector pos3((x3*cos(theta*math::Deg2Rad()) - y3*sin(theta*math::Deg2Rad())),
		       (x3*sin(theta*math::Deg2Rad()) + y3*cos(theta*math::Deg2Rad())),
		       0);
      ThreeVector pos4((x4*cos(theta*math::Deg2Rad()) - y4*sin(theta*math::Deg2Rad())),
		       (x4*sin(theta*math::Deg2Rad()) + y4*cos(theta*math::Deg2Rad())),
		       0);

      TLine *l1 = new TLine(pos1.x(), pos1.y(), pos2.x(), pos2.y());
      TLine *l2 = new TLine(pos2.x(), pos2.y(), pos3.x(), pos3.y());
      TLine *l3 = new TLine(pos3.x(), pos3.y(), pos4.x(), pos4.y());
      TLine *l4 = new TLine(pos4.x(), pos4.y(), pos1.x(), pos1.y());
      l1->Draw("same");
      l2->Draw("same");
      l3->Draw("same");
      l4->Draw("same");

      unit = i*NumOfPiIDInOneUnit + j + i;
      m_PiID_Line_cont[unit].push_back(l1);
      m_PiID_Line_cont[unit].push_back(l2);
      m_PiID_Line_cont[unit].push_back(l3);
      m_PiID_Line_cont[unit].push_back(l4);

      //std::cout << unit << std::endl;
    }
  }

  for (int i=0; i<NumOfPiIDUnit; i++) {
    double theta = 22.5 + (double)i*45.;

    for (int j=0; j<NumOfPiIDInOneUnit2; j++) {
      double x0 = RadiusOfPiID2Surface+PiID2_Y/2;
      double y0 = (double)(j)*PiID2_X;

      double x1 = x0+PiID2_Y/2;
      double y1 = y0+PiID2_X/2;

      double x2 = x0-PiID2_Y/2;
      double y2 = y0+PiID2_X/2;

      double x3 = x0-PiID2_Y/2;
      double y3 = y0-PiID2_X/2;

      double x4 = x0+PiID2_Y/2;
      double y4 = y0-PiID2_X/2;

      ThreeVector pos1((x1*cos(theta*math::Deg2Rad()) - y1*sin(theta*math::Deg2Rad())),
		       (x1*sin(theta*math::Deg2Rad()) + y1*cos(theta*math::Deg2Rad())),
		       0);
      ThreeVector pos2((x2*cos(theta*math::Deg2Rad()) - y2*sin(theta*math::Deg2Rad())),
		       (x2*sin(theta*math::Deg2Rad()) + y2*cos(theta*math::Deg2Rad())),
		       0);
      ThreeVector pos3((x3*cos(theta*math::Deg2Rad()) - y3*sin(theta*math::Deg2Rad())),
		       (x3*sin(theta*math::Deg2Rad()) + y3*cos(theta*math::Deg2Rad())),
		       0);
      ThreeVector pos4((x4*cos(theta*math::Deg2Rad()) - y4*sin(theta*math::Deg2Rad())),
		       (x4*sin(theta*math::Deg2Rad()) + y4*cos(theta*math::Deg2Rad())),
		       0);

      TLine *l1 = new TLine(pos1.x(), pos1.y(), pos2.x(), pos2.y());
      TLine *l2 = new TLine(pos2.x(), pos2.y(), pos3.x(), pos3.y());
      TLine *l3 = new TLine(pos3.x(), pos3.y(), pos4.x(), pos4.y());
      TLine *l4 = new TLine(pos4.x(), pos4.y(), pos1.x(), pos1.y());
      l1->Draw("same");
      l2->Draw("same");
      l3->Draw("same");
      l4->Draw("same");

      unit = (i+1)*NumOfPiIDInOneUnit + i*NumOfPiIDInOneUnit2 + j;
      m_PiID_Line_cont[unit].push_back(l1);
      m_PiID_Line_cont[unit].push_back(l2);
      m_PiID_Line_cont[unit].push_back(l3);
      m_PiID_Line_cont[unit].push_back(l4);

      //std::cout << unit << std::endl;
    }
  }

  return true;
}

//______________________________________________________________________________

bool EventDisplay::ConstructCATCH3d(void)
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  m_canvas_catch3d = new TCanvas( "canvas_catch3d", "CATCH Event Display",
				800, 800 );


  m_geometry_catch = new TGeometry( "evdisp_catch","CATCH Event Display" );

  ThreeVector worldSize( 200., 200., 400. ); /*mm*/
  new TBRIK( "world_catch", "world_catch", "void",
	     worldSize.x(), worldSize.y(), worldSize.z() );

  m_node_catch = new TNode( "node_catch", "node_catch", "world_catch", 0., 0., 0. );
  m_geometry_catch->GetNode("node_catch")->SetVisibility(0);


  double Rmin = 0.0;
  double Rmax = 0.75/2;
  double L    = 400./2.;

  new TTUBE( "CFTFiberTube", "CFTFiberTube", "void", Rmin, Rmax, L );

  for (int i=0; i<NumOfPlaneCFT/2; i++) {
    int layer = 2*i+1;
    m_CFT_node_cont[i].reserve(NumOfSegCFT[layer]);
    for (int seg=0; seg<NumOfSegCFT[layer]; seg++) {
      double x, y;
      FiberPosPhi(layer, seg, &x, &y);

      m_CFT_node_cont[i].push_back( new TNode( Form( "CFT%d_Node_%d", i, seg ),
					       Form( "CFT%d_Node_%d", i, seg ),
					       "CFTFiberTube",
					       x, y, L));
    }
  }


  new TBRIK( "BGOBrik", "BGOBrik", "void", BGO_Y/2., BGO_Z/2., BGO_X/2. );
  for (int i=0; i<NumOfBGOUnit; i++) {
    double theta = (double)i*45.;

    double rotMat[9] = {};
    CalcRotMatrix( theta, 0., 0., rotMat );

    new TRotMatrix( Form( "rotBGO_%d", i ),
		    Form( "rotBGO_%d", i ),
		    rotMat);

    for (int j=0; j<NumOfBGOInOneUnit; j++) {
      double x0 = RadiusOfBGOSurface+BGO_Y/2;
      double y0 = (double)(j-0.5)*BGO_X;

      ThreeVector pos0((x0*cos(theta*math::Deg2Rad()) - y0*sin(theta*math::Deg2Rad())),
		       (x0*sin(theta*math::Deg2Rad()) + y0*cos(theta*math::Deg2Rad())),
		       L+offsetBGO);

      m_BGOseg_node_cont.push_back( new TNode( Form( "BGOseg_node_%d", i*NumOfBGOInOneUnit + j + i ),
					       Form( "BGOseg_node_%d", i*NumOfBGOInOneUnit + j + i ),
					       "BGOBrik",
					       pos0.x(),
					       pos0.y(),
					       pos0.z(),
					       Form( "rotBGO_%d", i ),
					       "void") );



      //std::cout << unit << std::endl;
    }
  }

  for (int i=0; i<NumOfBGOUnit; i++) {
    double theta = 22.5 + (double)i*45.;

    double rotMat[9] = {};
    CalcRotMatrix( theta, 0., 0., rotMat );
    new TRotMatrix( Form( "rotBGO2_%d", i ),
		    Form( "rotBGO2_%d", i ),
		    rotMat);

    for (int j=0; j<NumOfBGOInOneUnit2; j++) {
      double x0 = RadiusOfBGOSurface2+BGO_Y/2;
      double y0 = (double)(j)*BGO_X;

      ThreeVector pos0((x0*cos(theta*math::Deg2Rad()) - y0*sin(theta*math::Deg2Rad())),
		       (x0*sin(theta*math::Deg2Rad()) + y0*cos(theta*math::Deg2Rad())),
		       L+offsetBGO);

      m_BGOseg_node_cont.push_back( new TNode( Form( "BGOseg_node_%d", (i+1)*NumOfBGOInOneUnit + i*NumOfBGOInOneUnit2 + j ),
					       Form( "BGOseg_node_%d", (i+1)*NumOfBGOInOneUnit + i*NumOfBGOInOneUnit2 + j ),
					       "BGOBrik",
					       pos0.x(),
					       pos0.y(),
					       pos0.z(),
					       Form( "rotBGO2_%d", i ),
					       "void") );
      //std::cout << unit << std::endl;
    }
  }

  std::string node_name;
  for (int seg=0; seg<2; seg++) {
    node_name = Form( "BGOseg_node_%d", seg );
    
    TNode *node = m_geometry_catch->GetNode( node_name.c_str() );
    if( !node ){
      hddaq::cout << "#E " << func_name << " "
		  << "no such node : " << node_name << std::endl;
      return false;
    }
    node->SetVisibility(0);    
  }


  new TBRIK( "PiIDBrik", "PiIDBrik", "void", PiID_Y/2., PiID_Z/2., PiID_X/2. );
  for (int i=0; i<NumOfPiIDUnit; i++) {
    double theta = (double)i*45.;

    double rotMat[9] = {};
    CalcRotMatrix( theta, 0., 0., rotMat );

    new TRotMatrix( Form( "rotPiID_%d", i ),
		    Form( "rotPiID_%d", i ),
		    rotMat);

    for (int j=0; j<NumOfPiIDInOneUnit; j++) {
      double x0 = RadiusOfPiIDSurface+PiID_Y/2;
      double y0 = (double)(j-1)*PiID_X;

      ThreeVector pos0((x0*cos(theta*math::Deg2Rad()) - y0*sin(theta*math::Deg2Rad())),
		       (x0*sin(theta*math::Deg2Rad()) + y0*cos(theta*math::Deg2Rad())),
		       L+offsetBGO);

      m_PiIDseg_node_cont.push_back( new TNode( Form( "PiIDseg_node_%d", i*NumOfPiIDInOneUnit + j + i ),
					       Form( "PiIDseg_node_%d", i*NumOfPiIDInOneUnit + j + i ),
					       "PiIDBrik",
					       pos0.x(),
					       pos0.y(),
					       pos0.z(),
					       Form( "rotPiID_%d", i ),
					       "void") );



      //std::cout << unit << std::endl;
    }
  }

  new TBRIK( "PiIDBrik2", "PiIDBrik2", "void", PiID2_Y/2., PiID2_Z/2., PiID2_X/2. );
  for (int i=0; i<NumOfPiIDUnit; i++) {
    double theta = 22.5 + (double)i*45.;

    double rotMat[9] = {};
    CalcRotMatrix( theta, 0., 0., rotMat );
    new TRotMatrix( Form( "rotPiID2_%d", i ),
		    Form( "rotPiID2_%d", i ),
		    rotMat);

    for (int j=0; j<NumOfPiIDInOneUnit2; j++) {
      double x0 = RadiusOfPiID2Surface+PiID_Y/2;
      double y0 = (double)(j)*PiID_X;

      ThreeVector pos0((x0*cos(theta*math::Deg2Rad()) - y0*sin(theta*math::Deg2Rad())),
		       (x0*sin(theta*math::Deg2Rad()) + y0*cos(theta*math::Deg2Rad())),
		       L+offsetBGO);

      m_PiIDseg_node_cont.push_back( new TNode( Form( "PiIDseg_node_%d", (i+1)*NumOfPiIDInOneUnit + i*NumOfPiIDInOneUnit2 + j ),
					       Form( "PiIDseg_node_%d", (i+1)*NumOfPiIDInOneUnit + i*NumOfPiIDInOneUnit2 + j ),
					       "PiIDBrik2",
					       pos0.x(),
					       pos0.y(),
					       pos0.z(),
					       Form( "rotPiID2_%d", i ),
					       "void") );
      //std::cout << unit << std::endl;
    }
  }

  for (int seg=0; seg<3; seg++) {
    node_name = Form( "PiIDseg_node_%d", seg );
    
    TNode *node = m_geometry_catch->GetNode( node_name.c_str() );
    if( !node ){
      hddaq::cout << "#E " << func_name << " "
		  << "no such node : " << node_name << std::endl;
      return false;
    }
    node->SetVisibility(0);    
  }

  m_geometry_catch->Draw();
  m_canvas_catch3d->Update();

  return true;
}

//______________________________________________________________________________
void
EventDisplay::DrawInitTrack( void )
{
  m_canvas->cd(2);
  if( m_init_step_mark ) m_init_step_mark->Draw();

  m_canvas->Update();

}

//______________________________________________________________________________
void
EventDisplay::DrawHitWire( int lid, int hit_wire, bool range_check, bool tdc_check )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( hit_wire<=0 ) return;

  std::string node_name;

  const std::string bcout_node_name[NumOfLayersBcOut]
    = { Form( "BC3x1_Node_%d", hit_wire ),
	Form( "BC3x2_Node_%d", hit_wire ),
	Form( "BC3u1_Node_%d", hit_wire ),
	Form( "BC3u2_Node_%d", hit_wire ),
	Form( "BC3v1_Node_%d", hit_wire ),
	Form( "BC3v2_Node_%d", hit_wire ),
	Form( "BC4x1_Node_%d", hit_wire ),
	Form( "BC4x2_Node_%d", hit_wire ),
	Form( "BC4u1_Node_%d", hit_wire ),
	Form( "BC4u2_Node_%d", hit_wire ),
	Form( "BC4v1_Node_%d", hit_wire ),
	Form( "BC4v2_Node_%d", hit_wire ) };

  const std::string sdcin_node_name[NumOfLayersSdcIn]
    = { Form( "SDC1v1_Node_%d", hit_wire ),
	Form( "SDC1v2_Node_%d", hit_wire ),
	Form( "SDC1x1_Node_%d", hit_wire ),
	Form( "SDC1x2_Node_%d", hit_wire ),
	Form( "SDC1u1_Node_%d", hit_wire ),
	Form( "SDC1u2_Node_%d", hit_wire ),
	Form( "SFTu_Node_%d", hit_wire ),
	Form( "SFTv_Node_%d", hit_wire ),
	Form( "SFTx_Node_%d", hit_wire )};

  const std::string sdcout_node_name[NumOfLayersSdcOut]
    = { Form( "SDC2x1_Node_%d", hit_wire ),
	Form( "SDC2x2_Node_%d", hit_wire ),
	Form( "SDC2y1_Node_%d", hit_wire ),
	Form( "SDC2y2_Node_%d", hit_wire ),
	Form( "SDC3y1_Node_%d", hit_wire ),
	Form( "SDC3y2_Node_%d", hit_wire ),
	Form( "SDC3x1_Node_%d", hit_wire ),
	Form( "SDC3x2_Node_%d", hit_wire ) };

  switch ( lid ) {

    // SDC1
  case 1: case 2: case 3: case 4: case 5: case 6:
    if( hit_wire>MaxWireSDC1 ) return;
    node_name = sdcin_node_name[lid-1];
    break;

    // SFT-UV
  case 7: case 8:
    if( hit_wire>NumOfSegSFT_UV ) return;
    node_name = sdcin_node_name[lid-1];
    break;

    // SFT-X
  case 9:
    if( hit_wire>NumOfSegSFT_X ) return;
    node_name = sdcin_node_name[lid-1];
    break;

    // SDC2
  case 31: case 32: case 33: case 34:
    if( hit_wire>MaxWireSDC2 ) return;
    node_name = sdcout_node_name[lid-31];
    break;

    // SDC3Y
  case 35: case 36:
    if( hit_wire>MaxWireSDC3Y ) return;
    node_name = sdcout_node_name[lid-31];
    break;

    // SDC3X
  case 37: case 38:
    if( hit_wire>MaxWireSDC3X ) return;
    node_name = sdcout_node_name[lid-31];
    break;

    // BC3
  case 113: case 114: case 115: case 116: case 117: case 118:
    if( hit_wire>MaxWireBC4 ) return;
    node_name = bcout_node_name[lid-113];
    break;

    // BC4
  case 119: case 120: case 121: case 122: case 123: case 124:
    if( hit_wire>MaxWireBC3 ) return;
    node_name = bcout_node_name[lid-113];
    break;

  default:
    hddaq::cout << "#E " << func_name << " "
		<< "no such plane : " << lid << std::endl;
    return;
  }

  TNode *node = m_geometry->GetNode( node_name.c_str() );
  if( !node ){
    hddaq::cout << "#E " << func_name << " "
		<< "no such node : " << node_name << std::endl;
    return;
  }

  node->SetVisibility(1);
  if( range_check && tdc_check )
    node->SetLineColor(kBlue);
  else if( range_check && !tdc_check )
    node->SetLineColor(28);
  else
    node->SetLineColor(kBlack);

  m_canvas->cd(2);
  m_geometry->Draw();
  m_canvas->Update();

}

//______________________________________________________________________________
void
EventDisplay::DrawHitHodoscope( int lid, int seg, int Tu, int Td )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  TString node_name;
  if( seg<0 ) return;

  if( lid == IdBH2 ){
    if( seg>=NumOfSegBH2 ) return;
    node_name = Form( "BH2seg_node_%d", seg );
  } else if( lid == IdSFT_U ){
    if( seg>=NumOfSegSFT_UV ) return;
    node_name = Form( "SFTu_Node_%d", seg );
  } else if( lid == IdSFT_V ){
    if( seg>=NumOfSegSFT_UV ) return;
    node_name = Form( "SFTv_Node_%d", seg );
  }  else if( lid == IdSFT_X ){
    if( seg>=NumOfSegSFT_X ) return;
    node_name = Form( "SFTx_Node_%d", seg );
  } else if( lid == IdSCH ){
    if( seg>=NumOfSegSCH ) return;
    node_name = Form( "SCHseg_node_%d", seg );
  } else if( lid == IdFBT1_D1 ){
    if( seg>=MaxSegFBT1 ) return;
    node_name = Form( "FBT1D1seg_node_%d", seg );
  } else if( lid == IdFBT1_U1 ){
    if( seg>=MaxSegFBT1 ) return;
    node_name = Form( "FBT1U1seg_node_%d", seg );
  } else if( lid == IdFBT1_D2 ){
    if( seg>=MaxSegFBT1 ) return;
    node_name = Form( "FBT1D2seg_node_%d", seg );
  } else if( lid == IdFBT1_U2 ){
    if( seg>=MaxSegFBT1 ) return;
    node_name = Form( "FBT1U2seg_node_%d", seg );
  } else if( lid == IdFBT2_D1 ){
    if( seg>=MaxSegFBT2 ) return;
    node_name = Form( "FBT2D1seg_node_%d", seg );
  } else if( lid == IdFBT2_U1 ){
    if( seg>=MaxSegFBT2 ) return;
    node_name = Form( "FBT2U1seg_node_%d", seg );
  } else if( lid == IdFBT2_D2 ){
    if( seg>=MaxSegFBT2 ) return;
    node_name = Form( "FBT2D2seg_node_%d", seg );
  } else if( lid == IdFBT2_U2 ){
    if( seg>=MaxSegFBT2 ) return;
    node_name = Form( "FBT2U2seg_node_%d", seg );
  } else if( lid == IdTOF ){
    if( seg>=NumOfSegTOF ) return;
    node_name = Form( "TOFseg_node_%d", seg );
  } else {
    hddaq::cout << "#E " << func_name << " "
		<< "no such plane : " << lid << std::endl;
    return;
  }

  TNode *node = m_geometry->GetNode( node_name );
  if( !node ){
    hddaq::cout << "#E " << func_name << " "
		<< "no such node : " << node_name << std::endl;
    return;
  }

  node->SetVisibility(1);

  if( Tu>0 && Td>0 ){
    node->SetLineColor( kBlue );
  }
  else {
    node->SetLineColor( kGreen );
  }

  m_canvas->cd(2);
  m_geometry->Draw();
  m_canvas->Update();
}

//______________________________________________________________________________
void
EventDisplay::DrawBcOutLocalTrack( DCLocalTrack *tp )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  const double offsetZ = -3108.4;

  double z0 = gGeom.GetLocalZ("BC3-X1") - 100.;
  double x0 = tp->GetX( z0 ) + BeamAxis;
  double y0 = tp->GetY( z0 );
  z0 += offsetZ;

  double z1 = zK18Target + 100.;
  double x1 = tp->GetX( z1 ) + BeamAxis;
  double y1 = tp->GetY( z1 );
  z1 += offsetZ;

  ThreeVector gPos0( x0, y0, z0 );
  ThreeVector gPos1( x1, y1, z1 );

  TPolyLine3D *p = new TPolyLine3D(2);
  p->SetLineColor(kRed);
  p->SetLineWidth(1);
  p->SetPoint( 0, gPos0.x(), gPos0.y(), gPos0.z() );
  p->SetPoint( 1, gPos1.x(), gPos1.y(), gPos1.z() );
  m_BcOutTrack.push_back(p);
  m_canvas->cd(2);
  p->Draw();

#if Vertex
  z0 = zK18Target + MinZ;
  z1 = zK18Target;
  x0 = tp->GetX( z0 ); y0 = tp->GetY( z0 );
  x1 = tp->GetX( z1 ); y1 = tp->GetY( z1 );
  z0 -= zK18Target;
  z1 -= zK18Target;
  {
    m_canvas_vertex->cd(1);
    TPolyLine *line = new TPolyLine(2);
    line->SetPoint( 0, z0, x0 );
    line->SetPoint( 1, z1, x1 );
    line->SetLineColor(kRed);
    line->SetLineWidth(1);
    line->Draw();
    m_BcOutXZ_line.push_back( line );
  }
  {
    m_canvas_vertex->cd(2);
    TPolyLine *line = new TPolyLine(2);
    line->SetPoint( 0, z0, y0 );
    line->SetPoint( 1, z1, y1 );
    line->SetLineColor(kRed);
    line->SetLineWidth(1);
    line->Draw();
    m_BcOutYZ_line.push_back( line );
  }
#endif

#if CATCH
  double zc0 = -500.;
  double xc0 = tp->GetX( zc0 + zK18Target ); 
  double yc0 = tp->GetY( zc0 + zK18Target );
  double zc1 = 250.;
  double xc1 = tp->GetX( zc1 + zK18Target ); 
  double yc1 = tp->GetY( zc1 + zK18Target );
  {
    TPolyLine3D *p = new TPolyLine3D(2);
    p->SetLineColor(kRed);
    p->SetLineWidth(1);
    p->SetPoint( 0, xc0, yc0, zc0+offsetCATCH );
    p->SetPoint( 1, xc1, yc1, zc1+offsetCATCH );
    m_BcOutTrack_Catch_cont.push_back(p);

    TPolyLine *lxy = new TPolyLine(2);
    lxy->SetPoint( 0, xc0, yc0 );
    lxy->SetPoint( 1, xc1, yc1 );
    lxy->SetLineColor(kRed);
    lxy->SetLineWidth(1);
    m_BcOutTrack_Catch_xy_cont.push_back(lxy);

    TPolyLine *lzx = new TPolyLine(2);
    lzx->SetPoint( 0, zc0+offsetCATCH, xc0 );
    lzx->SetPoint( 1, zc1+offsetCATCH, xc1 );
    lzx->SetLineColor(kRed);
    lzx->SetLineWidth(1);
    m_BcOutTrack_Catch_zx_cont.push_back(lzx);

    TPolyLine *lzy = new TPolyLine(2);
    lzy->SetPoint( 0, zc0+offsetCATCH, xc0 );
    lzy->SetPoint( 1, zc1+offsetCATCH, xc1 );
    lzy->SetLineColor(kRed);
    lzy->SetLineWidth(1);
    m_BcOutTrack_Catch_zy_cont.push_back(lzy);


  }
#endif

}

//______________________________________________________________________________
void
EventDisplay::DrawBcOutLocalTrack( double x0, double y0, double u0, double v0 )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  const double offsetZ = -3108.4;

  double z0 = offsetZ;

  double z1 = zK18Target + 100.;
  double x1 = x0+u0*z1;
  double y1 = y0+v0*z1;
  z1 += offsetZ;

  ThreeVector gPos0( x0+BeamAxis, y0, z0 );
  ThreeVector gPos1( x1+BeamAxis, y1, z1 );

  TPolyLine3D *p = new TPolyLine3D(2);
  p->SetLineColor(kRed);
  p->SetLineWidth(1);
  p->SetPoint( 0, gPos0.x(), gPos0.y(), gPos0.z() );
  p->SetPoint( 1, gPos1.x(), gPos1.y(), gPos1.z() );
  m_BcOutTrack.push_back(p);
  m_canvas->cd(2);
  p->Draw();
  gPad->Update();

#if Vertex
  double z2 = zK18Target + MinZ;
  z1 = zK18Target;
  double x2 = x0+u0*z0, y2 = y0+v0*z0;
  x1 = x0+u0*z1; y1 = y0+v0*z1;
  z2 -= zK18Target;
  z1 -= zK18Target;
  {
    TPolyLine *line = new TPolyLine(2);
    line->SetPoint( 0, z2, x2 );
    line->SetPoint( 1, z1, x1 );
    line->SetLineColor(kRed);
    line->SetLineWidth(1);
    m_BcOutXZ_line.push_back( line );
    m_canvas_vertex->cd(1);
    line->Draw();
  }
  {
    m_canvas_vertex->cd(2);
    TPolyLine *line = new TPolyLine(2);
    line->SetPoint( 0, z2, y2 );
    line->SetPoint( 1, z1, y1 );
    line->SetLineColor(kRed);
    line->SetLineWidth(1);
    line->Draw();
    m_BcOutYZ_line.push_back( line );
  }
  m_canvas_vertex->Update();
#endif

}

//______________________________________________________________________________
void
EventDisplay::DrawSdcInLocalTrack( DCLocalTrack *tp )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

#if SdcIn
  double x0 = tp->GetX0(), y0 = tp->GetY0();
  static const int lid = gGeom.GetDetectorId("SFT-U");
  static const double zSdc1x1 = gGeom.GetLocalZ( lid ) - 100.;
  double x1 = tp->GetX( zSdc1x1 ), y1 = tp->GetY( zSdc1x1 );

  //ThreeVector gPos0( x0+BeamAxis, y0, 0. );
  //ThreeVector gPos1( x1+BeamAxis, y1, zSdc1x1 );
  ThreeVector gPos0( x0, y0, 0. );
  ThreeVector gPos1( x1, y1, zSdc1x1 );

  TPolyLine3D *p = new TPolyLine3D(2);
  p->SetLineColor(kRed);
  p->SetLineWidth(1);
  p->SetPoint( 0, gPos0.x(), gPos0.y(), gPos0.z() );
  p->SetPoint( 1, gPos1.x(), gPos1.y(), gPos1.z() );
  m_SdcInTrack.push_back(p);
  m_canvas->cd(2);
  p->Draw();
  gPad->Update();
#endif

#if Vertex
  double z0 = zTarget;
  double z1 = zTarget + MaxZ;
  x0 = tp->GetX( z0 ); y0 = tp->GetY( z0 );
  x1 = tp->GetX( z1 ); y1 = tp->GetY( z1 );
  z0 -= zTarget;
  z1 -= zTarget;
  {
    m_canvas_vertex->cd(1);
    TPolyLine *line = new TPolyLine(2);
    line->SetPoint( 0, z0, x0 );
    line->SetPoint( 1, z1, x1 );
    line->SetLineColor(kRed);
    line->SetLineWidth(1);
    line->Draw();
    m_SdcInXZ_line.push_back( line );
  }
  {
    m_canvas_vertex->cd(2);
    TPolyLine *line = new TPolyLine(2);
    line->SetPoint( 0, z0, y0 );
    line->SetPoint( 1, z1, y1 );
    line->SetLineColor(kRed);
    line->SetLineWidth(1);
    line->Draw();
    m_SdcInYZ_line.push_back( line );
  }
  m_canvas_vertex->Update();
#endif

#if CATCH
  double zc0 = -150, zc1 = 400;
  double xc0 = tp->GetX( zc0 + gzK18Target ) - BeamAxis, yc0 = tp->GetY( zc0 + gzK18Target );
  double xc1 = tp->GetX( zc1 + gzK18Target ) - BeamAxis, yc1 = tp->GetY( zc1 + gzK18Target );
  {
    TPolyLine3D *p = new TPolyLine3D(2);
    p->SetLineColor(kRed);
    p->SetLineWidth(1);
    p->SetPoint( 0, xc0, yc0, zc0 + offsetCATCH );
    p->SetPoint( 1, xc1, yc1, zc1 + offsetCATCH );
    m_SdcInTrack_Catch_cont.push_back(p);

    TPolyLine *lxy = new TPolyLine(2);
    lxy->SetPoint( 0, xc0, yc0 );
    lxy->SetPoint( 1, xc1, yc1 );
    lxy->SetLineColor(kRed);
    lxy->SetLineWidth(1);
    m_SdcInTrack_Catch_xy_cont.push_back(lxy);

    TPolyLine *lzx = new TPolyLine(2);
    lzx->SetPoint( 0, zc0+offsetCATCH, xc0 );
    lzx->SetPoint( 1, zc1+offsetCATCH, xc1 );
    lzx->SetLineColor(kRed);
    lzx->SetLineWidth(1);
    m_SdcInTrack_Catch_zx_cont.push_back(lzx);

    TPolyLine *lzy = new TPolyLine(2);
    lzy->SetPoint( 0, zc0+offsetCATCH, yc0 );
    lzy->SetPoint( 1, zc1+offsetCATCH, yc1 );
    lzy->SetLineColor(kRed);
    lzy->SetLineWidth(1);
    m_SdcInTrack_Catch_zy_cont.push_back(lzy);


  }
#endif

}

//______________________________________________________________________________
void
EventDisplay::DrawSdcOutLocalTrack( DCLocalTrack *tp )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

#if SdcOut
  double x0 = tp->GetX0(), y0 = tp->GetY0();
  double zSdcOut = gGeom.GetLocalZ( "VTOF" );
  double x1 = tp->GetX( zSdcOut ), y1 = tp->GetY( zSdcOut );

  ThreeVector gPos0( x0, y0, 0. );
  ThreeVector gPos1( x1, y1, zSdcOut );

  TPolyLine3D *p = new TPolyLine3D(2);
  p->SetLineColor(kRed);
  p->SetLineWidth(1);
  p->SetPoint( 0, gPos0.x(), gPos0.y(), gPos0.z() );
  p->SetPoint( 1, gPos1.x(), gPos1.y(), gPos1.z() );
  m_SdcOutTrack.push_back(p);
  m_canvas->cd(2);
  p->Draw();
  gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawVertex( const ThreeVector& vertex )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

#if Vertex
  del::DeleteObject( m_VertexPointXZ );
  del::DeleteObject( m_VertexPointYZ );

  double x = vertex.x();
  double y = vertex.y();
  double z = vertex.z();

  m_VertexPointXZ = new TMarker( z, x, 34 );
  m_VertexPointXZ->SetMarkerSize(1);
  m_VertexPointXZ->SetMarkerColor(kBlue);
  m_VertexPointYZ = new TMarker( z, y, 34 );
  m_VertexPointYZ->SetMarkerSize(1);
  m_VertexPointYZ->SetMarkerColor(kBlue);

  m_canvas_vertex->cd(1);
  m_VertexPointXZ->Draw();
  m_canvas_vertex->cd(2);
  m_VertexPointYZ->Draw();
  m_canvas_vertex->Update();
  hddaq::cout <<"Draw Vertex! " << x << " " << y << " " << z << std::endl;
#endif

}

//______________________________________________________________________________
void
EventDisplay::DrawMissingMomentum( const ThreeVector& mom, const ThreeVector& pos )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

#if Vertex
  if( std::abs( pos.x() )>40. ) return;
  if( std::abs( pos.y() )>20. ) return;
  if( std::abs( pos.z() )>20. ) return;

  double z0 = pos.z();
  double x0 = pos.x();
  double y0 = pos.y();
  double u0 = mom.x()/mom.z();
  double v0 = mom.y()/mom.z();
  double z1 = MaxZ;
  double x1 = x0 + u0*z1;
  double y1 = y0 + v0*z1;
  m_canvas_vertex->cd(1);
  m_MissMomXZ_line = new TPolyLine(2);
  m_MissMomXZ_line->SetPoint( 0, z0, x0 );
  m_MissMomXZ_line->SetPoint( 1, z1, x1 );
  m_MissMomXZ_line->SetLineColor(kBlue);
  m_MissMomXZ_line->SetLineWidth(1);
  m_MissMomXZ_line->Draw();
  m_canvas_vertex->cd(2);
  m_MissMomYZ_line = new TPolyLine(2);
  m_MissMomYZ_line->SetPoint( 0, z0, y0 );
  m_MissMomYZ_line->SetPoint( 1, z1, y1 );
  m_MissMomYZ_line->SetLineColor(kBlue);
  m_MissMomYZ_line->SetLineWidth(1);
  m_MissMomYZ_line->Draw();
  m_canvas_vertex->Update();
  ::sleep(3);
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawKuramaTrack( int nStep, ThreeVector *StepPoint, int Polarity )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  del::DeleteObject( m_kurama_step_mark );

  m_kurama_step_mark = new TPolyMarker3D( nStep );
  for( int i=0; i<nStep; ++i ){
    m_kurama_step_mark->SetPoint( i,
				  StepPoint[i].x(),
				  StepPoint[i].y(),
				  StepPoint[i].z() );
  }

  Color_t color = kBlack;
  if( Polarity>0 )
    color = kBlue;
  else
    color = kRed;

  m_kurama_step_mark->SetMarkerSize(10);
  m_kurama_step_mark->SetMarkerColor(color);
  m_kurama_step_mark->SetMarkerStyle(6);

  m_canvas->cd(2);
  m_kurama_step_mark->Draw();
  m_canvas->Update();

#if Vertex
  del::DeleteObject( m_KuramaMarkVertexX );
  m_KuramaMarkVertexX = new TPolyMarker( nStep );
  for( int i=0; i<nStep; ++i ){
    double x = StepPoint[i].x()-BeamAxis;
    double z = StepPoint[i].z()-zTarget;
    m_KuramaMarkVertexX->SetPoint( i, z, x );
  }
  m_KuramaMarkVertexX->SetMarkerSize(0.4);
  m_KuramaMarkVertexX->SetMarkerColor(color);
  m_KuramaMarkVertexX->SetMarkerStyle(6);
  m_canvas_vertex->cd(1);
  m_KuramaMarkVertexX->Draw();
  del::DeleteObject( m_KuramaMarkVertexY );
  m_KuramaMarkVertexY = new TPolyMarker( nStep );
  for( int i=0; i<nStep; ++i ){
    double y = StepPoint[i].y();
    double z = StepPoint[i].z()-zTarget;
    m_KuramaMarkVertexY->SetPoint( i, z, y );
  }
  m_KuramaMarkVertexY->SetMarkerSize(0.4);
  m_KuramaMarkVertexY->SetMarkerColor(color);
  m_KuramaMarkVertexY->SetMarkerStyle(6);
  m_canvas_vertex->cd(2);
  m_KuramaMarkVertexY->Draw();
  m_canvas_vertex->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawTarget( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( !m_target_node ){
    hddaq::cout << "#E " << func_name << " "
		<< "node is null" << std::endl;
    return;
  }

  m_target_node->SetLineColor( kMagenta );
  m_canvas->cd(2);
  m_canvas->Update();

}

//______________________________________________________________________________
void
EventDisplay::DrawText( double xpos, double ypos, const std::string& arg )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if (arg.find("Run") != std::string::npos) {
    std::cout << arg << " find " << std::endl;
    m_canvas->cd(1)->Clear();
  }
  m_canvas->cd(1);
  TLatex tex;
  //tex.SetTextSize(0.5);
  tex.SetTextSize(0.15);
  tex.SetNDC();
  tex.DrawLatex( xpos, ypos, arg.c_str() );
  m_canvas->Update();
}

//______________________________________________________________________________
void
EventDisplay::EndOfEvent( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  del::DeleteObject( m_init_step_mark );
  del::DeleteObject( m_BcOutXZ_line );
  del::DeleteObject( m_BcOutYZ_line );
  del::DeleteObject( m_SdcInXZ_line );
  del::DeleteObject( m_SdcInYZ_line );
  del::DeleteObject( m_BcInTrack );
  del::DeleteObject( m_BcOutTrack );
  del::DeleteObject( m_BcOutTrack2 );
  del::DeleteObject( m_BcOutTrack3 );
  del::DeleteObject( m_SdcInTrack );
  del::DeleteObject( m_SdcInTrack2 );
  del::DeleteObject( m_SdcOutTrack );
  del::DeleteObject( m_kurama_step_mark );
  del::DeleteObject( m_VertexPointXZ );
  del::DeleteObject( m_VertexPointYZ );
  del::DeleteObject( m_MissMomXZ_line );
  del::DeleteObject( m_MissMomYZ_line );
  del::DeleteObject( m_KuramaMarkVertexX );
  del::DeleteObject( m_KuramaMarkVertexY );

  ResetVisibility();

  ResetHist();
  ResetCATCH();
}

//______________________________________________________________________________
void
EventDisplay::ResetVisibility( TNode *& node, Color_t c )
{
  if( !node ) return;
  if( c==kWhite )
    node->SetVisibility(kFALSE);
  else
    node->SetLineColor(c);
}

//______________________________________________________________________________
void
EventDisplay::ResetVisibility( std::vector<TNode*>& node, Color_t c )
{
  const std::size_t n = node.size();
  for( std::size_t i=0; i<n; ++i ){
    ResetVisibility( node[i], c );
  }
}

//______________________________________________________________________________
void
EventDisplay::ResetVisibility( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  ResetVisibility( m_BC3x1_node );
  ResetVisibility( m_BC3x2_node );
  ResetVisibility( m_BC3v1_node );
  ResetVisibility( m_BC3v2_node );
  ResetVisibility( m_BC3u1_node );
  ResetVisibility( m_BC3u2_node );
  ResetVisibility( m_BC4u1_node );
  ResetVisibility( m_BC4u2_node );
  ResetVisibility( m_BC4v1_node );
  ResetVisibility( m_BC4v2_node );
  ResetVisibility( m_BC4x1_node );
  ResetVisibility( m_BC4x2_node );
  ResetVisibility( m_SFTu_node );
  ResetVisibility( m_SFTv_node );
  ResetVisibility( m_SFTx_node );
  ResetVisibility( m_SDC1v1_node );
  ResetVisibility( m_SDC1v2_node );
  ResetVisibility( m_SDC1x1_node );
  ResetVisibility( m_SDC1x2_node );
  ResetVisibility( m_SDC1u1_node );
  ResetVisibility( m_SDC1u2_node );
  ResetVisibility( m_SDC2x1_node );
  ResetVisibility( m_SDC2x2_node );
  ResetVisibility( m_SDC2y1_node );
  ResetVisibility( m_SDC2y2_node );
  ResetVisibility( m_SDC3y1_node );
  ResetVisibility( m_SDC3y2_node );
  ResetVisibility( m_SDC3x1_node );
  ResetVisibility( m_SDC3x2_node );
  ResetVisibility( m_BH2seg_node, kBlack );
  ResetVisibility( m_SCHseg_node, kBlack );
  ResetVisibility( m_FBT1d1seg_node);
  ResetVisibility( m_FBT1u1seg_node);
  ResetVisibility( m_FBT1d2seg_node);
  ResetVisibility( m_FBT1u2seg_node);
  ResetVisibility( m_FBT2d1seg_node);
  ResetVisibility( m_FBT2u1seg_node);
  ResetVisibility( m_FBT2d2seg_node);
  ResetVisibility( m_FBT2u2seg_node);
  ResetVisibility( m_TOFseg_node, kBlack );
  ResetVisibility( m_target_node, kBlack );

  for (int layer=0; layer<NumOfPlaneCFT/2; layer++)
    ResetVisibility( m_CFT_node_cont[layer] );

  ResetVisibility( m_BGOseg_node_cont );
  ResetVisibility( m_PiIDseg_node_cont );

}

//______________________________________________________________________________
void
EventDisplay::ResetHist(  )
{
#if Hist_Timing
  m_hist_bh2->Reset();
  m_hist_sft_x->Reset();
  m_hist_sft_u->Reset();
  m_hist_sft_v->Reset();
  m_hist_sch->Reset();
  m_hist_tof->Reset();
  m_hist_sdc1->Reset();
  m_hist_sdc1p->Reset();

  m_hist_bc3->Reset();
  m_hist_bc3p->Reset();
  m_hist_bc3u->Reset();
  m_hist_bc3up->Reset();
  m_hist_bc3v->Reset();
  m_hist_bc3vp->Reset();

  m_hist_bc4->Reset();
  m_hist_bc4p->Reset();
  m_hist_bc4u->Reset();
  m_hist_bc4up->Reset();
  m_hist_bc4v->Reset();
  m_hist_bc4vp->Reset();

  m_hist_bc3_time->Reset();
  m_hist_bc3_time->SetMaximum(-1111);
  m_hist_bc3p_time->Reset();

  m_hist_bc4_time->Reset();
  m_hist_bc4_time->SetMaximum(-1111);
  m_hist_bc4p_time->Reset();
#endif 

#if Hist_SdcOut
  m_hist_fbt1u->Reset();
  m_hist_fbt1up->Reset();
  m_hist_fbt1d->Reset();
  m_hist_fbt1dp->Reset();
  m_hist_fbt2u->Reset();
  m_hist_fbt2up->Reset();
  m_hist_fbt2d->Reset();
  m_hist_fbt2dp->Reset();

  m_hist_sdc2_l->Reset();
  m_hist_sdc2_t->Reset();
  m_hist_sdc2p_l->Reset();
  m_hist_sdc2p_t->Reset();

  m_hist_sdc2y_l->Reset();
  m_hist_sdc2y_t->Reset();
  m_hist_sdc2yp_l->Reset();
  m_hist_sdc2yp_t->Reset();

  m_hist_sdc3_l->Reset();
  m_hist_sdc3_t->Reset();
  m_hist_sdc3p_l->Reset();
  m_hist_sdc3p_t->Reset();

  m_hist_sdc3y_l->Reset();
  m_hist_sdc3y_t->Reset();
  m_hist_sdc3yp_l->Reset();
  m_hist_sdc3yp_t->Reset();

#endif

#if Hist_Timing
  m_hist_bh1->Reset();
  m_hist_bft->Reset();
  m_hist_bft_p->Reset();
  m_hist_bcIn->Reset();
  m_hist_bcOut->Reset();

  int nc=m_BH1box_cont.size();
  for (int i=0; i<nc; i++)
    m_BH1box_cont[i]->SetFillColor(kWhite);

  int ncBh2=m_BH2box_cont.size();
  for (int i=0; i<ncBh2; i++) {
    m_BH2box_cont[i]->SetFillColor(kWhite);
    m_BH2box_cont2[i]->SetFillColor(kWhite);
  }

  int ncSch=m_SCHbox_cont.size();
  for (int i=0; i<ncSch; i++) 
    m_SCHbox_cont[i]->SetFillColor(kWhite);

  m_hist_bcOut_sdcIn->Reset();
  m_hist_sdcIn_predict->Reset();
  m_hist_sdcIn_predict2->Reset();
#endif

#if CATCH
  TH2 *hp_l=0, *hp_t=0;
  for (int layer=0; layer<8; layer++) {
    if (layer == 0) {
      hp_l = m_hist_cft1_l;
      hp_t = m_hist_cft1_t;
    } else if (layer == 1) {
      hp_l = m_hist_cft2_l;
      hp_t = m_hist_cft2_t;
    } else if (layer == 2) {
      hp_l = m_hist_cft3_l;
      hp_t = m_hist_cft3_t;
    } else if (layer == 3) {
      hp_l = m_hist_cft4_l;
      hp_t = m_hist_cft4_t;
    } else if (layer == 4) {
      hp_l = m_hist_cft5_l;
      hp_t = m_hist_cft5_t;
    } else if (layer == 5) {
      hp_l = m_hist_cft6_l;
      hp_t = m_hist_cft6_t;
    } else if (layer == 6) {
      hp_l = m_hist_cft7_l;
      hp_t = m_hist_cft7_t;
    } else if (layer == 7) {
      hp_l = m_hist_cft8_l;
      hp_t = m_hist_cft8_t;
    }
    hp_l->Reset();
    hp_t->Reset();
  }

  m_hist_piid_l->Reset();
  m_hist_piid_t->Reset();
  m_hist_bgo->Reset();

  TH2 *hp_hi=0, *hp_lo=0;
  for (int layer=0; layer<8; layer++) {
    if (layer == 0) {
      hp_hi = m_hist_cft1_hi;
      hp_lo = m_hist_cft1_lo;
    } else if (layer == 1) {
      hp_hi = m_hist_cft2_hi;
      hp_lo = m_hist_cft2_lo;
    } else if (layer == 2) {
      hp_hi = m_hist_cft3_hi;
      hp_lo = m_hist_cft3_lo;
    } else if (layer == 3) {
      hp_hi = m_hist_cft4_hi;
      hp_lo = m_hist_cft4_lo;
    } else if (layer == 4) {
      hp_hi = m_hist_cft5_hi;
      hp_lo = m_hist_cft5_lo;
    } else if (layer == 5) {
      hp_hi = m_hist_cft6_hi;
      hp_lo = m_hist_cft6_lo;
    } else if (layer == 6) {
      hp_hi = m_hist_cft7_hi;
      hp_lo = m_hist_cft7_lo;
    } else if (layer == 7) {
      hp_hi = m_hist_cft8_hi;
      hp_lo = m_hist_cft8_lo;
    }
    hp_hi->Reset();
    hp_lo->Reset();
  }

#endif

}

//______________________________________________________________________________

void EventDisplay::ResetCATCH( void )
{
  for (int layer=0; layer<NumOfPlaneCFT; layer++) {
    for (int seg=0; seg<NumOfSegCFT[layer]; seg++) {
      m_CFT_Arc_cont[layer][seg]->SetLineColor(kBlack);
    }
  }

  for (int seg=0; seg<NumOfSegBGO; seg++) {
    int size = m_BGO_Line_cont[seg].size();
    for (int i=0; i<size; i++)
      m_BGO_Line_cont[seg][i]->SetLineColor(kBlack);
  }

  for (int seg=0; seg<NumOfSegPiID; seg++) {
    int size = m_PiID_Line_cont[seg].size();
    for (int i=0; i<size; i++)
      m_PiID_Line_cont[seg][i]->SetLineColor(kBlack);
  }

  m_hbase_catch->Reset("ICES");

  del::DeleteObject( m_CFT_UV_cont);
  del::DeleteObject( m_BcOutTrack_Catch_cont);
  del::DeleteObject( m_SdcInTrack_Catch_cont);
  del::DeleteObject( m_CFTTrack_cont);
  del::DeleteObject( m_CFTTrack_xy_cont);
  del::DeleteObject( m_CFTTrack_zx_cont);
  del::DeleteObject( m_CFTTrack_zy_cont);
  del::DeleteObject( m_BcOutTrack_Catch_xy_cont);
  del::DeleteObject( m_BcOutTrack_Catch_zx_cont);
  del::DeleteObject( m_BcOutTrack_Catch_zy_cont);
  del::DeleteObject( m_SdcInTrack_Catch_xy_cont);
  del::DeleteObject( m_SdcInTrack_Catch_zx_cont);
  del::DeleteObject( m_SdcInTrack_Catch_zy_cont);
}

//______________________________________________________________________________
void
EventDisplay::DrawMomentum( double momentum )
{
#if Hist
  m_hist_p->Fill( momentum );
  m_canvas_hist->cd(1);
  gPad->Modified();
  gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawMassSquare( double mass_square )
{
#if Hist
  m_hist_m2->Fill( mass_square );
  m_canvas_hist->cd(2);
  gPad->Modified();
  gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawMissMass( double missmass )
{
#if Hist
  m_hist_missmass->Fill( missmass );
  m_canvas_hist->cd(3);
  gPad->Modified();
  gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawBH1( int seg, int tdc )
{
#if Hist_Timing
  double p0 = gHodo.GetOffset(DetIdBH1, 0, seg, 0);
  double p1 = gHodo.GetGain(DetIdBH1, 0, seg, 0);

  m_hist_bh1->Fill( seg, p1*((double)tdc-p0) );
#endif
}

//______________________________________________________________________________
void
EventDisplay::SetCorrectTimeBH1( int seg, double de )
{
#if Hist_BcIn
  int color;
  if (de<0.5)
    color = kBlue;
  else if (de >= 0.5 && de <= 1.5)
    color = kOrange;
  else
    color = kRed;
  m_BH1box_cont[seg]->SetFillColor(color);
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawBFT( int layer, int seg, int tdc )
{
#if Hist_Timing
  double p0 = gHodo.GetOffset(DetIdBFT, layer, seg, 0);
  double p1 = gHodo.GetGain(DetIdBFT, layer, seg, 0);

  TH2 *hp=0;

  if (layer == 0)
    hp = m_hist_bft;
  else
    hp = m_hist_bft_p;

  hp->Fill( seg, p1*((double)tdc-p0) );
  //m_canvas_hist2->cd(1);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawCFT( int layer, int seg, int LorT, int tdc )
{
#if CATCH
  //std::cout << "DrawCFT , layer = " << layer << ", seg = " << seg << ", LorT = " << LorT 
  //<< ", tdc = " << tdc << std::endl;

  double p0 = gHodo.GetOffset(DetIdCFT, layer, seg, 0);
  double p1 = gHodo.GetGain(DetIdCFT, layer, seg, 0);

  //std::cout << "p0 = " << p0 << ", p1 = " << p1 << std::endl;
  TH2 *hp=0;

  if (layer == 0 && LorT == 0) {
    hp = m_hist_cft1_l;
  } else if (layer == 0 && LorT == 1) {
    hp = m_hist_cft1_t;
  } else if (layer == 1 && LorT == 0) {
    hp = m_hist_cft2_l;
  } else if (layer == 1 && LorT == 1) {
    hp = m_hist_cft2_t;
  } else if (layer == 2 && LorT == 0) {
    hp = m_hist_cft3_l;
  } else if (layer == 2 && LorT == 1) {
    hp = m_hist_cft3_t;
  } else if (layer == 3 && LorT == 0) {
    hp = m_hist_cft4_l;
  } else if (layer == 3 && LorT == 1) {
    hp = m_hist_cft4_t;
  } else if (layer == 4 && LorT == 0) {
    hp = m_hist_cft5_l;
  } else if (layer == 4 && LorT == 1) {
    hp = m_hist_cft5_t;
  } else if (layer == 5 && LorT == 0) {
    hp = m_hist_cft6_l;
  } else if (layer == 5 && LorT == 1) {
    hp = m_hist_cft6_t;
  } else if (layer == 6 && LorT == 0) {
    hp = m_hist_cft7_l;
  } else if (layer == 6 && LorT == 1) {
    hp = m_hist_cft7_t;
  } else if (layer == 7 && LorT == 0) {
    hp = m_hist_cft8_l;
  } else if (layer == 7 && LorT == 1) {
    hp = m_hist_cft8_t;
  }


  //hp->Fill( seg, p1*((double)tdc-p0) );
  hp->Fill( seg, (double)tdc );
  //m_canvas_hist2->cd(1);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawCFT_Adc( int layer, int seg, int HorL, int adc )
{
#if CATCH
  //std::cout << "DrawCFT , layer = " << layer << ", seg = " << seg << ", HorL = " << HorL 
  //<< ", adc = " << adc << std::endl;

  double p0 = gHodo.GetP0(DetIdCFT, layer, seg, 0);
  double p1 = gHodo.GetP1(DetIdCFT, layer, seg, 0);

  //std::cout << "p0 = " << p0 << ", p1 = " << p1 << std::endl;
  TH2 *hp=0;

  if (layer == 0 && HorL == 0) {
    hp = m_hist_cft1_hi;
  } else if (layer == 0 && HorL == 1) {
    hp = m_hist_cft1_lo;
  } else if (layer == 1 && HorL == 0) {
    hp = m_hist_cft2_hi;
  } else if (layer == 1 && HorL == 1) {
    hp = m_hist_cft2_lo;
  } else if (layer == 2 && HorL == 0) {
    hp = m_hist_cft3_hi;
  } else if (layer == 2 && HorL == 1) {
    hp = m_hist_cft3_lo;
  } else if (layer == 3 && HorL == 0) {
    hp = m_hist_cft4_hi;
  } else if (layer == 3 && HorL == 1) {
    hp = m_hist_cft4_lo;
  } else if (layer == 4 && HorL == 0) {
    hp = m_hist_cft5_hi;
  } else if (layer == 4 && HorL == 1) {
    hp = m_hist_cft5_lo;
  } else if (layer == 5 && HorL == 0) {
    hp = m_hist_cft6_hi;
  } else if (layer == 5 && HorL == 1) {
    hp = m_hist_cft6_lo;
  } else if (layer == 6 && HorL == 0) {
    hp = m_hist_cft7_hi;
  } else if (layer == 6 && HorL == 1) {
    hp = m_hist_cft7_lo;
  } else if (layer == 7 && HorL == 0) {
    hp = m_hist_cft8_hi;
  } else if (layer == 7 && HorL == 1) {
    hp = m_hist_cft8_lo;
  }


  hp->Fill( seg, ((double)adc-p0) );
  //m_canvas_hist2->cd(1);
  //gPad->Modified();
  //gPad->Update();
#endif
}

void
EventDisplay::DrawCFT_AdcCor( int layer, int seg, int HorL, int adc )
{
#if CATCH
  //std::cout << "DrawCFT , layer = " << layer << ", seg = " << seg << ", HorL = " << HorL 
  //<< ", adc = " << adc << std::endl;

  //std::cout << "p0 = " << p0 << ", p1 = " << p1 << std::endl;
  TH2 *hp=0;

  if (layer == 0 && HorL == 0) {
    hp = m_hist_cft1_hi;
  } else if (layer == 0 && HorL == 1) {
    hp = m_hist_cft1_lo;
  } else if (layer == 1 && HorL == 0) {
    hp = m_hist_cft2_hi;
  } else if (layer == 1 && HorL == 1) {
    hp = m_hist_cft2_lo;
  } else if (layer == 2 && HorL == 0) {
    hp = m_hist_cft3_hi;
  } else if (layer == 2 && HorL == 1) {
    hp = m_hist_cft3_lo;
  } else if (layer == 3 && HorL == 0) {
    hp = m_hist_cft4_hi;
  } else if (layer == 3 && HorL == 1) {
    hp = m_hist_cft4_lo;
  } else if (layer == 4 && HorL == 0) {
    hp = m_hist_cft5_hi;
  } else if (layer == 4 && HorL == 1) {
    hp = m_hist_cft5_lo;
  } else if (layer == 5 && HorL == 0) {
    hp = m_hist_cft6_hi;
  } else if (layer == 5 && HorL == 1) {
    hp = m_hist_cft6_lo;
  } else if (layer == 6 && HorL == 0) {
    hp = m_hist_cft7_hi;
  } else if (layer == 6 && HorL == 1) {
    hp = m_hist_cft7_lo;
  } else if (layer == 7 && HorL == 0) {
    hp = m_hist_cft8_hi;
  } else if (layer == 7 && HorL == 1) {
    hp = m_hist_cft8_lo;
  }


  hp->Fill( seg, adc );
  //m_canvas_hist2->cd(1);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawPiID(int seg, int LorT, int tdc )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

#if CATCH
  std::cout << "DrawPiID, seg = " << seg << ", LorT = " << LorT 
	    << ", tdc = " << tdc << std::endl;

  double p0 = gHodo.GetOffset(DetIdPiID, 0, seg, 0);
  double p1 = gHodo.GetGain(DetIdPiID, 0, seg, 0);

  std::cout << "p0 = " << p0 << ", p1 = " << p1 << std::endl;
  TH2 *hp=0;

  if ( LorT == 0) {
    hp = m_hist_piid_l;
  } else if (LorT == 1) {
    hp = m_hist_piid_t;
  }

  hp->Fill( seg, p1*((double)tdc-p0) );

  //m_canvas_hist2->cd(1);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawBGO(int seg, int tdc )
{
#if CATCH
  std::cout << "DrawBGO, seg = " << seg 
	    << ", tdc = " << tdc << std::endl;

  double p0 = gHodo.GetOffset(DetIdBGO, 0, seg, 0);
  double p1 = gHodo.GetGain(DetIdBGO, 0, seg, 0);

  std::cout << "p0 = " << p0 << ", p1 = " << p1 << std::endl;

  m_hist_bgo->Fill( seg, p1*((double)tdc-p0) );

  //m_canvas_hist2->cd(1);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::SetBGOWaveformCanvas(int nhit )
{
#if CATCH
  m_canvas_hist9->Clear();

  if (nhit == 2)
    m_canvas_hist9->Divide(2, 1);
  else if (nhit>= 3 && nhit<=4)
    m_canvas_hist9->Divide(2, 2);
  else if (nhit>= 5 && nhit<=6)
    m_canvas_hist9->Divide(2, 3);
  else if (nhit>= 7 && nhit<=8)
    m_canvas_hist9->Divide(2, 4);
  else if (nhit>= 9 && nhit<=12)
    m_canvas_hist9->Divide(3, 4);
  else if (nhit>= 13 && nhit<=16)
    m_canvas_hist9->Divide(4, 4);
  else if (nhit>= 17 && nhit<=20)
    m_canvas_hist9->Divide(5, 4);
  else if (nhit>= 21 && nhit<=24)
    m_canvas_hist9->Divide(6, 4);


#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawBGOWaveform(int nc, int ngraph, int seg, TGraphErrors* gr )
{
#if CATCH
  m_canvas_hist9->cd(nc);

  gr->SetNameTitle(Form("gra_%d_%d", seg, ngraph),
		  Form("Segment %d", seg));

  if (ngraph == 0) {
    gr->SetMarkerSize(1.);
    gr->SetMarkerStyle(20);
    gr->Draw("ap");
  } else {
    gr->SetMarkerSize(1.);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kRed);
    gr->Draw("p");

  }

#endif
}


//______________________________________________________________________________
void
EventDisplay::DrawBGOFitFunc(int nc, int ngraph, int seg, TF1* func )
{
#if CATCH
  m_canvas_hist9->cd(nc);

  func->Draw("same");

#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawBcInTrack( double x0, double u0 )
{
#if Hist_Timing

  double z1 = -150, z2 = 50;
  TLine *l = new TLine(x0+u0*(z1-zBFT), z1, x0+u0*(z2-zBFT), z2);
  m_BcInTrack.push_back(l);
#endif
}

//______________________________________________________________________________
void
EventDisplay::SetCorrectTimeBFT( double pos )
{
#if Hist_BcIn
  m_hist_bcIn->Fill(pos, zBFT);
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawBH2( int seg, int tdc )
{
#if Hist_Timing
  double p0 = gHodo.GetOffset(DetIdBH2, 0, seg, 0);
  double p1 = gHodo.GetGain(DetIdBH2, 0, seg, 0);

  m_hist_bh2->Fill( seg, p1*((double)tdc-p0) );
  //m_canvas_hist2->cd(1);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::SetCorrectTimeBH2( int seg, double de )
{
#if Hist_BcIn
  int color;
  if (de<0.5)
    color = kBlue;
  else if (de >= 0.5 && de <= 1.5)
    color = kOrange;
  else
    color = kRed;
  m_BH2box_cont[seg]->SetFillColor(color);
  m_BH2box_cont2[seg]->SetFillColor(color);
#endif
}

//______________________________________________________________________________
void
EventDisplay::SetCorrectTimeBcOut( int layer, double pos )
{
#if Hist_BcIn
  double z = gGeom.GetLocalZ(layer);
  m_hist_bcOut->Fill(pos, z);

  m_hist_bcOut_sdcIn->Fill(pos + gGeom.GetGlobalPosition(layer).x(), 
			   gGeom.GetGlobalPosition(layer).z());
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawBcOutTrack( double x0, double u0, double y0, double v0, bool flagGoodForTracking  )
{
#if Hist_BcIn

  double z1 = 0, z2 = 600;
  TLine *l = new TLine(x0+u0*z1, z1, x0+u0*z2, z2);
  if (! flagGoodForTracking )
    l->SetLineColor(kOrange);

  m_BcOutTrack2.push_back(l);

  double z3 = -50, z4 = 2950;
  TLine *l2 = new TLine(x0+u0*z3 + gxK18Target, z3 - (zK18Target-gzK18Target), 
			x0+u0*z4 + gxK18Target, z4 - (zK18Target-gzK18Target));
  if (! flagGoodForTracking )
    l2->SetLineColor(kOrange);

  m_BcOutTrack3.push_back(l2);

  TLine *l3 = new TLine(y0+v0*z3, z3 - (zK18Target-gzK18Target), 
			y0+v0*z4, z4 - (zK18Target-gzK18Target));
  if (! flagGoodForTracking )
    l3->SetLineColor(kOrange);

  m_BcOutTrack3.push_back(l3);

  for (int layer=1; layer<=9; layer++) {
    double z_bcOut = gGeom.GetLocalZ(layer) - gzK18Target + zK18Target;
    double x = x0 + u0*z_bcOut + gxK18Target;
    double y = y0 + v0*z_bcOut;
    double z = gGeom.GetGlobalPosition(layer).z();

    ThreeVector gloPos(x, y, z);
    ThreeVector localPos = gGeom.Global2LocalPos(layer, gloPos);
    m_hist_sdcIn_predict->Fill(localPos.x(), z);

  }



#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawSdcInTrack( double x0, double u0, double y0, double v0, bool flagKurama, bool flagBeam )
{
#if Hist_BcIn

  double z1 = -1500, z2 = 0;
  TLine *l = new TLine(x0+u0*z1, z1, x0+u0*z2, z2);
  if (flagKurama)
    l->SetLineColor(kRed);
  if (flagBeam)
    l->SetLineColor(kOrange);
  m_SdcInTrack2.push_back(l);

  TLine *l2 = new TLine(y0+v0*z1, z1, y0+v0*z2, z2);
  if (flagKurama)
    l2->SetLineColor(kRed);
  if (flagBeam)
    l2->SetLineColor(kOrange);

  m_SdcInTrack2.push_back(l2);


  for (int layer=1; layer<=9; layer++) {
    double z = gGeom.GetLocalZ(layer);
    double x = x0 + u0*z;
    double y = y0 + v0*z;
    double gz = gGeom.GetGlobalPosition(layer).z();

    ThreeVector gloPos(x, y, gz);
    ThreeVector localPos = gGeom.Global2LocalPos(layer, gloPos);
    m_hist_sdcIn_predict2->Fill(localPos.x(), z);

  }

#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawSFT( int layer, int seg, int tdc )
{
#if Hist_Timing
  double p0 = gHodo.GetOffset(DetIdSFT, layer, seg, 0);
  double p1 = gHodo.GetGain(DetIdSFT, layer, seg, 0);

  TH2 *hp=0;

  if (layer == 0)
    hp = m_hist_sft_u;
  else if (layer == 1)
    hp = m_hist_sft_v;
  else 
    hp = m_hist_sft_x;

  hp->Fill( seg, p1*((double)tdc-p0) );
  //m_canvas_hist2->cd(1);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawSFT_X( int seg, int tdc )
{
#if Hist_Timing
  double p0 = gHodo.GetOffset(DetIdSFT, 2, seg, 0);
  double p1 = gHodo.GetGain(DetIdSFT, 2, seg, 0);

  m_hist_sft_x->Fill( seg, p1*((double)tdc-p0) );
  //m_canvas_hist2->cd(1);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawSFT_U( int seg, int tdc )
{
#if Hist_Timing
  double p0 = gHodo.GetOffset(DetIdSFT, 1, seg, 0);
  double p1 = gHodo.GetGain(DetIdSFT, 1, seg, 0);

  m_hist_sft_u->Fill( seg, p1*((double)tdc-p0) );
  //m_canvas_hist2->cd(1);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawSFT_V( int seg, int tdc )
{
#if Hist_Timing
  double p0 = gHodo.GetOffset(DetIdSFT, 0, seg, 0);
  double p1 = gHodo.GetGain(DetIdSFT, 0, seg, 0);

  m_hist_sft_v->Fill( seg, p1*((double)tdc-p0) );
  //m_canvas_hist2->cd(1);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawSCH( int seg, int tdc )
{
#if Hist_Timing
  double p0 = gHodo.GetOffset(DetIdSCH, 0, seg, 0);
  double p1 = gHodo.GetGain(DetIdSCH, 0, seg, 0);

  m_hist_sch->Fill( seg, p1*((double)tdc-p0) );
  //m_canvas_hist2->cd(2);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::SetCorrectTimeSCH( int seg, double de )
{
#if Hist_BcIn
  int color;
  if (de<30)
    color = kBlue;
  else if (de >= 30 && de <= 60)
    color = kOrange;
  else
    color = kRed;
  m_SCHbox_cont[seg]->SetFillColor(color);
#endif
}

//______________________________________________________________________________
void
EventDisplay::SetCorrectTimeSdcIn( int layer, double pos )
{
#if Hist_BcIn
  double z = gGeom.GetLocalZ(layer);
  m_hist_bcOut_sdcIn->Fill(pos, z);

#endif
}

//______________________________________________________________________________
void
  EventDisplay::DrawFBT( int detId, int layer, int UorD, int seg, int tdc )
{
#if Hist_SdcOut
  TH2 *hp=0;

  if (detId == 0 && layer == 0 && UorD == 0)
    hp = m_hist_fbt1u;
  else if (detId == 0 && layer == 1 && UorD == 0)
    hp = m_hist_fbt1up;
  else if (detId == 0 && layer == 0 && UorD == 1)
    hp = m_hist_fbt1d;
  else if (detId == 0 && layer == 1 && UorD == 1)
    hp = m_hist_fbt1dp;
  else if (detId == 1 && layer == 0 && UorD == 0)
    hp = m_hist_fbt2u;
  else if (detId == 1 && layer == 1 && UorD == 0)
    hp = m_hist_fbt2up;
  else if (detId == 1 && layer == 0 && UorD == 1)
    hp = m_hist_fbt2d;
  else if (detId == 1 && layer == 1 && UorD == 1)
    hp = m_hist_fbt2dp;

  double p0_FBT = 0;
  double p1_FBT = 0;
  if (detId == 0) {
    p0_FBT = gHodo.GetOffset(DetIdFBT1, layer, seg, UorD);
    p1_FBT = gHodo.GetGain(DetIdFBT1, layer, seg, UorD);
  } else {
    p0_FBT = gHodo.GetOffset(DetIdFBT2, layer, seg, UorD);
    p1_FBT = gHodo.GetGain(DetIdFBT2, layer, seg, UorD);
  }
  //std::cout << "p0 = " << p0_FBT << ", p1 = " << p1_FBT << std::endl;

  hp->Fill( seg, p1_FBT*(tdc-p0_FBT));
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawTOF( int seg, int tdc )
{
#if Hist_Timing
  double p0 = gHodo.GetOffset(DetIdTOF, 0, seg, 0);
  double p1 = gHodo.GetGain(DetIdTOF, 0, seg, 0);

  m_hist_tof->Fill( seg, p1*((double)tdc-p0 ));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawBcOutHit(int layer,  int wire, int tdc )
{
#if Hist_Timing
  TH2 *hp=0;

  if (layer == 1)
    hp = m_hist_bc3;
  else if (layer == 2)
    hp = m_hist_bc3p;
  else if (layer == 3)
    hp = m_hist_bc3v;
  else if (layer == 4)
    hp = m_hist_bc3vp;
  else if (layer == 5)
    hp = m_hist_bc3u;
  else if (layer == 6)
    hp = m_hist_bc3up;
  else if (layer == 7)
    hp = m_hist_bc4u;
  else if (layer == 8)
    hp = m_hist_bc4up;
  else if (layer == 9)
    hp = m_hist_bc4v;
  else if (layer == 10)
    hp = m_hist_bc4vp;
  else if (layer == 11)
    hp = m_hist_bc4;
  else if (layer == 12)
    hp = m_hist_bc4p;

  double p0=0.0, p1=-0.;
  gTdc.GetParameter(layer + 112, wire, p0, p1);

  hp->Fill( wire, -p1*(tdc+p0));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawBC3( int wire, int tdc )
{
#if Hist_Timing
  m_hist_bc3->Fill( wire, -1.*(tdc-351));
  m_hist_bc3_time->Fill( -1.*(tdc-351));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawBC3p( int wire, int tdc )
{
#if Hist_Timing
  m_hist_bc3p->Fill( wire, -1.*(tdc-351));
  m_hist_bc3p_time->Fill( -1.*(tdc-351));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}


//______________________________________________________________________________
void
EventDisplay::DrawBC4( int wire, int tdc )
{
#if Hist_Timing
  m_hist_bc4->Fill( wire, -1.*(tdc-351));
  m_hist_bc4_time->Fill( -1.*(tdc-351));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawBC4p( int wire, int tdc )
{
#if Hist_Timing
  m_hist_bc4p->Fill( wire, -1.*(tdc-351));
  m_hist_bc4p_time->Fill(-1.*(tdc-351));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawSDC1( int wire, int tdc )
{
#if Hist_Timing
  m_hist_sdc1->Fill( wire, -1.*(tdc-351));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}


//______________________________________________________________________________
void
EventDisplay::DrawSDC1p( int wire, int tdc )
{
#if Hist_Timing
  m_hist_sdc1p->Fill( wire, -1.*(tdc-351));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}


//______________________________________________________________________________
void
EventDisplay::DrawSdcOutHit(int layer,  int wire, int LorT, int tdc )
{
#if Hist_SdcOut
  TH2 *hp=0;

  if (layer == 1 && LorT == 0)
    hp = m_hist_sdc2_l;
  else if (layer == 1 && LorT == 1)
    hp = m_hist_sdc2_t;
  else if (layer == 2 && LorT == 0)
    hp = m_hist_sdc2p_l;
  else if (layer == 2 && LorT == 1)
    hp = m_hist_sdc2p_t;
  else if (layer == 3 && LorT == 0)
    hp = m_hist_sdc2y_l;
  else if (layer == 3 && LorT == 1)
    hp = m_hist_sdc2y_t;
  else if (layer == 4 && LorT == 0)
    hp = m_hist_sdc2yp_l;
  else if (layer == 4 && LorT == 1)
    hp = m_hist_sdc2yp_t;
  else if (layer == 5 && LorT == 0)
    hp = m_hist_sdc3y_l;
  else if (layer == 5 && LorT == 1)
    hp = m_hist_sdc3y_t;
  else if (layer == 6 && LorT == 0)
    hp = m_hist_sdc3yp_l;
  else if (layer == 6 && LorT == 1)
    hp = m_hist_sdc3yp_t;
  else if (layer == 7 && LorT == 0)
    hp = m_hist_sdc3_l;
  else if (layer == 7 && LorT == 1)
    hp = m_hist_sdc3_t;
  else if (layer == 8 && LorT == 0)
    hp = m_hist_sdc3p_l;
  else if (layer == 8 && LorT == 1)
    hp = m_hist_sdc3p_t;

  double p0=0.0, p1=-0.;
  gTdc.GetParameter(layer + 30, wire, p0, p1);

  hp->Fill( wire, -p1*(tdc+p0));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}


//______________________________________________________________________________
void
EventDisplay::DrawSDC2_Leading( int wire, int tdc )
{
#if Hist_Timing
  m_hist_sdc2_l->Fill( wire, -0.833*(tdc-890));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}


//______________________________________________________________________________
void
EventDisplay::DrawSDC2_Trailing( int wire, int tdc )
{
#if Hist_Timing
  m_hist_sdc2_t->Fill( wire, -0.833*(tdc-890));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawSDC2p_Leading( int wire, int tdc )
{
#if Hist_Timing
  m_hist_sdc2p_l->Fill( wire, -0.833*(tdc-890));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawSDC2p_Trailing( int wire, int tdc )
{
#if Hist_Timing
  m_hist_sdc2p_t->Fill( wire, -0.833*(tdc-890));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}


//______________________________________________________________________________
void
EventDisplay::DrawSDC3_Leading( int wire, int tdc )
{
#if Hist_Timing
  m_hist_sdc3_l->Fill( wire, -0.833*(tdc-885));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawSDC3_Trailing( int wire, int tdc )
{
#if Hist_Timing
  m_hist_sdc3_t->Fill( wire, -0.833*(tdc-885));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawSDC3p_Leading( int wire, int tdc )
{
#if Hist_Timing
  m_hist_sdc3p_l->Fill( wire, -0.833*(tdc-885));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawSDC3p_Trailing( int wire, int tdc )
{
#if Hist_Timing
  m_hist_sdc3p_t->Fill( wire, -0.833*(tdc-885));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//______________________________________________________________________________

void EventDisplay::ShowHitFiber(int layer, int segment, double pe)// const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  //printf("ShowHitFiber : layer %d, seg %d, pe %f\n", layer, segment, pe);

  Color_t colorPallet[5] = {kAzure, kTeal, kSpring, kOrange, kPink};
  Color_t color = kBlack;

  if (pe <= 0)
    color = kGray;
  else {
    int color_unit = 200;
    int color_index = ((int)pe/color_unit);
    int sub_color = ((int)pe%color_unit)/(color_unit/10);
    if (color_index>=5) {
      color_index = 4;
      sub_color = 10;
    } 
    color = colorPallet[color_index] + sub_color ;

  }


#if CATCH
  if (segment>=0 && segment<NumOfSegCFT[layer]) {
    double x, y;
    FiberPosPhi(layer, segment, &x, &y);
    m_hbase_catch->Fill(x, y, pe);



    m_CFT_Arc_cont[layer][segment]->SetLineColor(kRed);

    if (layer%2 == 1) {
      // phi layer
      std::string node_name;
      node_name = Form( "CFT%d_Node_%d", layer/2, segment );

      TNode *node = m_geometry_catch->GetNode( node_name.c_str() );
      if( !node ){
	hddaq::cout << "#E " << func_name << " "
		    << "no such node : " << node_name << std::endl;
	return;
      }

      node->SetVisibility(1);
      node->SetLineColor(color);
      
    } else {
      // UV layer
      const DCGeomMan & geomMan=DCGeomMan::GetInstance();

      int lnum=301+layer;

      double R     = geomMan.CalcCFTPositionR(lnum, segment);
      double Phi0 = geomMan.CalcCFTPositionPhi(lnum, segment) ;
      double slope = 0.;
      double d_phi = 5.0;
      int nStep = (int)(360/d_phi);
      if(layer==0||layer==4){
	slope = 400. /360.;
      }else if(layer==2||layer==6){
	slope = -400. /360.;
	d_phi  *= -1.;
      }

      TPolyMarker3D *uv_fiber = new TPolyMarker3D( nStep );
      for (int i=0; i<nStep; ++i) {
	double phi = Phi0 + d_phi*(double)i;
	double z = slope*(phi-Phi0);
	double x = R * cos(phi * math::Deg2Rad());
	double y = R * sin(phi * math::Deg2Rad());
	uv_fiber->SetPoint( i, x, y, z);
      }
      uv_fiber->SetMarkerSize(10);
      uv_fiber->SetMarkerColor(color);
      uv_fiber->SetMarkerStyle(6);
      m_CFT_UV_cont.push_back(uv_fiber);
    }
  }
#endif
}

//______________________________________________________________________________
void EventDisplay::ShowHitBGO(int segment, double de) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  printf("ShowHitBGO : seg %d, de %f\n", segment, de);

  Color_t colorPallet[5] = {kAzure, kTeal, kSpring, kOrange, kPink};
  Color_t color = kBlack;

  if (de <= 0)
    color = kGray;
  else {
    int color_unit = 5000;
    int color_index = ((int)de/color_unit);
    int sub_color = ((int)de%color_unit)/(color_unit/10);
    if (color_index>=5) {
      color_index = 4;
      sub_color = 10;
    } 
    color = colorPallet[color_index] + sub_color ;

  }

#if CATCH
  int size = m_BGO_Line_cont[segment].size();
  for (int i=0; i<size; i++)
    m_BGO_Line_cont[segment][i]->SetLineColor(kRed);

  double x, y;
  BGOPos(segment, &x, &y);
  m_hbase_catch->Fill(x, y, de);

  std::string node_name;
  node_name = Form( "BGOseg_node_%d", segment );
  
  TNode *node = m_geometry_catch->GetNode( node_name.c_str() );
  if( !node ){
    hddaq::cout << "#E " << func_name << " "
		<< "no such node : " << node_name << std::endl;
    return;
  }

  node->SetVisibility(1);
  node->SetLineColor(color);

#endif
}

//______________________________________________________________________________
void EventDisplay::ShowHitPiID(int segment)
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

#if CATCH
  int size = m_PiID_Line_cont[segment].size();
  for (int i=0; i<size; i++)
    m_PiID_Line_cont[segment][i]->SetLineColor(kRed);

  std::string node_name;
  node_name = Form( "PiIDseg_node_%d", segment );
  
  TNode *node = m_geometry_catch->GetNode( node_name.c_str() );
  if( !node ){
    hddaq::cout << "#E " << func_name << " "
		<< "no such node : " << node_name << std::endl;
    return;
  }

  node->SetVisibility(1);
  node->SetLineColor(kRed);

#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawCFTLocalTrack( DCLocalTrack *tp )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

#if CATCH
  ThreeVector Pos0 = tp->GetPos0();
  ThreeVector Dir = tp->GetDir();

  {
    TPolyLine3D *p = new TPolyLine3D(2);
    p->SetLineColor(kGreen);
    p->SetLineWidth(2);

    
    ThreeVector pos1 = Pos0 - 3.0*Dir; 
    p->SetPoint( 0, pos1.x(), pos1.y(), pos1.z() );
    ThreeVector pos2 = Pos0 + 2.5*Dir; 
    p->SetPoint( 1, pos2.x(), pos2.y(), pos2.z() );
    m_CFTTrack_cont.push_back(p);

    TPolyLine *lxy = new TPolyLine(2);
    lxy->SetPoint( 0, pos1.x(), pos1.y() );
    lxy->SetPoint( 1, pos2.x(), pos2.y() );
    lxy->SetLineColor(kGreen);
    lxy->SetLineWidth(1);
    m_CFTTrack_xy_cont.push_back(lxy);

    TPolyLine *lzx = new TPolyLine(2);
    lzx->SetPoint( 0, pos1.z(), pos1.x() );
    lzx->SetPoint( 1, pos2.z(), pos2.x() );
    lzx->SetLineColor(kGreen);
    lzx->SetLineWidth(1);
    m_CFTTrack_zx_cont.push_back(lzx);

    TPolyLine *lzy = new TPolyLine(2);
    lzy->SetPoint( 0, pos1.z(), pos1.y() );
    lzy->SetPoint( 1, pos2.z(), pos2.y() );
    lzy->SetLineColor(kGreen);
    lzy->SetLineWidth(1);
    m_CFTTrack_zy_cont.push_back(lzy);

  }
#endif

}


//______________________________________________________________________________
void EventDisplay::UpdateCATCH( void ) 
{
#if CATCH
  m_canvas_catch->cd(1);
  m_hbase_catch->Draw("colz");

  for (int layer=0; layer<NumOfPlaneCFT; layer++) {
    for (int seg=0; seg<NumOfSegCFT[layer]; seg++) {
      m_CFT_Arc_cont[layer][seg]->Draw("same");
    }
  }

  for (int seg=0; seg<NumOfSegBGO; seg++) {
    int size = m_BGO_Line_cont[seg].size();
    for (int i=0; i<size; i++)
      m_BGO_Line_cont[seg][i]->Draw("same");
  }

  for (int seg=0; seg<NumOfSegPiID; seg++) {
    int size = m_PiID_Line_cont[seg].size();
    for (int i=0; i<size; i++)
      m_PiID_Line_cont[seg][i]->Draw("same");
  }

  
  for (int i=0; i<m_CFTTrack_xy_cont.size(); i++)
    m_CFTTrack_xy_cont[i]->Draw();

  for (int i=0; i<m_BcOutTrack_Catch_xy_cont.size(); i++)
    m_BcOutTrack_Catch_xy_cont[i]->Draw();

  for (int i=0; i<m_SdcInTrack_Catch_xy_cont.size(); i++)
    m_SdcInTrack_Catch_xy_cont[i]->Draw();

  m_canvas_catch->cd(2);

  for (int i=0; i<m_CFTTrack_zx_cont.size(); i++)
    m_CFTTrack_zx_cont[i]->Draw();

  for (int i=0; i<m_BcOutTrack_Catch_zx_cont.size(); i++)
    m_BcOutTrack_Catch_zx_cont[i]->Draw();

  for (int i=0; i<m_SdcInTrack_Catch_zx_cont.size(); i++)
    m_SdcInTrack_Catch_zx_cont[i]->Draw();

  m_canvas_catch->cd(3);

  for (int i=0; i<m_CFTTrack_zy_cont.size(); i++)
    m_CFTTrack_zy_cont[i]->Draw();

  for (int i=0; i<m_BcOutTrack_Catch_zy_cont.size(); i++)
    m_BcOutTrack_Catch_zy_cont[i]->Draw();

  for (int i=0; i<m_SdcInTrack_Catch_zy_cont.size(); i++)
    m_SdcInTrack_Catch_zy_cont[i]->Draw();

  m_canvas_catch->cd();

  m_canvas_catch->Update();
  m_canvas_catch->Modified();

  m_canvas_catch3d->cd();
  m_geometry_catch->Draw();

  for (int i=0; i<m_CFT_UV_cont.size(); i++)
    m_CFT_UV_cont[i]->Draw();

  for (int i=0; i<m_BcOutTrack_Catch_cont.size(); i++)
    m_BcOutTrack_Catch_cont[i]->Draw();

  for (int i=0; i<m_SdcInTrack_Catch_cont.size(); i++)
    m_SdcInTrack_Catch_cont[i]->Draw();

  for (int i=0; i<m_CFTTrack_cont.size(); i++)
    m_CFTTrack_cont[i]->Draw();

  gPad->GetView()->ZoomIn();
  gPad->GetView()->ZoomIn();
  gPad->GetView()->ZoomIn();

  m_canvas_catch3d->Update();
  m_canvas_catch3d->Modified();


#endif
}

//______________________________________________________________________________
void
EventDisplay::UpdateHist( )
{
#if Hist_Timing
  for (int i=0; i<9; i++) {
    m_canvas_hist2->cd(i+1);
    gPad->Modified();
    gPad->Update();
  }

  int max1 = m_hist_bc3_time->GetMaximum();
  int max2 = m_hist_bc3p_time->GetMaximum();
  if (max2 > max1) 
    m_hist_bc3_time->SetMaximum(max2*1.1);

  max1 = m_hist_bc4_time->GetMaximum();
  max2 = m_hist_bc4p_time->GetMaximum();
  if (max2 > max1) 
    m_hist_bc4_time->SetMaximum(max2*1.1);

  for (int i=0; i<6; i++) {
    m_canvas_hist3->cd(i+1);
    gPad->Modified();
    gPad->Update();
  }
#endif

#if Hist_SdcOut
  for (int i=0; i<8; i++) {
    m_canvas_hist4->cd(i+1);
    gPad->Modified();
    gPad->Update();
  }
#endif

#if Hist_BcIn
  for (int i=0; i<4; i++) {
    m_canvas_hist5->cd(i+1);


    if (i==2) {
      int nc=m_BcInTrack.size();
      for (int i=0; i<nc; i++)
	m_BcInTrack[i]->Draw("same");
    } else if (i==3) {
      int nc=m_BcOutTrack2.size();
      for (int i=0; i<nc; i++)
	m_BcOutTrack2[i]->Draw("same");
    }


    gPad->Modified();
    gPad->Update();
  }

  {
    m_canvas_hist6->cd();
    
    int nc=m_BcOutTrack3.size();
    for (int i=0; i<nc; i++)
      m_BcOutTrack3[i]->Draw("same");

    int ncSdcIn=m_SdcInTrack2.size();
    for (int i=0; i<ncSdcIn; i++)
	m_SdcInTrack2[i]->Draw("same");

    gPad->Modified();
    gPad->Update();
  }

#endif

#if CATCH
  for (int i=0; i<12; i++) {
    m_canvas_hist7->cd(i+1);
    gPad->Modified();
    gPad->Update();
  }
  for (int i=0; i<12; i++) {
    m_canvas_hist8->cd(i+1);
    gPad->Modified();
    gPad->Update();
  }

  m_canvas_hist9->cd();
  m_canvas_hist9->Update();
  /*
  for (int i=0; i<24; i++) {
    m_canvas_hist9->cd(i+1);
    gPad->Modified();
    gPad->Update();
  }
  */

  UpdateCATCH();
#endif
}

//______________________________________________________________________________
void EventDisplay::FiberPosPhi(int layer, int seg, double *x, double *y) const
{
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  int lnum=301+layer;

  double R     = geomMan.CalcCFTPositionR(lnum, seg);
  double Phi = geomMan.CalcCFTPositionPhi(lnum, seg) ;

  *x = R * cos(Phi*math::Deg2Rad());
  *y = R * sin(Phi*math::Deg2Rad());
  /*
  std::cout << "FiberPos layer" << layer << ", seg" << seg 
	    << "(" << R * cos(Theta*Deg2Rad) 
	    << "," << R * sin(Theta*Deg2Rad)
	    << ")" << std::endl;
  */
}

//______________________________________________________________________________

void EventDisplay::BGOPos(int seg, double *x, double *y) const
{

  int UnitNum = seg/(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);
  int SegInUnit = seg%(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);

  double theta;
  double x0;
  double y0;

  if (SegInUnit==0 || SegInUnit==1 ) {
    theta = (double)UnitNum*45.;
    x0 = RadiusOfBGOSurface+BGO_Y/2;
    y0 = (double)((double)SegInUnit-0.5)*BGO_X;
  } else {
    theta = 22.5+(double)UnitNum*45.;
    x0 = RadiusOfBGOSurface2+BGO_Y/2;
    y0 = 0.;

  }


  *x = x0*cos(theta*math::Deg2Rad()) - y0*sin(theta*math::Deg2Rad());
  *y = x0*sin(theta*math::Deg2Rad()) + y0*cos(theta*math::Deg2Rad());

}

//______________________________________________________________________________
void
EventDisplay::CalcRotMatrix( double TA, double RA1, double RA2, double *rotMat )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  double ct0 = std::cos( TA*math::Deg2Rad()  );
  double st0 = std::sin( TA*math::Deg2Rad()  );
  double ct1 = std::cos( RA1*math::Deg2Rad() );
  double st1 = std::sin( RA1*math::Deg2Rad() );
  double ct2 = std::cos( RA2*math::Deg2Rad() );
  double st2 = std::sin( RA2*math::Deg2Rad() );

  double rotMat1[3][3], rotMat2[3][3];

  /* rotation matrix which is same as DCGeomRecord.cc */

#if 1
  /* new definition of RA2 */
  rotMat1[0][0] =  ct0*ct2+st0*st1*st2;
  rotMat1[0][1] = -st0*ct2+ct0*st1*st2;
  rotMat1[0][2] =  ct1*st2;

  rotMat1[1][0] =  st0*ct1;
  rotMat1[1][1] =  ct0*ct1;
  rotMat1[1][2] = -st1;

  rotMat1[2][0] = -ct0*st2+st0*st1*ct2;
  rotMat1[2][1] =  st0*st2+ct0*st1*ct2;
  rotMat1[2][2] =  ct1*ct2;

#else
  /* old definition of RA2 */
  rotMat1[0][0] =  ct0*ct2-st0*ct1*st2;
  rotMat1[0][1] = -st0*ct2-ct0*ct1*st2;
  rotMat1[0][2] =  st1*st2;

  rotMat1[1][0] =  ct0*st2+st0*ct1*ct2;
  rotMat1[1][1] = -st0*st2+ct0*ct1*ct2;
  rotMat1[1][2] = -st1*ct2;

  rotMat1[2][0] =  st0*st1;
  rotMat1[2][1] =  ct0*st1;
  rotMat1[2][2] =  ct1;
#endif

  /* rotation matrix which rotate -90 deg at x axis*/
  rotMat2[0][0] =  1.0;
  rotMat2[0][1] =  0.0;
  rotMat2[0][2] =  0.0;

  rotMat2[1][0] =  0.0;
  rotMat2[1][1] =  0.0;
  rotMat2[1][2] =  1.0;

  rotMat2[2][0] =  0.0;
  rotMat2[2][1] = -1.0;
  rotMat2[2][2] =  0.0;

  for (int i=0; i<9; i++)
    rotMat[i]=0.0;

  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
  	//rotMat[3*i+j] += rotMat1[i][k]*rotMat2[k][j];
  	rotMat[i+3*j] += rotMat1[i][k]*rotMat2[k][j];
      }
    }
  }

}

//______________________________________________________________________________
int
EventDisplay::GetCommand( void ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

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
    hddaq::cout << "q|n|p>" << std::endl;

    scanf("%c",&ch);
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
EventDisplay::Run( bool flag )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  hddaq::cout << func_name << " TApplication is running" << std::endl;

  m_theApp->Run(flag);
}
