/**
 *  file: EventDisplay.cc
 *  date: 2017.04.10
 *  note: Unit is mm.
 *
 */

#include "EventDisplayDST.hh"

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

#define BH2        0
#define BcOut      0
#define KURAMA     0
#define SdcIn      0
#define SdcOut     0
#define SCH        0
#define FBT        0
#define TOF        0
#define Vertex     0
#define Hist       0
#define Hist_Timing 0
#define Hist_SdcOut 0
#define Hist_BcIn   0
#define CATCH      1

namespace
{
  const std::string& class_name("EventDisplayDST");
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
  const double offsetBGO   = 61.5; // offset from CFT
}

//______________________________________________________________________________
EventDisplayDST::EventDisplayDST( void )
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
    m_hbase_catch(0),
    m_Tgt_zx(0),
    m_Tgt_zy(0)
{
}

//______________________________________________________________________________
EventDisplayDST::~EventDisplayDST( void )
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
EventDisplayDST::Initialize( void )
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
  /*
  m_geometry = new TGeometry( "evdisp","K1.8 Event Display" );

  ThreeVector worldSize( 1000., 1000., 1000. ); //mm
  new TBRIK( "world", "world", "void",
	     worldSize.x(), worldSize.y(), worldSize.z() );

  m_node = new TNode( "node", "node", "world", 0., 0., 0. );
  m_geometry->GetNode("node")->SetVisibility(0);
  */

  ConstructCATCH();


  //ResetVisibility();

  m_is_ready = true;
  return m_is_ready;
}


//______________________________________________________________________________

bool EventDisplayDST::ConstructCATCH(void)
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  m_canvas_catch = new TCanvas( "canvas_catch", "CATCH Event Display",
				1600, 900 );
  m_canvas_catch->Divide(1,4);
  m_canvas_catch->cd(1)->SetPad( 0.001,0.001,0.499,0.888);
  m_canvas_catch->cd(2)->SetPad( 0.501,0.445,0.999,0.888);
  m_canvas_catch->cd(3)->SetPad( 0.501,0.001,0.999,0.443);
  m_canvas_catch->cd(4)->SetPad( 0.001,0.890,0.999,0.999);

  m_canvas_catch->cd(1)->SetGrid(); 

  m_hbase_catch = new TH2F("hbase_catch","Event Display XY plane", 180, -180, 180, 180, -180, 180);
  m_hbase_catch->SetMaximum(200);
  m_hbase_catch->SetMinimum(-1);
  m_hbase_catch->Draw();
  m_Tgt_Arc = new TArc(0, 0, 40./2);
  m_Tgt_Arc->SetLineColor(kCyan);  
  //m_Tgt_Arc->SetFillStyle(0);
  m_Tgt_Arc->SetLineWidth(2);
  m_Tgt_Arc->Draw("same");

  m_CFRP_Arc = new TArc(0, 0, 80./2);
  m_CFRP_Arc->SetLineColor(kGray);  
  m_CFRP_Arc->SetFillStyle(0);
  m_CFRP_Arc->SetLineWidth(2);
  m_CFRP_Arc->Draw("same");

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();


  ConstructCFT();
  ConstructBGO();
  ConstructPiID();

  m_canvas_catch->cd(2)->SetGrid(); 
  m_hbase_catch_zx = new TH2F("hbase_catch_zx","Event Display ZX plane", 600, -200, 400, 180, -180, 180);
  m_hbase_catch_zx->Draw();
  m_Tgt_zx = new TBox(-150,-20,150,20);
  m_Tgt_zx->SetLineColor(kCyan);  
  m_Tgt_zx->SetFillStyle(0);
  m_Tgt_zx->SetLineWidth(1);
  m_Tgt_zx->Draw("same");
  m_CFT_zx1 = new TBox(-155,-84,-155+400,-49);
  m_CFT_zx1->SetLineColor(kGray+1);
  m_CFT_zx1->SetFillStyle(0);
  m_CFT_zx1->SetLineWidth(1);
  m_CFT_zx1->Draw("same");
  m_CFT_zx2 = new TBox(-155,49,-155+400,84);
  m_CFT_zx2->SetLineColor(kGray+1);  
  m_CFT_zx2->SetFillStyle(0);
  m_CFT_zx2->SetLineWidth(1);
  m_CFT_zx2->Draw("same");
  m_BGO_zx1 = new TBox(-155+offsetBGO,-125,-155+offsetBGO+400,-100);
  m_BGO_zx1->SetLineColor(kGray+3);
  m_BGO_zx1->SetFillStyle(0);
  m_BGO_zx1->SetLineWidth(1);
  m_BGO_zx1->Draw("same");
  m_BGO_zx2 = new TBox(-155+offsetBGO,100,-155+offsetBGO+400,125);
  m_BGO_zx2->SetLineColor(kGray+3);  
  m_BGO_zx2->SetFillStyle(0);
  m_BGO_zx2->SetLineWidth(1);
  m_BGO_zx2->Draw("same");
  

  m_canvas_catch->cd(3)->SetGrid(); 
  m_hbase_catch_zy = new TH2F("hbase_catch_zy","Event Display ZY plane", 600, -200, 400, 180, -180, 180);
  m_hbase_catch_zy->Draw();
  m_Tgt_zy = new TBox(-150,-20,150,20);
  m_Tgt_zy->SetLineColor(kCyan);  
  m_Tgt_zy->SetFillStyle(0);
  m_Tgt_zy->SetLineWidth(1);
  m_Tgt_zy->Draw("same");
  m_CFT_zy1 = new TBox(-155,-84,-155+400,-49);
  m_CFT_zy1->SetLineColor(kGray+1);  
  m_CFT_zy1->SetFillStyle(0);
  m_CFT_zy1->SetLineWidth(1);
  m_CFT_zy1->Draw("same");
  m_CFT_zy2 = new TBox(-155,49,-155+400,84);
  m_CFT_zy2->SetLineColor(kGray+1);  
  m_CFT_zy2->SetFillStyle(0);
  m_CFT_zy2->SetLineWidth(1);
  m_CFT_zy2->Draw("same");

  m_BGO_zy1 = new TBox(-155+offsetBGO,-125,-155+offsetBGO+400,-100);
  m_BGO_zy1->SetLineColor(kGray+3);
  m_BGO_zy1->SetFillStyle(0);
  m_BGO_zy1->SetLineWidth(1);
  m_BGO_zy1->Draw("same");
  m_BGO_zy2 = new TBox(-155+offsetBGO,100,-155+offsetBGO+400,125);
  m_BGO_zy2->SetLineColor(kGray+3);  
  m_BGO_zy2->SetFillStyle(0);
  m_BGO_zy2->SetLineWidth(1);
  m_BGO_zy2->Draw("same");

  return true;
}

//______________________________________________________________________________
bool EventDisplayDST::ConstructCFT(void)
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
bool EventDisplayDST::ConstructBGO(void)
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
bool EventDisplayDST::ConstructPiID(void)
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

bool EventDisplayDST::ConstructCATCH3d(void)
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
EventDisplayDST::SetVertex( const ThreeVector& vertex )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

#if CATCH

  double x = vertex.x();
  double y = vertex.y();
  double z = vertex.z();

  TMarker *vert_zx = new TMarker( z, x, 34 );
  vert_zx->SetMarkerSize(1);
  vert_zx->SetMarkerColor(kBlack);
  m_vert_CATCH_zx_cont.push_back(vert_zx);

  TMarker *vert_zy = new TMarker( z, y, 34 );
  vert_zy->SetMarkerSize(1);
  vert_zy->SetMarkerColor(kBlack);
  m_vert_CATCH_zy_cont.push_back(vert_zy);

  TMarker *vert_xy = new TMarker( x, y, 34 );
  vert_xy->SetMarkerSize(1);
  vert_xy->SetMarkerColor(kBlack);
  m_vert_CATCH_xy_cont.push_back(vert_xy);

  hddaq::cout <<"Set Vertex! " << x << " " << y << " " << z << std::endl;

#endif

}

//______________________________________________________________________________


//______________________________________________________________________________
void
EventDisplayDST::DrawText( double xpos, double ypos, const std::string& arg )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if (arg.find("Run") != std::string::npos) {
    std::cout << arg << " find " << std::endl;
    m_canvas_catch->cd(4)->Clear();
  }
  m_canvas_catch->cd(4);
  TLatex tex;
  tex.SetTextSize(0.3);
  //tex.SetTextSize(0.15);
  tex.SetNDC();
  tex.DrawLatex( xpos, ypos, arg.c_str() );
  m_canvas_catch->Update();
}

//______________________________________________________________________________
void
EventDisplayDST::EndOfEvent( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  /*
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
  */
  ResetCATCH();
}

//______________________________________________________________________________
void
EventDisplayDST::ResetVisibility( TNode *& node, Color_t c )
{
  if( !node ) return;
  if( c==kWhite )
    node->SetVisibility(kFALSE);
  else
    node->SetLineColor(c);
}

//______________________________________________________________________________
void
EventDisplayDST::ResetVisibility( std::vector<TNode*>& node, Color_t c )
{
  const std::size_t n = node.size();
  for( std::size_t i=0; i<n; ++i ){
    ResetVisibility( node[i], c );
  }
}

//______________________________________________________________________________
void
EventDisplayDST::ResetVisibility( void )
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

void EventDisplayDST::ResetCATCH( void )
{
  for (int layer=0; layer<NumOfPlaneCFT; layer++) {
    for (int seg=0; seg<NumOfSegCFT[layer]; seg++) {
      m_CFT_Arc_cont[layer][seg]->SetLineColor(kBlack);
    }
  }

  for (int seg=0; seg<NumOfSegBGO; seg++) {
    int size = m_BGO_Line_cont[seg].size();
    for (int i=0; i<size; i++){
      m_BGO_Line_cont[seg][i]->SetLineColor(kBlack);
      if(seg<2){
	m_BGO_Line_cont[seg][i]->SetLineColor(kWhite);
      }
    }
  }

  for (int seg=0; seg<NumOfSegPiID; seg++) {
    int size = m_PiID_Line_cont[seg].size();
    for (int i=0; i<size; i++){
      m_PiID_Line_cont[seg][i]->SetLineColor(kBlack);
      if(seg<3){
	m_PiID_Line_cont[seg][i]->SetLineColor(kWhite);
      }
    }
  }

  m_hbase_catch->Reset("ICES");

  del::DeleteObject( m_CFT_UV_cont);
  del::DeleteObject( m_BcOutTrack_Catch_cont);
  del::DeleteObject( m_SdcInTrack_Catch_cont);
  del::DeleteObject( m_CFTTrack_cont);
  del::DeleteObject( m_CFTTrack_xy_cont);
  del::DeleteObject( m_CFTTrack_zx_cont);
  del::DeleteObject( m_CFTTrack_zy_cont);
  del::DeleteObject( m_K18Track_Catch_xy_cont);
  del::DeleteObject( m_K18Track_Catch_zx_cont);
  del::DeleteObject( m_K18Track_Catch_zy_cont);
  del::DeleteObject( m_KuramaTrack_Catch_xy_cont);
  del::DeleteObject( m_KuramaTrack_Catch_zx_cont);
  del::DeleteObject( m_KuramaTrack_Catch_zy_cont);
  del::DeleteObject( m_vert_CATCH_xy_cont);
  del::DeleteObject( m_vert_CATCH_zx_cont);
  del::DeleteObject( m_vert_CATCH_zy_cont);
  del::DeleteObject( m_Sigma_xy_cont);
  del::DeleteObject( m_Sigma_zx_cont);
  del::DeleteObject( m_Sigma_zy_cont);

}

//______________________________________________________________________________
void EventDisplayDST::ShowHitBGO(int segment, double de) const
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
  //m_hbase_catch->Fill(x, y, de);

  /*
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
  */

#endif
}

//______________________________________________________________________________
void EventDisplayDST::ShowHitPiID(int segment)
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

#if CATCH
  int size = m_PiID_Line_cont[segment].size();
  for (int i=0; i<size; i++)
    m_PiID_Line_cont[segment][i]->SetLineColor(kRed);
  /*
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
  */
#endif
}

//______________________________________________________________________________
void
EventDisplayDST::DrawCFTLocalTrack( DCLocalTrack *tp )
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
void
EventDisplayDST::SetCFTTrack(ThreeVector& X0, ThreeVector& X1, int particle, int selected )
{
  //particle 0:others 1:proton
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

#if CATCH

  {
    /*
    TPolyLine3D *p = new TPolyLine3D(2);
    p->SetLineColor(kGreen);
    p->SetLineWidth(2);

    
    ThreeVector pos1 = Pos0 - 3.0*Dir; 
    p->SetPoint( 0, pos1.x(), pos1.y(), pos1.z() );
    ThreeVector pos2 = Pos0 + 2.5*Dir; 
    p->SetPoint( 1, pos2.x(), pos2.y(), pos2.z() );
    m_CFTTrack_cont.push_back(p);
    */

    int color;

    if(particle==0){
      color=kBlue;
    }

    if(particle==1){
      color=kRed;
    }

    if(selected==0){
      color-=10;
    }
    

    TPolyLine *lxy = new TPolyLine(2);
    lxy->SetPoint( 0, X0.x(), X0.y() );
    lxy->SetPoint( 1, X1.x(), X1.y() );
    lxy->SetLineColor(color);
    lxy->SetLineWidth(2);
    m_CFTTrack_xy_cont.push_back(lxy);

    TPolyLine *lzx = new TPolyLine(2);
    lzx->SetPoint( 0, X0.z(), X0.x() );
    lzx->SetPoint( 1, X1.z(), X1.x() );
    lzx->SetLineColor(color);
    lzx->SetLineWidth(2);
    m_CFTTrack_zx_cont.push_back(lzx);

    TPolyLine *lzy = new TPolyLine(2);
    lzy->SetPoint( 0, X0.z(), X0.y() );
    lzy->SetPoint( 1, X1.z(), X1.y() );
    lzy->SetLineColor(color);
    lzy->SetLineWidth(2);
    m_CFTTrack_zy_cont.push_back(lzy);

  }
#endif

}


//______________________________________________________________________________
void
EventDisplayDST::SetSigmaTrack(ThreeVector& X0, ThreeVector& X1, int status )
{
  //particle 0:sigma beam 1:scattered sigma
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

#if CATCH

  {
    /*
    TPolyLine3D *p = new TPolyLine3D(2);
    p->SetLineColor(kGreen);
    p->SetLineWidth(2);

    
    ThreeVector pos1 = Pos0 - 3.0*Dir; 
    p->SetPoint( 0, pos1.x(), pos1.y(), pos1.z() );
    ThreeVector pos2 = Pos0 + 2.5*Dir; 
    p->SetPoint( 1, pos2.x(), pos2.y(), pos2.z() );
    m_CFTTrack_cont.push_back(p);
    */

    int color;

    color=kGreen;

    if(status==1){
      color-=5;
    }
    

    TPolyLine *lxy = new TPolyLine(2);
    lxy->SetPoint( 0, X0.x(), X0.y() );
    lxy->SetPoint( 1, X1.x(), X1.y() );
    lxy->SetLineColor(color);
    lxy->SetLineWidth(2);
    m_Sigma_xy_cont.push_back(lxy);

    TPolyLine *lzx = new TPolyLine(2);
    lzx->SetPoint( 0, X0.z(), X0.x() );
    lzx->SetPoint( 1, X1.z(), X1.x() );
    lzx->SetLineColor(color);
    lzx->SetLineWidth(2);
    m_Sigma_zx_cont.push_back(lzx);

    TPolyLine *lzy = new TPolyLine(2);
    lzy->SetPoint( 0, X0.z(), X0.y() );
    lzy->SetPoint( 1, X1.z(), X1.y() );
    lzy->SetLineColor(color);
    lzy->SetLineWidth(2);
    m_Sigma_zy_cont.push_back(lzy);

  }
#endif

}

//______________________________________________________________________________
void
EventDisplayDST::SetK18Track(ThreeVector& X0, ThreeVector& X1, int selected )
{
  //particle 0:sigma beam 1:scattered sigma
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

#if CATCH

  {
    /*
    TPolyLine3D *p = new TPolyLine3D(2);
    p->SetLineColor(kGreen);
    p->SetLineWidth(2);

    
    ThreeVector pos1 = Pos0 - 3.0*Dir; 
    p->SetPoint( 0, pos1.x(), pos1.y(), pos1.z() );
    ThreeVector pos2 = Pos0 + 2.5*Dir; 
    p->SetPoint( 1, pos2.x(), pos2.y(), pos2.z() );
    m_CFTTrack_cont.push_back(p);
    */

    int color;

    color=kGray+2;

    if(selected==0){
      color=kGray;
    }
    
    /*    
    TPolyLine *lxy = new TPolyLine(2);
    lxy->SetPoint( 0, X0.x(), X0.y() );
    lxy->SetPoint( 1, X1.x(), X1.y() );
    lxy->SetLineColor(color);
    lxy->SetLineWidth(1);
    m_K18Track_Catch_xy_cont.push_back(lxy);
    */

    TPolyLine *lzx = new TPolyLine(2);
    lzx->SetPoint( 0, X0.z(), X0.x() );
    lzx->SetPoint( 1, X1.z(), X1.x() );
    lzx->SetLineColor(color);
    lzx->SetLineWidth(1);
    m_K18Track_Catch_zx_cont.push_back(lzx);

    TPolyLine *lzy = new TPolyLine(2);
    lzy->SetPoint( 0, X0.z(), X0.y() );
    lzy->SetPoint( 1, X1.z(), X1.y() );
    lzy->SetLineColor(color);
    lzy->SetLineWidth(1);
    m_K18Track_Catch_zy_cont.push_back(lzy);

  }
#endif

}

//______________________________________________________________________________
void
EventDisplayDST::SetKuramaTrack(ThreeVector& X0, ThreeVector& X1, int selected )
{
  //particle 0:sigma beam 1:scattered sigma
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

#if CATCH

  {
    /*
    TPolyLine3D *p = new TPolyLine3D(2);
    p->SetLineColor(kGreen);
    p->SetLineWidth(2);

    
    ThreeVector pos1 = Pos0 - 3.0*Dir; 
    p->SetPoint( 0, pos1.x(), pos1.y(), pos1.z() );
    ThreeVector pos2 = Pos0 + 2.5*Dir; 
    p->SetPoint( 1, pos2.x(), pos2.y(), pos2.z() );
    m_CFTTrack_cont.push_back(p);
    */

    int color;

    color=kGray+2;

    if(selected==0){
      color=kGray;
    }
    
    
    TPolyLine *lxy = new TPolyLine(2);
    lxy->SetPoint( 0, X0.x(), X0.y() );
    lxy->SetPoint( 1, X1.x(), X1.y() );
    lxy->SetLineColor(color);
    lxy->SetLineWidth(1);
    m_KuramaTrack_Catch_xy_cont.push_back(lxy);
    

    TPolyLine *lzx = new TPolyLine(2);
    lzx->SetPoint( 0, X0.z(), X0.x() );
    lzx->SetPoint( 1, X1.z(), X1.x() );
    lzx->SetLineColor(color);
    lzx->SetLineWidth(1);
    m_KuramaTrack_Catch_zx_cont.push_back(lzx);

    TPolyLine *lzy = new TPolyLine(2);
    lzy->SetPoint( 0, X0.z(), X0.y() );
    lzy->SetPoint( 1, X1.z(), X1.y() );
    lzy->SetLineColor(color);
    lzy->SetLineWidth(1);
    m_KuramaTrack_Catch_zy_cont.push_back(lzy);

  }
#endif

}


//______________________________________________________________________________
void EventDisplayDST::UpdateCATCH( void ) 
{
#if CATCH
  m_canvas_catch->cd(1);
  m_hbase_catch->Draw();
  m_Tgt_Arc->Draw("same");
  m_CFRP_Arc->Draw("same");

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
  /*
  for (int i=0; i<m_K18Track_Catch_xy_cont.size(); i++)
    m_K18Track_Catch_xy_cont[i]->Draw();
  */
  for (int i=0; i<m_KuramaTrack_Catch_xy_cont.size(); i++)
    m_KuramaTrack_Catch_xy_cont[i]->Draw();

  for (int i=0; i<m_vert_CATCH_xy_cont.size(); i++)
    m_vert_CATCH_xy_cont[i]->Draw();

  for (int i=0; i<m_Sigma_xy_cont.size(); i++)
    m_Sigma_xy_cont[i]->Draw();

  m_canvas_catch->cd(2);

  for (int i=0; i<m_CFTTrack_zx_cont.size(); i++)
    m_CFTTrack_zx_cont[i]->Draw();

  for (int i=0; i<m_K18Track_Catch_zx_cont.size(); i++)
    m_K18Track_Catch_zx_cont[i]->Draw();

  for (int i=0; i<m_KuramaTrack_Catch_zx_cont.size(); i++)
    m_KuramaTrack_Catch_zx_cont[i]->Draw();

  for (int i=0; i<m_vert_CATCH_zx_cont.size(); i++)
    m_vert_CATCH_zx_cont[i]->Draw();

  for (int i=0; i<m_Sigma_zx_cont.size(); i++)
    m_Sigma_zx_cont[i]->Draw();

  m_canvas_catch->cd(3);

  for (int i=0; i<m_CFTTrack_zy_cont.size(); i++)
    m_CFTTrack_zy_cont[i]->Draw();

  for (int i=0; i<m_K18Track_Catch_zy_cont.size(); i++)
    m_K18Track_Catch_zy_cont[i]->Draw();

  for (int i=0; i<m_KuramaTrack_Catch_zy_cont.size(); i++)
    m_KuramaTrack_Catch_zy_cont[i]->Draw();

  for (int i=0; i<m_vert_CATCH_zy_cont.size(); i++)
    m_vert_CATCH_zy_cont[i]->Draw();

  for (int i=0; i<m_Sigma_zy_cont.size(); i++)
    m_Sigma_zy_cont[i]->Draw();

  m_canvas_catch->cd();

  m_canvas_catch->Update();
  m_canvas_catch->Modified();
  /*
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
  */

#endif
}

//______________________________________________________________________________
void EventDisplayDST::FiberPosPhi(int layer, int seg, double *x, double *y) const
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

void EventDisplayDST::BGOPos(int seg, double *x, double *y) const
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
EventDisplayDST::CalcRotMatrix( double TA, double RA1, double RA2, double *rotMat )
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
EventDisplayDST::GetCommand( void ) const
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
EventDisplayDST::Run( bool flag )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  hddaq::cout << func_name << " TApplication is running" << std::endl;

  m_theApp->Run(flag);
}
