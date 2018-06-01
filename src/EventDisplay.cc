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
#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TEnv.h>
#include <TFile.h>
#include <TGeometry.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TMarker3DBox.h>
#include <TMixture.h>
#include <TNode.h>
#include <TPad.h>
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

#include <std_ostream.hh>

#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"
#include "MathTools.hh"
#include "DeleteUtility.hh"

#define BH2        1
#define BcOut      1
#define KURAMA     1
#define SdcIn      1
#define SdcOut     1
#define SCH        1
#define TOF        1
#define Vertex     1
#define Hist       0

namespace
{
  const std::string& class_name("EventDisplay");
  const DCGeomMan& gGeom = DCGeomMan::GetInstance();
  const int& IdBH2 = gGeom.DetectorId("BH2");
  const int& IdSCH = gGeom.DetectorId("SCH");
  const int& IdTOF = gGeom.DetectorId("TOF");
  const int& IdTarget    = gGeom.DetectorId("Target");
  const double& zTarget    = gGeom.LocalZ("Target");
  const double& zK18Target = gGeom.LocalZ("K18Target");

  const double BeamAxis = -150.;
#if Vertex
  const double MinX = -50.;
  const double MaxX =  50.;
  const double MinY = -50.;
  const double MaxY =  50.;
  const double MinZ = -25.;
#endif
  const double MaxZ =  50.;
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
    m_hist_vertex_x(0),
    m_hist_vertex_y(0),
    m_hist_p(0),
    m_hist_m2(0),
    m_target_node(0),
    m_kurama_inner_node(0),
    m_kurama_outer_node(0),
    m_BH2wall_node(0),
    m_SCHwall_node(0),
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
    m_MissMomYZ_line(0)
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

#if TOF
  ConstructTOF();
#endif

  m_canvas = new TCanvas( "canvas", "K1.8 Event Display",
			  1000, 1000 );
  m_canvas->Divide(1,2);
  m_canvas->cd(1)->SetPad( 0.00, 0.92, 1.00, 1.00 );
  m_canvas->cd(2)->SetPad( 0.00, 0.00, 1.00, 0.92 );
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
  m_canvas_hist->Divide(1,2);
  m_hist_p  = new TH1F( "hist1", "Momentum", 100, 0., 3. );
  m_hist_m2 = new TH1F( "hist2", "Mass Square", 100, -0.8, 2.2 );
  m_canvas_hist->cd(1)->SetGrid();
  m_hist_p->Draw();
  m_canvas_hist->cd(2)->SetGrid();
  m_hist_m2->Draw();
  // gStyle->SetOptTitle(0);
  gStyle->SetOptStat(1111110);
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
  new TBRIK( "SCHwall_brik", "SCHwall_brik", "void",
	     SCHwallX, SCHwallY, SCHwallZ);
  m_SCHwall_node = new TNode( "SCHwall_node", "SCHwall_node", "SCHwall_brik",
			      SCHwallPos.x(),
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
  new TBRIK( "TOFwall_brik", "TOFwall_brik", "void",
	     TOFwallX, TOFwallY, TOFwallZ);
  m_TOFwall_node = new TNode( "TOFwall_node", "TOFwall_node", "TOFwall_brik",
			      TOFwallPos.x(),
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
	Form( "SDC1u2_Node_%d", hit_wire ) };

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
  } else if( lid == IdSCH ){
    if( seg>=NumOfSegSCH ) return;
    node_name = Form( "SCHseg_node_%d", seg );
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
  static const int lid = gGeom.GetDetectorId("SDC1-X1");
  static const double zSdc1x1 = gGeom.GetLocalZ( lid ) - 100.;
  double x1 = tp->GetX( zSdc1x1 ), y1 = tp->GetY( zSdc1x1 );

  ThreeVector gPos0( x0+BeamAxis, y0, 0. );
  ThreeVector gPos1( x1+BeamAxis, y1, zSdc1x1 );

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

  m_canvas->cd(1)->Clear();
  TLatex tex;
  tex.SetTextSize(0.5);
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
  del::DeleteObject( m_BcOutTrack );
  del::DeleteObject( m_SdcInTrack );
  del::DeleteObject( m_SdcOutTrack );
  del::DeleteObject( m_kurama_step_mark );
  del::DeleteObject( m_VertexPointXZ );
  del::DeleteObject( m_VertexPointYZ );
  del::DeleteObject( m_MissMomXZ_line );
  del::DeleteObject( m_MissMomYZ_line );
  del::DeleteObject( m_KuramaMarkVertexX );
  del::DeleteObject( m_KuramaMarkVertexY );

  ResetVisibility();

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
  ResetVisibility( m_SCHseg_node );
  ResetVisibility( m_TOFseg_node, kBlack );
  ResetVisibility( m_target_node, kBlack );

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
void
EventDisplay::GetCommand( void ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  char ch;
  char data[100];
  static int stat   = 0;
  static int Nevent = 0;
  static int ev     = 0;

  if (stat == 1 && Nevent > 0 && ev<Nevent) {
    ev++;
    return;
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
      hddaq::cout << "Continue " << Nevent << "event" << std::endl;
      break;
    case 'p':
      m_theApp->Run(kTRUE);
      break;
    }
  }
}

//______________________________________________________________________________
void
EventDisplay::Run( bool flag )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  hddaq::cout << func_name << " TApplication is running" << std::endl;

  m_theApp->Run(flag);
}
