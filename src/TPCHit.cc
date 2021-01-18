/**
 *  file: TPCHit.cc
 *  date: 2020.12.21
 *
 */

#include "DCHit.hh"
#include "TPCHit.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>

#include <std_ostream.hh>

//#include "DCDriftParamMan.hh"
#include "DCGeomMan.hh"
//#include "DCParameters.hh"
//#include "DCTdcCalibMan.hh"
//#include "DCLTrackHit.hh"
#include "TPCLTrackHit.hh"
#include "DebugCounter.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "TPCPadHelper.hh"
#include "UserParamMan.hh"

namespace
{
  const std::string& class_name("TPCHit");
  const DCGeomMan&       gGeom  = DCGeomMan::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();
  const double zTgtTPC = -143.;
}

//__for single hit____________________________________________________________
TPCHit::TPCHit( int padid, double y, double charge )
  : m_pad(padid),
    m_charge(charge)
{
  m_pos = tpc::getPosition(padid);
  m_pos.SetY(y);
  m_layer = tpc::getLayerID(padid);
  m_row   = tpc::getRowID(padid);
  m_is_good = true;
  m_is_calculated = false;
  m_hit_xz = new DCHit(m_layer, m_row);
  m_hit_xz->SetWirePosition(m_pos.x());
  m_hit_xz->SetZ(m_pos.z());
  m_hit_xz->SetTiltAngle(0.);
  m_hit_xz->SetDummyPair();

  m_hit_yz = new DCHit(m_layer, m_row);
  m_hit_yz->SetWirePosition(m_pos.y());
  m_hit_yz->SetZ(m_pos.z());
  m_hit_yz->SetTiltAngle(90.);
  m_hit_yz->SetDummyPair();

  debug::ObjectCounter::increase(class_name);
}

//__for cluster hit____________________________________________________________
TPCHit::TPCHit( int padid, TVector3 pos, double charge )
  : m_pad(padid),
    m_charge(charge)
{
  m_pos=pos;
  m_layer = tpc::getLayerID(padid);
  m_row   = tpc::getRowID(padid);
  m_is_good = true;
  m_is_calculated = false;
  m_hit_xz = new DCHit(m_layer, m_row);
  m_hit_xz->SetWirePosition(m_pos.x());
  m_hit_xz->SetZ(m_pos.z());
  m_hit_xz->SetTiltAngle(0.);
  m_hit_xz->SetDummyPair();

  m_hit_yz = new DCHit(m_layer, m_row);
  m_hit_yz->SetWirePosition(m_pos.y());
  m_hit_yz->SetZ(m_pos.z());
  m_hit_yz->SetTiltAngle(90.);
  m_hit_yz->SetDummyPair();

  debug::ObjectCounter::increase(class_name);
}


//______________________________________________________________________________
TPCHit::~TPCHit( void )
{
  ClearRegisteredHits();
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
void
TPCHit::ClearRegisteredHits( void )
{
  int n = m_register_container.size();
  for(int i=0; i<n; ++i){
    delete m_register_container[i];
  }
  if(m_hit_xz) delete m_hit_xz;
  if(m_hit_yz) delete m_hit_yz;
}

//______________________________________________________________________________                                                              
bool
TPCHit::CalcTPCObservables( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_calculated ){
    hddaq::cout << func_name << " already calculated" << std::endl;
    return false;
  }

  return true;
}


//______________________________________________________________________________
double 
TPCHit::GetResolutionX( void )
{
  //calculated by using NIM paper
  //To do:change the resolution by checking cluster size
  double y_pos= m_pos.Y();
//double s0 = 0.204;// mm HIMAC result //To do parameter
  double s0 = gUser.GetParameter("TPC_sigma0");
  //s0 is considered for common resolution
//  double Dt = 0.18;//mm/sqrt(cm) at 1T //To do parameter
  double Dt = gUser.GetParameter("TPC_Dt");
  double L_D = 30.+(y_pos*0.1);//cm
  double N_eff = 42.8;
  double A = 0.0582*0.01;//m-1 -> cm-1     
  double e_ALD = exp(-1.*A*L_D);
  //double sT2 = s0*s0 + (Dt*Dt*L_D/(N_eff*e_ALD));
  double sT2 = (Dt*Dt*L_D/(N_eff*e_ALD));
  double sT = sqrt(sT2);

  double x_pos= m_pos.X();
  double z_pos= m_pos.Z() - zTgtTPC;
  double alpha =  atan2(x_pos,z_pos);
  double rho =  sqrt(pow(z_pos,2)+pow(x_pos,2));
  double dalpha =sT/rho;
  double smear_alpha = alpha + dalpha;
  double res_xdiff = fabs(rho*(sin(smear_alpha)-sin(alpha)));
  double res_x = sqrt(s0*s0+res_xdiff*res_xdiff);

  //std::cout<<"res_x: "<<res_x<<std::endl;
  
  return res_x;
  //return 0.2;
}

//______________________________________________________________________________
double 
TPCHit::GetResolutionZ( void )
{
  //calculated by using NIM paper
  //To do:change the resolution by checking cluster size
  double y_pos= m_pos.Y();
  double s0 = gUser.GetParameter("TPC_sigma0");
  //s0 is considered for common resolution
//  double Dt = 0.18;//mm/sqrt(cm) at 1T //To do parameter
  double Dt = gUser.GetParameter("TPC_Dt");
  double L_D = 30.+(y_pos*0.1);//cm
  double N_eff = 42.8;
  double A = 0.0582*0.01;//m-1 -> cm-1     
  double e_ALD = exp(-1.*A*L_D);
  //  double sT2 = s0*s0 + (Dt*Dt*L_D/(N_eff*e_ALD));
  double sT2 = (Dt*Dt*L_D/(N_eff*e_ALD));
  double sT = sqrt(sT2);


  double x_pos= m_pos.X();
  double z_pos= m_pos.Z() - zTgtTPC;
  double alpha =  atan2(x_pos,z_pos);
  double rho =  sqrt(pow(z_pos,2)+pow(x_pos,2));
  double dalpha =sT/rho;
  double smear_alpha = alpha + dalpha;
  double res_zdiff = fabs(rho*(cos(smear_alpha)-cos(alpha)));
  double res_z = sqrt(s0*s0 + res_zdiff*res_zdiff);

  //std::cout<<"res_z: "<<res_z<<std::endl;

  return res_z;
  //return 0.2;
}

//______________________________________________________________________________
double 
TPCHit::GetResolutionY( void )
{
  // temporary 
  return 0.5;
}

//______________________________________________________________________________
double 
TPCHit::GetResolution( void )
{
  //calculated by using NIM paper
  //To do:change the resolution by checking cluster size
  double y_pos= m_pos.Y();
  double s0 = gUser.GetParameter("TPC_sigma0");
  //s0 is considered for common resolution
//  double Dt = 0.18;//mm/sqrt(cm) at 1T //To do parameter
  double Dt = gUser.GetParameter("TPC_Dt");
  double L_D = 30.+(y_pos*0.1);//cm
  double N_eff = 42.8;
  double A = 0.0582*0.01;//m-1 -> cm-1     
  double e_ALD = exp(-1.*A*L_D);
  //  double sT2 = s0*s0 + (Dt*Dt*L_D/(N_eff*e_ALD));
  double sT2 = (Dt*Dt*L_D/(N_eff*e_ALD));
  double sT = sqrt(sT2);


  double x_pos= m_pos.X();
  double z_pos= m_pos.Z() - zTgtTPC;
  double alpha =  atan2(x_pos,z_pos);
  double rho =  sqrt(pow(z_pos,2)+pow(x_pos,2));
  double dalpha =sT/rho;
  double smear_alpha = alpha + dalpha;
  double res_x_diff = fabs(rho*(sin(smear_alpha)-sin(alpha)));
  double res_z_diff = fabs(rho*(cos(smear_alpha)-cos(alpha)));

  double res_x = sqrt(s0*s0+res_x_diff*res_x_diff);
  double res_y = 0.5;
  double res_z = sqrt(s0*s0+res_z_diff*res_z_diff);

  double tot_res = sqrt(res_x*res_x + res_y*res_y + res_z*res_z);
  //double tot_res = sqrt(0.5*0.5 + 0.5*0.5 + 0.5*0.5);
  //double tot_res = sqrt(0.2*0.2 + 0.5*0.5 + 0.2*0.2);

  return tot_res;
}


//______________________________________________________________________________
void
TPCHit::Print( const std::string& arg, std::ostream& ost ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  const int w = 16;
  ost << "#D " << func_name << " " << arg << std::endl
      << std::setw(w) << std::left << "layer" << m_layer << std::endl
      << std::setw(w) << std::left << "pad"  << m_pad  << std::endl
      << std::setw(w) << std::left << "row"  << m_row  << std::endl
      << std::setw(w) << std::left << "charge" << m_charge << std::endl
      << std::setw(w) << std::left << "posx" << m_pos.x() << std::endl
      << std::setw(w) << std::left << "posy" << m_pos.y() << std::endl
      << std::setw(w) << std::left << "posz" << m_pos.z() << std::endl;  
}
