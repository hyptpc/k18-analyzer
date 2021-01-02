/**
 *  file: TPCCluster.cc
 *  date: 2020.04.11
 *
 */

#include "TPCCluster.hh"

#include <iostream>
#include <iterator>

#include "std_ostream.hh"

#include "DebugCounter.hh"
#include "TPCPadHelper.hh"

namespace
{
  const std::string& class_name("TPCCluster");
}

//______________________________________________________________________________
TPCCluster::TPCCluster( int layer, TPCHitContainer& HitCont )
  : m_layer_id( layer ),
    m_pos_calculated( false )
{
  std::copy( HitCont.begin(), HitCont.end(), m_tpchits.begin() );
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
TPCCluster::TPCCluster( double x, double y, double z, double de )
  : m_pos_calculated( true ) // for MC data
{
  m_pos.SetXYZ(x,y,z);
  m_charge = de;
  m_layer_id = tpc::getLayerID( tpc::findPadID(z,x) );
  m_pad_id = tpc::findPadID(z,x);
}

//______________________________________________________________________________
TPCCluster::~TPCCluster( void )
{
  ClearTPCHits();
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
void 
TPCCluster::ClearTPCHits( void )
{
  m_tpchits.clear();
  m_charge=0;
}

//______________________________________________________________________________
void 
TPCCluster::AddTPCHit( TPCHit* hit )
{
  if( hit ){
    m_tpchits.push_back(hit);
    m_charge+=hit->GetCharge();
    m_pos_calculated = false;
  }
}

//______________________________________________________________________________
void
TPCCluster::CalculateWeightedMean( void )
{
  if( m_pos_calculated ) return;
  double x=0, y=0, z=0, charge=0, dummy_padid=0, mrow=0;
  for(int hiti=0; hiti<m_tpchits.size(); hiti++) {
    x+=m_tpchits[hiti]->GetX()*m_tpchits[hiti]->GetCharge();
    y+=m_tpchits[hiti]->GetY()*m_tpchits[hiti]->GetCharge();
    z+=m_tpchits[hiti]->GetZ()*m_tpchits[hiti]->GetCharge();

    dummy_padid+=(double)(m_tpchits[hiti]->GetPad())*m_tpchits[hiti]->GetCharge();
    mrow+=(double)(m_tpchits[hiti]->GetRow())*m_tpchits[hiti]->GetCharge();
    charge+=m_tpchits[hiti]->GetCharge();
  }
  if( charge ){
    m_pos.SetXYZ( x/charge, y/charge, z/charge );
    m_pad_id = (int)(dummy_padid/charge);
    m_mrow = mrow/charge;
  }
  else {
    m_pos.SetXYZ( 0, 0, 0 ); 
    m_pad_id = 0;
    m_mrow = 0;
  }
  m_pos_calculated = true;
}

//______________________________________________________________________________
TVector3 
TPCCluster::Position( void )
{
  if(!m_pos_calculated) CalculateWeightedMean();
  return m_pos;
}

//______________________________________________________________________________
double 
TPCCluster::X( void )
{
  if(!m_pos_calculated) CalculateWeightedMean();
  return m_pos.X();
}

//______________________________________________________________________________
double 
TPCCluster::Y( void )
{
  if(!m_pos_calculated) CalculateWeightedMean();
  return m_pos.Y();
}

//______________________________________________________________________________
double 
TPCCluster::Z( void )
{
  if(!m_pos_calculated) CalculateWeightedMean();
  return m_pos.Z();
}

//______________________________________________________________________________
double 
TPCCluster::ResX( void )
{
  if(!m_pos_calculated) CalculateWeightedMean();
  //calculated by using NIM paper
  //To do:change the resolution by checking cluster size
  double y_pos= m_pos.Y();
  double s0 = 0.204;// mm HIMAC result //To do parameter
  double Dt = 0.18;//mm/sqrt(cm) at 1T //To do parameter
  double L_D = 30.+(y_pos*0.1);//cm
  double N_eff = 42.8;
  double A = 0.0582*0.01;//m-1 -> cm-1                                                                                                                                                  
  double e_ALD = exp(-1.*A*L_D);
  double sT2 = s0*s0 + (Dt*Dt*L_D/(N_eff*e_ALD));
  double sT = sqrt(sT2);

  return sT;
}

//______________________________________________________________________________
double 
TPCCluster::ResY( void )
{
  //Temp value (need to change)
  double y_res = 0.5;
  return y_res;
}

//______________________________________________________________________________
double 
TPCCluster::ResZ( void )
{
  if(!m_pos_calculated) CalculateWeightedMean();
  //calculated by using NIM paper
  //To do:change the resolution by checking cluster size
  double y_pos= m_pos.Y();
  double s0 = 0.204;// mm HIMAC result //To do parameter
  double Dt = 0.18;//mm/sqrt(cm) at 1T //To do parameter
  double L_D = 30.+(y_pos*0.1);//cm
  double N_eff = 42.8;
  double A = 0.0582*0.01;//m-1 -> cm-1                                                                                                                                                  
  double e_ALD = exp(-1.*A*L_D);
  double sT2 = s0*s0 + (Dt*Dt*L_D/(N_eff*e_ALD));
  double sT = sqrt(sT2);

  return sT;
}


//______________________________________________________________________________
void
TPCCluster::Print( const std::string& arg ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  std::cerr << func_name << " " << arg << std::endl
	      << "layerID= "  << m_layer_id  << std::endl
	      << "(x,y,z)= (" << m_pos.X() <<
	      		   ","<< m_pos.Y() <<
			   ","<< m_pos.Z()<<")" << std::endl
	      << "charge = "  << m_charge  << std::endl;
  std::cerr << std::endl;
}

