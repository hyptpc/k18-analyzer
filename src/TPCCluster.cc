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
#include "DeleteUtility.hh"

namespace
{
  const std::string& class_name("TPCCluster");
}

#define WeightedMean 0
#define WeightedMeanTheta 1



//______________________________________________________________________________
TPCCluster::TPCCluster( int layer, const TPCHitContainer& HitCont )
  : m_layer_id( layer ),
    m_pad_id(),
    m_charge(),
    m_pos(),
    m_tpchits(),
    m_pos_calculated(false),
    m_mrow()
{
  m_tpchits.resize(HitCont.size());
  std::copy(HitCont.begin(), HitCont.end(), m_tpchits.begin());
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
TPCCluster::TPCCluster( double x, double y, double z, double de )
  : m_layer_id(tpc::getLayerID(tpc::findPadID(z, x))),
    m_pad_id(tpc::findPadID(z, x)),
    m_charge(de),
    m_pos(x, y, z),
    m_tpchits(),
    m_pos_calculated(true), // for MC data
    m_mrow()
{
  debug::ObjectCounter::increase(class_name);
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
  //  m_tpchits.clear();
  //del::ClearContainer( m_tpchits );
  //del::DeleteObject( m_tpchits );
   // int n = m_tpchits.size();
   // for(int i=0; i<n; ++i){
   //    delete m_tpchits[i];
   // }
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
TPCCluster::Calculate( void )
{
#if  WeightedMean
  CalculateWeightedMean();
#endif
#if WeightedMeanTheta
  CalculateWeightedMeanTheta();
#endif
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
    m_charge = charge;
  }
  else {
    m_pos.SetXYZ( 0, 0, 0 );
    m_pad_id = 0;
    m_mrow = 0;
  }
  m_pos_calculated = true;
}

//______________________________________________________________________________
void
TPCCluster::CalculateWeightedMeanTheta( void )
{
  if( m_pos_calculated ) return;
  [[maybe_unused]] double x=0, y=0, z=0, charge=0, dummy_padid=0, mrow=0;
  for(int hiti=0; hiti<m_tpchits.size(); hiti++) {
    // x+=m_tpchits[hiti]->GetX()*m_tpchits[hiti]->GetCharge();
    y+=m_tpchits[hiti]->GetY()*m_tpchits[hiti]->GetCharge();
    // z+=m_tpchits[hiti]->GetZ()*m_tpchits[hiti]->GetCharge();

    TVector3 pos_check = tpc::getPosition(m_tpchits[hiti]->GetLayer(), (double)m_tpchits[hiti]->GetRow());
    // std::cout<<"(x,z)=("<<m_tpchits[hiti]->GetX()
    //    	     <<","<<m_tpchits[hiti]->GetZ()<<")"<<std::endl;
    // std::cout<<"theta, (x,z)=("<<pos_check.x()
    //  	     <<","<<pos_check.z()<<")"<<std::endl;

    dummy_padid+=(double)(m_tpchits[hiti]->GetPad())*m_tpchits[hiti]->GetCharge();
    mrow+=(double)(m_tpchits[hiti]->GetRow())*m_tpchits[hiti]->GetCharge();
    charge+=m_tpchits[hiti]->GetCharge();
  }
  if( charge ){
    //m_pos.SetXYZ( x/charge, y/charge, z/charge );
    //TVector3 m_pos_dummy( x/charge, y/charge, z/charge );
    m_pad_id = (int)(dummy_padid/charge);
    m_mrow = mrow/charge;
    TVector3 pos_xz = tpc::getPosition(tpc::getLayerID(m_pad_id), m_mrow);
    m_pos.SetXYZ( pos_xz.x(), y/charge, pos_xz.z() );
    m_charge = charge;
    // std::cout<<"Weight (x,z) = ("<<m_pos_dummy.x()<<", "
    // 	     <<m_pos_dummy.z()<<")"<<std::endl;
    // std::cout<<"Theta (x,z) = ("<<m_pos.x()<<", "
    // 	     <<m_pos.z()<<")"<<std::endl;

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
  if(!m_pos_calculated) Calculate();
  return m_pos;
}

//______________________________________________________________________________
double
TPCCluster::X( void )
{
  if(!m_pos_calculated) Calculate();
  return m_pos.X();
}

//______________________________________________________________________________
double
TPCCluster::Y( void )
{
  if(!m_pos_calculated) Calculate();
  return m_pos.Y();
}

//______________________________________________________________________________
double
TPCCluster::Z( void )
{
  if(!m_pos_calculated) Calculate();
  return m_pos.Z();
}

//______________________________________________________________________________
double
TPCCluster::ResX( void )
{
  if(!m_pos_calculated) Calculate();
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
  if(!m_pos_calculated) Calculate();
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
