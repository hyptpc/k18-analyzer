#if TKIN
#include "KinematicalFit.hh"
#include <DatabasePDG.hh>

#include <std_ostream.hh>
#define SAKUMA 0
namespace{
  const std::string& class_name("KinematicalFit");
  const double covVal[6][16] = {
    // ### obtained from (p_meas[j]-p_gene[j])*(p_meas[k]-p_gene[k])
    // ###  using G4-data with TH1F(Form("cov_%d_%d_%d", i, j, k), 100, -cov_MAX, cov_MAX);
    // {beam kaon}, {Lambda}, {neutron},
    // {proton}, {proton from Lambda}, {pi- from Lambda}
    { 5.89783e-05, 3.466e-05, 7.87064e-05, 7.0336e-05,
      3.466e-05, 5.65764e-05, 7.58042e-05, 6.77393e-05,
      7.87064e-05, 7.58042e-05, 5.28858e-05, 4.73357e-05,
      7.0336e-05, 6.77393e-05, 4.73357e-05, 4.23708e-05 },
    { 0.000207257, 0.000179973, 0.000148674, 0.000199125,
      0.000179973, 0.000213038, 0.000155228, 0.000208214,
      0.000148674, 0.000155228, 0.000188014, 0.000144255,
      0.000199125, 0.000208214, 0.000144255, 0.000195619 },
    { 0.000310575, 0.000310411, 0.000313857, 0.000302319,
      0.000310411, 0.000309396, 0.000326896, 0.000312685,
      0.000313857, 0.000326896, 0.000302884, 0.000307791,
      0.000302319, 0.000312685, 0.000307791, 0.000326694 },
    { 0.000249203, 0.000218144, 0.000188326, 0.000237145,
      0.000218144, 0.000258923, 0.000189825, 0.000245481,
      0.000188326, 0.000189825, 0.000213973, 0.000186421,
      0.000237145, 0.000245481, 0.000186421, 0.000232521 },
    { 0.000211854, 0.000185869, 0.000139497, 0.000190273,
      0.000185869, 0.000205783, 0.000155463, 0.000208168,
      0.000139497, 0.000155463, 0.000167076, 0.000140599,
      0.000190273, 0.000208168, 0.000140599, 0.000192585 },
    { 2.38346e-05, 2.27798e-05, 1.05223e-05, 1.19754e-05,
      2.27798e-05, 1.52574e-05, 1.0513e-05, 1.84231e-05,
      1.05223e-05, 1.0513e-05, 1.79304e-05, 8.66373e-06,
      1.19754e-05, 1.84231e-05, 8.66373e-06, 1.05972e-05 }
  };
  const double covariance[6][4] = {
    {0.0050873,0.0051207,0.00200588,0.00179655},
    {0.0208179,0.0208186,0.0159654,0.0205903},
    {0.0308008,0.0308272,0.0239514,0.0319291},
    {0.0241764,0.0242385,0.0180298,0.0256994},
    {0.0223017,0.0223124,0.0151864,0.0202703},
    {0.00872987,0.00874526,0.00512054,0.00334967}
  };

  TMatrixD *covZero;
  TMatrixD *covParticle[6];
  const int PDG[6] = {321, 3122, 2112, pdg::kDeuteron, 2212, -211};
  TString str_particle[6] = {"Lbeam", "Llam", "Lmn", "Ld", "Lp", "Lpim"};
  const int PDG_lpn[6] = {321, 3122, 2112, 2212, 2212, -211};
  TString str_particle_lpn[6] = {"Lbeam", "Llam", "Lmn", "Lp", "Llp", "Lpim"};
  const int  target_pdg=pdg::kHe4;
  const int  target_pdg_lpn=pdg::kHe3;
  TVector3* TV_target;
}
//______________________________________________________________________________
KinematicalFit::KinematicalFit( void )
  : m_is_ready(false)
{
}
//______________________________________________________________________________
KinematicalFit::~KinematicalFit( void )
{
  if(kinfitter) delete kinfitter;
}
//______________________________________________________________________________
bool
KinematicalFit::Initialize( bool yamaga )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_ready ){
    hddaq::cerr << "#W " << func_name
                << " already initialied" << std::endl;
    return false;
  }

  if(kinfitter) delete kinfitter;
  kinfitter=new TKinFitter();
  TV_target=new TVector3(0,0,0);
  covZero = new TMatrixD(4, 4);
  covZero->Zero();
  covZero->ResizeTo(3, 3); // resize from 4x4 to 3x3
  for( int i=0; i<6; i++ ){
    covParticle[i] = new TMatrixD(4, 4);
    int n = 0;
    for( int j=0; j<4; j++ ){
      for( int k=0; k<4; k++ ){
	if( j==k ){
	  if(yamaga){
	    (*covParticle[i])[j][k] = covariance[i][j]*covariance[i][j]; // only diagonal elements
	  }else{
	    (*covParticle[i])[j][k] = covVal[i][n]; // only diagonal elements
	  }
	}
	else  (*covParticle[i])[j][k] = 0;
	n++;
      }
    }
    covParticle[i]->ResizeTo(3, 3); // resize from 4x4 to 3x3
    hddaq::cout<<"Particle: "<<str_particle[i]<<"  "<<PDG[i]<<std::endl;
    covParticle[i]->Print(); // Print all
  }

  m_is_ready = true;
  return true;
}
//______________________________________________________________________________
//--- KinFitter :: execution ---//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Naively,
//  measured values are momenta (3-vectors) of K-, p, p, pi- -> 3*4=12 
//  constraints for kinematical fit are masses of lambda and missing-neutron -> 1*2=2
//    where energy and momentum of missing-neutron is obtained from 4-momentum conservation of K- 3He -> L p n
//   => number of parameters is 12-2=10
//   => DOF is 12-10=2
//
// In the kinematical fit routine, KinFitter,
//  fitting values are momenta (3-vectors) of K-, p, p, pi-, n -> 3*5=15
//    where mass of neutron is fixed to the PDG value
//  constraints for kinematical fit are mass of lambda and 4-momentum conservation -> 1+4=5
//   => number of parameters is 15-5=10
//  measured values are momenta (3-vectors) of K-, p, p, pi- -> 3*4=12 
//   => DOF is 12-10=2
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// "Status" obtaind with kinfitter.getStatus() is as follows:
//  (see /w/e15/common/KinFitter/KinFitter/TKinFitter.cxx)
//  case -1:  statusstring = "NO FIT PERFORMED";
//  case 10:  statusstring = "RUNNING";
//  case 0:   statusstring = "CONVERGED";
//  case 1:   statusstring = "NOT CONVERGED";
//  case -10: statusstring = "ABORTED";
//  default:  statusstring = "NOT DEFINED";
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void
KinematicalFit::Ldn( TVector3 TV_meas[6] )
{
   static const std::string func_name("["+class_name+"::"+__func__+"()]");
   //--- KinFitter :: initialization ---//
   //*** definition of fit particles in cartesian coordinates ***//
   TFitParticlePxPyPz ParticleTgt = TFitParticlePxPyPz("target", "target", TV_target,
						       //						       pdg->GetParticle(target_name)->Mass(), covZero);
						       pdg::Mass(target_pdg), covZero);
   TFitParticlePxPyPz Particle[6];
   for( int i=0; i<6; i++ ){
     Particle[i] = TFitParticlePxPyPz(str_particle[i], str_particle[i], &TV_meas[i],
				      pdg::Mass(PDG[i]), covParticle[i]);
   } 
   //*** definition of constraints ***//
   // constraint :: mass of Lambda
   TFitConstraintM ConstML = TFitConstraintM("M_L", "M_L", 0, 0, pdg::Mass(PDG[1]));
   ConstML.addParticles1(&Particle[4], &Particle[5]);
   // constraint :: 4-momentum conservation
   TFitConstraintEp ConstEp[4];
   TString str_constEp[4]  = {"Px", "Py", "Pz", "E"};
   for( int i=0; i<4; i++ ){
     ConstEp[i] = TFitConstraintEp(str_constEp[i], str_constEp[i], 0, TFitConstraintEp::component(i), 0);
     ConstEp[i].addParticles1(&ParticleTgt, &Particle[0]);
     ConstEp[i].addParticles2(&Particle[2], &Particle[3], &Particle[4], &Particle[5]);
   }
   // add measured particles
   kinfitter->reset();
   kinfitter->addMeasParticles(&Particle[0], &Particle[3], &Particle[4], &Particle[5]); // K, d, p, pi-
   kinfitter->addUnmeasParticles(&Particle[2]); // n
   // add constraints
   kinfitter->addConstraint(&ConstML); // mass of Lambda
   for( int i=0; i<4; i++ ){
     kinfitter->addConstraint(&ConstEp[i]); // 4-momentum conservation
   }
   //*** perform the fit ***//
   kinfitter->setMaxNbIter(50);       // max number of iterations
   kinfitter->setMaxDeltaS(5e-5);     // max delta chi2
   kinfitter->setMaxF(1e-4);          // max sum of constraints
   kinfitter->setVerbosity(0);
   //kinfitter->setVerbosity(KFDEBUG);  // verbosity level
   kinfitter->fit();

   for( int i=0; i<6; i++ ){
     TL_kfit[i] = (*Particle[i].getCurr4Vec());
   }
   TL_kfit[1] = TL_kfit[4]+TL_kfit[5];
}

void
KinematicalFit::Lpn( TVector3 TV_meas[6] )
{
   static const std::string func_name("["+class_name+"::"+__func__+"()]");
   //--- KinFitter :: initialization ---//
   //*** definition of fit particles in cartesian coordinates ***//
   TFitParticlePxPyPz ParticleTgt = TFitParticlePxPyPz("target", "target", TV_target,
						       //						       pdg->GetParticle(target_name)->Mass(), covZero);
						       pdg::Mass(target_pdg_lpn), covZero);
   TFitParticlePxPyPz Particle[6];
   for( int i=0; i<6; i++ ){
     Particle[i] = TFitParticlePxPyPz(str_particle_lpn[i], str_particle_lpn[i], &TV_meas[i],
				      pdg::Mass(PDG_lpn[i]), covParticle[i]);
   } 
   //*** definition of constraints ***//
   // constraint :: mass of Lambda
   TFitConstraintM ConstML = TFitConstraintM("M_L", "M_L", 0, 0, pdg::Mass(PDG_lpn[1]));
   ConstML.addParticles1(&Particle[4], &Particle[5]);
   // constraint :: 4-momentum conservation
   TFitConstraintEp ConstEp[4];
   TString str_constEp[4]  = {"Px", "Py", "Pz", "E"};
   for( int i=0; i<4; i++ ){
     ConstEp[i] = TFitConstraintEp(str_constEp[i], str_constEp[i], 0, TFitConstraintEp::component(i), 0);
     ConstEp[i].addParticles1(&ParticleTgt, &Particle[0]);
     ConstEp[i].addParticles2(&Particle[2], &Particle[3], &Particle[4], &Particle[5]);
   }
   // add measured particles
   kinfitter->reset();
   kinfitter->addMeasParticles(&Particle[0], &Particle[3], &Particle[4], &Particle[5]); // K, d, p, pi-
   kinfitter->addUnmeasParticles(&Particle[2]); // n
   // add constraints
   kinfitter->addConstraint(&ConstML); // mass of Lambda
   for( int i=0; i<4; i++ ){
     kinfitter->addConstraint(&ConstEp[i]); // 4-momentum conservation
   }
   //*** perform the fit ***//
   kinfitter->setMaxNbIter(50);       // max number of iterations
   kinfitter->setMaxDeltaS(5e-5);     // max delta chi2
   kinfitter->setMaxF(1e-4);          // max sum of constraints
   kinfitter->setVerbosity(0);
   //kinfitter->setVerbosity(KFDEBUG);  // verbosity level
   kinfitter->fit();

   for( int i=0; i<6; i++ ){
     TL_kfit[i] = (*Particle[i].getCurr4Vec());
   }
   TL_kfit[1] = TL_kfit[4]+TL_kfit[5];
}
#endif
