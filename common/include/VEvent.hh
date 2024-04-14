// -*- C++ -*-

#ifndef VEVENT_HH
#define VEVENT_HH

#include <iostream>

class DetectorData;
class MCData;
class ReactionData;

//_____________________________________________________________________
class VEvent
{
public:
                VEvent( void );
  virtual      ~VEvent( void )           = 0;
  virtual bool  ProcessingBeginMC( DetectorData*, MCData*, ReactionData*, const int& )
  {
    std::cout<<std::endl;
    std::cout<<"#W VEvent::ProcessingBeginMC() not implemented !!!"<<std::endl;
    return false;
  }
  virtual bool  ProcessingBegin( void )  = 0;
  virtual bool  ProcessingEnd( void )    = 0;
  virtual bool  ProcessingNormal( void ) = 0;
  virtual bool  ProcessingScaler( int nsca, unsigned int *sca ) { return true; } // for unidaq scaler
};

#endif
