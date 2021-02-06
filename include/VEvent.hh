// -*- C++ -*-

#ifndef VEVENT_HH
#define VEVENT_HH

#include <TString.h>

//_____________________________________________________________________________
class VEvent
{
public:
  static TString ClassName( void );
                 VEvent( void );
  virtual       ~VEvent( void ) = 0;
  virtual bool   ProcessingBegin( void )  = 0;
  virtual bool   ProcessingEnd( void )    = 0;
  virtual bool   ProcessingNormal( void ) = 0;
};

//_____________________________________________________________________________
inline TString
VEvent::ClassName( void )
{
  static TString s_name( "UserEvent" );
  return s_name;
}

#endif
