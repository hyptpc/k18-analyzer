// -*- C++ -*-

#ifndef VEVENT_HH
#define VEVENT_HH

#include <TString.h>

//_____________________________________________________________________________
class VEvent
{
public:
  static TString ClassName();
  VEvent();
  virtual ~VEvent() = 0;
  virtual Bool_t ProcessingBegin() = 0;
  virtual Bool_t ProcessingEnd() = 0;
  virtual Bool_t ProcessingNormal() = 0;
};

//_____________________________________________________________________________
inline TString
VEvent::ClassName()
{
  static TString s_name("UserEvent");
  return s_name;
}

#endif
