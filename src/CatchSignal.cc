// -*- C++ -*-

#include "CatchSignal.hh"

#include <signal.h>
#include <sstream>

#include <spdlog/spdlog.h>

namespace CatchSignal
{
namespace
{
bool user_stop = false;
}

//_____________________________________________________________________________
Bool_t
Stop()
{
  return user_stop;
}

//_____________________________________________________________________________
void
Catch(Int_t sig)
{
  user_stop = true;
  std::ostringstream oss;
  oss << "[CatchSignal::Catch()] exit process by signal " << sig;
  spdlog::info(oss.str());
}

//______________________________________________________________________________
void
Set(Int_t sig)
{
  ::signal(sig, Catch);
}

}
