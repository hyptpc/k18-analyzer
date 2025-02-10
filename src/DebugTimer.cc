// -*- C++ -*-

#include "DebugTimer.hh"

#include <ctime>

#include <spdlog/spdlog.h>

namespace debug
{
//_____________________________________________________________________________
Timer::Timer(const std::string& msg, bool verbose)
  : m_start(new ::timespec),
    m_stop(0),
    m_cstart(),
    m_cstop(),
    m_msg(msg),
    m_verbose(verbose)
{
  ::clock_gettime(CLOCK_REALTIME, m_start);
  const std::time_t& start = std::time(0);
  m_cstart = std::string(std::ctime(&start));
  m_cstart.pop_back();
}

//______________________________________________________________________________
Timer::~Timer()
{
  if (!m_stop){
    stop();
    if(m_verbose)
      print_verbose();
    else
      print();
  }

  delete m_stop;  m_stop = 0;
  delete m_start; m_start = 0;
}

//______________________________________________________________________________
double
Timer::sec() const
{
  if (m_start && m_stop)
    return m_stop->tv_sec - m_start->tv_sec;
  else
    return 0.;
}

//______________________________________________________________________________
double
Timer::nsec() const
{
  if (m_start && m_stop)
    return m_stop->tv_nsec - m_start->tv_nsec;
  else
    return 0.;
}

//______________________________________________________________________________
void
Timer::print(const std::string& arg) const
{
  const double sec  = m_stop->tv_sec  - m_start->tv_sec;
  const double nsec = m_stop->tv_nsec - m_start->tv_nsec;
  spdlog::info("{} {}", m_msg, arg);
  spdlog::info("{} nsec ({} sec)",
               (sec*.1e9 + nsec),
               (sec + nsec*1.e-9));
}

//______________________________________________________________________________
void
Timer::print_verbose(const std::string& arg) const
{
  const std::clock_t& time = m_stop->tv_sec - m_start->tv_sec;
  spdlog::info("{} {}", m_msg, arg);
  spdlog::info("   Process Start   : {}\r", m_cstart);
  spdlog::info("   Process Stop    : {}", m_cstop);
  spdlog::info("   Processing Time : {}:{:02}:{:02}",
               time/3600, (time/60)%60, time%60);
}

//______________________________________________________________________________
void
Timer::stop()
{
  m_stop = new ::timespec;
  ::clock_gettime(CLOCK_REALTIME, m_stop);
  const std::time_t& stop = std::time(0);
  m_cstop = std::string(std::ctime(&stop));
  m_cstop.pop_back();
}

}
