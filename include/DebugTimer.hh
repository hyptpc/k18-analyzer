// -*- C++ -*-

#ifndef DEBUG_TIMER_HH
#define DEBUG_TIMER_HH

#include <string>

//_____________________________________________________________________________
namespace debug
{
class Timer
{
public:
  explicit Timer(const std::string& msg="",
                 bool verbose=true);
  ~Timer();

private:
  Timer(const Timer&);
  Timer& operator=(const Timer&);

private:
  ::timespec* m_start;
  ::timespec* m_stop;
  std::string m_cstart;
  std::string m_cstop;
  std::string m_msg;
  bool        m_verbose;

public:
  double sec() const;
  double nsec() const;
  void stop();
  void print(const std::string& arg="") const;
  void print_verbose(const std::string& arg="") const;
};
}
#endif
