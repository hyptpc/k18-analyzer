// -*- C++ -*-

#include "BH1Match.hh"

#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>

#include <escape_sequence.hh>
#include <std_ostream.hh>
#include <UnpackerConfig.hh>
#include <UnpackerXMLReadDigit.hh>

#include "DetectorID.hh"
#include "FuncName.hh"
#include "Exception.hh"
#include "PrintHelper.hh"

//_____________________________________________________________________________
BH1Match::Param::Param()
  : m_xmin(0.),
    m_xmax(0.)
{
}

//_____________________________________________________________________________
BH1Match::Param::~Param()
{
}

//_____________________________________________________________________________
void
BH1Match::Param::Print() const
{
  PrintHelper helper(1, std::ios::fixed);
  static const Int_t w = 5;
  hddaq::cout << " BH1 seg " << std::right << std::setw(w) << m_seg << ":"
              << " (" << std::right
              << std::setw(w) << m_xmin << ", "
              << std::setw(w) << m_xmax << ")"
              << std::endl;
}

//_____________________________________________________________________________
BH1Match::BH1Match()
  : m_status(),
    m_param()
{
  m_status.reset();
}

//_____________________________________________________________________________
BH1Match::~BH1Match()
{
}

//_____________________________________________________________________________
Bool_t
BH1Match::Initialize(const TString& file_name)
{
  std::ifstream ifs(file_name);
  if(!ifs.is_open()){
    std::cerr << FUNC_NAME << " No such parameter file("
              << file_name << ")" << std::endl;
    return false;
  }

  TString line;
  while(ifs.good() && line.ReadLine(ifs)){
    if(line.IsNull() || line[0] == '#') continue;
    std::istringstream iss(line.Data());
    std::istream_iterator<Double_t> itr_begin(iss);
    std::istream_iterator<Double_t> itr_end;
    std::vector<Double_t> cont(itr_begin, itr_end);
    if(cont.size() != kNParam){
      std::cerr << FUNC_NAME << std::endl
		<< " Number of parameters: " << cont.size() << std::endl
		<< " Required: " << kNParam << std::endl;
      throw Exception(FUNC_NAME + " invalid parameter.");
    }

    const Double_t bh1seg = cont[kBH1Segment];
    const Double_t xmin   = cont[kXMin];
    const Double_t xmax   = cont[kXMax];

    m_param[bh1seg].m_seg = bh1seg;
    m_param[bh1seg].m_xmin = xmin;
    m_param[bh1seg].m_xmax = xmax;
  }// read line

  if(m_status[kVerbose]) Print();
  m_status.set(kReady);

  return true;
}

//_____________________________________________________________________________
Bool_t
BH1Match::Judge(Double_t bft_xpos, Double_t bh1seg)
{
  if(!m_status[kReady])
    throw Exception(FUNC_NAME + " is not initialized.");
  if(m_param.find(bh1seg) == m_param.end())
    return false;

  const Double_t xmin  = m_param.at(bh1seg).m_xmin;
  const Double_t xmax  = m_param.at(bh1seg).m_xmax;

  if(m_status[kVerbose]){
    hddaq::cout << FUNC_NAME << std::endl;
    m_param.at(bh1seg).Print();
    PrintHelper helper(2, std::ios::fixed);
    hddaq::cout << " BFT pos " << std::setw(7) << bft_xpos << " -> ";
    if(xmin < bft_xpos && bft_xpos < xmax){
      hddaq::cout << hddaq::unpacker::esc::k_green
                  << "accept"
                  << hddaq::unpacker::esc::k_default_color << std::endl;
    } else {
      hddaq::cout << hddaq::unpacker::esc::k_red
                  << "reject"
                  << hddaq::unpacker::esc::k_default_color << std::endl;
    }
  }

  return (xmin < bft_xpos && bft_xpos < xmax);
}

//_____________________________________________________________________________
void
BH1Match::Print() const
{
  hddaq::cout << FUNC_NAME << std::endl
              << "  seg: (xmin, xmax)" << std::endl;
  for(const auto& p : m_param) p.second.Print();
}
