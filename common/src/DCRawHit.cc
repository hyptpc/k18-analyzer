// -*- C++ -*-

#include "DCRawHit.hh"

#include <iostream>
#include <iterator>

#include <UnpackerConfig.hh>
#include <UnpackerManager.hh>
#include <UnpackerXMLReadDigit.hh>
#include <std_ostream.hh>

#include "DCGeomMan.hh"
#include "DebugCounter.hh"
#include "FuncName.hh"

namespace
{
const auto& gUConf = hddaq::unpacker::GConfig::get_instance();
const auto& gGeom = DCGeomMan::GetInstance();
}

//_____________________________________________________________________________
DCRawHit::DCRawHit(const TString& detector_name, Int_t plane_id, Int_t wire_id)
  : m_detector_name(detector_name),
    m_detector_id(),
    m_plane_name(),
    m_plane_id(plane_id),
    m_dcgeom_layer(),
    m_wire_id(wire_id),
    m_tdc(),
    m_trailing(),
    m_oftdc(false)
{
  static const auto& digit_info = gUConf.get_digit_info();
  m_detector_id = digit_info.get_device_id(detector_name.Data());
  const auto& plane_names = digit_info.get_name_list(m_detector_id);
  m_plane_name = plane_names.at(plane_id);
  m_dcgeom_layer = plane_id; // gGeom.GetLayerId(m_detector_name+"-"+m_plane_name);
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
DCRawHit::~DCRawHit()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
DCRawHit::Print(const TString& arg) const
{
  hddaq::cerr << FUNC_NAME << " " << arg << std::endl
	      << "plane_id = " << m_plane_id    << std::endl
	      << "wire_id  = " << m_wire_id     << std::endl;
  std::vector<Int_t>::const_iterator itr, end;
  hddaq::cout << "tdc      = " << m_tdc.size() << " ";
  std::copy(m_tdc.begin(), m_tdc.end(),
	     std::ostream_iterator<Int_t>(hddaq::cout," "));
  hddaq::cout << std::endl
	      << "trailing = " << m_trailing.size() << " ";
  std::copy(m_trailing.begin(), m_trailing.end(),
	     std::ostream_iterator<Int_t>(hddaq::cout," "));
  hddaq::cout << std::endl;
}
