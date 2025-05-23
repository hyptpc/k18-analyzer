// -*- C++ -*-

#include "ConfMan.hh"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>

#include <TFile.h>
#include <TMacro.h>
#include <TNamed.h>
#include <TSystem.h>

#include <lexical_cast.hh>
#include <filesystem_util.hh>
#include <replace_string.hh>

#include <spdlog/spdlog.h>

#include "FuncName.hh"
#include "UnpackerManager.hh"
#include "UserParamMan.hh"

namespace
{
const TString kConfFile("CONF");
TString sConfDir;
auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
// const auto& gMatrix = MatrixParamMan::GetInstance();
auto& gUser = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
ConfMan::ConfMan()
  : m_is_ready(false),
    m_file(),
    m_string(),
    m_double(),
    m_int(),
    m_bool(),
    m_buf(),
    m_object()
{
}

//_____________________________________________________________________________
ConfMan::~ConfMan()
{
}

//_____________________________________________________________________________
Bool_t
ConfMan::Initialize()
{
  if(m_is_ready){
    spdlog::error("{} already initialied", FUNC_NAME.Data());
    return false;
  }

  std::ifstream ifs(m_file[kConfFile]);
  if(!ifs.is_open()){
    spdlog::error("{} cannot open file : {}", FUNC_NAME.Data(),
                  m_file[kConfFile].Data());
    return false;
  }

  sConfDir = hddaq::dirname(m_file[kConfFile].Data());
  m_buf = "\n";

  std::ostringstream oss;
  oss << FUNC_NAME << " " << m_file[kConfFile] << std::endl;

  TString line;
  while(ifs.good() && line.ReadLine(ifs)){
    m_buf += line + "\n";
    if(line.IsNull() || line[0]=='#') continue;

    line.ReplaceAll(",",  ""); // remove ,
    line.ReplaceAll(":",  ""); // remove :
    line.ReplaceAll("\"",  ""); // remove "

    std::istringstream iss(line.Data());
    std::istream_iterator<std::string> begin(iss);
    std::istream_iterator<std::string> end;
    std::vector<TString> v(begin, end);
    if(v.size()<2) continue;

    TString key = v[0];
    TString val = v[1];
    oss << " key = "   << std::setw(10) << std::left << key
        << " value = " << std::setw(30) << std::left << val << std::endl;

    m_file[key] = FilePath(val);
    m_string[key] = val;
    m_double[key] = val.Atof();
    m_int[key] = val.Atoi();
    m_bool[key] = (val.Atoi() == 1);
  }
  spdlog::info(oss.str());

  if(!InitializeParameterFiles())
    return false;
  // if(gMatrix.IsReady())
  //   gMatrix.Print2D();
  // if(gUser.IsReady())
  //   gUser.Print();

  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::Initialize(const TString& file_name)
{
  m_file[kConfFile] = file_name;
  return Initialize();
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeUnpacker()
{
  gUnpacker.set_config_file(m_file["UNPACK"].Data(),
                            m_file["DIGIT"].Data(),
                            m_file["CMAP"].Data());
  if(gUnpacker.get_skip() == 0){
    TNamed git("git", ("\n"+gSystem->GetFromPipe("git log -1")).Data());
    git.Write();
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::Finalize()
{
  return FinalizeProcess();
}

//_____________________________________________________________________________
TString
ConfMan::FilePath(const TString& src) const
{
  std::ifstream tmp(src);
  if(tmp.good()) return src;
  else           return sConfDir + "/" + src;
}

//_____________________________________________________________________________
void
ConfMan::WriteParameters()
{
  if(gUnpacker.get_skip() == 0){
    gFile->mkdir("param");
    gFile->cd("param");
    for(const auto& itr: m_file){
      TMacro paramfile;
      paramfile.SetName(itr.first);
      paramfile.SetTitle(itr.second);
      paramfile.ReadFile(itr.second);
      paramfile.Write();
    }
    gFile->cd();
  }
}
