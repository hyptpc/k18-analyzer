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
#include <TMD5.h>
#include <TNamed.h>
#include <TSystem.h>

#include "FuncName.hh"
#include "UserParamMan.hh"

#include <filesystem_util.hh>
#include <lexical_cast.hh>
#include <replace_string.hh>
#include <spdlog/spdlog.h>
#include <UnpackerManager.hh>

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
    spdlog::error("{} already initialized", FUNC_NAME.Data());
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

  spdlog::info("{} Read {}", FUNC_NAME.Data(), m_file[kConfFile].Data());

  TString line;
  while(ifs.good() && line.ReadLine(ifs)){
    m_buf += line + "\n";
    if(line.IsNull() || line[0]=='#') continue;
    spdlog::debug(" [Load] : \"{}\"", line.Data());

    line.ReplaceAll(",",  ""); 
    line.ReplaceAll(":",  ""); 
    line.ReplaceAll("\"",  ""); 

    std::istringstream iss(line.Data());
    std::istream_iterator<std::string> begin(iss);
    std::istream_iterator<std::string> end;
    std::vector<TString> v(begin, end);
    if(v.size()<2) continue;

    TString key = v[0];
    TString val = v[1];
    spdlog::info(" key = {:<15} value = {:<30}", key.Data(), val.Data());

    m_file[key]   = FilePath(val);
    m_string[key] = val;
    m_double[key] = val.Atof();
    m_int[key]    = val.Atoi();
    m_bool[key]   = (val.Atoi() == 1);
  }

  if(!InitializeParameterFiles())
    return false;

  //if(!InitializeHistograms())
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
// Content fingerprint only; does not restore the file. Keep large maps by path.
static TString
FileMd5Hex(const TString& path)
{
  std::ifstream in(path.Data(), std::ios::binary);
  if (!in) return "unavailable";

  TMD5 md5;
  constexpr std::size_t kBufSize = 1 << 20; // 1 MiB
  std::vector<char> buf(kBufSize);
  while (in) {
    in.read(buf.data(), static_cast<std::streamsize>(buf.size()));
    const auto n = in.gcount();
    if (n > 0)
      md5.Update(reinterpret_cast<const UChar_t*>(buf.data()),
                 static_cast<UInt_t>(n));
  }
  if (in.bad()) return "unavailable";
  md5.Final();
  return md5.AsString();
}

//_____________________________________________________________________________
void
ConfMan::WriteParameters()
{
  if(gUnpacker.get_skip() == 0){
    gFile->mkdir("param");
    gFile->cd("param");
    for(const auto& itr: m_file){
      // Numeric configuration entries (e.g. HSFLDCALC) are not files.
      if(gSystem->AccessPathName(itr.second))
        continue;
      TMacro paramfile;
      paramfile.SetName(itr.first);
      paramfile.SetTitle(itr.second);
      // Large files (e.g. the field map) cost GBs as a TMacro, one TObjString
      // per line: store path, size, and md5 fingerprint (not the contents).
      static const Long64_t MaxParamFileSize = 10*1024*1024; // [byte]
      FileStat_t stat;
      if(gSystem->GetPathInfo(itr.second, stat) == 0 &&
         stat.fSize > MaxParamFileSize){
        const TString md5 = FileMd5Hex(itr.second);
        paramfile.AddLine(Form("# not stored (%lld bytes, md5=%s): %s",
                               stat.fSize, md5.Data(), itr.second.Data()));
      }else{
        paramfile.ReadFile(itr.second);
      }
      paramfile.Write();
    }
    gFile->cd();
  }
}
