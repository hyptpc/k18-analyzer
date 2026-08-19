// -*- C++ -*-

#include "ConfMan.hh"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>

#include <TSystem.h>
#include <TFile.h>
#include <TNamed.h>
#include <TMacro.h>
#include <TObject.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TObjArray.h>
#include <TDirectory.h>

#include <lexical_cast.hh>
#include <filesystem_util.hh>
#include <replace_string.hh>
#include <std_ostream.hh>

// #include "BH1Filter.hh"
#include "BH1Match.hh"
#include "BH2Filter.hh"
#include "DCGeomMan.hh"
#include "DCTdcCalibMan.hh"
#include "DCDriftParamMan.hh"
#include "EventDisplay.hh"
#include "FieldMan.hh"
#include "FuncName.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "K18TransMatrix.hh"
#include "MatrixParamMan.hh"
#include "MsTParamMan.hh"
#include "UnpackerManager.hh"
#include "UserParamMan.hh"

namespace
{
using hddaq::unpacker::GUnpacker;
const TString kConfFile("CONF");
TString sConfDir;
auto& gUnpacker = GUnpacker::get_instance();
const auto& gMatrix = MatrixParamMan::GetInstance();
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
  fParamObjectsCache.clear(); 

  // close TFile and delete it 
  for(auto& pair : fParamTFiles){
    if(pair.second){
      std::cout << "--- [ConfMan::~ConfMan] Closing parameter file: " 
                << pair.second->GetName() << std::endl;
      if(pair.second->IsOpen()) pair.second->Close();
      delete pair.second;
      pair.second = nullptr;
    }
  }
  fParamTFiles.clear();
  fParamRootPaths.clear();
}

//_____________________________________________________________________________
void
ConfMan::AddObject()
{
  if(m_object) delete m_object;
  m_object = new TNamed("conf", m_buf.Data());
  m_object->Write();
}

//_____________________________________________________________________________
Bool_t
ConfMan::Initialize()
{
  if(m_is_ready){
    hddaq::cerr << FUNC_NAME << " already initialied" << std::endl;
    return false;
  }

  std::ifstream ifs(m_file[kConfFile]);
  if(!ifs.is_open()){
    hddaq::cerr << FUNC_NAME << " cannot open file : "
                << m_file[kConfFile] << std::endl;
    return false;
  }

  hddaq::cout << FUNC_NAME << " " << m_file[kConfFile] << std::endl;
  sConfDir = hddaq::dirname(m_file[kConfFile].Data());

  m_buf = "\n";

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

    if (key.BeginsWith("ParamFile.")) {
      TString fileKey = key(10, key.Length() - 10); 
      val = val.Strip(TString::kBoth);

      if (fileKey.IsNull()) {
        hddaq::cerr << FUNC_NAME << " Warning: Empty key for ParamFile (Key=" << key << ")" << std::endl;
        continue;
      }
      if (val.IsNull()) {
        hddaq::cerr << FUNC_NAME << " Warning: Empty path for ParamFile." << fileKey << std::endl;
        continue;
      }

      fParamRootPaths[fileKey] = FilePath(val);
      
      hddaq::cout << " " << "key = " << std::setw(10) << std::left << key
                  << " " << "value = " << std::setw(30) << std::left << val
                  << std::endl;
                  
      continue;
    }
    hddaq::cout << " key = "   << std::setw(10) << std::left << key
		<< " value = " << std::setw(30) << std::left << val
		<< std::endl;    

    m_file[key] = FilePath(val);
    m_string[key] = val;
    m_double[key] = val.Atof();
    m_int[key] = val.Atoi();
    m_bool[key] = (val.Atoi() == 1);
  }

  for(const auto& pair : fParamRootPaths){
    const TString& key = pair.first;
    const TString& path = pair.second;
    
    TFile* file = TFile::Open(path.Data(), "READ");
    if(!file || !file->IsOpen()){
      hddaq::cerr << FUNC_NAME << " !!! Failed to open parameter file (Key=" 
                  << key << ", Path=" << path << ")" << std::endl;
      delete file;
    } else {
      hddaq::cout << FUNC_NAME << " --- Loaded parameter file (Key=" 
                  << key << ", Path=" << path << ")" << std::endl;
      fParamTFiles[key] = file;
    }
  }
  
  AddObject();

  // For E42
  gUnpacker.enable_istream_bookmark();
  //

  if(!InitializeParameterFiles() || !InitializeHistograms()){
    return false;
  }

  if(gMatrix.IsReady()){
    gMatrix.Print2D();
  }
  if(gUser.IsReady()){
    gUser.Print();
  }

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


TObject* ConfMan::GetParamObject(const TString& fileKey, const TString& objName)
{
  auto it_file_cache = fParamObjectsCache.find(fileKey);
  if (it_file_cache != fParamObjectsCache.end()) {
    auto it_obj_cache = it_file_cache->second.find(objName);
    if (it_obj_cache != it_file_cache->second.end()) {
      return it_obj_cache->second;
    }
  }
  auto it_tfile = fParamTFiles.find(fileKey);
  if (it_tfile == fParamTFiles.end() || !it_tfile->second || !it_tfile->second->IsOpen()) {
    std::cerr << "!!! [ConfMan::GetParamObject] Error: Parameter file key '" << fileKey 
              << "' not found, not loaded, or not open." << std::endl;
    return nullptr;
  }

  TFile* file = it_tfile->second;
  TDirectory* g_dir_save = gDirectory;
  
  TObject* obj = file->Get(objName.Data());
  if (g_dir_save) {
      g_dir_save->cd(); // 保存したディレクトリ (出力ファイルのはず) に戻す
  }
  if (!obj) {
    std::cerr << "!!! [ConfMan::GetParamObject] Error: Object '" << objName 
              << "' not found in parameter file (Key=" << fileKey 
              << ", Path=" << file->GetName() << ")." << std::endl;
    return nullptr;
  }

  fParamObjectsCache[fileKey][objName] = obj;
  return obj;
}

TF1* ConfMan::GetParamTF1(const TString& fileKey, const TString& objName)
{
  TObject* obj = GetParamObject(fileKey, objName);
  TF1* func = dynamic_cast<TF1*>(obj);
  if (!func && obj) {
    std::cerr << "!!! [ConfMan::GetParamTF1] Error: Object '" << objName 
              << "' (Key=" << fileKey << ") was found but is not a TF1 (it is a " 
              << obj->ClassName() << ")." << std::endl;
  }
  return func;
}

TH1* ConfMan::GetParamTH1(const TString& fileKey, const TString& objName)
{
  TObject* obj = GetParamObject(fileKey, objName);
  TH1* hist = dynamic_cast<TH1*>(obj);
  if (!hist && obj) {
    std::cerr << "!!! [ConfMan::GetParamTH1] Error: Object '" << objName 
              << "' (Key=" << fileKey << ") was found but is not a TH1 (it is a " 
              << obj->ClassName() << ")." << std::endl;
  }
  return hist;
}

TH2* ConfMan::GetParamTH2(const TString& fileKey, const TString& objName)
{
  TObject* obj = GetParamObject(fileKey, objName);
  TH2* hist = dynamic_cast<TH2*>(obj);
  if (!hist && obj) {
    std::cerr << "!!! [ConfMan::GetParamTH2] Error: Object '" << objName 
              << "' (Key=" << fileKey << ") was found but is not a TH2 (it is a " 
              << obj->ClassName() << ")." << std::endl;
  }
  return hist;
}
