#ifndef TPCPID_HH
#define TPCPID_HH

#define NUMTPCPID 1

#include "TObject.h"
#include <map>
#include "TString.h"

class TPCPID : public TObject {

public:
  enum PID {
    kNotDetermined = -1,
    kPion = 0,
    kKaon = 1,
    kProton = 2,
    kDeuteron = 3,
    kTriton = 4
  };
  
  friend std::ostream& operator<<(std::ostream &out, const PID value) {
    static std::map<PID, std::string> particleNames;
    if (particleNames.size() == 0) {
      particleNames[kNotDetermined] = "NotDetermined";
      particleNames[kPion]          = "Pion";
      particleNames[kKaon]          = "Kaon";
      particleNames[kProton]        = "Proton";
      particleNames[kDeuteron]      = "Deuteron";
      particleNames[kTriton]        = "Triton";
    }
    return out << particleNames[value];
  }

  static TString GetName(PID pid) {
    static std::map<PID, std::string> particleNames;
    if (particleNames.size() == 0) {
      particleNames[kNotDetermined] = "NotDetermined";
      particleNames[kPion]          = "Pion";
      particleNames[kKaon]          = "Kaon";
      particleNames[kProton]        = "Proton";
      particleNames[kDeuteron]      = "Deuteron";
      particleNames[kTriton]        = "Triton";
    }
    return TString(particleNames[pid].c_str());
  }

  static Int_t GetPDG(PID pid) {
    static std::map<PID, Int_t> particlePDG;
    if (particlePDG.size() == 0) {
      particlePDG[kNotDetermined] = 0;
      particlePDG[kPion]          = 211;
      particlePDG[kKaon]          = 321;
      particlePDG[kProton]        = 2212;
      particlePDG[kDeuteron]      = 1000010020;
      particlePDG[kTriton]        = 1000010030;
    }
    return particlePDG[pid];
  }

  static Int_t GetNumTPCPID() { return NUMTPCPID; }

private:

  //  const Int_t NUMTPCPID = 1;
  
  ClassDef(TPCPID, 1)
};

#endif

                                                                                           
                                                          
    
    
  
