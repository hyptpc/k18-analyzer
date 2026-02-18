// -*- C++ -*-

#include "CounterMapMan.hh"

#include <cstdio>
#include <iomanip>
#include <iostream>
#include <new>
#include <string>

#ifndef FUNC_NAME
#define FUNC_NAME "[" << className << "::" << __func__ << "()]"
#endif

namespace
{
const TString className = "CounterMapMan";
const UInt_t KEYMASK  = 0x0007;
const UInt_t AMASK    = 0x003F;      /* A Mask 6 Bits (0-63) */
const UInt_t NMASK    = 0x001F;      /* N Mask 5 Bits (0-31) */
const UInt_t CMASK    = 0x000F;      /* C Mask 4 Bits (0-15) */
const Int_t  ASHIFT   =  4;
const Int_t  NSHIFT   = 12;
const Int_t  CSHIFT   = 20;
const UInt_t KEYFLAG  = 0x0003;

const UInt_t SEGMASK  = 0x00FF;    /* Segment  (0-255) */
const UInt_t CIDMASK  = 0x007F;    /* CId      (0-127) */
const UInt_t UDMASK   = 0x0007;    /* UD       (0-7)   */
const UInt_t ATMASK   = 0x0003;    /* AT       (0-1)   */
const UInt_t LAYMASK  = 0x007F;    /* Layer    (0-127)  */
const Int_t  SEGSHIFT =  3;
const Int_t  CIDSHIFT = 12;
const Int_t  UDSHIFT  = 20;
const Int_t  ATSHIFT  = 23;
const Int_t  LAYSHIFT = 25;
const UInt_t RKEYFLAG = 0x0004;

const Int_t MAXCHAR = 144;

UInt_t KEY(Int_t c, Int_t n, Int_t a)
{
  return ((((c) & CMASK) << CSHIFT) |
          (((n) & NMASK) << NSHIFT) |
          (((a) & AMASK) << ASHIFT) | KEYFLAG);
}

UInt_t RKEY(Int_t at, Int_t ud, Int_t cid, Int_t seg, Int_t lay)
{
  return ((((at) & ATMASK) << ATSHIFT) |
          (((ud) & UDMASK) << UDSHIFT) |
          (((cid) & CIDMASK) << CIDSHIFT) |
          (((seg) & SEGMASK) << SEGSHIFT) |
          (((lay) & LAYMASK) << LAYSHIFT) | RKEYFLAG);
}
}

//_____________________________________________________________________________
CounterMapMan::CounterMapMan()
  : NSMP(0),
    NCH_SCA(0),
    m_isready(false)
{
}

//_____________________________________________________________________________
CounterMapMan::~CounterMapMan()
{
}

//_____________________________________________________________________________
void
CounterMapMan::Clear()
{
  fContainer.clear();
  bContainer.clear();
  fCrateDef.clear();
  bCrateDef.clear();
  nameCNAContainer.clear();
  nameCounterContainer.clear();
  NSMP = NCH_SCA = 0;
  for(Int_t i=0; i<MMAXSMP; ++i) {
    for(Int_t j=0; j<MAXSLOT; ++j) {
      CRATE_TYPE[i][j] = DRT;
    }
  }
}

//_____________________________________________________________________________
Bool_t
CounterMapMan::Initialize()
{
  Clear();
  for(Int_t cr=0; cr<15; ++cr) {
    for(Int_t sl=1; sl<24; ++sl) {
      for(Int_t ch=0; ch<32; ++ch) {
        fContainer[KEY(cr, sl, ch)] = RKEY(0, 0, 127, 0, 0);
      }
    }
  }

  for(const auto& file : FileName) {
    if(!ReadFile(file)) return false;
  }

  m_isready = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
CounterMapMan::Initialize(const char* file_name)
{
  FileName.clear();
  FileName.push_back(file_name);
  return Initialize();
}

//_____________________________________________________________________________
Bool_t
CounterMapMan::Initialize(const std::string& file_name)
{
  FileName.clear();
  FileName.push_back(file_name);
  return Initialize();
}

//_____________________________________________________________________________
Bool_t
CounterMapMan::Initialize(const std::string& file_name1,
                          const std::string& file_name2)
{
  FileName.clear();
  FileName.push_back(file_name1);
  FileName.push_back(file_name2);
  return Initialize();
}

//_____________________________________________________________________________
Bool_t
CounterMapMan::ReadFile(const TString& filename)
{
  FILE *fp = fopen(filename.Data(), "r");
  if(!fp) {
    std::cerr << FUNC_NAME << " file open fail: " << filename << std::endl;
    return false;
  }

  Char_t str[MAXCHAR];
  while(fgets(str, MAXCHAR, fp)) {
    if(str[0] == '#') continue;

    Int_t c, n, a, cid, lay, wire, at, ud, seg, address;
    Int_t drtmin, drtmax, normmin, normmax;
    Char_t name[MAXCHAR];

    if(sscanf(str, "CrateDef: %d %x", &c, &address) == 2) {
      fCrateDef[c] = address;
      bCrateDef[address] = c;
    } else if(sscanf(str, "CrateType: %d MIX DRT %d %d NORMAL %d %d",
                     &n, &drtmin, &drtmax, &normmin, &normmax) == 5) {
      for(Int_t i=drtmin; i<=drtmax; ++i) CRATE_TYPE[n][i-1] = DRT;
      for(Int_t i=normmin; i<=normmax; ++i) CRATE_TYPE[n][i-1] = NORMAL;
    } else if(sscanf(str, "CrateType: %d %s", &n, name) == 2) {
      Int_t type = (TString(name) == "DRT") ? DRT : NORMAL;
      for(Int_t i=0; i<MAXSLOT; ++i) CRATE_TYPE[n][i] = type;
    } else if(sscanf(str, "NumSMP: %d", &n) == 1) {
      NSMP = n;
    } else if(sscanf(str, "NumScaler: %d", &n) == 1) {
      NCH_SCA = n;
    } else if(sscanf(str, "%d %d %d %d %d %d %d %d %s",
                     &c, &n, &a, &cid, &lay, &wire, &at, &ud, name) == 9) {
      UInt_t key = KEY(c, n, a);
      UInt_t rkey = RKEY(at, ud, cid, wire, lay);
      fContainer[key] = rkey;
      bContainer[rkey] = key;
      nameCNAContainer[key] = name;
      nameCounterContainer[rkey] = name;
    } else if(sscanf(str, "%d %d %d %d %d %d %d %s",
                     &c, &n, &a, &cid, &seg, &at, &ud, name) == 8) {
      UInt_t key = KEY(c, n, a);
      UInt_t rkey = RKEY(at, ud, cid, seg, 0);
      fContainer[key] = rkey;
      bContainer[rkey] = key;
      nameCNAContainer[key] = name;
      nameCounterContainer[rkey] = name;
    } else if(sscanf(str, "%d %d %d %d %d %d %d",
                     &c, &n, &a, &cid, &seg, &at, &ud) == 7) {
      UInt_t key = KEY(c, n, a);
      UInt_t rkey = RKEY(at, ud, cid, seg, 0);
      fContainer[key] = rkey;
      bContainer[rkey] = key;
    }
  }
  fclose(fp);
  return true;
}

//_____________________________________________________________________________
TString
CounterMapMan::GetName(Int_t c, Int_t n, Int_t a)
{
  auto it = nameCNAContainer.find(KEY(c, n, a));
  return (it != nameCNAContainer.end()) ? it->second : "";
}

//_____________________________________________________________________________
TString
CounterMapMan::GetName(Int_t cid, Int_t seg, Int_t at, Int_t ud)
{
  auto it = nameCounterContainer.find(RKEY(at, ud, cid, seg, 0));
  return (it != nameCounterContainer.end()) ? it->second : "";
}

//_____________________________________________________________________________
TString
CounterMapMan::GetName(Int_t cid, Int_t lay, Int_t wire, Int_t at, Int_t ud)
{
  auto it = nameCounterContainer.find(RKEY(at, ud, cid, wire, lay));
  return (it != nameCounterContainer.end()) ? it->second : "";
}

//_____________________________________________________________________________
Int_t
CounterMapMan::GetCrateNum(Int_t address)
{
  auto it = bCrateDef.find(address);
  return (it != bCrateDef.end()) ? it->second : -1;
}

//_____________________________________________________________________________
Int_t
CounterMapMan::GetSMPAddress(Int_t c)
{
  auto it = fCrateDef.find(c);
  return (it != fCrateDef.end()) ? it->second : -1;
}

//_____________________________________________________________________________
Bool_t
CounterMapMan::GetCNA(Int_t cid, Int_t seg, Int_t at, Int_t ud,
                      Int_t& c, Int_t& n, Int_t& a)
{
  UInt_t rkey = RKEY(at, ud, cid, seg, 0);
  auto it = bContainer.find(rkey);
  if(it != bContainer.end() && (it->second & KEYMASK) == KEYFLAG) {
    c = (it->second >> CSHIFT) & CMASK;
    n = (it->second >> NSHIFT) & NMASK;
    a = (it->second >> ASHIFT) & AMASK;
    return true;
  }
  return false;
}

//_____________________________________________________________________________
Bool_t
CounterMapMan::GetInfo(Int_t c, Int_t n, Int_t a,
                       Int_t& cid, Int_t& lay, Int_t& seg, Int_t& at, Int_t& ud)
{
  auto it = fContainer.find(KEY(c, n, a));
  if(it != fContainer.end() && (it->second & KEYMASK) == RKEYFLAG) {
    UInt_t rkey = it->second;
    cid = (rkey >> CIDSHIFT) & CIDMASK;
    seg = (rkey >> SEGSHIFT) & SEGMASK;
    lay = (rkey >> LAYSHIFT) & LAYMASK;
    at  = (rkey >> ATSHIFT ) & ATMASK;
    ud  = (rkey >> UDSHIFT ) & UDMASK;
    return true;
  }
  return false;
}

//_____________________________________________________________________________
Int_t
CounterMapMan::GetCID(Int_t c, Int_t n, Int_t a)
{
  auto it = fContainer.find(KEY(c, n, a));
  if(it != fContainer.end() && (it->second & KEYMASK) == RKEYFLAG) {
    return (it->second >> CIDSHIFT) & CIDMASK;
  }
  return -1;
}

//_____________________________________________________________________________
Bool_t
CounterMapMan::GetCNA(Int_t cid, Int_t lay, Int_t wire, Int_t at, Int_t ud,
                      Int_t& c, Int_t& n, Int_t& a)
{
  UInt_t rkey = RKEY(at, ud, cid, wire, lay);
  auto it = bContainer.find(rkey);
  if(it != bContainer.end() && (it->second & KEYMASK) == KEYFLAG) {
    c = (it->second >> CSHIFT) & CMASK;
    n = (it->second >> NSHIFT) & NMASK;
    a = (it->second >> ASHIFT) & AMASK;
    return true;
  }
  c = n = a = -1;
  return false;
}

//_____________________________________________________________________________
void
CounterMapMan::PrintSimpleMap(std::ostream& p_out)
{
  p_out << FUNC_NAME << std::endl;
  for(const auto& it : fContainer) {
    UInt_t key = it.first;
    UInt_t rkey = it.second;
    Int_t cr = (key >> CSHIFT) & CMASK;
    if(cr == 0 || cr == 1 || cr == 2 || cr == 6) {
      p_out << std::setw(5) << cr
            << std::setw(5) << ((key >> NSHIFT) & NMASK)
            << std::setw(5) << ((key >> ASHIFT) & AMASK)
            << std::setw(5) << ((rkey >> CIDSHIFT) & CIDMASK)
            << std::setw(5) << ((rkey >> LAYSHIFT) & LAYMASK)
            << std::setw(5) << ((rkey >> SEGSHIFT) & SEGMASK)
            << std::setw(5) << ((rkey >> ATSHIFT ) & ATMASK)
            << std::setw(5) << ((rkey >> UDSHIFT ) & UDMASK)
            << std::endl;
    }
  }
}

//_____________________________________________________________________________
void
CounterMapMan::PrintMap()
{
  std::cout << FUNC_NAME << std::endl;
  for(const auto& it : fContainer) {
    UInt_t key = it.first;
    UInt_t rkey = it.second;
    std::cout << std::setw(10) << key << " ("
              << std::setw(2) << ((key >> CSHIFT) & CMASK) << ","
              << std::setw(2) << ((key >> NSHIFT) & NMASK) << ","
              << std::setw(2) << ((key >> ASHIFT) & AMASK) << ") -> "
              << std::setw(10) << rkey << " ("
              << std::setw(3) << ((rkey >> CIDSHIFT) & CIDMASK) << ","
              << std::setw(3) << ((rkey >> LAYSHIFT) & LAYMASK) << ","
              << std::setw(3) << ((rkey >> SEGSHIFT) & SEGMASK) << ","
              << std::setw(1) << ((rkey >> ATSHIFT ) & ATMASK) << ","
              << std::setw(1) << ((rkey >> UDSHIFT ) & UDMASK) << ")"
              << std::endl;
  }
}
