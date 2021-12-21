// -*- C++ -*-

#ifndef TPC_PARAM_MAN_HH
#define TPC_PARAM_MAN_HH

#include <map>
#include <TMath.h>
#include <TString.h>

//_____________________________________________________________________________
class TPCAParam
{
public:
  TPCAParam(Double_t offset, Double_t gain)
    : m_offset(offset),
      m_gain(gain)
    {}
  ~TPCAParam()
    {}

private:
  TPCAParam();
  TPCAParam(const TPCAParam&);
  TPCAParam& operator =(const TPCAParam&);

private:
  Double_t m_offset;
  Double_t m_gain;

public:
  Double_t Offset() const { return m_offset; }
  Double_t Gain() const { return m_gain; }
  Double_t CDeltaE(Double_t de) const { return (de - m_offset) * m_gain; }
  Double_t Offset(Double_t cde) const { return cde / m_gain + m_offset; }
  Double_t P0() const { return m_offset; }
  Double_t P1() const { return m_gain; }
};

//_____________________________________________________________________________
class TPCTParam
{
public:
  TPCTParam(Double_t offset, Double_t gain)
    : m_offset(offset),
      m_gain(gain)
    {}
  ~TPCTParam()
    {}

private:
  TPCTParam();
  TPCTParam(const TPCTParam&);
  TPCTParam& operator =(const TPCTParam&);

private:
  Double_t m_offset;
  Double_t m_gain;

public:
  Double_t Offset() const { return m_offset; }
  Double_t Gain() const { return m_gain; }
  Double_t P0() const { return m_offset; }
  Double_t P1() const { return m_gain; }
  Double_t Time(Double_t tdc) const { return (tdc - m_offset) * m_gain; }
  Double_t Tdc(Double_t time) const { return time / m_gain + m_offset; }
};

//_____________________________________________________________________________
class TPCYParam
{
public:
  TPCYParam(Double_t offset, Double_t drift_velocity)
    : m_offset(offset),
      m_drift_velocity(drift_velocity)
    {}
  ~TPCYParam()
    {}

private:
  TPCYParam();
  TPCYParam(const TPCYParam&);
  TPCYParam& operator =(const TPCYParam&);

private:
  Double_t m_offset;
  Double_t m_drift_velocity;

public:
  Double_t Offset() const { return m_offset; }
  Double_t DriftLength(Double_t time) const { return Y(time); }
  Double_t DriftVelocity() const { return m_drift_velocity; }
  Double_t P0() const { return m_offset; }
  Double_t P1() const { return m_drift_velocity; }
  Double_t Time(Double_t dl) const
    { return dl / m_drift_velocity + m_offset; }
  Double_t Y(Double_t time) const
    { return (time - m_offset) * m_drift_velocity; }
};

//_____________________________________________________________________________
class TPCCoboParam
{
public:
  TPCCoboParam(const std::vector<Double_t> params)
    : m_params(params)
    {}
  // TPCCoboParam(Double_t phase_shift)
  //   : m_phase_shift(phase_shift)
  //   {}
  ~TPCCoboParam()
    {}

private:
  TPCCoboParam();
  TPCCoboParam(const TPCCoboParam&);
  TPCCoboParam& operator =(const TPCCoboParam&);

private:
  std::vector<Double_t> m_params;
  // Double_t              m_phase_shift;

public:
  Double_t PhaseShift(Double_t clk) const
  { return clk + m_params[0]*TMath::Freq((clk-m_params[1])/m_params[2]); }
  // Double_t PhaseShift() const { return m_phase_shift; }
  const std::vector<Double_t>& Params() const { return m_params; }
};



//_____________________________________________________________________________
class TPCParamMan
{
public:
  static const TString& ClassName();
  static TPCParamMan&   GetInstance();
  ~TPCParamMan();

private:
  TPCParamMan();
  TPCParamMan(const TPCParamMan&);
  TPCParamMan& operator =(const TPCParamMan&);

private:
  //  enum eAorT { kAdc, kTdc, kY, kATY };
  enum eAorT { kAdc, kTdc, kY, kCobo };
  typedef std::map<Int_t, TPCAParam*> AContainer;
  typedef std::map<Int_t, TPCTParam*> TContainer;
  typedef std::map<Int_t, TPCYParam*> YContainer;
  typedef std::map<Int_t, TPCCoboParam*> CoboContainer;
  typedef AContainer::const_iterator AIterator;
  typedef TContainer::const_iterator TIterator;
  typedef YContainer::const_iterator YIterator;
  typedef CoboContainer::const_iterator CoboIterator;
  Bool_t        m_is_ready;
  TString       m_file_name;
  AContainer    m_APContainer;
  TContainer    m_TPContainer;
  YContainer    m_YPContainer;
  CoboContainer m_CoboContainer;

public:
  Bool_t GetCDe(Int_t layer, Int_t row, Double_t de, Double_t &cde) const;
  Bool_t GetCTime(Int_t layer, Int_t row, Double_t time,
                  Double_t &ctime) const;
  Bool_t GetCClock(Int_t layer, Int_t row, Double_t time,
                   Double_t& cclk) const;
  Bool_t GetDriftLength(Int_t layer, Int_t row, Double_t time,
                        Double_t& y) const
    { return GetY(layer, row, time, y); }
  Bool_t GetY(Int_t layer, Int_t row, Double_t time, Double_t& y) const;
  Bool_t Initialize();
  Bool_t Initialize(const TString& file_name);
  Bool_t IsReady() const { return m_is_ready; }
  void   SetFileName(const TString& file_name) { m_file_name = file_name; }

private:
  void          ClearACont();
  void          ClearTCont();
  void          ClearYCont();
  void          ClearCoboCont();
  TPCAParam*    GetAmap(Int_t layer, Int_t row) const;
  TPCTParam*    GetTmap(Int_t layer, Int_t row) const;
  TPCYParam*    GetYmap(Int_t layer, Int_t row) const;
  TPCCoboParam* GetCobomap(Int_t Cobo) const;
};

//_____________________________________________________________________________
inline const TString&
TPCParamMan::ClassName()
{
  static TString g_name("TPCParamMan");
  return g_name;
}

//_____________________________________________________________________________
inline TPCParamMan&
TPCParamMan::GetInstance()
{
  static TPCParamMan s_instance;
  return s_instance;
}

#endif
