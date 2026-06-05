// -*- C++ -*-

#ifndef MC_PARAM_MAN_HH
#define MC_PARAM_MAN_HH

#include <map>
#include <vector>
#include <TMath.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TRandom.h>
//_____________________________________________________________________________
class MCParamMan{
private:
    std::map<TString, TH1*> hMap;
    std::map<TString, TH2*> hMap2D;

public:
    static const TString& ClassName();
    static MCParamMan& GetInstance();
    ~MCParamMan();

private:
    MCParamMan(){};
    //MCParamMan(const MCParamMan&);
    //MCParamMan& operator =(const MCParamMan&);

public:
    void SetHisto(const TString& name, TH1* histo){
        hMap[name] = histo;
    };
    void SetHisto2D(const TString& name, TH2* histo){
        hMap2D[name] = histo;
    };
    TH1* GetHisto(const TString& name) const{
        auto it = hMap.find(name);
        return (it != hMap.end()) ? it->second : nullptr;
    };
    TH2* GetHisto2D(const TString& name) const{
        auto it = hMap2D.find(name);
        return (it != hMap2D.end()) ? it->second : nullptr;
    };

    //From here, User-specific functions are defined.
private:
    double dth = 2;
    int nbin_th = 12; // ~24 deg
    vector<double> resZ = {
        85,42,26,19,14,11,10,9,8,8,7,6
    };
    int nbin_dist = 300;
    double d0 = 0, d1 = 30; 
    struct HistChi2{
        int id; 
        double d_min,d_max;
        int nbin_chi2;
        double chi2_min,chi2_max;
    };
    std::vector<HistChi2> HistChi2Lists = { 
        {0,0.,2.,1000,0,3},
        {1,2.,4.,100,0,7},
        {2,4.,6.,100,1,15},
        {3,6.,8.,100,4,20},
        {4,8.,10.,100,9,30},
        {5,10.,14.,100,14,40},
        {6,14.,18.,100,25,60},
        {7,18.,25.,100,45,120},
        {8,25.,30.,100,85,140}
    };
    TFile* kkfile = nullptr;

public:
    TString DistHistName(double th);
    TString XDistHistName(double th);
    TString YDistHistName(double th);
    TString XYDistHistName(double th);
    TString Chi2HistName(double dist);
    void LoadKKResHistograms(const TString& file_name);
    double GetResZ(double th);
    void GetDistChi2(double th, double& dist, double& chi2);
    void GetRandomParams(double th, double& z, double& dist, double& chi2);
    void GetRandomdXY(double th, double& z, double& dist_x, double& dist_y);
};
inline const TString&
MCParamMan::ClassName(){
    static TString g_name("MCParamMan");
    return g_name;
}
inline MCParamMan&
MCParamMan::GetInstance(){
    static MCParamMan s_instance;
    //We need a really precise timing to set random seed for a multiple parallel jub, submitted at nearly the same time.
    auto time_seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    gRandom->SetSeed(time_seed);
    return s_instance;
}
inline MCParamMan::~MCParamMan(){
    for(auto& pair : hMap){
        //if(pair.second) delete pair.second;
    }
    for(auto& pair : hMap2D){
        //if(pair.second) delete pair.second;
    }
    //delete kkfile;
    //Deleting these by hand is not a good idea. It causes segfaults when the program ends.
}
//Following histograms are expected to be loaded from the file.
inline TString
MCParamMan::DistHistName(double th){
    int ith = (int)(th / dth);
    double th_min = ith * dth;
    double th_max = (ith+1) * dth;
    return Form("h_dist_th_%g_%g", th_min, th_max);
}
inline TString
MCParamMan::XDistHistName(double th){
    int ith = (int)(th / dth);
    double th_min = ith * dth;
    double th_max = (ith+1) * dth;
    return Form("h_Xdist_th_%g_%g", th_min, th_max);
}
inline TString
MCParamMan::YDistHistName(double th){
    int ith = (int)(th / dth);
    double th_min = ith * dth;
    double th_max = (ith+1) * dth;
    return Form("h_Ydist_th_%g_%g", th_min, th_max);
}
inline TString
MCParamMan::XYDistHistName(double th){
    int ith = (int)(th / dth);
    double th_min = ith * dth;
    double th_max = (ith+1) * dth;
    return Form("h_XYdist_th_%g_%g", th_min, th_max);
}
inline TString
MCParamMan::Chi2HistName(double dist){
    for(auto& h: HistChi2Lists){
        if(dist >= h.d_min && dist < h.d_max){
        return Form("h_chi2_d_%g_%g", h.d_min, h.d_max);
        }
    }
    return "";
}
inline void
MCParamMan::LoadKKResHistograms(const TString& file_name){
    kkfile = new TFile(file_name);
    if(!kkfile || kkfile->IsZombie()){
        std::cerr << "Error: Cannot open file " << file_name << std::endl;
        return;
    }
    for(int ith = 0; ith < nbin_th; ++ith){
        double th_min = ith * dth;
        double th_max = (ith+1) * dth;
        TH1D* hist = (TH1D*)kkfile->Get(Form("h_dist_th_%g_%g", th_min, th_max));
        if(hist){
            SetHisto(Form("h_dist_th_%g_%g", th_min, th_max), hist);
        } else {
            std::cerr << "Warning: Cannot find histogram " << Form("h_dist_th_%g_%g", th_min, th_max) << " in file " << file_name << std::endl;
        }
        hist = (TH1D*)kkfile->Get(XDistHistName((th_min+th_max)/2));
        if(hist){
            SetHisto(XDistHistName((th_min+th_max)/2), hist);
        }
        else {
            std::cerr << "Warning: Cannot find histogram " << XDistHistName((th_min+th_max)/2) << " in file " << file_name << std::endl;
        }
        hist = (TH1D*)kkfile->Get(YDistHistName((th_min+th_max)/2));
        if(hist){
            SetHisto(YDistHistName((th_min+th_max)/2), hist);
        }
        else {
            std::cerr << "Warning: Cannot find histogram " << YDistHistName((th_min+th_max)/2) << " in file " << file_name << std::endl;
        }
        TH2D* hist2d = (TH2D*)kkfile->Get(XYDistHistName((th_min+th_max)/2));
        if(hist2d){
            SetHisto2D(XYDistHistName((th_min+th_max)/2), hist2d);
        }
        else {
            std::cerr << "Warning: Cannot find histogram " << XYDistHistName((th_min+th_max)/2) << " in file " << file_name << std::endl; 
        }
    }
    for(auto& h: HistChi2Lists){
        TH1D* hist = (TH1D*)kkfile->Get(Form("h_chi2_d_%g_%g", h.d_min, h.d_max));
        if(hist){
            SetHisto(Form("h_chi2_d_%g_%g", h.d_min, h.d_max), hist);
        } else {
            std::cerr << "Warning: Cannot find histogram " << Form("h_chi2_d_%g_%g", h.d_min, h.d_max) << " in file " << file_name << std::endl;
        }
    }
}
inline double
MCParamMan::GetResZ(double th){
    int ith = (int)(th / dth);
    if(ith >= 0 and ith < (resZ.size())){
        return resZ[ith];
    }
    return resZ[resZ.size()-1];
}
inline void
MCParamMan::GetDistChi2(double th, double& dist, double& chi2){
    TH1* h_dist = GetHisto(DistHistName(th));
    if(h_dist){
        dist = h_dist->GetRandom();
    } else {
        std::cerr << "Error: Cannot find histogram for dist with th = " << th << std::endl;
        dist = 0;
    }
    TH1* h_chi2 = GetHisto(Chi2HistName(dist));
    if(h_chi2){
        chi2 = h_chi2->GetRandom();
    } else {
        std::cerr << "Error: Cannot find histogram for chi2 with dist = " << dist << std::endl;
        chi2 = 0;
    }
}
inline
void 
MCParamMan::GetRandomParams(double th, double& z, double& dist, double& chi2){
    double res_z = GetResZ(th);
    z = gRandom->Gaus(0, res_z);
    GetDistChi2(th, dist, chi2);
}
inline
void 
MCParamMan::GetRandomdXY(double th, double& z, double& dist_x, double& dist_y){
    double res_z = GetResZ(th);
    z = gRandom->Gaus(0, res_z);
    TH2* h_xy_dist = GetHisto2D(XYDistHistName(th));
    h_xy_dist->GetRandom2(dist_x, dist_y);
}
#endif 