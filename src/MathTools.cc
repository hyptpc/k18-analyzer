// -*- C++ -*-

#include "MathTools.hh"

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <string>

#include <spdlog/spdlog.h>

#include "FuncName.hh"
#include "PrintHelper.hh"

#define ERROROUT 1

namespace MathTools
{
//______________________________________________________________________________
Bool_t SolveGaussJordan(const std::vector<Double_t>& z,
                        const std::vector<Double_t>& w,
                        const std::vector<Double_t>& s,
                        const std::vector<Double_t>& ct,
                        const std::vector<Double_t>& st,
                        Double_t& x0,
                        Double_t& u0,
                        Double_t& y0,
                        Double_t& v0)
{
  const Size_t n = z.size();
  if(n != w.size()  ||
     n != s.size()  ||
     n != ct.size() ||
     n != st.size()){
    std::ostringstream oss;
    oss << FUNC_NAME << " wrong vector size" << std::endl
        << " z.size()="  << n
        << " w.size()="  << w.size()
        << " s.size()="  << s.size()
        << " ct.size()=" << ct.size()
        << " st.size()=" << st.size();
    spdlog::warn(oss.str());
    return kFALSE;
  }

  Double_t matrx[16] = {};
  Double_t *mtp[4]   = {};
  Double_t fitp[4]   = {};
  mtp[0] = &matrx[ 0];
  mtp[1] = &matrx[ 4];
  mtp[2] = &matrx[ 8];
  mtp[3] = &matrx[12];

  for(Int_t i=0; i<4; ++i)
    for(Int_t j=0; j<4; ++j)
      mtp[i][j] = 0.0;

  for(Size_t i=0; i<n; ++i){
    Double_t ww=w[i], zz=z[i], ss=s[i], ctt=ct[i], stt=st[i];
    mtp[0][0] += ww*ctt*ctt;
    mtp[0][1] += ww*zz*ctt*ctt;
    mtp[0][2] += ww*ctt*stt;
    mtp[0][3] += ww*zz*ctt*stt;
    mtp[1][1] += ww*zz*zz*ctt*ctt;
    mtp[1][2] += ww*zz*ctt*stt;
    mtp[1][3] += ww*zz*zz*ctt*stt;
    mtp[2][2] += ww*stt*stt;
    mtp[2][3] += ww*zz*stt*stt;
    mtp[3][3] += ww*zz*zz*stt*stt;

    fitp[0] += ww*ss*ctt;
    fitp[1] += ww*zz*ss*ctt;
    fitp[2] += ww*ss*stt;
    fitp[3] += ww*zz*ss*stt;
  }
  mtp[1][0] = mtp[0][1];
  mtp[2][0] = mtp[0][2];
  mtp[3][0] = mtp[0][3];
  mtp[2][1] = mtp[1][2];
  mtp[3][1] = mtp[1][3];
  mtp[3][2] = mtp[2][3];

#if 0
  PrintMatrix(mtp, "Original Matrix: M_ij");
  PrintVector(fitp, "Original Vector: Q_i");
#endif

  Double_t Org[4][4]  = {{}};
  Double_t Org_vec[4] = {};
  for(Int_t i=0; i<4; ++i){
    Org_vec[i] = fitp[i];
    for(Int_t j=0; j<4; ++j){
      Org[i][j] = mtp[i][j];
    }
  }

  Int_t indxr[4] = {};
  Int_t indxc[4] = {};
  Int_t ipiv[4]  = {};
  if(!GaussJordan(mtp, 4, fitp, indxr, indxc, ipiv)){
    hddaq::cerr << FUNC_NAME << " Fitting failed" << std::endl;
    return kFALSE;
  }

  x0 = fitp[0];
  u0 = fitp[1];
  y0 = fitp[2];
  v0 = fitp[3];

  Double_t Inv[4][4] = {{}};
  Double_t Sol[4]    = {};
  for(Int_t i=0; i<4; ++i){
    for(Int_t j=0; j<4; ++j){
      Inv[i][j] = mtp[i][j];
      Sol[i] += Inv[i][j]*Org_vec[j];
    }
  }

  Double_t Final[4][4] = {{}};
  for(Int_t i=0; i<4; ++i){
    for(Int_t j=0; j<4; ++j){
      for(Int_t k=0; k<4; ++k){
        Final[i][j] += Inv[i][k]*Org[k][j];
      }
      if(Final[i][j]<1.e-10) Final[i][j]=0.;
    }
  }

#if 0
  PrintVector(indxr, "indxr[4] Gauss-Jordan");
  PrintVector(indxc, "indxc[4] Gauss-Jordan");
  PrintVector(ipiv,  "ipiv[4] Gauss-Jordan");
  PrintMatrix(Inv,   "Inverse Matrix");
  PrintVector(fitp,  "Solution from Gauss-Jordan");
  PrintVector(Sol,   "Solution from (M^-1)*Q_i");
  PrintMatrix(Final, "(M^-1)*M");
#endif

  return kTRUE;
}

//______________________________________________________________________________
Bool_t
GaussElim(Double_t **a, Int_t n, Double_t *b, Int_t *indx, Int_t *ipiv)
{
  Double_t big, c, pivinv, sum, c2;
  Int_t js   = 0;
  Int_t irow = 0;

  for(Int_t i=0; i<n; ++i) ipiv[i]=0;
  for(Int_t i=0; i<n; ++i){
    big = 0.;
    for(Int_t j=0; j<n; ++j){
      if(!ipiv[j]){
        if((c=TMath::Abs(a[j][i]))>=big){ big=c; irow=j; }
      }
      else if(ipiv[j]>1){
#ifdef ERROROUT
        std::ostringstream oss;
        oss << FUNC_NAME << " Singular Matrix" << std::endl;
        spdlog::warn(oss.str());
#endif
        return kFALSE;
      }
    }
    ++(ipiv[irow]); indx[i]=irow;
    if(a[irow][i]==0.0){
#ifdef ERROROUT
      std::ostringstream oss;
      oss << FUNC_NAME << " Singular Matrix" << std::endl;
      spdlog::warn(oss.str());
#endif
      return kFALSE;
    }
    pivinv=1.0/a[irow][i];

    for(Int_t j=0; j<n; ++j){
      if(ipiv[j]==0){
        c=a[j][i]; a[j][i]=0.0;
        for(Int_t k=i+1; k<n; ++k)
          a[j][k]-=a[irow][k]*pivinv*c;
        b[j]-=b[irow]*pivinv*c;
      }
    }
  }

  b[indx[n-1]]/=a[indx[n-1]][n-1];
  for(Int_t i=n-2; i>=0; --i){
    sum=b[indx[i]];
    for(Int_t j=i+1; j<n; ++j)
      sum-=a[indx[i]][j]*b[indx[j]];
    b[indx[i]]=sum/a[indx[i]][i];
  }
  for(Int_t i=0; i<n; ++i){
    if(indx[i]!=i){
      c2=b[i];
      for(Int_t j=indx[i]; indx[j]!=j; j=js){
        c=b[j]; b[j]=c2; c2=c; js=indx[j]; indx[j]=j;
      }
    }
  }
  return kTRUE;
}

//______________________________________________________________________________
Bool_t
GaussJordan(Double_t **a, // matrix a[n][n]
            Int_t n,      // dimension
            Double_t *b,  // vector b[n]
            Int_t *indxr, Int_t *indxc, Int_t *ipiv)
{
  for(Int_t i=0; i<n; ++i)
    ipiv[i]=0;

  for(Int_t i=0; i<n; ++i){
    Double_t big=0.0;
    Int_t irow=-1, icol=-1;
    for(Int_t j=0; j<n; ++j)
      if(ipiv[j]!=1)
        for(Int_t k=0; k<n; ++k){
          if(ipiv[k]==0){
            if(TMath::Abs(a[j][k])>=big){
              big=TMath::Abs(a[j][k]);
              irow=j; icol=k;
            }
          }
          else if(ipiv[k]>1){
#ifdef ERROROUT
            std::ostringstream oss;
            oss << FUNC_NAME << " Singular Matrix" << std::endl;
            spdlog::warn(oss.str());
#endif
            return kFALSE;
          }
        }
    ++(ipiv[icol]);

    if(irow!=icol){
      for(Int_t k=0; k<n; ++k){
        Double_t ta=a[irow][k];
        a[irow][k]=a[icol][k];
        a[icol][k]=ta;
      }
      Double_t tb=b[irow];
      b[irow]=b[icol];
      b[icol]=tb;
    }

    indxr[i]=irow; indxc[i]=icol;

    if(a[icol][icol]==0.0){
#ifdef ERROROUT
      std::ostringstream oss;
      oss << FUNC_NAME << " Singular Matrix" << std::endl;
      spdlog::warn(oss.str());
#endif
      return kFALSE;
    }
    Double_t pivinv=1./a[icol][icol];
    a[icol][icol]=1.;
    for(Int_t k=0; k<n; ++k) a[icol][k]*=pivinv;
    b[icol]*=pivinv;
    for(Int_t k=0; k<n; ++k){
      if(k!=icol){
        Double_t d=a[k][icol];
        a[k][icol]=0.;
        for(Int_t l=0; l<n; ++l) a[k][l] -= a[icol][l]*d;
        b[k] -= b[icol]*d;
      }
    }
  }

  for(Int_t l=n-1; l>=0; --l){
    if(indxr[l]!=indxc[l]){
      for(Int_t k=0; k<n; ++k){
        Double_t t=a[k][indxr[l]];
        a[k][indxr[l]]=a[k][indxc[l]];
        a[k][indxc[l]]=t;
      }
    }
  }
  return kTRUE;
}

//______________________________________________________________________________
Bool_t
InterpolateRatio(Int_t n, const Double_t *xa, const Double_t *ya,
                 Double_t *w1, Double_t *w2,
                 Double_t x, Double_t &y, Double_t &dy)
{
  Int_t i, m, ns=1;
  Double_t w, t, hh, h, dd;

  hh=TMath::Abs(x-xa[0]);
  for(i=1; i<=n; ++i){
    h=TMath::Abs(x-xa[i-1]);
    if(h==0.0) { y=ya[i-1]; dy=0.0; return kTRUE; }
    else if(h<hh){ ns=i; hh=h; }
    w1[i-1]=ya[i-1]; w2[i-1]=ya[i-1]*(1.+Epsilon());
  }
  y=ya[ns-1]; ns--;
  for(m=1; m<n; ++m){
    for(i=1; i<=n-m; ++i){
      w=w1[i]-w2[i-1]; h=xa[i+m-1]-x;
      t=(xa[i-1]-x)*w2[i-1]/h; dd=t-w1[i];
      if(dd==0.0){
#ifdef ERROROUT
        std::ostringstream oss;
        oss << FUNC_NAME << " : Error" << std::endl;
        spdlog::warn(oss.str());
#endif
        y=Infinity(); dy=Infinity(); return kFALSE;
      }
      dd=w/dd; w2[i-1]=w1[i]*dd; w1[i-1]=t*dd;
    }
    if(2*ns < (n-m)) dy=w1[ns];
    else { dy=w2[ns-1]; ns--; }
    y+=dy;
  }
#if 0
  hddaq::cout << FUNC_NAME << ": x=" << std::setw(10) << x
              << " y=" << std::setw(10) << y << std::endl;
#endif
  return kTRUE;
}

//______________________________________________________________________________
Bool_t
InterpolatePol(Int_t n, const Double_t *xa, const Double_t *ya,
               Double_t *w1, Double_t *w2,
               Double_t x, Double_t &y, Double_t &dy)
{
  Int_t i, m, ns=1;
  Double_t den, dif, dift, ho, hp, w;

  dif=TMath::Abs(x-xa[0]);
  for(i=1; i<=n; ++i){
    if((dift=TMath::Abs(x-xa[i-1]))<dif){ ns=i; dif=dift; }
    w1[i-1]=w2[i-1]=ya[i-1];
  }
  y=ya[ns-1]; --ns;
  for(m=1; m<n; ++m){
    for(i=1; i<=n-m; ++i){
      ho=xa[i-1]-x; hp=xa[i+m-1]-x;
      w=w1[i]-w2[i-1];
      den=ho-hp;
      if(den==0.0){
#ifdef ERROROUT
        std::ostringstream oss;
        oss << FUNC_NAME << " : Error" << std::endl;
        spdlog::warn(oss.str());
#endif
        y=Infinity(); dy=Infinity(); return kFALSE;
      }
      den=w/den;
      w2[i-1]=hp*den; w1[i-1]=ho*den;
    }
    if(2*ns<(n-m)) dy=w1[ns];
    else           { dy=w2[ns-1]; --ns; }
    y+=dy;
  }
  return kFALSE;
}

//______________________________________________________________________________
Bool_t
SVDksb(Double_t **u, const Double_t *w, Double_t **v,
       Int_t m, Int_t n, const Double_t *b, Double_t *x, Double_t *wv)
{
  for(Int_t j=0; j<n; ++j){
    Double_t s=0.0;
    if(w[j]!=0.0){
      for(Int_t i=0; i<m; ++i)
        s += u[i][j]*b[i];
      s /= w[j];
    }
    wv[j]=s;
  }
  for(Int_t i=0; i<n; ++i){
    Double_t s=0.0;
    for(Int_t j=0; j<n; ++j)
      s += v[i][j]*wv[j];
    x[i]=s;
  }
  return kTRUE;
}

//______________________________________________________________________________
inline Double_t
pythag(Double_t a, Double_t b)
{
  Double_t aa=TMath::Abs(a), ab=TMath::Abs(b);
  if(aa>ab)
    return aa*TMath::Sqrt(1.+(ab/aa)*(ab/aa));
  else if(ab!=0.)
    return ab*TMath::Sqrt(1.+(aa/ab)*(aa/ab));
  else
    return 0.0;
}

//______________________________________________________________________________
Bool_t
SVDcmp(Double_t **a, Int_t m, Int_t n, Double_t *w, Double_t **v, Double_t *wv)
{
  Double_t g     = 0.;
  Double_t scale = 0.;
  Double_t anorm = 0.;
  Double_t s, f, h, c;
  Int_t nm;

#ifdef DebugPrint
  for(Int_t i=0; i<n; ++i) {
    w[i]=wv[i]=0.0;
    for(Int_t j=0; j<n; ++j) v[j][i]=0.0;
  }

  {
    PrintHelper helper(3, std::ios::scientific);
    hddaq::cout << FUNC_NAME << ": A in SVDcmp 1" <<  std::endl;
    for(Int_t ii=0; ii<m; ++ii){
      for(Int_t ij=0; ij<n; ++ij){
        hddaq::cout << std::setw(12) << a[ii][ij];
        if(ij!=n-1) hddaq::cout << ",";
      }
      hddaq::cout << std::endl;
    }
    hddaq::cout << FUNC_NAME << ": V in SVDcmp 1" << std::endl;
    for(Int_t ii=0; ii<n; ++ii){
      for(Int_t ij=0; ij<n; ++ij){
        hddaq::cout << std::setw(12) << v[ii][ij];
        if(ij!=n-1) hddaq::cout << ",";
      }
      hddaq::cout << std::endl;
    }
    hddaq::cout << FUNC_NAME << ": W in SVDcmp 1" << std::endl;
    for(Int_t ij=0; ij<n; ++ij){
      hddaq::cout << std::setw(12) << w[ij];
      if(ij!=n-1) hddaq::cout << ",";
    }
    hddaq::cout << std::endl;

    hddaq::cout << FUNC_NAME << ": WV in SVDcmp 1" << std::endl;
    for(Int_t ij=0; ij<n; ++ij){
      hddaq::cout << std::setw(12) << wv[ij];
      if(ij!=n-1) hddaq::cout << ",";
    }
    hddaq::cout << std::endl;
    hddaq::cout << std::endl;
  }
#endif

  // Householder method
  for(Int_t i=0; i<n; ++i){
    wv[i]=scale*g;
    g = scale = 0.0;
    if(i<m){
      for(Int_t k=i; k<m; ++k) scale += TMath::Abs(a[k][i]);
      if(scale!=0.){
        s = 0;
        for(Int_t k=i; k<m; ++k){
          a[k][i] /= scale;
          s += a[k][i]*a[k][i];
        }
        f = a[i][i];
        g = ((f>0.0) ? -TMath::Sqrt(s) : TMath::Sqrt(s));
        h = f*g-s;
        a[i][i] = f-g;
        for(Int_t j=i+1; j<n; ++j){
          s = 0.0;
          for(Int_t k=i; k<m; ++k) s += a[k][i]*a[k][j];
          f = s/h;
          for(Int_t k=i; k<m; ++k)  a[k][j] += f*a[k][i];
        }
        for(Int_t k=i; k<m; ++k) a[k][i] *= scale;
      }
    }     /* if(i<m) */
    w[i] = scale*g;
    g = s = scale = 0.0;

    if(i<m && i!=n-1){
      for(Int_t k=i+1; k<n; ++k) scale += TMath::Abs(a[i][k]);
      if(scale!=0.0){
        for(Int_t k=i+1; k<n; ++k){
          a[i][k] /= scale;
          s += a[i][k]*a[i][k];
        }
        f = a[i][i+1];
        g = ((f>0.0) ? -TMath::Sqrt(s) : TMath::Sqrt(s));
        h = f*g-s;
        a[i][i+1] = f-g;
        for(Int_t k=i+1; k<n; ++k) wv[k] = a[i][k]/h;
        for(Int_t j=i+1; j<m; ++j){
          s = 0.0;
          for(Int_t k=i+1; k<n; ++k) s += a[j][k]*a[i][k];
          for(Int_t k=i+1; k<n; ++k) a[j][k] += s*wv[k];
        }
        for(Int_t k=i+1; k<n; ++k) a[i][k] *= scale;
      }
    }   /* if(i<m && i!=n-1) */
    Double_t tmp=TMath::Abs(w[i])+TMath::Abs(wv[i]);
    if(tmp>anorm) anorm = tmp;
  }     /* for(Int_t i ...) */

#if DebugPrint
  {
    PrintHelper helper(3, std::ios::scientific);
    hddaq::cout << FUNC_NAME << ": A in SVDcmp 2" <<  std::endl;
    for(Int_t ii=0; ii<m; ++ii){
      for(Int_t ij=0; ij<n; ++ij){
        hddaq::cout << std::setw(12) << a[ii][ij];
        if(ij!=n-1) hddaq::cout << ",";
      }
      hddaq::cout << std::endl;
    }
    hddaq::cout << FUNC_NAME << ": V in SVDcmp 2" << std::endl;
    for(Int_t ii=0; ii<n; ++ii){
      for(Int_t ij=0; ij<n; ++ij){
        hddaq::cout << std::setw(12) << v[ii][ij];
        if(ij!=n-1) hddaq::cout << ",";
      }
      hddaq::cout << std::endl;
    }
    hddaq::cout << FUNC_NAME << ": W in SVDcmp 2" << std::endl;
    for(Int_t ij=0; ij<n; ++ij){
      hddaq::cout << std::setw(12) << w[ij];
      if(ij!=n-1) hddaq::cout << ",";
    }
    hddaq::cout << std::endl;

    hddaq::cout << FUNC_NAME << ": WV in SVDcmp 2" << std::endl;
    for(Int_t ij=0; ij<n; ++ij){
      hddaq::cout << std::setw(12) << wv[ij];
      if(ij!=n-1) hddaq::cout << ",";
    }
    hddaq::cout << std::endl;
    hddaq::cout << std::endl;
  }
#endif

  for(Int_t i=n-1; i>=0; --i){
    if(i<n-1){
      if(g!=0.0){
        for(Int_t j=i+1; j<n; ++j) v[j][i] = (a[i][j]/a[i][i+1])/g;
        for(Int_t j=i+1; j<n; ++j){
          s = 0.0;
          for(Int_t k=i+1; k<n; ++k) s += a[i][k]*v[k][j];
          for(Int_t k=i+1; k<n; ++k) v[k][j] += s*v[k][i];
        }
      }
      for(Int_t j=i+1; j<n; ++j)
        v[i][j] = v[j][i] = 0.0;
    }
    v[i][i]=1.0;  g=wv[i];
  }   /* for(Int_t i= ...) */

#ifdef DebugPrint
  {
    PrintHelper helper(3, std::ios::scientific);
    hddaq::cout << FUNC_NAME << ": A in SVDcmp 3" <<  std::endl;
    for(Int_t ii=0; ii<m; ++ii){
      for(Int_t ij=0; ij<n; ++ij){
        hddaq::cout << std::setw(12) << a[ii][ij];
        if(ij!=n-1) hddaq::cout << ",";
      }
      hddaq::cout << std::endl;
    }
    hddaq::cout << FUNC_NAME << ": V in SVDcmp 3" << std::endl;
    for(Int_t ii=0; ii<n; ++ii){
      for(Int_t ij=0; ij<n; ++ij){
        hddaq::cout << std::setw(12) << v[ii][ij];
        if(ij!=n-1) hddaq::cout << ",";
      }
      hddaq::cout << std::endl;
    }
    hddaq::cout << FUNC_NAME << ": W in SVDcmp 3" << std::endl;
    for(Int_t ij=0; ij<n; ++ij){
      hddaq::cout << std::setw(12) << w[ij];
      if(ij!=n-1) hddaq::cout << ",";
    }
    hddaq::cout << std::endl;

    hddaq::cout << FUNC_NAME << ": WV in SVDcmp 3" << std::endl;
    for(Int_t ij=0; ij<n; ++ij){
      hddaq::cout << std::setw(12) << wv[ij];
      if(ij!=n-1) hddaq::cout << ",";
    }
    hddaq::cout << std::endl;
    hddaq::cout << std::endl;
  }
#endif

  Int_t mn = ((m<n) ? m : n);

  for(Int_t i=mn-1; i>=0; --i){
    g=w[i];
    for(Int_t j=i+1; j<n; ++j) a[i][j]=0.0;
    if(g!=0.0){
      g = 1./g;
      for(Int_t j=i+1; j<n; ++j){
        s = 0.0;
        for(Int_t k=i+1; k<m; ++k) s += a[k][i]*a[k][j];
        f = (s/a[i][i])*g;
        for(Int_t k=i; k<m; ++k) a[k][j] += f*a[k][i];
      }
      for(Int_t j=i; j<m; ++j) a[j][i] *= g;
    }
    else
      for(Int_t j=i; j<m; ++j) a[j][i] = 0.0;

    a[i][i] += 1.0;
  }   /* for(Int_t i= ...) */

#ifdef DebugPrint
  {
    PrintHelper helper(3, std::ios::scientific);
    hddaq::cout << FUNC_NAME << ": A in SVDcmp 4" <<  std::endl;
    for(Int_t ii=0; ii<m; ++ii){
      for(Int_t ij=0; ij<n; ++ij){
        hddaq::cout << std::setw(12) << a[ii][ij];
        if(ij!=n-1) hddaq::cout << ",";
      }
      hddaq::cout << std::endl;
    }
    hddaq::cout << FUNC_NAME << ": V in SVDcmp 4" << std::endl;
    for(Int_t ii=0; ii<n; ++ii){
      for(Int_t ij=0; ij<n; ++ij){
        hddaq::cout << std::setw(12) << v[ii][ij];
        if(ij!=n-1) hddaq::cout << ",";
      }
      hddaq::cout << std::endl;
    }
    hddaq::cout << FUNC_NAME << ": W in SVDcmp 4" << std::endl;
    for(Int_t ij=0; ij<n; ++ij){
      hddaq::cout << std::setw(12) << w[ij];
      if(ij!=n-1) hddaq::cout << ",";
    }
    hddaq::cout << std::endl;

    hddaq::cout << FUNC_NAME << ": WV in SVDcmp 4" << std::endl;
    for(Int_t ij=0; ij<n; ++ij){
      hddaq::cout << std::setw(12) << wv[ij];
      if(ij!=n-1) hddaq::cout << ",";
    }
    hddaq::cout << std::endl;
    hddaq::cout << std::endl;
  }
#endif


  Int_t ll=1;

  for(Int_t k=n-1; k>=0; --k){
    for(Int_t its=1; its<=30; ++its){
      Int_t flag=1; nm=ll;
      for(ll=k; ll>=0; --ll){
        nm=ll-1;
        if(TMath::Abs(wv[ll])+anorm == anorm){
          flag=0; break;
        }
        if(TMath::Abs(w[nm])+anorm == anorm)
          break;
      }

      if(flag){
        c=0.0; s=1.0;
        for(Int_t i=ll; i<=k; ++i){
          f = s*wv[i]; wv[i] *= c;
          if(TMath::Abs(f)+anorm == anorm)
            break;
          g=w[i]; h=pythag(f,g); w[i]=h;
          h=1./h; c=g*h; s=-f*h;
          for(Int_t j=0; j<m; ++j){
            Double_t y=a[j][nm], z=a[j][i];
            a[j][nm]=y*c+z*s; a[j][i]=z*c-y*s;
          }
        }
      }   /* if(flag) */

      Double_t z = w[k];
      if(ll==k){
        if(z<0.){
          w[k]=-z;
          for(Int_t j=0; j<n; ++j) v[j][k]=-v[j][k];
        }
        break;
      }
#ifdef ERROROUT
      if(its==30){
        return kFALSE;
      }
#endif
      nm=k-1;
      Double_t x=w[ll], y=w[nm];
      g=wv[nm]; h=wv[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.*h*y);
      g=pythag(f,1.0);
      Double_t gtmp = ((f>0.) ? g : -g);
      f=((x-z)*(x+z)+h*((y/(f+gtmp))-h))/x;
      c=s=1.0;
      for(Int_t j=ll; j<=nm; ++j){
        g=wv[j+1]; y=w[j+1]; h=s*g; g=c*g;
        z=pythag(f,h);
        wv[j]=z; c=f/z; s=h/z;
        f=x*c+g*s; g=g*c-x*s;
        h=y*s; y=y*c;
        for(Int_t jj=0; jj<n; ++jj){
          x=v[jj][j]; z=v[jj][j+1];
          v[jj][j]=x*c+z*s; v[jj][j+1]=z*c-x*s;
        }
        z=pythag(f,h);
        w[j]=z;
        if(z!=0.0){ z=1./z; c=f*z; s=h*z; }
        f=c*g+s*y; x=c*y-s*g;
        for(Int_t jj=0; jj<m; ++jj){
          y=a[jj][j]; z=a[jj][j+1];
          a[jj][j]=y*c+z*s; a[jj][j+1]=z*c-y*s;
        }
      }
      wv[ll]=0.0; wv[k]=f; w[k]=x;
    }   /* for(Int_t its ...) */
  }     /* for(Int_t k= ...) */


#ifdef DebugPrint
  {
    PrintHelper helper(3, std::ios::scientific);
    hddaq::cout << FUNC_NAME << ": A in SVDcmp 5" <<  std::endl;
    for(Int_t ii=0; ii<m; ++ii){
      for(Int_t ij=0; ij<n; ++ij){
        hddaq::cout << std::setw(12) << a[ii][ij];
        if(ij!=n-1) hddaq::cout << ",";
      }
      hddaq::cout << std::endl;
    }
    hddaq::cout << FUNC_NAME << ": V in SVDcmp 5" << std::endl;
    for(Int_t ii=0; ii<n; ++ii){
      for(Int_t ij=0; ij<n; ++ij){
        hddaq::cout << std::setw(12) << v[ii][ij];
        if(ij!=n-1) hddaq::cout << ",";
      }
      hddaq::cout << std::endl;
    }
    hddaq::cout << FUNC_NAME << ": W in SVDcmp 5" << std::endl;
    for(Int_t ij=0; ij<n; ++ij){
      hddaq::cout << std::setw(12) << w[ij];
      if(ij!=n-1) hddaq::cout << ",";
    }
    hddaq::cout << std::endl;

    hddaq::cout << FUNC_NAME << ": WV in SVDcmp 5" << std::endl;
    for(Int_t ij=0; ij<n; ++ij){
      hddaq::cout << std::setw(12) << wv[ij];
      if(ij!=n-1) hddaq::cout << ",";
    }
    hddaq::cout << std::endl;
    hddaq::cout << std::endl;
  }
#endif

  return kTRUE;
}

//______________________________________________________________________________
template <typename T>
void
PrintMatrix(T *mat, const std::string& arg,
            const Int_t column, const Int_t line)
{
  PrintHelper helper(5, std::ios::scientific);

  hddaq::cout << FUNC_NAME << " " << arg << std::endl;
  for(Int_t l=0; l<line; ++l){
    for(Int_t c=0; c<line; ++c){
      hddaq::cout << "  " << std::setw(12) << mat[l][c];
    }
    hddaq::cout << std::endl;
  }
}

//______________________________________________________________________________
template <typename T>
void
PrintVector(T *vec, const std::string& arg, const Size_t size)
{
  PrintHelper helper(5, std::ios::scientific);

  hddaq::cout << FUNC_NAME << " " << arg << std::endl;
  for(Size_t i=0; i<size; ++i){
    hddaq::cout << "  " << std::setw(12) << vec[i];
  }
  hddaq::cout << std::endl;
}
}
