#include <stdio.h>
#include <ctime>
#include <cstdlib>
#include <string>
#include "TMath.h"
#include "TString.h"
#include "TMultiGraph.h"

//#include "Definitions_sn2.h"

using namespace std;

class Functions
{
 public: 


  inline double cross_section(){

    double d_ene=0.1; // MeV
    iN_el=iA_el-iZ_el;
    //cout << iN_el << endl;
    iQw = iN_el - (1-4*sin2wnb)*iZ_el;
    irn = TMath::Sqrt(1.5129*TMath::Power(iA_el,2./3.) -1.4760*TMath::Power(iA_el,1./3.) + 2.5371);

    q = TMath::Sqrt(2*iM_el*(E_nu*E_nu*TMath::Power(10,-6)*(1-cos_theta))/(E_nu*TMath::Power(10,-3)*(1-cos_theta)+iM_el));
    qrn = q*irn*fm_in_GeV;
    qs = q*s*fm_in_GeV;
    Fq = (3*(TMath::Sin(qrn)-qrn*TMath::Cos(qrn))/TMath::Power(qrn,3))*TMath::Exp(-(qs*qs)/2.);  // diverso rispetto a plot di esclusione
    //cout << q << " " << qrn << " " << qs << " " << Fq << endl;
    crsect = ((GF*GF*iQw*iQw*E_nu*E_nu*Fq*Fq*(1+cos_theta)*(2./d_theta))/(8*pi))*TMath::Power(iMeV_in_fm,2)*TMath::Power(10,-26)*d_ene;

    return crsect;
  }


  inline double num_eventi_bkg(){
    
    int_const_bkg = (M_riv*imass_fract)/(iA_el*uma_in_ton);
    dN_vs_dE_bkg = int_const_bkg*crsect*iflux_boron_vs_E;    
    return dN_vs_dE_bkg;
  }

  
 protected:

 private:
  
};



