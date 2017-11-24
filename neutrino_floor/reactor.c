#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include <iomanip>      // std::setprecision
#include <fstream>
#include <vector>
#include <time.h>

#include "TH1F.h"
#include "TMath.h"
#include "TComplex.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TString.h"
#include "TLine.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TPaveStats.h"

#include "Definitions_sn2.h"
//#include "Functions.h"
#include "EnergyRangeCorrelation.C"
#include "TRandom3.h"
#include <TSystem.h>



using namespace std;

void reactor() {

  Int_t resto=0;
  Double_t reactor_antinu_en[15];
  Double_t reactor_antinu_rate[15];
  Int_t iline=0;
  Double_t reactor_antinu_flux=TMath::Power(10,18);
  
  TString dir = gSystem->UnixPathName(__FILE__);
  dir.ReplaceAll("reactor.c","");
  dir.ReplaceAll("/./","/");
  ifstream in;

  in.open(Form("%santinu_reactor_spectrum.txt",dir.Data()));
  while (1) {
    in >> reactor_antinu_en[iline] >> reactor_antinu_rate[iline];
    reactor_antinu_rate[iline]=reactor_antinu_rate[iline]*reactor_antinu_flux;
    //cout << reactor_antinu_en[iline] << " " << reactor_antinu_rate[iline] << endl;
    if (!in.good()) break;
    iline++;
  }


  TH1D * h_wimp_ev_Erec = new TH1D("wimp_ene","",100000,0.001,100);
  TH1D * h_Erec = new TH1D("ene","",100000,0.001,100);
  TH1D * h_Ltrue = new TH1D("len2","",1000,0,1000);
  TH1D * h_antinu = new TH1D("len","",100,1,10);
  Int_t ipoint=0;

  
  for(int k=0;k<20000;k++){
    num_ev_bkg=0;
    num_wimp_ev=0;
    E_rec=k/1000.;
    //sum_crsect=0;
    
    for(int i=0; i<n_el;i++){   
      iel = el[i];
      iZ_el = Z_el[i];
      iA_el = A_el[i];
      ifract_el = atomic_fract_el[i]; // uma
      imass_fract = mass_fract[i];
      iM_el=iA_el*uma;

      L_true = GetRangeFromEnergy(iZ_el,E_rec);
      //cout << E_rec << " " << L_true << endl;
      
      iN_el=iA_el-iZ_el;
      iQw = iN_el - (1-4*sin2wnb)*iZ_el;
      irn = TMath::Sqrt(1.5129*TMath::Power(iA_el,2./3.) -1.4760*TMath::Power(iA_el,1./3.) + 2.5371);
      q=TMath::Sqrt((2*iM_el*E_rec*TMath::Power(10,-6))); //GeV
      qrn = q*irn*fm_in_GeV;    
      qs = q*s*fm_in_GeV;       
      Fq = (3*(TMath::Sin(qrn)-qrn*TMath::Cos(qrn))/TMath::Power(qrn,3))*TMath::Exp(-(qs*qs)/2.);  // diverso rispetto a plot di esclusione
      int_const_bkg = (M_riv*imass_fract)/(iA_el*uma_in_ton);
      for(int j=0;j<15;j++){
	E_nu=reactor_antinu_en[j];
	iflux_boron_vs_E=reactor_antinu_rate[j];
	if(k==0 && i==0)h_antinu->Fill(E_nu,iflux_boron_vs_E);
	crsect = ((GF*GF*iQw*iQw*Fq*Fq*iM_el*1000)/(4*pi))*(1-(iM_el*E_rec)/(2*E_nu*E_nu))*TMath::Power(iMeV_in_fm,2)*TMath::Power(10,-26)/TMath::Power(10,6);	
	if(E_nu >(TMath::Sqrt(iM_el*E_rec/2.)) && L_true>100){	  
	  num_ev_bkg += crsect*iflux_boron_vs_E*int_const_bkg;
	}
      }
    }
    cout << E_rec << " " << L_true << endl;
    //num_ev_bkg = num_ev_bkg*86400*365*1000; //TMath::Power(10,33)/1.66054;
    num_ev_bkg = num_ev_bkg*86400*7*10; //10Kg per week
    h_Erec->Fill(E_rec,num_ev_bkg); //num_ev_bkg
    h_Ltrue->Fill(L_true,num_ev_bkg); //num_ev_bkg
    sum_ev_bkg+=num_ev_bkg;
  }


  TCanvas *c3 = new TCanvas("c3","c3",1200,1200);
  c3->Divide(2,2);
  c3->cd(1);
  h_antinu->Draw("");
  c3->cd(2);
  h_Erec->Draw("");
  c3->cd(3);
  h_Ltrue->Draw("");
}
