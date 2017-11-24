#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include <iomanip>      // std::setprecision
#include <fstream>
#include <vector>
#include <time.h>

#include "TH1F.h"
#include "TH2D.h"
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
#include "EnergyRangeCorrelation.C"
#include "TRandom3.h"
#include <TSystem.h>



using namespace std;

void nu_floor_v2() {

  Double_t B_int_flux = 5.58*TMath::Power(10,6);
  Int_t B_iline=0;
  Double_t B_nu_en[820]={};
  Double_t B_nu_rate[820]={};

  Double_t sum_flux=0;
  
  TString dir = gSystem->UnixPathName(__FILE__);
  dir.ReplaceAll("nu_floor_v2.c","");
  dir.ReplaceAll("/./","/");
  ifstream in;
  
  in.open(Form("%sboron_spectrum.txt",dir.Data()));
  while (1) {
    in >> B_nu_en[B_iline] >> B_nu_rate[B_iline];
    B_nu_rate[B_iline]=B_nu_rate[B_iline]*B_int_flux;
    if (!in.good()) break;
    cout << B_iline << " " << B_nu_en[B_iline] << " " << B_nu_rate[B_iline] << endl;
    B_iline++;
    
  }

  in.close();
  
  TH1D * h_ev_Erec = new TH1D("ene","",1000,0,5);
  TH2D * h_q = new TH2D("q","q",1000,0,100,1000,0,1);
  
  
  for(int k=0;k<10000;k++){
    num_ev_bkg=0;
    E_rec=k/1000.;
    sum_crsect=0;
    sum_flux=0;
    for(int i=2; i<3;i++){   
      iel = el[i];
      iZ_el = Z_el[i];
      iA_el = A_el[i];
      ifract_el = atomic_fract_el[i]; // uma
      imass_fract = mass_fract[i];
      iM_el=iA_el*uma;
     
      iN_el=iA_el-iZ_el;
      iQw = iN_el - (1-4*sin2wnb)*iZ_el;
      irn = TMath::Sqrt(1.5129*TMath::Power(iA_el,2./3.) -1.4760*TMath::Power(iA_el,1./3.) + 2.5371);
      q=TMath::Sqrt((2*iM_el*E_rec*TMath::Power(10,-6))); //GeV
      qrn = q*irn*fm_in_GeV;    
      qs = q*s*fm_in_GeV;       
      Fq = (3*(TMath::Sin(qrn)-qrn*TMath::Cos(qrn))/TMath::Power(qrn,3))*TMath::Exp(-(qs*qs)/2.);  // diverso rispetto a plot di esclusione
      
      sum_flux=0;
      //h_q->Fill(q*1000,Fq);
      for(int j=0; j<818;j++){
	E_nu=B_nu_en[j];
	iflux_boron_vs_E = B_nu_rate[j]*0.02;
	sum_flux+=iflux_boron_vs_E;
	crsect = ((GF*GF)/4*pi)*iQw*iM_el*1000*(1-(iM_el*E_rec)/(2*E_nu*E_nu))*Fq*Fq*TMath::Power(iMeV_in_fm,2)*TMath::Power(10,-26)/TMath::Power(10,6);
	//cout << TMath::Sqrt((iM_el*E_rec)/2.) << endl;
	if(crsect>0 && E_nu>=TMath::Sqrt((iM_el*E_rec)/2.)){
	  int_const_bkg = (M_riv*imass_fract)/(iA_el*uma_in_ton);
	  sum_crsect = sum_crsect+crsect;
	  num_ev_bkg += int_const_bkg*crsect*iflux_boron_vs_E;
	  //cout << E_nu << " " << " " << E_rec << " " << sum_crsect << endl;
	}
	//	else cout << E_nu << " " << " " << E_rec << " " << endl;
      }
      // cout << iel << " " << iZ_el << " " << iA_el << " " << iN_el << " " << ifract_el << " " << imass_fract << " " << iM_el << endl;
      //cout << iQw << " " << irn << " " << q << " " <<  qrn << " " << qs << " " << Fq << endl;
    } 
    num_ev_bkg = num_ev_bkg*86400*365;
    //h_ev_Erec->Fill(E_rec,num_ev_bkg); //num_ev_bkg
    //sum_ev_bkg+=num_ev_bkg;
    //cout << sum_flux << endl;
    h_ev_Erec->Fill(E_rec,num_ev_bkg);
  }
  
  TCanvas *c3 = new TCanvas("c3","c3",1200,600);
  h_ev_Erec->Draw("");
  //h_q->Draw();
}
