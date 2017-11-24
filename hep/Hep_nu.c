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
#include "Functions.h"
#include "EnergyRangeCorrelation.C"
#include "TRandom3.h"
#include <TSystem.h>



using namespace std;

void Hep_nu() {

  const int tot_line=1000;
  Int_t iline=0;
  Int_t resto=-1;
  Double_t sum_flux=0;
  Double_t tmp_E_nu=0;
  Double_t tmp_nu_rate=0;
  Double_t nu_en[tot_line]={};
  Double_t nu_rate[tot_line]={};
  
  TString dir = gSystem->UnixPathName(__FILE__);
  dir.ReplaceAll("Hep_nu.c","");
  dir.ReplaceAll("/./","/");
  ifstream in;
  in.open(Form("%shep_spectrum2.txt",dir.Data()));
  
  // TGraph * gr = new TGraph();
  //gr->SetName("gr");
  while (1) {
    in >> nu_en[iline] >> nu_rate[iline];
    if (!in.good()) break;
    cout << iline << " " << nu_en[iline] << " " << nu_rate[iline] << endl;
    //gr->SetPoint(iline,nu_en[iline],nu_rate[iline]);
    iline++;
  }

  TGraph *gr = new TGraph(iline,nu_en,nu_rate);
  //TCanvas c1 = new TCanvas("c1","c1",600,600);
  gr->Draw("A*");
  
  
  
  ofstream log_file;
  log_file.open ("num_eventi.txt");
  log_file << "HEP induced recoil rate";
  log_file << endl;
  
  Functions* functions;

  functions = new Functions();
  

  //for(int i=2; i<3;i++){   // for Xenon cross-check    (SEE ALSO DEFINITIONS2.H)
  for(int i=0; i<n_el;i++){    
    iel = el[i];
    iZ_el = Z_el[i];
    iA_el = A_el[i];
    ifract_el = atomic_fract_el[i]; // uma
    imass_fract = mass_fract[i];
    iM_el=iA_el*uma;
    iN_el = N_el[i];
    num_ev_bkg=0;
    tmp_E_nu=0;
    sum_flux=0;
    tmp_nu_rate=0;
    //cout << iel << " " << iZ_el << " " << iA_el << " " << ifract_el << " " << imass_fract << " " << iM_el << " " << iN_el << endl;
    
    
    for(int j=0; j<tot_line;j++){
      resto = j % 5;
      if(resto==0){
	E_nu=nu_en[j];
	sum_crsect=0;
	//if(j==0)iflux_hep_vs_E=nu_rate[j]*integral_flux_hep*nu_en[iline];
	iflux_hep_vs_E=nu_rate[j]*integral_flux_hep*(E_nu-tmp_E_nu);
	if(j>0)sum_flux += iflux_hep_vs_E;
	//if(j==0)cout << j << " " << E_nu  << " " << tmp_E_nu <<  endl;

      
	for(int k=0;k<d_theta;k++){
	  
	  cos_theta = -1 + 1./d_theta + k*(2./d_theta); // devi moltiplicare per delta theta		
	  E_rec = (E_nu*E_nu*(1-cos_theta))/(E_nu*(1-cos_theta)+iM_el*TMath::Power(10,3)); // MeV
	  crsect = functions->cross_section();
	  //cout << E_rec << endl;
	  /// XENON CROSS-CHECK
	  //if(E_rec>0.002){ 
	  //  num_ev_bkg = functions->num_eventi_bkg();
	  // sum_ev_bkg = sum_ev_bkg + num_ev_bkg;
	  //  }
	  /////
	
	  L_true = GetRangeFromEnergy(iZ_el,E_rec*1000);
	  L_proj = L_true*cos_theta_rec;
	  //cout << j << " " << k << " " << crsect << " " << dN_vs_dE_bkg << " " << num_ev_bkg << " " << sum_ev_bkg << endl;
	  
	  if(L_true > 100 && L_true<1500 && E_rec>0.002){
	    num_ev_bkg = functions->num_eventi_bkg();
	  sum_ev_bkg = sum_ev_bkg + num_ev_bkg;
	  h_len->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);
	  }
	}
	tmp_E_nu=nu_en[j];
	tmp_nu_rate=nu_rate[j];
      }
      
    }
    cout << iel << " " << sum_ev_bkg << " " << sum_flux << endl;
  }
  num_ev_bkg=sum_ev_bkg*86400*365;
  cout << num_ev_bkg << endl;
  h_len->Draw("");
  h_len->Scale(num_ev_bkg/h_len->Integral());
  delete functions;
}
