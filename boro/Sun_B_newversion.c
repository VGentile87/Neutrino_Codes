#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include <iomanip>      // std::setprecision
#include <fstream>
#include <vector>
#include <time.h>

#include "TF1.h"
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

void Sun_B_newversion() {


  //// LEWIN-SMITH GALACTIC COORDINATES ////////////////////////////  
  Double_t bx=5.536*TMath::Pi()/180.;
  Double_t by=-59.574*TMath::Pi()/180.;
  Double_t bz=-29.811*TMath::Pi()/180.;
  Double_t lx=266.840*TMath::Pi()/180.;
  Double_t ly=347.340*TMath::Pi()/180.;
  Double_t lz=180.023*TMath::Pi()/180.;

  //cout << bx << " " << lx << endl;
  
  TF1 *fx = new TF1("fx","[0]*sin(x - [1] + [2])",0,2*TMath::Pi());
  TF1 *fy = new TF1("fy","[0]*sin(x - [1] + [2])",0,2*TMath::Pi());
  TF1 *fz = new TF1("fz","[0]*sin(x - [1] + [2])",0,2*TMath::Pi());
  fx->SetParameter(0,TMath::Cos(bx));
  fx->SetParameter(1,lx);
  fx->SetParameter(2,ly);
  fy->SetParameter(0,TMath::Cos(by));
  fy->SetParameter(1,ly);
  fy->SetParameter(2,ly);
  fz->SetParameter(0,TMath::Cos(bz));
  fz->SetParameter(1,lz);
  fz->SetParameter(2,ly);

  TF1 *fphi = new TF1("fphi","atan2(fx,fy)",0,2*TMath::Pi());
  TF1 *fthe = new TF1("fthe","atan2(fz,sqrt(fx*fx+fy*fy))",0,2*TMath::Pi());
  TH1F *hr = new TH1F("hr","hr",100,0,2*TMath::Pi());
  TH1F *hphi = new TH1F("hphi","hphi",100,-TMath::Pi(),TMath::Pi());
  //TH1F *hthe = new TH1F("hthe","hthe",100,-TMath::Pi()/2.,TMath::Pi()/2.);
  TH1F *hthe = new TH1F("hthe","hthe",200,0,1);

  TRandom3 *r = new TRandom3();
  Double_t t = 0;
  Double_t px = 0;
  Double_t py = 0;
  Double_t pz = 0;
  Double_t pphi = 0;
  Double_t pthe = 0;
  for(int i=0;i<180000;i++){
    t = r->Uniform(0,2*TMath::Pi());
    px = TMath::Cos(bx)*TMath::Sin(t-lx+ly);
    py = TMath::Cos(by)*TMath::Sin(t-ly+ly);
    pz = TMath::Cos(bz)*TMath::Sin(t-lz+ly);
    pphi = TMath::ATan2(px,py);
    pthe = TMath::ATan(pz/TMath::Sqrt(px*px + py*py));
    /*if(pthe > TMath::Pi()/2.) {
      cout << pthe << endl;
      pthe = pthe - TMath::Pi();
      cout << pthe << endl;
      }*/
    hr->Fill(t);
    hphi->Fill(pphi);
    hthe->Fill(TMath::Cos(pthe));
  }
  

  
  
  TCanvas *cc1 = new TCanvas("cc1","cc1",1200,600);
  cc1->Divide(2,1);
  cc1->cd(1);
  fphi->Draw("P");
  cc1->cd(2);
  fthe->Draw("P");

  TCanvas *cc2 = new TCanvas("cc2","cc2",1200,600);
  //hthe->Draw();
  cc2->Divide(2,1);
  cc2->cd(1);
  //7hr->Draw("");
  //cc2->cd(2);
  hphi->Draw("");
  cc2->cd(2);
  hthe->Draw("");

  Double_t whist[900]={};
  
  for(int ib=0;ib<900;ib++){
    whist[ib] = hthe->GetBinContent(ib+1)/180000.;
    //cout << ib << " " << whist[ib] << endl;
  }

  ///////////////////////////////////////////
  

  
  
  Double_t sum_flux=0;
  Double_t tmp_B_rate=0;
  Int_t iline=0;
  Int_t resto=-1;
  Double_t B_en[818]={};
  Double_t B_rate[818]={};
  
  TString dir = gSystem->UnixPathName(__FILE__);
  dir.ReplaceAll("Sun_B_newversion.c","");
  dir.ReplaceAll("/./","/");
  ifstream in;
  in.open(Form("%sboron_spectrum.txt",dir.Data()));

  while (1) {
    in >> B_en[iline] >> B_rate[iline];
    if (!in.good()) break;
    // cout << B_en[iline] << " " << B_rate[iline];
    iline++;
  }

  ofstream log_file;
  log_file.open ("num_eventi.txt");
  log_file << "Sun Boro induced recoil rate";
  log_file << endl;
  
  Functions* functions;

  functions = new Functions();


  TH1F* h_norm_cos_theta = new TH1F("norm_theta","nt",200,-1,1);
  for(int k=0;k<d_theta;k++){
    h_norm_cos_theta->Fill(-1 + 1./d_theta + k*(2./d_theta),1./d_theta + k*(2./d_theta));
  }
  TCanvas *cnt = new TCanvas("cnt","cnt",600,600);
  h_norm_cos_theta->Draw("");
  h_norm_cos_theta->Scale(1./h_norm_cos_theta->Integral());
  

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
    tmp_B_rate=0;
    sum_flux=0;
    //cout << iel << " " << iZ_el << " " << iA_el << " " << ifract_el << " " << imass_fract << " " << iM_el << " " << iN_el << endl;
    
    
    for(int j=5; j<818;j++){
      resto = j % 5;
      if(resto==0){
	E_nu=(j/50.);
	sum_crsect=0;
	//iflux_boron_vs_E=B_rate[j]*integral_flux_boron*(1./10.);
	iflux_boron_vs_E=integral_flux_boron*((TMath::Abs(B_rate[j]+tmp_B_rate))*(1./10.)/2.);
	sum_flux += iflux_boron_vs_E;
	//cout << j << " " << E_nu  << " " << B_rate[j] << " " << tmp_B_rate <<  " " << sum_flux << endl;

      
	for(int k=0;k<d_theta;k++){
	cos_theta = h_norm_cos_theta->GetRandom();	
	//cos_theta = -1 + 1./d_theta + k*(2./d_theta); // moltiplicato per delta theta (vedi crsect in Functions.h)		
	E_rec = (E_nu*E_nu*(1-cos_theta))/(E_nu*(1-cos_theta)+iM_el*TMath::Power(10,3)); // MeV
	crsect = functions->cross_section();
	sum_crsect+=crsect;
	//cout << E_rec << endl;
	/// XENON CROSS-CHECK
	//if(E_rec>0.002){ 
	 // num_ev_bkg = functions->num_eventi_bkg();
	 // sum_ev_bkg = sum_ev_bkg + num_ev_bkg;
	  //}
	/////
	
	L_true = GetRangeFromEnergy(iZ_el,E_rec*1000);
	//L_proj = L_true*cos_theta_rec;
	//cout << j << " " << k << " " << crsect << " " << dN_vs_dE_bkg << " " << num_ev_bkg << " " << sum_ev_bkg << endl;
	
	if(L_true > 100 && L_true<1000){
	  //num_ev_bkg = functions->num_eventi_bkg()*whist[k]*(d_theta/2.);
	  num_ev_bkg = functions->num_eventi_bkg()*whist[k]*(d_theta/2.);
	  sum_ev_bkg = sum_ev_bkg + num_ev_bkg;
	  h_len->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]*whist[k]);
	}
      }
      tmp_B_rate=B_rate[j];
    }
      //if(E_nu==10)cout << E_nu << " " << sum_crsect << " " << int_const_bkg << endl;
    }
    cout << iel << " " << sum_ev_bkg << endl;
  }
  num_ev_bkg=sum_ev_bkg*86400*365;
  cout << num_ev_bkg << endl;
  h_len->Draw("");
  h_len->Scale(num_ev_bkg/h_len->Integral());
  delete functions;

}
