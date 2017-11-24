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

#include "Definitions_boro.h"
#include "Functions.h"
#include "EnergyRangeCorrelation.C"
#include "TRandom3.h"
#include <TSystem.h>



using namespace std;

void Sun_B_newfluence() {


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
  TH1F *hphiB = new TH1F("hphiB","hphiB",100,-TMath::Pi(),TMath::Pi());
  //TH1F *htheB = new TH1F("htheB","htheB",100,-TMath::Pi()/2.,TMath::Pi()/2.);
  TH1F *htheB = new TH1F("htheB","htheB",200,0,1);

  TRandom3 *r = new TRandom3();
  Double_t t = 0;
  Double_t px = 0;
  Double_t py = 0;
  Double_t pz = 0;
  Double_t pphi = 0;
  Double_t pthe = 0;
  for(int i=0;i<180000;i++){
    t = r->Uniform(-TMath::Pi(),TMath::Pi());
    px = TMath::Cos(bx)*TMath::Sin(t-lx+ly);
    py = TMath::Cos(by)*TMath::Sin(t-ly+ly);
    pz = TMath::Cos(bz)*TMath::Sin(t-lz+ly);
    pphi = TMath::ATan2(px,py);
    pthe = TMath::ATan(pz/TMath::Sqrt(px*px + py*py));
    hr->Fill(t);
    hphiB->Fill(pphi);
    htheB->Fill(TMath::Cos(pthe));
    //htheB->Fill((pthe));
  }

  /////////////////////////////////////////////////////////////////


    
  
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
  hphiB->Draw("");
  cc2->cd(2);
  htheB->Draw("");

  Double_t whist[900]={};
  
  for(int ib=0;ib<900;ib++){
    whist[ib] = htheB->GetBinContent(ib+1)/180000.;
    //cout << ib << " " << whist[ib] << endl;
  }

  ///////////////////////////////////////////



  TRandom3 *tt = new TRandom3();
  Double_t sum_flux=0;
  Double_t tmp_B_rate=0;
  Int_t iline=0;
  Int_t resto=-1;
  Double_t B_en[818]={};
  Double_t B_rate[818]={};
  
  TString dir = gSystem->UnixPathName(__FILE__);
  dir.ReplaceAll("Sun_B_newfluence.c","");
  dir.ReplaceAll("/./","/");
  ifstream in;
  in.open(Form("%sboron_spectrum.txt",dir.Data()));

  while (1) {
    in >> B_en[iline] >> B_rate[iline];
    if (!in.good()) break;
    //cout << B_en[iline] << " " << B_rate[iline] << endl;
    iline++;
  }

  //////////////////////////////////////////////////////


  ofstream log_file;
  log_file.open ("num_eventi.txt");
  log_file << "Sun Boro induced recoil rate";
  log_file << endl;
  
  
  Functions* functions;

  functions = new Functions();
  
  gStyle->SetOptStat("eMR");
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.15);

  gr_spect[0]=new TGraph();
  gr_spect[1]=new TGraph();
  gr_spect[2]=new TGraph();

  h4->Sumw2();
  h44->Sumw2();
  h2->Sumw2();
  hthe->Sumw2();

  //// COSTHETA  //////////////////////////////////////////////////////
  TH1F* h_norm_cos_theta = new TH1F("norm_theta","nt",200,-1,1);
  for(int k=0;k<d_theta;k++){
    h_norm_cos_theta->Fill(-1 + 1./d_theta + k*(2./d_theta),1./d_theta + k*(2./d_theta));
  }
  TCanvas *cnt = new TCanvas("cnt","cnt",600,600);
  h_norm_cos_theta->Draw("");
  h_norm_cos_theta->Scale(1./h_norm_cos_theta->Integral());
  /////////////////////////////////////////////////////////////////

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
    
    for(int j=5; j<818;j++){
      resto = j % 5;
      if(resto==0){
	E_nu=(j/50.);
	sum_crsect=0;
	iflux_boron_vs_E=integral_flux_boron*((TMath::Abs(B_rate[j]+tmp_B_rate))*(1./10.)/2.);
	sum_flux += iflux_boron_vs_E;
	//cout << j << " " << E_nu  << " " << B_rate[j] << " " << tmp_B_rate <<  " " << sum_flux << endl;
	
	h1->Reset(0);
	
	for(int k=0;k<d_theta;k++){
	  
	  ////// CROSS SECTION  //////////////////////////////
	  //cos_theta = h_norm_cos_theta->GetRandom();
	  cos_theta = -1 + 1./d_theta + k*(2./d_theta); // devi moltiplicare per delta theta
	  crsect = functions->cross_section();
	  h1->Fill(cos_theta,crsect);
	  ///////////////////////////////////	
	}
	
	h1->Scale(1./h1->Integral());
	
	for(int k=0;k<d_theta;k++){
	  cos_theta = -1 + 1./d_theta + k*(2./d_theta); // devi moltiplicare per delta theta
	  //cos_theta_norm = h1->GetRandom();

	  ///////////////////////////////////////////////////
	  t = r->Uniform(-TMath::Pi(),TMath::Pi());
	  px = TMath::Cos(bx)*TMath::Sin(t-lx+ly);
	  py = TMath::Cos(by)*TMath::Sin(t-ly+ly);
	  pz = TMath::Cos(bz)*TMath::Sin(t-lz+ly);
	  pphi = TMath::ATan2(px,py);
	  pthe = TMath::ATan(pz/TMath::Sqrt(px*px + py*py));
	  ////////////////////////////////////////////////////
	  h0->Fill(t);
	  
	  cos_theta_norm = cos_theta;
	  h11->Fill(cos_theta_norm);
	  
	  crsect = functions->cross_section();
	  sum_crsect = crsect + sum_crsect;
	  
	  E_rec = (E_nu*E_nu*(1-cos_theta_norm))/(E_nu*(1-cos_theta_norm)+iM_el*TMath::Power(10,3)); // MeV
	  cos_theta_rec = ((E_nu + iM_el*TMath::Power(10,3))/(E_nu))*TMath::Sqrt(E_rec/(2*iM_el*TMath::Power(10,3)));
	  theta_rec = TMath::ACos(cos_theta_rec)-pthe;
	  if(theta_rec>pi/2.)theta_rec = theta_rec -pi/2.;
	  if(theta_rec<-pi/2.)theta_rec = pi/2. - theta_rec;
	  theta_nu = TMath::ASin((TMath::Sqrt(2*(iM_el/1000.)*E_rec)*TMath::Sin(theta_rec))/(E_nu-TMath::Sqrt(2*(iM_el/1000.)*E_rec)));
	  L_true = GetRangeFromEnergy(iZ_el,E_rec*1000);
	  L_proj = L_true*cos_theta_rec;

	  //if(crsect>0) cout << E_rec << " " << theta_rec << " " << L_true << endl;

	  ////////////////////  BOOST LORENTZ  /////////////////////////////////////////////////////////////
	  E_nu_cm = E_nu*TMath::Sqrt((iM_el*1000)/(iM_el*1000+2*E_nu)); //MeV    
	  E_rec_cm = TMath::Sqrt(E_nu_cm*E_nu_cm+iM_el*iM_el*1000000);     //MeV
	  
	  rnd_cos_theta_cm = tt->Uniform(-1,1);
	  rnd_phi_cm = tt->Uniform(0,2*TMath::Pi());
	  
	  px_cm = E_nu_cm*TMath::Sqrt((1-TMath::Power(rnd_cos_theta_cm,2)))*TMath::Sin(rnd_phi_cm);
	  py_cm = E_nu_cm*TMath::Sqrt((1-TMath::Power(rnd_cos_theta_cm,2)))*TMath::Cos(rnd_phi_cm);
	  pz_cm = E_nu_cm*rnd_cos_theta_cm;
	  p_mod_cm = TMath::Sqrt(px_cm*px_cm + py_cm*py_cm + pz_cm*pz_cm);
	  
	  beta_cm = E_nu/(E_nu + iM_el*1000);
	  gamma_cm = (E_nu + iM_el*1000)/(TMath::Sqrt(iM_el*iM_el*1000000 + 2*E_nu*iM_el*1000));
	  py_lab = beta_cm*gamma_cm*E_rec_cm + gamma_cm*py_cm;
	  E_tot_lab = gamma_cm*E_rec_cm + beta_cm*gamma_cm*py_cm;  //MeV
	  E_rec_lab = gamma_cm*E_rec_cm + beta_cm*gamma_cm*py_cm - iM_el*1000; //MeV
	  p_mod_lab = TMath::Sqrt(E_tot_lab*E_tot_lab - iM_el*iM_el*1000000);
	  the_rec = TMath::ATan2(TMath::Sqrt((1-TMath::Power(rnd_cos_theta_cm,2))),(gamma_cm*(beta_cm*(E_rec_cm/p_mod_cm)+rnd_cos_theta_cm)))-pthe;  //buono
	  phi_rec = TMath::ATan2(py_lab,px_cm)-pphi;
	  //if(phi_rec)
	  
	  hpx_cm->Fill(px_cm);
	  hpy_cm->Fill(py_cm);
	  hpz_cm->Fill(pz_cm);
	  hpy_lab->Fill(py_lab);
	  ///////////////////////////////////////////////////////////////////////////////////////////////
	 
	  	  
	  
	  //hphi->Fill(phi_rec,atomic_fract_el[i]*A_el[i]*A_el[i]*crsect);
	  //hthe->Fill(the_rec,atomic_fract_el[i]*A_el[i]*A_el[i]*crsect);
	  
	  hphi->Fill(phi_rec,atomic_fract_el[i]*A_el[i]*A_el[i]);
	  hthe->Fill(the_rec,atomic_fract_el[i]*A_el[i]*A_el[i]);
	  h2->Fill(theta_rec,atomic_fract_el[i]*A_el[i]*A_el[i]);
	  
	  h4->Fill(E_rec*1000,atomic_fract_el[i]*A_el[i]*A_el[i]);
	  h44->Fill(E_rec_lab*1000,atomic_fract_el[i]*A_el[i]*A_el[i]);
	  
	  if(/*L_true>=100 && L_true<=1000*/ L_true>0){
	    
	    num_ev_bkg = functions->num_eventi_bkg()*whist[k]*(d_theta/2.);
	    sum_ev_bkg = sum_ev_bkg + num_ev_bkg;

	    //cout << num_ev_bkg << " " << sum_ev_bkg << " " << whist[k] << endl;
	    
	    if(i==0)hh0->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);
	    if(i==1)hh1->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);
	    if(i==2)hh2->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);
	    if(i==3)hh3->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);
	    if(i==4)hh4->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);
	    if(i==5)hh5->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);
	    if(i==6)hh6->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);
	    
	    hcosthe_thr->Fill(TMath::ATan(theta_rec),atomic_fract_el[i]*A_el[i]*A_el[i]);
	    hthe_thr->Fill(theta_rec,atomic_fract_el[i]*A_el[i]*A_el[i]);
	    hphi_thr->Fill(phi_rec,atomic_fract_el[i]*A_el[i]*A_el[i]);
	    hlen_thr->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);
	    hene_thr->Fill(E_rec*1000,atomic_fract_el[i]*A_el[i]*A_el[i]);

	    // cout << L_true/1000. << " " << theta_rec << " " << phi_rec << " " << E_rec*1000 << " " << atomic_fract_el[i] << " " << A_el[i] << endl;
	  }
	}
	tmp_B_rate=B_rate[j];
      }
    }
    cout << iel << " " << sum_ev_bkg << endl;
  }
  num_ev_bkg=sum_ev_bkg*86400*365;
  cout << num_ev_bkg << endl;

  TCanvas *c0 = new TCanvas("c0","c0",600,600);
  h0->Draw("");
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  h1->Draw("");
  TCanvas *c11 = new TCanvas("c11","c11",600,600);
  h11->Draw("");

  /*
  TCanvas *cc4 = new TCanvas("cc4","cc4",600,600);
  hthe->Draw("");
  h2->Draw("same");
  hthe->SetLineColor(1);
  h2->SetLineColor(2);
  
  Double_t ks0 = hthe->KolmogorovTest(h2);
  cout << "KS_theta " << ks0 << endl;
  
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  h2->Draw("");
  
  
  TCanvas *c4 = new TCanvas("c4","c4",600,600);
  h4->Draw("");
  h44->Draw("same");
  h44->SetLineColor(1);
  h4->SetLineColor(2);

  Double_t ks = h4->KolmogorovTest(h44);
  cout << "KS_ene " << ks << endl;

  */
  
  TCanvas *clen = new TCanvas("clen","clen",600,600);
  hlen_thr->Draw("");
  hlen_thr->Scale(num_ev_bkg/hlen_thr->Integral());

  TCanvas *cene = new TCanvas("cene","cene",600,600);
  hene_thr->Draw("");
  hene_thr->Scale(num_ev_bkg/hene_thr->Integral());

  TCanvas *cphi = new TCanvas("cphi","cphi",600,600);
  hphi_thr->Draw("");
  hphi_thr->Scale(num_ev_bkg/hphi_thr->Integral());

  TCanvas *cthe = new TCanvas("cthe","cthe",600,600);
  hthe_thr->Draw("");
  hthe_thr->Scale(num_ev_bkg/hthe_thr->Integral());

  TCanvas *ccosthe = new TCanvas("ccosthe","ccosthe",600,600);
  hcosthe_thr->Draw("");
  hcosthe_thr->Scale(num_ev_bkg/hcosthe_thr->Integral());
 

  
  TCanvas *c15 = new TCanvas("c15","c15",1200,1200);
  c15->Divide(2,2);
  c15->cd(1);
  hpx_cm->Draw("");
  c15->cd(2);
  hpy_cm->Draw("");
  c15->cd(3);
  hpz_cm->Draw("");
  c15->cd(4);
  hpy_lab->Draw("");
  
  
  TCanvas *cc3 = new TCanvas("cc3","cc3",600,600);
  hphi->Draw("");
  hphi->Scale(num_ev/hphi->Integral());
  //hphi->Fit("fphi");
  
  
  delete functions;
}
