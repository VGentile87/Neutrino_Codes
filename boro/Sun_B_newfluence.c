//  AUTHOR: V. GENTILE  2017 // 

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
  //TH1F *hphiB = new TH1F("hphiB","hphiB",100,0,2*TMath::Pi());
  //TH1F *htheB = new TH1F("htheB","htheB",100,0,2*TMath::Pi());
  TH1F *hphiB = new TH1F("hphiB","hphiB",100,-pi,pi);
    
  
  TH2F *hphitheB = new TH2F("hphitheB","hphitheB",128,-pi,pi,128,-pi/2.,pi/2.);
  //TH2F *hphitheB = new TH2F("hphitheB","hphitheB",100,-pi,pi,100,0,pi);
  TH1F *htheB = new TH1F("htheB","htheB",100,-pi/2.,pi/2.);
  TH1F *hcostheB = new TH1F("hcostheB","hcostheB",900,-1,1);

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
    //pthe = TMath::ATan(TMath::Sqrt(px*px + py*py)/pz);
    // if(pthe<0)pthe=pthe+pi;
    hr->Fill(t);
    hphiB->Fill(pphi);
    hcostheB->Fill(TMath::Cos(pthe));
    htheB->Fill((pthe));
    hphitheB->Fill(pphi,pthe);
  }

  /*
  TCanvas *cc1 = new TCanvas("cc1","cc1",1200,600);
  cc1->Divide(2,1);
  cc1->cd(1);
  fphi->Draw("P");
  cc1->cd(2);
  fthe->Draw("P");
  */


  

  TCanvas *cd = new TCanvas("cd","cd",600,600);
  hphiB->Draw("");
  TCanvas *ce = new TCanvas("ce","ce",600,600);
  htheB->Draw("");
  TCanvas *cf = new TCanvas("cf","cf",600,600);
  hphitheB->Draw("");
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
  h22->Sumw2();

  //// COSTHETA  //////////////////////////////////////////////////////
  TH1F* h_norm_cos_theta = new TH1F("norm_theta","nt",200,-1,1);
  for(int k=0;k<d_theta;k++){
    h_norm_cos_theta->Fill(-1 + 1./d_theta + k*(2./d_theta),1./d_theta + k*(2./d_theta));
  }
  //TCanvas *cnt = new TCanvas("cnt","cnt",600,600);
  //h_norm_cos_theta->Draw("");
  //h_norm_cos_theta->Scale(1./h_norm_cos_theta->Integral());
  /////////////////////////////////////////////////////////////////

  for(int i=0; i<n_el;i++){
    h1[i] = new TH1D(Form("cos_el_%d",i),"",100,-1,1);
    h11[i] = new TH1D(Form("phi_el_%d",i),"",100,0,2*pi);
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
	
	//h1->Reset(0);
	
	for(int k=0;k<d_theta;k++){	  
	  //////  CROSS SECTION  /////////////////////////////
	  cos_theta = -1 + 1./d_theta + k*(2./d_theta);     // dtheta
	  crsect = functions->cross_section();              // sez. d'urto 
	  sum_crsect = crsect + sum_crsect;	            // somma sez. d'urto
	  E_rec = (E_nu*E_nu*(1-cos_theta))/(E_nu*(1-cos_theta)+iM_el*TMath::Power(10,3)); // energia nucluear recoil (MeV)
	  costheta_sc = ((E_nu + iM_el*TMath::Power(10,3))/(E_nu))*TMath::Sqrt(E_rec/(2*iM_el*TMath::Power(10,3)));  // angolo nuclear recoil rispetto neutrino incidente
	  L_true = GetRangeFromEnergy(iZ_el,E_rec*1000);    // lunghezza traccia rinculo
	  L_proj = L_true*cos_theta_rec;                    // lunghezza proiettata sul piano di emulsione
	  ////////////////////////////////////////////////////

	  //cout << k << " " << cos_theta << " " << E_rec << " " << L_true << " " << whist[k] << endl; 

	  ////////////////////  BOOST LORENTZ  /////////////////////////////////////////////////////////////
	  E_nu_cm = E_nu*TMath::Sqrt((iM_el*1000)/(iM_el*1000+2*E_nu)); //MeV    
	  E_rec_cm = TMath::Sqrt(E_nu_cm*E_nu_cm+iM_el*iM_el*1000000);  //MeV
	  
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
	  the_rec = TMath::ATan2(TMath::Sqrt((1-TMath::Power(rnd_cos_theta_cm,2))),(gamma_cm*(beta_cm*(E_rec_cm/p_mod_cm)+rnd_cos_theta_cm)));  //buono
	  phi_rec = TMath::ATan2(py_lab,px_cm);
	  
	  hpx_cm->Fill(px_cm);
	  hpy_cm->Fill(py_cm);
	  hpz_cm->Fill(pz_cm);
	  hpy_lab->Fill(py_lab);
	  ///////////////////////////////////////////////////////////////////////////////////////////////

	  h2->Fill(TMath::ACos(costheta_sc));
	  h22->Fill(the_rec);
	  h4->Fill(E_rec*1000);
	  h44->Fill(E_rec_lab*1000);

	  if(L_true>=50 && L_true<=1000){
	    num_ev_bkg = functions->num_eventi_bkg();
	    sum_ev_bkg = sum_ev_bkg + num_ev_bkg;
	    h1[i]->Fill(costheta_sc);
	    h11[i]->Fill(phi_rec);
	    hthe_thr->Fill(the_rec);
	    hcosthe_thr->Fill(TMath::Cos(the_rec));
	    hphi_thr->Fill(phi_rec);
	    hlen_thr->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);
	    hene_thr->Fill(E_rec*1000);
	  }
	  ///////////////////////////////////	
	}	
	tmp_B_rate=B_rate[j];
      }
    }
    cout << iel << " " << sum_ev_bkg << endl;


    
    for(int k=0;k<20000;k++){
      if(h1[i]->GetEntries()!=0){
	costheta_sc = h1[i]->GetRandom();
        phi_sc = h11[i]->GetRandom();
	h00->Fill(costheta_sc,atomic_fract_el[i]*A_el[i]*A_el[i]);
	theta_sc = TMath::ACos(costheta_sc);
	//hthe_thr->Fill(theta_sc,atomic_fract_el[i]*A_el[i]*A_el[i]);
	///////////////////////////////////////////////////
	t = r->Uniform(0,2*TMath::Pi());
	px = TMath::Cos(bx)*TMath::Sin(t-lx+ly);  // centro galattico
	py = TMath::Cos(by)*TMath::Sin(t-ly+ly);  // cigno
	pz = TMath::Cos(bz)*TMath::Sin(t-lz+ly);  // normale a centro galattico
	pphi = TMath::ATan2(px,py);
	pthe = TMath::ATan(pz/TMath::Sqrt(px*px + py*py));
	if(pthe<0)pthe = pthe + 2*pi;
	if(pphi<0)pphi = pphi + 2*pi;
	////////////////////////////////////////////////////
	//px_sum = TMath::Cos(phi_sc)*TMath::Cos(theta_sc) + TMath::Cos(pphi)*TMath::Cos(pthe);
	//py_sum = TMath::Sin(phi_sc)*TMath::Cos(theta_sc) + TMath::Sin(pphi)*TMath::Cos(pthe);
	//pz_sum = TMath::Sin(theta_sc) + TMath::Sin(pthe);
	//theta_sum  = TMath::ATan(pz_sum/TMath::Sqrt(px_sum*px_sum + py_sum*py_sum));
	//phi_sum = TMath::ATan2(py_sum,px_sum);
	/////////////////////////////////////////////////////
	theta_sum = theta_sc +  pthe;
	phi_sum = phi_sc + pphi;
	if(theta_sum>2*pi)theta_sum = theta_sum - 2*pi;
	if(theta_sum>pi)theta_sum = theta_sum - 2*pi;
	if(phi_sum>2*pi)phi_sum = phi_sum -2*pi;
	if(phi_sum>pi)phi_sum = phi_sum -2*pi;
	/////////////////////////////////////////////////////
	h0->Fill(pthe,atomic_fract_el[i]*A_el[i]*A_el[i]);
	cos_theta_sum = TMath::Cos(theta_sum); 
	h3->Fill(theta_sum,atomic_fract_el[i]*A_el[i]*A_el[i]);
	h33->Fill(cos_theta_sum,atomic_fract_el[i]*A_el[i]*A_el[i]);
	h5->Fill(phi_sum,atomic_fract_el[i]*A_el[i]*A_el[i]);
	h55->Fill(phi_sc,atomic_fract_el[i]*A_el[i]*A_el[i]);

	hphiB->Fill(pphi);
	htheB->Fill(pthe);
      }
    }
        
  }
  num_ev_bkg=sum_ev_bkg*86400*365;
  cout << num_ev_bkg << endl;
  
  TCanvas *cc2 = new TCanvas("cc2","cc2",1200,600);
  //hthe->Draw();
  cc2->Divide(2,1);
  cc2->cd(1);
  //7hr->Draw("");
  //cc2->cd(2);
  hphiB->Draw("");
  cc2->cd(2);
  htheB->Draw("");

  TCanvas *cthe3 = new TCanvas("cthe3","cthe3",600,600);
  h3->Draw("");
  
  
  /*
  TCanvas *c0 = new TCanvas("c0","c0",1200,600);
  c0->Divide(2,1);
  c0->cd(1);
  h0->Draw("");
  c0->cd(2);
  h00->Draw("");
  */

  /*
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(4,2);
  for(int s=0;s<7;s++){
    c1->cd(s+1);
    h1[s]->Draw("");
  }

  TCanvas *c11 = new TCanvas("c11","c11",1200,600);
  c11->Divide(4,2);
  for(int s=0;s<7;s++){
    c11->cd(s+1);
    h11[s]->Draw("");
  }
  */
  
  TCanvas *c3 = new TCanvas("c3","c3",1200,600);
  c3->Divide(2,1);
  c3->cd(1);
  h3->Draw("");
  c3->cd(2);
  h33->Draw("");

  TCanvas *c5 = new TCanvas("c5","c5",1200,600);
  c5->Divide(2,1);
  c5->cd(1);
  h5->Draw("");
  c5->cd(2);
  h55->Draw("");
  
  
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  h2->Draw("");
  h22->Draw("sames");
  h2->SetLineColor(1);
  h22->SetLineColor(2);
  
  Double_t ks0 = h2->KolmogorovTest(h22);
  cout << "KS_theta " << ks0 << endl;
  
  
  
  TCanvas *c4 = new TCanvas("c4","c4",600,600);
  h4->Draw("");
  h44->Draw("sames");
  h4->SetLineColor(1);
  h44->SetLineColor(2);

  Double_t ks = h4->KolmogorovTest(h44);
  cout << "KS_ene " << ks << endl;

  
  
  TCanvas *clen = new TCanvas("clen","clen",600,600);
  hlen_thr->Draw("");
  //hlen_thr->Scale(num_ev_bkg/hlen_thr->Integral());
  /*
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
*/ 

  /*
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
  */
  
  /*
  TH2F *skymap=new TH2F("Blank skymap","Blank skymap",40,-200,200,20,-100,100);
  gStyle->SetOptStat(0000);	
  float conv=3.1415926/180; // I am also aware of TMath::DegToRad() and TMath::Pi() which could be used there...
  const int Nl=19; // Number of drawn latitudes
  const int NL=19; // Number of drawn longitudes
  int M=30;
  TGraph *latitudes[Nl];
  TGraph *longitudes[NL];
  for(int j=0;j<Nl;++j)
    {
      latitudes[j]=new TGraph();
      float la=-90+180/(Nl-1)*j;
      
      for(int i=0;i<M+1;++i)
	{
	  float lo= -180+360/M*i;
	  float z = sqrt(1+cos(la*conv)*cos(lo*conv/2));
	  float x = 180*cos(la*conv)*sin(lo*conv/2)/z;
	  float y = 90*sin(la*conv)/z;
	  latitudes[j]->SetPoint(i,x,y);
	}
    }  
  
  for(int j=0;j<NL;++j)
    {
      longitudes[j]=new TGraph();
      float lo=-180+360/(NL-1)*j;
      
      for(int i=0;i<M+1;++i)
	{
	  float la= -90+180/M*i;
	  float z = sqrt(1+cos(la*conv)*cos(lo*conv/2));
	  float x = 180*cos(la*conv)*sin(lo*conv/2)/z;
	  float y = 90*sin(la*conv)/z;
	  longitudes[j]->SetPoint(i,x,y);
	}
    }
  skymap->Draw();
  for(int j=0;j<Nl;++j) latitudes[j]->Draw("c");
  for(int j=0;j<NL;++j) longitudes[j]->Draw("c");
  */
  
  
  // delete functions;
}
