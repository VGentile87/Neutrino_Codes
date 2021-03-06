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

#include "Definitions_News.h"
#include "Functions_bkg_ev.h"
#include "EnergyRangeCorrelation.C"
#include "TRandom3.h"
#include <TSystem.h>



using namespace std;

void Sun_nu() {


 Int_t resto=0;
 Int_t num_ev_bkg=0;
 Double_t B_int_flux = 5.58*TMath::Power(10,6);
 Double_t atm_int_flux = 10.5;
 Double_t dsnb_int_flux = 85.5;
 //Double_t pp_flux = 5.98*TMath::Power(10,10);
 //Double_t Be7_flux = 5*TMath::Power(10,9);
 //Double_t pep_flux = 1.44*TMath::Power(10,8);
 //Double_t O15_flux = 2.23*TMath::Power(10,8);
 Double_t hep_int_flux = 8.04*TMath::Power(10,3);
 Int_t dsnb_iline=0;
 Int_t atm_iline=0;
 Int_t B_iline=0;
 Int_t hep_iline=0;
 Double_t dsnb_nu_en[62]={};
 Double_t dsnb_nu_rate[62]={};
 Double_t atm_nu_en[48]={};
 Double_t atm_nu_rate[48]={};
 Double_t B_nu_en[819]={};
 Double_t B_nu_rate[819]={};
 Double_t hep_nu_en[1000]={};
 Double_t hep_nu_rate[1000]={};

 Double_t sum_flux=0;
  
  TString dir = gSystem->UnixPathName(__FILE__);
  dir.ReplaceAll("Sun_nu.c","");
  dir.ReplaceAll("/./","/");
  ifstream in;

  in.open(Form("%sdsnb_spectrum3.txt",dir.Data()));
  while (1) {
    in >> dsnb_nu_en[dsnb_iline] >> dsnb_nu_rate[dsnb_iline];
    if (!in.good()) break;
    dsnb_iline++;
  }

  Double_t xpoint[196]={};
  Double_t tot_rate[196]={};

  for(int i=0;i<196;i++){
    if(i<10)xpoint[i]=(i+1)/10.;
    if(i>=10 && i<110)xpoint[i]=(i+1)-10;
    if(i>=110)xpoint[i]=xpoint[i-1]+10;
    //cout << xpoint[i] << " " << tot_rate[i] << endl;
  }
  
  TGraph *dsnb_gr = new TGraph(dsnb_iline,dsnb_nu_en,dsnb_nu_rate);
  for(int i=0;i<100;i++){
    if(i<10)xpoint[i]=(i+1)/10.;
    if(i>=10)xpoint[i]=(i+1)-10;
    tot_rate[i]=tot_rate[i]+dsnb_gr->Eval(xpoint[i]);
    //cout << xpoint[i] << " " << dsnb_gr->Eval(xpoint[i]) << " " << tot_rate[i] << endl; 
  }

  in.close();

  in.open(Form("%satm_spectrum3.txt",dir.Data()));
  while (1) {
    in >> atm_nu_en[atm_iline] >> atm_nu_rate[atm_iline];
    //cout << atm_nu_en[atm_iline] << " " <<  atm_nu_rate[atm_iline] << endl;
    atm_nu_rate[atm_iline]=atm_nu_rate[atm_iline];//*atm_int_flux;
    if (!in.good()) break;
    atm_iline++;
  }

  TGraph *atm_gr = new TGraph(atm_iline,atm_nu_en,atm_nu_rate);
  for(int i=0;i<196;i++){
   
    tot_rate[i]=tot_rate[i]+atm_gr->Eval(xpoint[i]);
    //cout << xpoint[i] << " " << atm_gr->Eval(xpoint[i]) << " " << tot_rate[i] << endl; 
  }

  in.close();

  in.open(Form("%sboron_spectrum.txt",dir.Data()));
  while (1) {
    in >> B_nu_en[B_iline] >> B_nu_rate[B_iline];
    B_nu_rate[B_iline]=B_nu_rate[B_iline]*B_int_flux;
    if (!in.good()) break;
    B_iline++;
  }

  TGraph *B_gr = new TGraph(B_iline,B_nu_en,B_nu_rate);
  for(int i=0;i<26;i++){
    if(i<10)xpoint[i]=(i+1)/10.;
    if(i>=10)xpoint[i]=(i+1)-10;
    tot_rate[i]=tot_rate[i]+B_gr->Eval(xpoint[i]);
    cout << xpoint[i] << " " << B_gr->Eval(xpoint[i]) << " " << tot_rate[i] << endl; 
  }

  in.close();

  in.open(Form("%shep_spectrum2.txt",dir.Data()));
  while (1) {
    in >> hep_nu_en[hep_iline] >> hep_nu_rate[hep_iline];
    hep_nu_rate[hep_iline]=hep_nu_rate[hep_iline]*hep_int_flux;
    if (!in.good()) break;
    hep_iline++;
  }

  TGraph *hep_gr = new TGraph(hep_iline,hep_nu_en,hep_nu_rate);
  for(int i=0;i<28;i++){
    if(i<10)xpoint[i]=(i+1)/10.;
    if(i>=10)xpoint[i]=(i+1)-10;
    tot_rate[i]=tot_rate[i]+hep_gr->Eval(xpoint[i]);
    //cout << xpoint[i] << " " << hep_gr->Eval(xpoint[i]) << " " << tot_rate[i] << endl; 
  }

  in.close();

  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  TGraph *tot_rate_gr = new TGraph(196,xpoint,tot_rate);
  tot_rate_gr->Draw("AL");

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(dsnb_gr);
  mg->Add(atm_gr);
  mg->Add(B_gr);
  mg->Add(hep_gr);
  dsnb_gr->SetLineWidth(3);
  dsnb_gr->SetLineColor(kOrange+7);
  atm_gr->SetLineWidth(3);
  atm_gr->SetLineColor(kViolet);
  B_gr->SetLineWidth(3);
  B_gr->SetLineColor(kRed+1);
  hep_gr->SetLineWidth(3);
  hep_gr->SetLineColor(kBlue+1);
  mg->Draw("AL");
  c1->SetLogy();
  c1->SetLogx();
  mg->GetXaxis()->SetTitle("Neutrino Energy [MeV]");
  mg->GetYaxis()->SetTitle("Neutrino Flux [cm^{-2} s^{-1} MeV^{-1}]");
  TLegend *legend= new TLegend(0.6,0.65,0.88,0.85);
  //legend->SetTextFont(72);
  legend->SetTextSize(0.04);
  legend->AddEntry(dsnb_gr,"DSNB","l");
  legend->AddEntry(B_gr,"^{8}B","l");
  legend->AddEntry(atm_gr,"Atm","l");
  legend->AddEntry(hep_gr,"hep","l");
  legend->Draw();



  ofstream log_file;
  log_file.open ("num_eventi.txt");
  log_file << "Sun  induced recoil rate";
  log_file << endl;
  



  Functions* functions;

  functions = new Functions();
  

  for(int i=0; i<n_el;i++){
    
    iel = el[i];
    iZ_el = Z_el[i];
    iA_el = A_el[i];
    ifract_el = atomic_fract_el[i]; // uma
    imass_fract = mass_fract[i];
    iM_el=iA_el*uma;
    iN_el=iA_el-iZ_el;
    num_ev_bkg=0;
    sum_flux=0;
    cout << iel << " " << iZ_el << " " << iA_el << " " << ifract_el << " " << imass_fract << " " << iM_el << " " << iN_el << endl;
    
    
    for(int j=0; j<196;j++){
      sum_crsect=0;
      E_nu=xpoint[j];
      if(j<10)iflux_boron_vs_E=tot_rate[j]*0.1;
      if(j>=10 && j<110)iflux_boron_vs_E=tot_rate[j];
      if(j>=110)iflux_boron_vs_E=tot_rate[j]*10;
      sum_flux += iflux_boron_vs_E;
      //cout << j << " " << E_nu  << " " << iflux_boron_vs_E <<  " " << sum_flux << endl;
      
      
      for(int k=0;k<d_theta;k++){	
	cos_theta = -1 + 1./d_theta + k*(2./d_theta); // devi moltiplicare per delta theta		
	E_rec = (E_nu*E_nu*(1-cos_theta))/(E_nu*(1-cos_theta)+iM_el*TMath::Power(10,3)); // MeV
	crsect = functions->cross_section();
	sum_crsect+=crsect;	
	L_true = GetRangeFromEnergy(iZ_el,E_rec*1000);
	//L_proj = L_true*cos_theta_rec;
	//cout << j << " " << k << " " << crsect << " " << dN_vs_dE_bkg << " " << num_ev_bkg << " " << sum_ev_bkg << endl;
	//cout << L_true << endl;
	
	if(L_true >= 300 && L_true<1000 /*&& E_rec>0.002*/){
	  //num_ev_bkg = functions->num_eventi_bkg();
	  //sum_ev_bkg = sum_ev_bkg + num_ev_bkg;
	  sum_ev_bkg += functions->num_eventi_bkg();
	  h_len->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);
	  //cout << j << " " << k << " " << crsect << " " << dN_vs_dE_bkg << " " << sum_ev_bkg << endl;
	}
      }
    }
    //cout << iel << " " << sum_ev_bkg << endl;
  }
  sum_ev_bkg=sum_ev_bkg*86400*365;
  cout << sum_ev_bkg << endl;
  //h_len->Draw("");
  //h_len->Scale(num_ev_bkg/h_len->Integral());
  delete functions;
}
