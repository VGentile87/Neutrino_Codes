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
//#include "Functions.h"
#include "EnergyRangeCorrelation.C"
#include "TRandom3.h"
#include <TSystem.h>



using namespace std;

void nu_floor_News() {

  Int_t resto=0;
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
  dir.ReplaceAll("nu_floor_News.c","");
  dir.ReplaceAll("/./","/");
  ifstream in;

  in.open(Form("%sdsnb_spectrum3.txt",dir.Data()));
  while (1) {
    in >> dsnb_nu_en[dsnb_iline] >> dsnb_nu_rate[dsnb_iline];
    dsnb_nu_rate[dsnb_iline]=dsnb_nu_rate[dsnb_iline];//*dsnb_int_flux;
    //cout << dsnb_iline << " " << dsnb_nu_rate[dsnb_iline] << endl;
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
    //if(i<10)xpoint[i]=(i+1)/10.;
    //if(i>=10)xpoint[i]=(i+1)-10;
    //tot_rate[i]=tot_rate[i]+dsnb_gr->Eval(xpoint[i]);
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
   
    // tot_rate[i]=tot_rate[i]+atm_gr->Eval(xpoint[i]);
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
    //if(i<10)xpoint[i]=(i+1)/10.;
    //if(i>=10)xpoint[i]=(i+1)-10;
    tot_rate[i]=tot_rate[i]+B_gr->Eval(xpoint[i]);
    //cout << xpoint[i] << " " << B_gr->Eval(xpoint[i]) << " " << tot_rate[i] << endl; 
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
    //if(i<10)xpoint[i]=(i+1)/10.;
    //if(i>=10)xpoint[i]=(i+1)-10;
    //tot_rate[i]=tot_rate[i]+hep_gr->Eval(xpoint[i]);
    //cout << xpoint[i] << " " << hep_gr->Eval(xpoint[i]) << " " << tot_rate[i] << endl; 
  }

  in.close();

  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  TGraph *tot_rate_gr = new TGraph(196,xpoint,tot_rate);
  tot_rate_gr->Draw("A*");

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
  log_file << "Total neutrino induced recoil rate";
  log_file << endl;
  
  //Functions* functions;

  //functions = new Functions();
  

  //TGraph *gr_ev_Erec = new TGraph();
  TH1D * h_wimp_ev_Erec = new TH1D("wimp_ene","",100000,0.001,100);
  TH1D * h_ev_Erec = new TH1D("ene","",100000,0.001,100);
  TH1D * h_ev_Ltrue = new TH1D("len","",1000,0.01,10);
  Int_t ipoint=0;

  
  for(int k=0;k<50000;k++){
    num_ev_bkg=0;
    num_wimp_ev=0;
    E_rec=k/1000.;
    //sum_crsect=0;
    
    for(int i=3; i<n_el;i++){   
      iel = el[i];
      iZ_el = Z_el[i];
      iA_el = A_el[i];
      ifract_el = atomic_fract_el[i]; // uma
      imass_fract = mass_fract[i];
      iM_el=iA_el*uma;
      //iN_el = N_el[i];
      
      iN_el=iA_el-iZ_el;
      iQw = iN_el - (1-4*sin2wnb)*iZ_el;
      irn = TMath::Sqrt(1.5129*TMath::Power(iA_el,2./3.) -1.4760*TMath::Power(iA_el,1./3.) + 2.5371);
      q=TMath::Sqrt((2*iM_el*E_rec*TMath::Power(10,-6))); //GeV
      qrn = q*irn*fm_in_GeV;    
      qs = q*s*fm_in_GeV;       
      Fq = (3*(TMath::Sin(qrn)-qrn*TMath::Cos(qrn))/TMath::Power(qrn,3))*TMath::Exp(-(qs*qs)/2.);  // diverso rispetto a plot di esclusione
      int_const_bkg = (M_riv*imass_fract)/(iA_el*uma_in_ton);
      L_true = GetRangeFromEnergy(iZ_el,E_rec);
      sum_flux=0;
      //for(int j=0; j<195;j++){
      for(int j=0; j<818;j++){
	E_nu=B_nu_en[j];
	iflux_boron_vs_E=B_nu_rate[j]/50.;
	sum_flux+=iflux_boron_vs_E;
        crsect = ((GF*GF*iQw*iQw*Fq*Fq*iM_el*1000)/(4*pi))*(1-(iM_el*E_rec)/(2*E_nu*E_nu))*TMath::Power(iMeV_in_fm,2)*TMath::Power(10,-26)/TMath::Power(10,6);
	
	if(E_nu >(TMath::Sqrt(iM_el*E_rec/2.)) && L_true>100){
	  if(E_nu==10 && crsect>0)sum_crsect+=crsect;	  
	  num_ev_bkg += crsect*iflux_boron_vs_E*int_const_bkg;
	  // cout << E_nu <<" " << E_rec <<  " " << sum_crsect << " " << int_const_bkg << endl;
	}
	tmp_B_rate=B_nu_rate[j];
      }
      
      //cout << sum_crsect << " " << sum_flux <<  " " << int_const_bkg << endl;
       
      /// WIMP ///
      Md = 10;
      MN = (Md*0.932*iA_el)/(Md+iA_el*0.932);
      Mp = (Md*0.938)/(Md+0.938);
      xsec = TMath::Power(iA_el*(MN/Mp),2);
      r = (4*Md*iA_el*0.932)/(TMath::Power(Md+iA_el*0.932,2));
      R0 = (2./TMath::Power(pi,0.5))*(numav/iA_el)*(densdm/Md)*v0*TMath::Power(10,5)*xsec;      
      E0 = 0.5*Md*TMath::Power(v0/c,2)*TMath::Power(10,6);  // keV
      vmin = v0*TMath::Power(E_rec/(E0*r),0.5);
      if(vmin > TMath::Sqrt(2*iM_el*E_rec)/(2*r) && L_true>100){
	//exp1 = (R0/(E0*r))*((TMath::Power(pi,0.5)/4.)*(v0/vE)*(TMath::Erf((vmin+vE)/v0)-TMath::Erf((vmin-vE)/v0)));
	exp1 = 0.9965*(R0/(E0*r))*((TMath::Power(pi,0.5)/4.)*(v0/vE)*(TMath::Erf((vmin+vE)/v0)-TMath::Erf((vmin-vE)/v0))-TMath::Exp(-(vesc*vesc)/(v0*v0)));
	num_wimp_ev += exp1*Fq*Fq*imass_fract;
      }
      ////////
      
    }
    //num_ev_bkg = int_const_bkg*sum_crsect;
    num_ev_bkg = num_ev_bkg*86400*365*400; //TMath::Power(10,33)/1.66054;
    h_ev_Erec->Fill(E_rec,num_ev_bkg); //num_ev_bkg
    sum_ev_bkg+=num_ev_bkg;

    num_wimp_ev = num_wimp_ev*86400*365*400*0.21; //1000 sta per ton
    //cout << k << " " << num_ev_bkg << endl;
    h_wimp_ev_Erec->Fill(E_rec,num_wimp_ev); //num_ev_bkg
  }

 
  TCanvas *c3 = new TCanvas("c3","c3",1200,600);
  h_ev_Erec->Draw("");

  
  TCanvas *c4 = new TCanvas("c4","c4",1200,600);
  h_wimp_ev_Erec->Draw();
  h_wimp_ev_Erec->Scale(crsect_wimp);
  h_ev_Erec->Draw("same");
  h_ev_Erec->SetLineColor(kRed+1);
}
