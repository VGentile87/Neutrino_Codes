#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "RooStats/ModelConfig.h"
#include "TH1.h"
#include "TMath.h"
#include <vector>
#include "TLine.h"
#include "TF1.h"

Int_t Nobs  = 9;
Int_t Ntoys = 10000;
// 100 nm
/*
Double_t nSn1Ton = 2.182;
Double_t nB1Ton = 0.257;
Double_t nN1Ton = 0.274;
Double_t bkg_frB=0.48;
Double_t bkg_frN=0.52;
Double_t mass_exposure=3;
*/
// 50 nm
Double_t nSn1Ton = 6.42;
Double_t nB1Ton =1.32;
Double_t nN1Ton = 0.33;
Double_t bkg_frB=0.8;
Double_t bkg_frN=0.2;
Double_t mass_exposure=1;

TRandom3 rand;
Double_t mean_the = 0.;
Double_t sigma_the = 0.;
Double_t mean_phi = 0.;
Double_t sigma_phi = 0.;
Double_t const_sig_len = 0.;
Double_t lambda_sig_len = 0.;
Double_t const_bkg_len_neutron = 0.;
Double_t lambda_bkg_len_neutron = 0.;
Double_t const_bkg_len_boron = 0.;
Double_t lambda_bkg_len_boron = 0.;

Double_t  sum_sn_sn_the = 0;
Double_t  sum_sn_b_the = 0; 
Double_t  sum_sn_n_the = 0; 
Double_t  sum_b_sn_the = 0; 
Double_t  sum_b_b_the = 0;  
Double_t  sum_b_n_the = 0;  
Double_t  sum_n_sn_the = 0; 
Double_t  sum_n_b_the = 0;  
Double_t  sum_n_n_the = 0;  
Double_t  sum_sn_sn_phi = 0;
Double_t  sum_sn_b_phi = 0; 
Double_t  sum_sn_n_phi = 0; 
Double_t  sum_b_sn_phi = 0; 
Double_t  sum_b_b_phi = 0;  
Double_t  sum_b_n_phi = 0;  
Double_t  sum_n_sn_phi = 0; 
Double_t  sum_n_b_phi = 0;  
Double_t  sum_n_n_phi = 0;  
Double_t  sum_sn_sn_len = 0;
Double_t  sum_sn_b_len = 0; 
Double_t  sum_sn_n_len = 0; 
Double_t  sum_b_sn_len = 0; 
Double_t  sum_b_b_len = 0;  
Double_t  sum_b_n_len = 0;  
Double_t  sum_n_sn_len = 0; 
Double_t  sum_n_b_len = 0;  
Double_t  sum_n_n_len = 0;

Double_t LR_sig = 0;	  
Double_t LR_bkg = 0;	  
Double_t rsn_the=0;
Double_t rsn_phi=0;
Double_t rsn_len=0;
Double_t rb_the=0;
Double_t rb_phi=0;
Double_t rb_len=0;
Double_t rn_the=0;
Double_t rn_phi=0;
Double_t rn_len=0;

const Int_t Npoints=1;
Double_t time_exposure[Npoints]={};
Double_t Median_sig[Npoints];
Double_t Minus1s_sig[Npoints];
Double_t Plus1s_sig[Npoints];
Double_t Minus2s_sig[Npoints];
Double_t Plus2s_sig[Npoints];
Double_t Median_bkg[Npoints];
Double_t Minus1s_bkg[Npoints];
Double_t Plus1s_bkg[Npoints];
Double_t Minus2s_bkg[Npoints];
Double_t Plus2s_bkg[Npoints];

//TH1F *Lsig = new TH1F("Lsig","",256,-20,100);
//TH1F *Lbkg = new TH1F("Lbkg","",256,-20,100);

TH1F *Lsig = new TH1F("Lsig","",512,-200,600);
TH1F *Lbkg = new TH1F("Lbkg","",512,-200,600);

void Likelihood_SN_v3(){ 
   
  // 100 nm
  /*
  TFile *file = new TFile("allhisto_v5b.root");
  TH1F* hbthe = (TH1F*)file->Get("h00");
  TH1F* hbphi = (TH1F*)file->Get("phi_sum");
  TH1F* hblen = (TH1F*)file->Get("L_thr");
  TH1F* hsnthe = (TH1F*)file->Get("costhSN");
  TH1F* hsnphi = (TH1F*)file->Get("phiSN_v2");
  TH1F* hsnlen = (TH1F*)file->Get("lenSN");
  TH1F* hnthe = (TH1F*)file->Get("hthe_n");
  TH1F* hnphi = (TH1F*)file->Get("nphi");
  TH1F* hnlen = (TH1F*)file->Get("nlen");
  */

  
  // 50 nm
  /*TFile *file = new TFile("allhisto_v9_50nm.root");
  TH1F* hbthe = (TH1F*)file->Get("theta_B");
  TH1F* hbphi = (TH1F*)file->Get("phi_boro");
  TH1F* hblen = (TH1F*)file->Get("L_thr");
  TH1F* hsnthe = (TH1F*)file->Get("theta_v4");
  TH1F* hsnphi = (TH1F*)file->Get("phiSN_v2");
  //TH1F* hsnphi = (TH1F*)file->Get("nphi");
  TH1F* hsnlen = (TH1F*)file->Get("lenSN");
  TH1F* hnthe = (TH1F*)file->Get("theta_N");
  TH1F* hnphi = (TH1F*)file->Get("nphi_v2");
  TH1F* hnlen = (TH1F*)file->Get("nlen");*/

  TFile *file = new TFile("allhisto_v11_50nm.root");
  TH1F* hbthe = (TH1F*)file->Get("theB");
  TH1F* hbphi = (TH1F*)file->Get("phiB");
  TH1F* hblen = (TH1F*)file->Get("lenB");
  TH1F* hsnthe = (TH1F*)file->Get("theSN");
  TH1F* hsnphi = (TH1F*)file->Get("phiSN");
  TH1F* hsnlen = (TH1F*)file->Get("lenSN");
  TH1F* hnthe = (TH1F*)file->Get("theN");
  TH1F* hnphi = (TH1F*)file->Get("phiN");
  TH1F* hnlen = (TH1F*)file->Get("lenN");
  /*
  hbthe->Rebin(2);
  hbphi->Rebin(2);
  hblen->Rebin(2);
  hsnthe->Rebin(2);
  hsnphi->Rebin(2);
  hsnlen->Rebin(2);
  hnthe->Rebin(2);
  hnphi->Rebin(2);
  hnlen->Rebin(2);
  */
  Double_t mean_sign=0;

  for(int s=0;s<10;s++){
    
    vector<double> quantiles_sig;
    vector<double> quantiles_bkg;
    
    Lsig->Reset();
    Lbkg->Reset();

    
    EvaluateQuantiles(hsnthe,hsnphi,hsnlen,hbthe,hbphi,hblen,hnthe,hnphi,hnlen,quantiles_sig,quantiles_bkg);
    
    Double_t median = median(Lsig);
    Double_t p_value = get_fraction_above(Lbkg,median);///Lbkg->Integral();
    Double_t significance = TMath::NormQuantile(p_value);//*-1;
    Double_t p_value1_r = get_fraction_above(Lbkg,quantiles_sig[1]);///Lbkg->Integral();
    Double_t sigma1_r = TMath::NormQuantile(p_value1_r);//*-1;
    Double_t p_value1_l = get_fraction_above(Lbkg,quantiles_sig[3]);///Lbkg->Integral();
    Double_t sigma1_l = TMath::NormQuantile(p_value1_l);//*-1;
    
    printf("MEDIAN:  p-value %e significance %4.2f\n",p_value,significance);
    printf("-1sigma: p-value %e significance %4.2f\n",p_value1_r,sigma1_r);
    printf("+1sigma: p-value %e significance %4.2f\n",p_value1_l,sigma1_l);

    if(significance==0)significance=5;
    mean_sign +=significance;
  }
  mean_sign/=10.;
  cout << mean_sign << endl;

}
//----------------------------------------------------------\\





Double_t EvaluateQuantiles(TH1F *hsnthe, TH1F *hsnphi, TH1F *hsnlen, TH1F *hbthe, TH1F *hbphi, TH1F *hblen, TH1F *hnthe, TH1F *hnphi, TH1F *hnlen, vector<double>& quantiles_sig, vector<double>& quantiles_bkg){
  //Double_t EvaluateQuantiles(Double_t mean, Double_t sigma, vector<double>& quantiles_sig,  vector<double>& quantiles_bkg){
  /*
  hbthe->Scale(nB1Ton/hbthe->Integral());
  hbphi->Scale(nB1Ton/hbphi->Integral());
  hblen->Scale(nB1Ton/hblen->Integral());
  
  hnthe->Scale(nN1Ton/hnthe->Integral());
  hnphi->Scale(nN1Ton/hnphi->Integral()); 
  hnlen->Scale(nN1Ton/hnlen->Integral());
  
  hsnthe->Scale(nSn1Ton/hsnthe->Integral());
  hsnphi->Scale(nSn1Ton/hsnphi->Integral());
  hsnlen->Scale(nSn1Ton/hsnlen->Integral());
  */
  
  hbthe->Scale(1./hbthe->Integral());
  hbphi->Scale(1./hbphi->Integral());
  hblen->Scale(1./hblen->Integral());
  
  hnthe->Scale(1./hnthe->Integral());
  hnphi->Scale(1./hnphi->Integral()); 
  hnlen->Scale(1./hnlen->Integral());
  
  hsnthe->Scale(1./hsnthe->Integral());
  hsnphi->Scale(1./hsnphi->Integral());
  hsnlen->Scale(1./hsnlen->Integral());
  

  for(Int_t itoy=0; itoy<Ntoys; itoy++){

    sum_sn_sn_the = 0;
    sum_sn_b_the = 0;
    sum_sn_n_the = 0;
    sum_b_sn_the = 0;
    sum_b_b_the = 0;
    sum_b_n_the = 0;
    sum_n_sn_the = 0;
    sum_n_b_the = 0;
    sum_n_n_the = 0;
    sum_sn_sn_phi = 0;
    sum_sn_b_phi = 0;
    sum_sn_n_phi = 0;
    sum_b_sn_phi = 0;
    sum_b_b_phi = 0;
    sum_b_n_phi = 0;
    sum_n_sn_phi = 0;
    sum_n_b_phi = 0;
    sum_n_n_phi = 0;
    sum_sn_sn_len = 0;
    sum_sn_b_len = 0;
    sum_sn_n_len = 0;
    sum_b_sn_len = 0;
    sum_b_b_len = 0;
    sum_b_n_len = 0;
    sum_n_sn_len = 0;
    sum_n_b_len = 0;
    sum_n_n_len = 0;

     
    
    for(Int_t i=0; i<Nobs; i++){

      //
      /*
      rsn_the = hsnthe->GetRandom();
      if(hsnthe->GetBinContent(hsnthe->FindBin(rsn_the))>0)sum_sn_sn_the = sum_sn_sn_the + TMath::Log(hsnthe->GetBinContent(hsnthe->FindBin(rsn_the)));
      if(hbthe->GetBinContent(hbthe->FindBin(rsn_the))>0)sum_sn_b_the = sum_sn_b_the + TMath::Log(hbthe->GetBinContent(hbthe->FindBin(rsn_the)));
      if(hnthe->GetBinContent(hnthe->FindBin(rsn_the))>0)sum_sn_n_the = sum_sn_n_the + TMath::Log(hnthe->GetBinContent(hnthe->FindBin(rsn_the)));
      
      //cout << rsn_the << " " << sum_sn_sn_the << " " << sum_sn_b_the << " " << sum_sn_n_the << endl;
      //cout << rsn_the << " " << hsnthe->GetBinContent(hsnthe->FindBin(rsn_the)) << " " << hsnthe->FindBin(rsn_the) << " " <<  TMath::Log(hsnthe->GetBinContent(hsnthe->FindBin(rsn_the))) << endl;
      */
      rsn_phi = hsnphi->GetRandom();
      if(hsnphi->GetBinContent(hsnphi->FindBin(rsn_phi))>0)sum_sn_sn_phi = sum_sn_sn_phi + TMath::Log(hsnphi->GetBinContent(hsnphi->FindBin(rsn_phi)));
      if(hbphi->GetBinContent(hbphi->FindBin(rsn_phi))>0)sum_sn_b_phi = sum_sn_b_phi + TMath::Log(hbphi->GetBinContent(hbphi->FindBin(rsn_phi)));
      if(hnphi->GetBinContent(hnphi->FindBin(rsn_phi))>0)sum_sn_n_phi = sum_sn_n_phi + TMath::Log(hnphi->GetBinContent(hnphi->FindBin(rsn_phi)));
      /*rsn_len = hsnlen->GetRandom();
      if(hsnlen->GetBinContent(hsnlen->FindBin(rsn_len))>0)sum_sn_sn_len = sum_sn_sn_len + TMath::Log(hsnlen->GetBinContent(hsnlen->FindBin(rsn_len)));
      if(hblen->GetBinContent(hblen->FindBin(rsn_len))>0)sum_sn_b_len = sum_sn_b_len + TMath::Log(hblen->GetBinContent(hblen->FindBin(rsn_len)));
      if(hnlen->GetBinContent(hnlen->FindBin(rsn_len))>0)sum_sn_n_len = sum_sn_n_len + TMath::Log(hnlen->GetBinContent(hnlen->FindBin(rsn_len)));
      
      
      rb_the = hbthe->GetRandom();
      if(hsnthe->GetBinContent(hsnthe->FindBin(rb_the))>0)sum_b_sn_the = sum_b_sn_the + TMath::Log(hsnthe->GetBinContent(hsnthe->FindBin(rb_the)));
      if(hbthe->GetBinContent(hbthe->FindBin(rb_the))>0)sum_b_b_the = sum_b_b_the + TMath::Log(hbthe->GetBinContent(hbthe->FindBin(rb_the)));
      if(hnthe->GetBinContent(hnthe->FindBin(rb_the))>0)sum_b_n_the = sum_b_n_the + TMath::Log(hnthe->GetBinContent(hnthe->FindBin(rb_the)));
      */rb_phi = hbphi->GetRandom();
      if(hsnphi->GetBinContent(hsnphi->FindBin(rb_phi))>0)sum_b_sn_phi = sum_b_sn_phi + TMath::Log(hsnphi->GetBinContent(hsnphi->FindBin(rb_phi)));
      if(hbphi->GetBinContent(hbphi->FindBin(rb_phi))>0)sum_b_b_phi = sum_b_b_phi + TMath::Log(hbphi->GetBinContent(hbphi->FindBin(rb_phi)));
      if(hnphi->GetBinContent(hnphi->FindBin(rb_phi))>0)sum_b_n_phi = sum_b_n_phi + TMath::Log(hnphi->GetBinContent(hnphi->FindBin(rb_phi)));
      /*rb_len = hblen->GetRandom();
      if(hsnlen->GetBinContent(hsnlen->FindBin(rb_len))>0)sum_b_sn_len = sum_b_sn_len + TMath::Log(hsnlen->GetBinContent(hsnlen->FindBin(rb_len)));
      if(hblen->GetBinContent(hblen->FindBin(rb_len))>0)sum_b_b_len = sum_b_b_len + TMath::Log(hblen->GetBinContent(hblen->FindBin(rb_len)));
      if(hnlen->GetBinContent(hnlen->FindBin(rb_len))>0)sum_b_n_len = sum_b_n_len + TMath::Log(hnlen->GetBinContent(hnlen->FindBin(rb_len)));
      
      rn_the = hnthe->GetRandom();
      if(hsnthe->GetBinContent(hsnthe->FindBin(rn_the))>0)sum_n_sn_the = sum_n_sn_the + TMath::Log(hsnthe->GetBinContent(hsnthe->FindBin(rn_the)));
      if(hbthe->GetBinContent(hbthe->FindBin(rn_the))>0)sum_n_b_the = sum_n_b_the + TMath::Log(hbthe->GetBinContent(hbthe->FindBin(rn_the)));
      if(hnthe->GetBinContent(hnthe->FindBin(rn_the))>0)sum_n_n_the = sum_n_n_the + TMath::Log(hnthe->GetBinContent(hnthe->FindBin(rn_the)));
      */rn_phi = hnphi->GetRandom();
      if(hsnphi->GetBinContent(hsnphi->FindBin(rn_phi))>0)sum_n_sn_phi = sum_n_sn_phi + TMath::Log(hsnphi->GetBinContent(hsnphi->FindBin(rn_phi)));
      if(hbphi->GetBinContent(hbphi->FindBin(rn_phi))>0)sum_n_b_phi = sum_n_b_phi + TMath::Log(hbphi->GetBinContent(hbphi->FindBin(rn_phi)));
      if(hnphi->GetBinContent(hnphi->FindBin(rn_phi))>0)sum_n_n_phi = sum_n_n_phi + TMath::Log(hnphi->GetBinContent(hnphi->FindBin(rn_phi)));
      /*rn_len = hnlen->GetRandom();
      if(hsnlen->GetBinContent(hsnlen->FindBin(rn_len))>0)sum_n_sn_len = sum_n_sn_len + TMath::Log(hsnlen->GetBinContent(hsnlen->FindBin(rn_len)));
      if(hblen->GetBinContent(hblen->FindBin(rn_len))>0)sum_n_b_len = sum_n_b_len + TMath::Log(hblen->GetBinContent(hblen->FindBin(rn_len)));
      if(hnlen->GetBinContent(hnlen->FindBin(rn_len))>0)sum_n_n_len = sum_n_n_len + TMath::Log(hnlen->GetBinContent(hnlen->FindBin(rn_len)));
      */
    }

    LR_sig = 2*(sum_sn_sn_the + sum_sn_sn_phi + sum_sn_sn_len - bkg_frB*(sum_sn_b_the + sum_sn_b_phi + sum_sn_b_len) - bkg_frN*(sum_sn_n_the + sum_sn_n_phi + sum_sn_n_len));
    
    LR_bkg = 2*((sum_b_sn_the + sum_b_sn_phi + sum_b_sn_len) + (sum_n_sn_the + sum_n_sn_phi + sum_n_sn_len) - bkg_frB*(sum_b_b_the + sum_b_b_phi + sum_b_b_len + sum_n_b_the + sum_n_b_phi + sum_n_b_len) - bkg_frN*(sum_n_n_the + sum_n_n_phi + sum_n_n_len + sum_b_n_the + sum_b_n_phi + sum_b_n_len));


    //cout << sum_sn_b_the << " " << sum_sn_b_phi << " " << sum_sn_b_len << endl;
    

    if(0){
      //printf("sum_sn_sn_the %4.2f sum_sig_const_the %4.2f sum_bkg_gaus_the %4.2f sum_bkg_const_the %4.2f\n",sum_sn_sn_the,sum_sig_const_the,sum_bkg_gaus_the,sum_bkg_const_the);
      //printf("LR_sig %4.2f LR_bkg %4.2f\n",LR_sig,LR_bkg);
    }
    Lsig->Fill(LR_sig);
    Lbkg->Fill(LR_bkg);
    //cout << LR_sig << " " << LR_bkg << endl;
    
  }
  
  quantiles_sig = quantiles(Lsig);
  quantiles_bkg = quantiles(Lbkg);
  
  //----------------------------------------
  
  gStyle->SetOptStat(0);

  //TCanvas *c0 = new TCanvas("c0","",400,400);
  //hh->Draw();
  
  TCanvas *c2 = new TCanvas("c2","",400,400);
  hsnthe->Draw();
  hsnthe->GetXaxis()->SetTitle("2D angle (rad)");
  hsnthe->GetYaxis()->SetTitle("normalized");
  hsnthe->GetYaxis()->SetTitleOffset(1.5);
  hbthe->Draw("same");
  hnthe->Draw("same");
  //hbthe_neutron->Smooth(6);
  //hbthe_boron->Smooth(6);
  //hsnthe->Smooth(6);
  hsnthe->SetLineWidth(2);
  hnthe->SetLineWidth(2);
  hbthe->SetLineWidth(2);
  hsnthe->SetLineColor(kBlue);
  hnthe->SetLineColor(kRed);
  hbthe->SetLineColor(kGreen-1);
  hsnthe->SetFillColor(38);
  hnthe->SetFillColor(kRed);
  hbthe->SetFillColor(kGreen-1);
  hnthe->SetFillStyle(3005);
  hbthe->SetFillStyle(3004);
  TLegend * legP = new TLegend();
  legP = new TLegend(0.15,0.60,0.4,0.80);
  legP->SetTextSize(0.03);
  legP->SetBorderSize(0);
  legP->SetLineStyle(0);
  legP->SetTextSize(0.03);
  legP->SetFillStyle(1001);
  legP->SetFillColor(kWhite);
  legP->AddEntry(hsnthe,"SN signal","f");
  legP->AddEntry(hnthe,"neutron bkg","f");
  legP->AddEntry(hbthe,"solar #nu (^{8}B)","f");
  legP->Draw();

  TCanvas *c3 = new TCanvas("c3","",400,400);
  hsnphi->Draw();
  hsnphi->GetXaxis()->SetTitle("2D angle (rad)");
  hsnphi->GetYaxis()->SetTitle("normalized");
  hsnphi->GetYaxis()->SetTitleOffset(1.5);
  hbphi->Draw("same");
  hnphi->Draw("same");
  //hbphi_neutron->Smooth(6);
  //hbphi_boron->Smooth(6);
  //hsnphi->Smooth(6);
  hsnphi->SetLineWidth(2);
  hnphi->SetLineWidth(2);
  hbphi->SetLineWidth(2);
  hsnphi->SetLineColor(kBlue);
  hnphi->SetLineColor(kRed);
  hbphi->SetLineColor(kGreen-1);
  hsnphi->SetFillColor(38);
  hnphi->SetFillColor(kRed);
  hbphi->SetFillColor(kGreen-1);
  hnphi->SetFillStyle(3005);
  hbphi->SetFillStyle(3004);
  TLegend * legP = new TLegend();
  legP = new TLegend(0.15,0.60,0.4,0.80);
  legP->SetTextSize(0.03);
  legP->SetBorderSize(0);
  legP->SetLineStyle(0);
  legP->SetTextSize(0.03);
  legP->SetFillStyle(1001);
  legP->SetFillColor(kWhite);
  legP->AddEntry(hsnphi,"SN signal","f");
  legP->AddEntry(hnphi,"neutron bkg","f");
  legP->AddEntry(hbphi,"solar #nu (^{8}B)","f");
  legP->Draw();


  TCanvas *c4 = new TCanvas("c4","",400,400);
  hsnlen->Draw();
  hsnlen->GetXaxis()->SetTitle("trk len (#mum)");
  hsnlen->GetYaxis()->SetTitle("normalized");
  hsnlen->GetYaxis()->SetTitleOffset(1.5);
  hblen->Draw("same");
  hnlen->Draw("same");
  //hblen_neutron->Smooth(6);
  //hblen_boron->Smooth(6);
  //hsnlen->Smooth(6);
  hsnlen->SetLineWidth(2);
  hnlen->SetLineWidth(2);
  hnlen->SetLineWidth(2);
  hsnlen->SetLineColor(kBlue);
  hnlen->SetLineColor(kRed);
  hblen->SetLineColor(kGreen-1);
  hsnlen->SetFillColor(38);
  hnlen->SetFillColor(kRed);
  hblen->SetFillColor(kGreen-1);
  hnlen->SetFillStyle(3005);
  hblen->SetFillStyle(3004);
  TLegend * legP = new TLegend();
  legP = new TLegend(0.15,0.60,0.4,0.80);
  legP->SetTextSize(0.03);
  legP->SetBorderSize(0);
  legP->SetLineStyle(0);
  legP->SetTextSize(0.03);
  legP->SetFillStyle(1001);
  legP->SetFillColor(kWhite);
  legP->AddEntry(hsnlen,"SN signal","f");
  legP->AddEntry(hnlen,"neutron bkg","f");
  legP->AddEntry(hblen,"solar #nu (^{8}B)","f");
  legP->Draw();
  
  
  TCanvas *csum = new TCanvas("csum","",400,400);
  //Lbkg->Smooth(6);
  //Lsig->Smooth(6);
  Lsig->Draw();
  //Lsig->Scale(1./hsnthe->Integral());
  //Lbkg->Scale(1./hnthe->Integral());
  Lsig->SetLineStyle(2);
  Lbkg->SetLineStyle(2);
  Lsig->SetLineColor(kRed);
  Lbkg->SetLineColor(kBlue);
  Lbkg->Draw("same");
  Lbkg->GetXaxis()->SetTitle("2#times ln(L_{s}/L_{b})");
  Lbkg->GetYaxis()->SetTitle("Pseudoexperiments");
  Lbkg->GetYaxis()->SetTitleOffset(1.5);
  Lbkg->GetYaxis()->SetRangeUser(0,Lsig->GetMaximum()*1.1);
  Lbkg->GetXaxis()->CenterTitle();
  Lbkg->GetYaxis()->CenterTitle();
  Lbkg->SetFillColorAlpha(kAzure+7,0.75);
  Lsig->SetFillColorAlpha(kOrange-3,0.75);

  Double_t median = median(Lsig);
  printf("mediana %f\n",median);
  Float_t ymax = Lsig->GetMaximum();
  TLine *line = new TLine(median,0,median,ymax);
  line->SetLineColor(kRed);
  line->Draw();

  pt = new TPaveText(0.30,0.92,0.7,0.98, "NDC"); // NDC sets coords
  // relative to pad dimensions
  pt->SetFillColor(0); // text is black on white
  pt->SetTextSize(0.04); 
  pt->SetTextAlign(12);
  pt->SetBorderSize(0);
  pt->SetLineStyle(0);
  char string[255];
  sprintf(string,"Nobs = %d",Nobs);
  text = pt->AddText(string);
  pt->Draw();       //to draw your text object
  
  TLegend * legE = new TLegend();
  legE = new TLegend(0.15,0.60,0.4,0.80);
  legE->SetTextSize(0.03);
  legE->SetBorderSize(0);
  legE->SetLineStyle(0);
  legE->SetTextSize(0.03);
  legE->SetFillStyle(1001);
  legE->SetFillColor(kWhite);
  //legE->AddEntry(Lbkg2,"boron","f");
  legE->AddEntry(Lbkg,"background","f");
  legE->AddEntry(Lsig,"SN signal","f");
  legE->AddEntry(line,"median","l");
  legE->Draw();

  /*
  Double_t p_value = get_fraction_above(Lbkg,median);///Lbkg->Integral();
  Double_t significance = TMath::NormQuantile(p_value);//*-1;
  Double_t p_value1_r = get_fraction_above(Lbkg,quantiles_sig[1]);///Lbkg->Integral();
  Double_t sigma1_r = TMath::NormQuantile(p_value1_r);//*-1;
  Double_t p_value1_l = get_fraction_above(Lbkg,quantiles_sig[3]);///Lbkg->Integral();
  Double_t sigma1_l = TMath::NormQuantile(p_value1_l);//*-1;
  
  printf("MEDIAN:  p-value %e significance %4.2f\n",p_value,significance);
  printf("-1sigma: p-value %e significance %4.2f\n",p_value1_r,sigma1_r);
  printf("+1sigma: p-value %e significance %4.2f\n",p_value1_l,sigma1_l);
  */
  //char name[255];
  //sprintf(name,"results/Likelihood_SN_%0.3f.pdf",sigma_the);
  //csum->SaveAs(name);
  //csum->SaveAs("results/Likelihood_SN_all.gif+50");
    
  return 1;
 }



//----------------------------------------

// Quantiles returns a vector<double> with 5 entries.
// Entries 0 and 4 are the values on the histogram x-axis
// so that 95% of the content lies between these values.
// Entries 1 and 3 bracket 68% in the same way.
// Entry 2 is the median of the histogram.
//==================================================
vector<double> quantiles( TH1F* h ){
//==================================================
  double q[5];
  double probs[5] = {0.025, 0.16, 0.5, 1 - 0.16, 0.975 };
  h->GetQuantiles( 5, q, probs );
  
  vector<double> r(5);
  for (int i=0; i<5; i++) 
    {
      cout << " quantile " << i << " " << q[i] << endl;
      r[i] = q[i];
    }
  return r;
}


// Return the median M of a histogram : 50% of the area under the histogram
// will be below M, and 50% above M.
//======================
double median( TH1F* h ){
//======================
  return quantiles(h)[2];
}

// This function plots a histogram and draws the central 68(95)% 
// in green(yellow).
//===========================
void plot_with_bands(TH1F* h){
//===========================
  TH1F* h_68 = (TH1F*) h->Clone( Form("%s_q68", h->GetName() ) );
  TH1F* h_95 = (TH1F*) h->Clone( Form("%s_q95", h->GetName() ) );

  vector<double> q = quantiles( h );

  for (int i=1 ; i < h ->GetNbinsX() ; i++)
    {
      double x = h->GetBinCenter( i );
      if ( x < q[0] || x > q[4] ) h_95 -> SetBinContent( i , 0 );
      if ( x < q[1] || x > q[3] ) h_68 -> SetBinContent( i , 0 );
    }

  h_95 -> SetFillColor(5);
  h_68 -> SetFillColor(3);
  
  h->SetFillStyle(0);
  h->SetLineWidth(2);
  h->Draw();
  h_95->Draw("same");
  h_68->Draw("same");

  TLine * l = new TLine( q[2],0, q[2] , h->GetBinContent( h->FindBin(q[2]) ) );
  l->SetLineColor(2);
  l->SetLineWidth(2);
  l->Draw();
}


// Return the fraction of the area below the histogram with x<x_value
// (no interpolation between bins, so expect binning effects).
// the bin that contains x_value is included in the fraction.
//==================================================
double get_fraction_below( TH1F* h, double x_value ){
//==================================================
  int bin = h->FindBin( x_value );
  return h->Integral( -100000, bin ) / h->Integral();
}

// Return the fraction of the area below the histogram with x>x_value
// the bin that contains x_value is not included in the fraction.
//===================================================
double get_fraction_above(  TH1F* h, double x_value ){
//===================================================
  return 1 - get_fraction_below( h, x_value);
}

