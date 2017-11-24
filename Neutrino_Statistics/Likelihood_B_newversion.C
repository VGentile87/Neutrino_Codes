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

Int_t Nobs  = 42;
Int_t Ntoys = 10000;


TRandom3 rand;

Double_t  sum_b_b_the = 0;  
Double_t  sum_b_n_the = 0;  
Double_t  sum_n_b_the = 0;  
Double_t  sum_n_n_the = 0;  
Double_t  sum_b_b_phi = 0;  
Double_t  sum_b_n_phi = 0;   
Double_t  sum_n_b_phi = 0;  
Double_t  sum_n_n_phi = 0;   
Double_t  sum_b_b_len = 0;  
Double_t  sum_b_n_len = 0;   
Double_t  sum_n_b_len = 0;  
Double_t  sum_n_n_len = 0;

Double_t LR_sig = 0;	  
Double_t LR_bkg = 0;	  
Double_t rb_the=0;
Double_t rb_phi=0;
Double_t rb_len=0;
Double_t rn_the=0;
Double_t rn_phi=0;
Double_t rn_len=0;

const Int_t Npoints=10;
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

Double_t time_exposure[Npoints]={1,2,3,4,5,6,7,8,9,10};

//TH1F *Lsig = new TH1F("Lsig","",128,-100,300);
//TH1F *Lbkg = new TH1F("Lbkg","",128,-100,300);

TH1F *Lsig = new TH1F("Lsig","",256,-200,200);
TH1F *Lbkg = new TH1F("Lbkg","",256,-200,200);

void Likelihood_B_newversion(){ 
   
   // 50 nm
  /*TFile *file = new TFile("allhisto_v8_50nm.root");
  TH1F* hbthe = (TH1F*)file->Get("costhetaB");
  TH1F* hbphi = (TH1F*)file->Get("phi_boro");
  TH1F* hblen = (TH1F*)file->Get("L_thr");
  TH1F* hnthe = (TH1F*)file->Get("nthe2");
  TH1F* hnphi = (TH1F*)file->Get("nphi_v2");
  TH1F* hnlen = (TH1F*)file->Get("nlen");
  */

  TFile *file = new TFile("allhisto_v11_50nm.root");
  TH1F* hbthe = (TH1F*)file->Get("theB");
  TH1F* hbphi = (TH1F*)file->Get("phiB");
  TH1F* hblen = (TH1F*)file->Get("lenB");
  TH1F* hnthe = (TH1F*)file->Get("theN");
  TH1F* hnphi = (TH1F*)file->Get("phiN");
  TH1F* hnlen = (TH1F*)file->Get("lenN");

  hbthe->Rebin(2);
  hbphi->Rebin(2);
  hblen->Rebin(2);
  hnthe->Rebin(2);
  hnphi->Rebin(2);
  hnlen->Rebin(2);
  
  Double_t mean_sign=0;

  for(int s=0;s<10;s++){
    
  vector<double> quantiles_sig;
  vector<double> quantiles_bkg;
  
  Lsig->Reset();
  Lbkg->Reset();
  EvaluateQuantiles(hbthe,hbphi,hblen,hnthe,hnphi,hnlen,quantiles_sig,quantiles_bkg);
  
    //cout << median(Lsig) << endl;
    Double_t median = median(Lsig);
    Double_t p_value = get_fraction_above(Lbkg,median);///Lbkg->Integral();
    Double_t significance = TMath::NormQuantile(p_value)*-1;
    Double_t p_value1_r = get_fraction_above(Lbkg,quantiles_sig[1]);///Lbkg->Integral();
    Double_t sigma1_r = TMath::NormQuantile(p_value1_r)*-1;
    Double_t p_value1_l = get_fraction_above(Lbkg,quantiles_sig[3]);///Lbkg->Integral();
    Double_t sigma1_l = TMath::NormQuantile(p_value1_l)*-1;
    
    printf("MEDIAN:  p-value %e significance %4.2f\n",p_value,significance);
    printf("-1sigma: p-value %e significance %4.2f\n",p_value1_r,sigma1_r);
    printf("+1sigma: p-value %e significance %4.2f\n",p_value1_l,sigma1_l);
    
    mean_sign +=significance;
     }
  mean_sign/=10.;
  cout << mean_sign << endl;
}
//----------------------------------------------------------\\


Double_t EvaluateQuantiles(TH1F *hbthe, TH1F *hbphi, TH1F *hblen, TH1F *hnthe, TH1F *hnphi, TH1F *hnlen, vector<double>& quantiles_sig, vector<double>& quantiles_bkg){

  
  hbthe->Scale(1./hbthe->Integral());
  hbphi->Scale(1./hbphi->Integral());
  hblen->Scale(1./hblen->Integral());
  
  hnthe->Scale(1./hnthe->Integral());
  hnphi->Scale(1./hnphi->Integral()); 
  hnlen->Scale(1./hnlen->Integral());
  
  

  for(Int_t itoy=0; itoy<Ntoys; itoy++){

    sum_b_b_the = 0;
    sum_b_n_the = 0;
    sum_n_b_the = 0;
    sum_n_n_the = 0;
    sum_b_b_phi = 0;
    sum_b_n_phi = 0;
    sum_n_b_phi = 0;
    sum_n_n_phi = 0;
    sum_b_b_len = 0;
    sum_b_n_len = 0;
    sum_n_b_len = 0;
    sum_n_n_len = 0;

     
    
    for(Int_t i=0; i<Nobs; i++){
      
      /*
      rb_the = hbthe->GetRandom();
      if(hbthe->GetBinContent(hbthe->FindBin(rb_the))>0)sum_b_b_the = sum_b_b_the + TMath::Log(hbthe->GetBinContent(hbthe->FindBin(rb_the)));
      if(hnthe->GetBinContent(hnthe->FindBin(rb_the))>0)sum_b_n_the = sum_b_n_the + TMath::Log(hnthe->GetBinContent(hnthe->FindBin(rb_the)));
      */rb_phi = hbphi->GetRandom();
      if(hbphi->GetBinContent(hbphi->FindBin(rb_phi))>0)sum_b_b_phi = sum_b_b_phi + TMath::Log(hbphi->GetBinContent(hbphi->FindBin(rb_phi)));
      if(hnphi->GetBinContent(hnphi->FindBin(rb_phi))>0)sum_b_n_phi = sum_b_n_phi + TMath::Log(hnphi->GetBinContent(hnphi->FindBin(rb_phi)));
      /*rb_len = hblen->GetRandom();
      if(hblen->GetBinContent(hblen->FindBin(rb_len))>0)sum_b_b_len = sum_b_b_len + TMath::Log(hblen->GetBinContent(hblen->FindBin(rb_len)));
      if(hnlen->GetBinContent(hnlen->FindBin(rb_len))>0)sum_b_n_len = sum_b_n_len + TMath::Log(hnlen->GetBinContent(hnlen->FindBin(rb_len)));      
      /*rn_the = hnthe->GetRandom();
      if(hbthe->GetBinContent(hbthe->FindBin(rn_the))>0)sum_n_b_the = sum_n_b_the + TMath::Log(hbthe->GetBinContent(hbthe->FindBin(rn_the)));
      if(hnthe->GetBinContent(hnthe->FindBin(rn_the))>0)sum_n_n_the = sum_n_n_the + TMath::Log(hnthe->GetBinContent(hnthe->FindBin(rn_the)));
      */rn_phi = hnphi->GetRandom();
      if(hbphi->GetBinContent(hbphi->FindBin(rn_phi))>0)sum_n_b_phi = sum_n_b_phi + TMath::Log(hbphi->GetBinContent(hbphi->FindBin(rn_phi)));
      if(hnphi->GetBinContent(hnphi->FindBin(rn_phi))>0)sum_n_n_phi = sum_n_n_phi + TMath::Log(hnphi->GetBinContent(hnphi->FindBin(rn_phi)));
      /*rn_len = hnlen->GetRandom();
      if(hblen->GetBinContent(hblen->FindBin(rn_len))>0)sum_n_b_len = sum_n_b_len + TMath::Log(hblen->GetBinContent(hblen->FindBin(rn_len)));
      if(hnlen->GetBinContent(hnlen->FindBin(rn_len))>0)sum_n_n_len = sum_n_n_len + TMath::Log(hnlen->GetBinContent(hnlen->FindBin(rn_len)));
      */
    }

    LR_sig = 2*(sum_b_b_the + sum_b_b_phi + sum_b_b_len - (sum_b_n_the + sum_b_n_phi + sum_b_n_len));    
    LR_bkg = 2*(sum_n_b_the + sum_n_b_phi + sum_n_b_len - (sum_n_n_the + sum_n_n_phi + sum_n_n_len));


    //cout << sum_n_b_the << " " << sum_n_b_phi << " " << sum_n_b_len << " " << sum_n_n_the << " " << sum_n_n_phi << " " << sum_n_n_len << endl;
    //cout << sum_b_n_the << " " << sum_b_n_phi << " " << sum_b_n_len << " " << sum_b_b_the << " " << sum_b_b_phi << " " << sum_b_b_len << endl;
    

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
  hbthe->Draw();
  hbthe->GetXaxis()->SetTitle("2D angle (rad)");
  hbthe->GetYaxis()->SetTitle("normalized");
  hbthe->GetYaxis()->SetTitleOffset(1.5);
  hnthe->Draw("same");
  hnthe->SetLineWidth(2);
  hbthe->SetLineWidth(2);
  hnthe->SetLineColor(kRed);
  hbthe->SetLineColor(kGreen-1);
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
  legP->AddEntry(hnthe,"neutron bkg","f");
  legP->AddEntry(hbthe,"solar #nu (^{8}B)","f");
  legP->Draw();

  TCanvas *c3 = new TCanvas("c3","",400,400);
  hbphi->Draw();
  hbphi->GetXaxis()->SetTitle("2D angle (rad)");
  hbphi->GetYaxis()->SetTitle("normalized");
  hbphi->GetYaxis()->SetTitleOffset(1.5);
  hnphi->Draw("same");
  hnphi->SetLineWidth(2);
  hbphi->SetLineWidth(2);
  hnphi->SetLineColor(kRed);
  hbphi->SetLineColor(kGreen-1);
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
  legP->AddEntry(hnphi,"neutron bkg","f");
  legP->AddEntry(hbphi,"solar #nu (^{8}B)","f");
  legP->Draw();


  TCanvas *c4 = new TCanvas("c4","",400,400);
  hblen->Draw();
  hblen->GetXaxis()->SetTitle("trk len (#mum)");
  hblen->GetYaxis()->SetTitle("normalized");
  hblen->GetYaxis()->SetTitleOffset(1.5);
  hnlen->Draw("same");
  hnlen->SetLineWidth(2);
  hnlen->SetLineWidth(2);
  hnlen->SetLineColor(kRed);
  hblen->SetLineColor(kGreen-1);
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
  legP->AddEntry(hnlen,"neutron bkg","f");
  legP->AddEntry(hblen,"solar #nu (^{8}B)","f");
  legP->Draw();
  
  
  TCanvas *csum = new TCanvas("csum","",400,400);
  Lsig->Draw();
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
  legE->AddEntry(Lbkg,"neutron bkg","f");
  legE->AddEntry(Lsig,"#nu from ^{8}B signal","f");
  legE->AddEntry(line,"median","l");
  legE->Draw();

  /*
  Double_t p_value = get_fraction_above(Lbkg,median);///Lbkg->Integral();
  Double_t significance = TMath::NormQuantile(p_value)*-1;
  Double_t p_value1_r = get_fraction_above(Lbkg,quantiles_sig[1]);///Lbkg->Integral();
  Double_t sigma1_r = TMath::NormQuantile(p_value1_r)*-1;
  Double_t p_value1_l = get_fraction_above(Lbkg,quantiles_sig[3]);///Lbkg->Integral();
  Double_t sigma1_l = TMath::NormQuantile(p_value1_l)*-1;
  
  printf("MEDIAN:  p-value %e significance %4.2f\n",p_value,significance);
  printf("-1sigma: p-value %e significance %4.2f\n",p_value1_r,sigma1_r);
  printf("+1sigma: p-value %e significance %4.2f\n",p_value1_l,sigma1_l);
  */
    
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

