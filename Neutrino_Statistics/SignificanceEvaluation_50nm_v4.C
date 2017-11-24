#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooHistPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "RooStats/ModelConfig.h"
#include "TH1.h"
#include "TMath.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <TGraph.h>
#include <TLine.h>
#include <TLegend.h>

using namespace RooFit;
using namespace RooStats;

gSystem->Load("libRooFit");
gSystem->Load("libRooStats");

RooWorkspace w("w",kTRUE);

TH1F *hSig    = new TH1F("hSig","",100,0,100);
TH1I *hN      = new TH1I("hN","",100,0,100);
//TCanvas *cSig = new TCanvas("cSig","",600,600);
//TCanvas *cN   = new TCanvas("cN","",600,600);


void SignificanceEvaluation_50nm_v4()
{
  //ofstream log_file("log_1ton_1y.txt",ios::app);
  const int Npoints=7;
  const double mass_exposure=2; // ton
  double nsig = 0;
  double nsig_1y = 6.42;
  double nbkg = 0;
  double nbkg_1y = 1.65;
  double sigmab = 0.2;
  double significance=0;
  double npois = 0;
  TRandom3 rand;

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

  vector<double> quantiles_sig;
  vector<double> quantiles_bkg;

  for(int ip=0; ip<Npoints; ip++){

    hSig->Reset(0);
    hN->Reset(0);

    time_exposure[ip]=ip+1;
    nbkg = nbkg_1y*time_exposure[ip]*mass_exposure;
    nsig = nsig_1y*mass_exposure;
    
    for(int ir=0; ir<100; ir++){
      RooWorkspace w("w");
      npois=rand.Poisson(TMath::Nint(nsig));
      MakeModel(npois,nbkg,sigmab);
      significance = MakeSignificance();
      if(significance!=0) hSig->Fill(significance);
      hN->Fill(npois);
  }
    
    cout << ip << " " << nsig << " " << nbkg << " " << sigmab <<  endl;
    
    //cSig->cd();
    //hSig->Draw();
    //cN->cd();
    //hN->Draw();
    
    Double_t median = median(hSig);
    cout << " MEDIANA " << median << endl;   
    TString fileName = "MyModel.root";

    quantiles_sig = quantiles(hSig);
    
    Median_sig[ip]  = quantiles_sig.at(2);
    Minus2s_sig[ip] = quantiles_sig.at(2) - quantiles_sig.at(0);
    Minus1s_sig[ip] = quantiles_sig.at(2) - quantiles_sig.at(1);
    Plus1s_sig[ip]  = quantiles_sig.at(3) - quantiles_sig.at(2);
    Plus2s_sig[ip]  = quantiles_sig.at(4) - quantiles_sig.at(2);
  }

  TLine *line = new TLine(1,3,Npoints,3);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->SetLineColor(kRed+2);
    
  Double_t exl[Npoints] = {0};
  Double_t exh[Npoints] = {0};
  
  TMultiGraph *mgExcl = new TMultiGraph();
  gr0_sig = new TGraph(Npoints,time_exposure,Median_sig);
  gr1_sig = new TGraphAsymmErrors(Npoints,time_exposure,Median_sig,exl,exh,Minus1s_sig,Plus1s_sig);
  gr2_sig = new TGraphAsymmErrors(Npoints,time_exposure,Median_sig,exl,exh,Minus2s_sig,Plus2s_sig);
  gr1_sig->SetFillColor(kGreen);
  gr2_sig->SetFillColor(kYellow);
  gr0_sig->SetMarkerStyle(20);
  gr0_sig->SetMarkerColor(kBlue);
  gr0_sig->SetMarkerSize(0.8);
  gr0_sig->SetLineStyle(2);
  gr0_sig->SetLineWidth(2.8);
  gr0_sig->SetLineColor(kBlue);
  
  
  TCanvas *c1 = new TCanvas("c1","",700,700);
  //c1->SetLogx();
  //c1->SetLogy();
  //c1->SetGridx();
  //c1->SetGridy();
  
  mgExcl->Add(gr2_sig);
  mgExcl->Add(gr1_sig);
  mgExcl->Add(gr0_sig);
  mgExcl->Draw("alp3");
  line->Draw();
  c1->Update();
  mgExcl->GetYaxis()->SetTitleOffset(1.1);
  mgExcl->GetXaxis()->SetTitleOffset(1.1);
  mgExcl->GetYaxis()->SetTitle("Significance [#sigma]");
  mgExcl->GetXaxis()->SetTitle("exposure time [y]");
  mgExcl->GetXaxis()->CenterTitle();
  mgExcl->GetYaxis()->CenterTitle();
  //mgExcl->GetYaxis()->SetRangeUser(1.0E-43,1.0E-37);
  c1->Modified();
  
  TLegend * legS = new TLegend();
  legS = new TLegend(0.50,0.60,0.7,0.80);
  legS->SetTextSize(0.02);
  legS->SetBorderSize(0);
  legS->SetLineStyle(0);
  legS->SetTextSize(0.03);
  legS->SetFillStyle(1001);
  legS->SetFillColor(kWhite);
  legS->AddEntry(gr0_sig,"Median","lp");
  legS->AddEntry(gr1_sig,"1#sigma band ","f");
  legS->AddEntry(gr2_sig,"2#sigma band ","f");
  legS->AddEntry(line,"3#sigma level ","l");
  legS->Draw();
  
  
  pt = new TPaveText(0.23,0.95,0.35,0.91, "NDC"); // NDC sets coords
  // relative to pad dimensions
  pt->SetFillColor(0); // text is black on white
  pt->SetTextSize(0.03); 
  pt->SetTextAlign(12);
  pt->SetBorderSize(0);
  pt->SetLineStyle(0);
  char string[255];
  sprintf(string,"%d ton mass detector", mass_exposure);
  text = pt->AddText(string);
  pt->Draw();       //to draw your text object
    
}

//-------------------------------------

void MakeModel(double nsig,                  // number of signal events 
	       double nbkg, double sigmab)   // number of background events
{

  Double_t pi = TMath::Pi();
  // Declare observable x
  RooRealVar the("the","the",-pi/2.,pi/2.);
  RooRealVar len("len","len",0.05,1);
  RooRealVar phi("phi","phi",-pi,pi);
  the.setBins(10);
  len.setBins(10);
  phi.setBins(10);
  /*
  TFile *file = new TFile("allhisto_v9_50nm.root");
  TH1F* hbthe = (TH1F*)file->Get("theta_B");
  TH1F* hbphi = (TH1F*)file->Get("phi_boro");
  TH1F* hblen = (TH1F*)file->Get("L_thr");
  TH1F* hsnthe = (TH1F*)file->Get("theta_v4");
  TH1F* hsnphi = (TH1F*)file->Get("phiSN_v2");
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

 
  // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
  RooDataHist dbthe("dbthe","dbthe",the,Import(*hbthe));
  RooDataHist dbphi("dbphi","dbphi",phi,Import(*hbphi));
  RooDataHist dblen("dblen","dblen",len,Import(*hblen));
  RooDataHist dsnthe("dsnthe","dsnthe",the,Import(*hsnthe));
  RooDataHist dsnphi("dsnphi","dsnphi",phi,Import(*hsnphi));
  RooDataHist dsnlen("dsnlen","dsnlen",len,Import(*hsnlen));
  RooDataHist dnthe("dnthe","dnthe",the,Import(*hnthe));
  RooDataHist dnphi("dnphi","dnphi",phi,Import(*hnphi));
  RooDataHist dnlen("dnlen","dnlen",len,Import(*hnlen));

  // Represent data in dh as pdf in x, apply 2nd order interpolation  
  RooHistPdf bthepdf("bthepdf","bthepdf",the,dbthe,2);
  RooHistPdf bphipdf("bphipdf","bphipdf",phi,dbphi,2);
  RooHistPdf blenpdf("blenpdf","blenpdf",len,dblen,2);
  RooHistPdf snthepdf("snthepdf","snthepdf",the,dsnthe,2);
  RooHistPdf snphipdf("snphipdf","snphipdf",phi,dsnphi,2);
  RooHistPdf snlenpdf("snlenpdf","snlenpdf",len,dsnlen,2);
  RooHistPdf nthepdf("nthepdf","nthepdf",the,dnthe,2);
  RooHistPdf nphipdf("nphipdf","nphipdf",phi,dnphi,2);
  RooHistPdf nlenpdf("nlenpdf","nlenpdf",len,dnlen,2);


  // Make plot of binned dataset showing Poisson error bars (RooFit default)
  //RooPlot* frame = the.frame(Title("Imported TH1 with Poisson error bars")) ;
  //dbthe.plotOn(frame); 
  //frame->Draw();

  w.import(nthepdf);
  w.import(bthepdf);
  w.import(snthepdf);
  w.import(nphipdf);
  w.import(bphipdf);
  w.import(snphipdf);
  w.import(nlenpdf);
  w.import(blenpdf);
  w.import(snlenpdf);

  //w.factory("SUM:bkg_the(0.21*nthepdf,0.79*bthepdf)");
  //w.factory("SUM:bkg_phi(0.21*nphipdf,0.79*bphipdf)");
  //w.factory("SUM:bkg_len(0.21*nlenpdf,0.79*blenpdf)");

  w.factory("SUM:bkg_the(frac[0.20]*nthepdf,bthepdf)");
  w.factory("SUM:bkg_phi(frac[0.20]*nphipdf,bphipdf)");
  w.factory("SUM:bkg_len(frac[0.20]*nlenpdf,blenpdf)");

  w.factory("PROD:sig_pdf(snthepdf,snphipdf,snlenpdf)");
  w.factory("PROD:bkg_pdf(bkg_the,bkg_phi,bkg_len)");

  //w.factory("PROD:sig_pdf(snphipdf)");
  //w.factory("PROD:bkg_pdf(bkg_phi)");

  //w.factory("PROD:sig_pdf(snlenpdf,snthepdf)");
  //w.factory("PROD:bkg_pdf(bkg_len,bkg_the)");
  
  w.factory("SUM:model(nsig[0,1000]*sig_pdf, nbkg[0,1000]*bkg_pdf)");  // for extended model

  // w.factory("SUM:model(nsig[0,1000]*snlenpdf, nbkg[0,1000]*bkg_len)");  // for extended model


  RooAbsPdf * pdf = w.pdf("model");
  w.var("nbkg")->setConstant(1);
  w.var("nsig")->setVal(nsig);
  w.var("nbkg")->setVal(nbkg);
  
  // use fixed random numbers for reproducibility (use 0 for changing every time)
  RooRandom::randomGenerator()->SetSeed(0);

  //RooDataSet *data = pdf->generate(RooArgSet(phi),NumEvents(TMath::Nint(nbkg+nsig)), Extended(1));  // will generate accordint to total S+B events
  RooDataSet *data = pdf->generate(RooArgSet(the,phi,len),NumEvents(TMath::Nint(nbkg+nsig)), Extended(1));  // will generate accordint to total S+B events
  //RooDataSet *data = pdf->generate(RooArgSet(len,the),NumEvents(TMath::Nint(nbkg+nsig)), Extended(1));  // will generate accordint to total S+B events

  
  data->SetName("data");
  w.import(*data);
  /* 
  RooPlot* the_frame = the.frame() ;
  data->plotOn(the_frame);
  pdf->plotOn(the_frame);
  pdf->paramOn(the_frame,Layout(0.55));
  data->statOn(the_frame,Layout(0.55,0.99,0.8));
  
  TCanvas *cthe = new TCanvas("cthe","cthe",600,600);
  the_frame->Draw();

  RooPlot* phi_frame = phi.frame() ;
  data->plotOn(phi_frame);
  pdf->plotOn(phi_frame);
  pdf->paramOn(phi_frame,Layout(0.55));
  data->statOn(phi_frame,Layout(0.55,0.99,0.8));
  
  TCanvas *cphi = new TCanvas("cphi","cphi",600,600);
  phi_frame->Draw();

  RooPlot* len_frame = len.frame() ;
  data->plotOn(len_frame);
  pdf->plotOn(len_frame);
  pdf->paramOn(len_frame,Layout(0.55));
  data->statOn(len_frame,Layout(0.55,0.99,0.8));
  
  TCanvas *clen = new TCanvas("clen","clen",600,600);
  len_frame->Draw();*/

  data->Print();
  
  RooStats::ModelConfig mc("ModelConfig",&w);
  mc.SetPdf(*pdf);
  mc.SetParametersOfInterest(*w.var("nsig"));
  mc.SetObservables(*w.var("the"));
  mc.SetObservables(*w.var("phi"));
  mc.SetObservables(*w.var("len"));
  /* mc.SetNuisanceParameters(*w.var("nbkg"));
     
  // these are needed for the hypothesis tests
  mc.SetSnapshot(*w.var("nsig"));
  mc.SetGlobalObservables(*w.var("b0"));*/
  
  // import model in the workspace 
  w.import(mc);
}

//-------------------------------------------------

double MakeSignificance()
{

  // get the data  out of the file
  RooAbsData* data = w.data("data");

  printf("ECCOCI 3\n");
    
  // get the modelConfig (S+B) out of the file
  // and create the B model from the S+B model
  ModelConfig*  sbModel = (RooStats::ModelConfig*) w.obj("ModelConfig");
  sbModel->SetName("S+B Model");      
  RooRealVar* poi = (RooRealVar*) sbModel->GetParametersOfInterest()->first();
  //poi->setVal(50);
  //poi->setRange(0,50.);  // set POI snapshot in S+B model for expected significance
  sbModel->SetSnapshot(*poi);
  ModelConfig * bModel = (ModelConfig*) sbModel->Clone();
  bModel->SetName("B Model");      
  poi->setVal(0);
  bModel->SetSnapshot( *poi  );

  // create the AsymptoticCalculator from data,alt model, null model
  AsymptoticCalculator  ac(*data, *sbModel, *bModel);
  ac.SetOneSidedDiscovery(true);  // for one-side discovery test
  ac.SetPrintLevel(-1);  // to suppress print level 

  // run the calculator
  HypoTestResult * asResult = ac.GetHypoTest();
  //cout << asResult->NullPValue() << endl;
  //cout << asResult->Significance() << endl;

  double Significance = asResult->Significance();
  return Significance;
}


//-------------------------------------------------------------







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
  return h->Integral( 1, bin ) / h->Integral();
}

// Return the fraction of the area below the histogram with x>x_value
// the bin that contains x_value is not included in the fraction.
//===================================================
double get_fraction_above(  TH1F* h, double x_value ){
//===================================================
  return 1 - get_fraction_below( h, x_value);
}
