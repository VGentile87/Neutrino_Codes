//
// A. Di Crescenzo (November 2015)
//
// Functions implemented:
//
// DrawCorrelation()
// GetEnergyFromRange(charge, range_nm)
// GetRangeFromEnergy(charge, ene_kev)
//

#include "TLegend.h"
#include "TSpline.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "TCanvas.h"
#include "TMath.h"
#include "TH1F.h"

using namespace std;


TCanvas *vC1,*vC2;
TGraph *gr;

//----------------------------------

TSpline3 *RangeVersusEnergy(Int_t charge){

  ifstream inputFile("DATA/H.txt");

  if(charge==1)
    std::ifstream inputFile("DATA/H.txt");
  if(charge==6)
    std::ifstream inputFile("DATA/C.txt");
  if(charge==7)
    std::ifstream inputFile("DATA/N.txt");
  if(charge==8)
    std::ifstream inputFile("DATA/O.txt");
  if(charge==16)
    std::ifstream inputFile("DATA/S.txt");
  if(charge==35)
    std::ifstream inputFile("DATA/Br.txt");
  if(charge==47)
    std::ifstream inputFile("DATA/Ag.txt");
  if(charge==53)
    std::ifstream inputFile("DATA/I.txt");
    
  int num_line=0;
  int index_vector=0;
  double mass_wimp_vec[500];
  double ene_vec[500];
  double len_vec[500];
  double mass_wimp=0;
  double vel_mean_wimp = 0.75*TMath::Power(10,-3); // [c]
  double amu_in_mevc2=931.494;
  double mass_element;
  
  std::string line;

  Double_t x[79];
  Double_t y[79];

  Int_t Narray = 0;

  while(getline(inputFile, line)) {
    num_line++;
    if (!line.length() || line[0] == '#')
      continue;
    std::istringstream iss(line);
    double ene = 0.;
    char ene_unit[3] = {};
    double fake=0.;
    double fake2=0;
    double len=0.;
    char len_unit[3] = {};

    if(num_line<=80){

      iss>>ene>>ene_unit>>fake>>fake2>>len>>len_unit;

      if(ene_unit[0]=='k')ene=ene; //GeV
      if(ene_unit[0]=='M')ene*=1000;  //GeV
      if(len_unit[0]=='A')len/=10; // nm
      if(len_unit[0]=='u')len*=1000; // nm
      //std::cout<<"point: "<<ene<< " keV " << len<<' '<< " nm " <<std::endl;
      //printf("i %d -> ENE %4.4f LEN %4.4f\n",Narray,ene,len);
      
      y[Narray] = len;
      x[Narray] = ene;
      Narray++;
    }
  }

  gr = new TGraph(Narray,x,y);
  //vC1 = new TCanvas("vC1","square",200,10,700,700);

  TSpline3 *s3 = new TSpline3("s3",gr->GetX(),gr->GetY(),gr->GetN());
  //s3->SetLineColor(kRed);
  //gr->SetMarkerColor(kBlue);
  //gr->SetMarkerStyle(21);
  //gr->SetMarkerSize(0.5);
  //gr->Draw("alp");
  //s3->Draw("l same");

  return s3;
}


//----------------------------------

TSpline3 *EnergyVersusRange(Int_t charge){

  ifstream inputFile("DATA/H.txt");
  
  if(charge==1)
    std::ifstream inputFile("DATA/H.txt");
  if(charge==6)
    std::ifstream inputFile("DATA/C.txt");
  if(charge==7)
    std::ifstream inputFile("DATA/N.txt");
  if(charge==8)
    std::ifstream inputFile("DATA/O.txt");
  if(charge==16)
    std::ifstream inputFile("DATA/S.txt");
  if(charge==35)
    std::ifstream inputFile("DATA/Br.txt");
  if(charge==47)
    std::ifstream inputFile("DATA/Ag.txt");
  if(charge==53)
    std::ifstream inputFile("DATA/I.txt");

  
  int num_line=0;
  int index_vector=0;
  double mass_wimp_vec[500];
  double ene_vec[500];
  double len_vec[500];
  double mass_wimp=0;
  double vel_mean_wimp = 0.75*TMath::Power(10,-3); // [c]
  double amu_in_mevc2=931.494;
  double mass_element;
  
  std::string line;

  Double_t x[79];
  Double_t y[79];

  Int_t Narray = 0;

  while(getline(inputFile, line)) {
    num_line++;
    if (!line.length() || line[0] == '#')
      continue;
    std::istringstream iss(line);
    double ene = 0.;
    char ene_unit[3] = {};
    double fake=0.;
    double fake2=0;
    double len=0.;
    char len_unit[3] = {};

    if(num_line<=80){
      
      iss>>ene>>ene_unit>>fake>>fake2>>len>>len_unit;

      if(ene_unit[0]=='k')ene=ene; //GeV
      if(ene_unit[0]=='M')ene*=1000;  //GeV
      if(len_unit[0]=='A')len/=10; // nm
      if(len_unit[0]=='u')len*=1000; // nm
      //std::cout<<"point: "<<ene<< " keV " << len<<' '<< " nm " <<std::endl;
      //printf("i %d -> ENE %4.4f LEN %4.4f\n",Narray,ene,len);
      
      x[Narray] = len;
      y[Narray] = ene;
      Narray++;
    }
  }

  gr = new TGraph(Narray,x,y);
  //vC1 = new TCanvas("vC1","square",200,10,700,700);

  TSpline3 *s3 = new TSpline3("s3",gr->GetX(),gr->GetY(),gr->GetN());
  //s3->SetLineColor(kRed);
  //gr->SetMarkerColor(kBlue);
  //gr->SetMarkerStyle(21);
  //gr->SetMarkerSize(0.5);
  //gr->Draw("alp");
  //s3->Draw("l same");

  return s3;
}


//--------------------------------

void DrawCorrelation(){
  TSpline3 *Line_H  = EnergyVersusRange(1);
  TSpline3 *Line_C  = EnergyVersusRange(6);
  TSpline3 *Line_N  = EnergyVersusRange(7);
  TSpline3 *Line_O  = EnergyVersusRange(8);
  TSpline3 *Line_Br = EnergyVersusRange(35);
  TSpline3 *Line_Ag = EnergyVersusRange(47);
  
  vC2 = new TCanvas("vC2","square",200,10,700,700);
  vC2->cd();

  TH1F *vFrame = gPad->DrawFrame(0,0,3000,4000);
  vFrame->Draw();
  vFrame->GetXaxis()->SetTitleOffset(1.3);
  vFrame->GetYaxis()->SetTitleOffset(1.6);
  vFrame->SetXTitle("Range (nm)");
  vFrame->SetYTitle("Energy (keV)");
  
  Line_H->SetLineColor(kCyan);
  Line_H->SetLineWidth(2.0);
  Line_C->SetLineColor(kBlue);
  Line_C->SetLineWidth(2.0);
  Line_N->SetLineColor(kMagenta);
  Line_N->SetLineWidth(2.0);
  Line_O->SetLineColor(kGreen);
  Line_O->SetLineWidth(2.0);
  Line_Br->SetLineColor(kRed);
  Line_Br->SetLineWidth(2.0);
  Line_Ag->SetLineColor(kBlack);
  Line_Ag->SetLineWidth(2.0);

  Line_O->Draw("l same");
  Line_Br->Draw("l same");
  Line_N->Draw("l same");
  Line_C->Draw("l same");
  Line_Ag->Draw("l same");
  Line_H->Draw("l same");
    
  TLegend * legE = new TLegend();
  legE = new TLegend(0.70,0.60,1.0,0.80);
  legE->SetTextSize(0.03);
  legE->SetBorderSize(0);
  legE->SetLineStyle(0);
  legE->SetTextSize(0.03);
  legE->SetFillStyle(0);
  legE->SetFillColor(0);
  legE->AddEntry(Line_Ag,"Ag","l");
  legE->AddEntry(Line_Br,"Br","l");
  legE->AddEntry(Line_N,"N","l");
  legE->AddEntry(Line_O,"O","l");
  legE->AddEntry(Line_C,"C","l");
  legE->AddEntry(Line_H,"H","l");
  legE->Draw();
  
}

//--------------------------------

Float_t GetEnergyFromRange(Int_t charge, Float_t range_nm){

  TSpline3 *EneVsRange = EnergyVersusRange(charge);

  Float_t ene_kev = 0;
  ene_kev = EneVsRange->Eval(range_nm);
  
  return ene_kev;
}

//-------------------------------

Float_t GetRangeFromEnergy(Int_t charge, Float_t ene_kev){

  TSpline3 *RangeVsEne = RangeVersusEnergy(charge);

  Float_t range_nm = 0;
  range_nm = RangeVsEne->Eval(ene_kev);
  //printf("CHARGE %d, range %4.1f, energy %4.1f\n",charge,range_nm,ene_kev);

  return range_nm;
}
