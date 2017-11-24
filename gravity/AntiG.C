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
#include "TRandom3.h"



using namespace std;

void AntiG() {

  //int opzioni = 0;

  //cout << "Test di Controllo (tesi Gallo Rosso): opzione '0'"<<endl;
  //cout << "SN neutrinos in NEWS (V. Gentile): opzione '1'" <<endl;
  //cout <<"\n Scegli opzione: ";
  //cin << opzioni;
  //start(opzioni);

  ofstream log_file;
  //log_file.open ("num_eventi.txt");
  //log_file << "SN Neutrino induced recoil rate";
  //log_file << endl;

  TGraph *gr_gamma = new TGraph();
  TGraph *gr_dx = new TGraph();
  TGraph *gr_dt = new TGraph();
  
  const double c=1;
  double v=0;
  double t=0;
  double E_kin=0;
  double sign=0;
  double m=1;
  double new_gamma=0;
  double dt0=1;
  double dx0=1;
  double dt=0;
  double dx=0;

  for(int i=0;i<100000;i++){
    if(v<=c)sign=1;
    else sign=-1;
    new_gamma = 1./TMath::Sqrt((1-(v*v)/(c*c))*sign);
    //E_kin=new_gamma*m*c*c;
    //x_pr = (x-v*t)*new_gamma;
    gr_gamma->SetPoint(i,v,new_gamma);
    dt = new_gamma*dt0;
    dx = dx0/new_gamma;
    gr_dx->SetPoint(i,v,dx);
    gr_dt->SetPoint(i,v,dt);
    v=v+0.0001;
    t=t+0.0001;
    //cout << i << " " << new_gamma << endl;
  }

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  gr_gamma->Draw("A*");
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  gr_dx->Draw("A*");
  TCanvas *c3 = new TCanvas("c3","c3",600,600);
  gr_dt->Draw("A*");

  
  /*
  TCanvas *cthe = new TCanvas("cthe","cthe",600,600);
  hthe->Draw("");
  h2->Draw("same");
  hthe->SetLineColor(1);
  h2->SetLineColor(2);
  */
  

}
