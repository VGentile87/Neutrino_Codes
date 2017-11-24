#include <stdio.h>
#include <ctime>
#include <cstdlib>
#include <string>
#include "TMath.h"
#include "TString.h"
#include "TMultiGraph.h"

//#include "Definitions_sn2.h"

using namespace std;

class Functions
{
 public: 

  //inline void Initialize(){}

  inline double cross_section(){

    iN_el=iA_el-iZ_el;
    //cout << iN_el << endl;
    iQw = iN_el - (1-4*sin2wnb)*iZ_el;
    irn = TMath::Sqrt(1.5129*TMath::Power(iA_el,2./3.) -1.4760*TMath::Power(iA_el,1./3.) + 2.5371);

    q = TMath::Sqrt(2*iM_el*(E_nu*E_nu*TMath::Power(10,-6)*(1-cos_theta))/(E_nu*TMath::Power(10,-3)*(1-cos_theta)+iM_el));
    qrn = q*irn*fm_in_GeV;
    qs = q*s*fm_in_GeV;
    Fq = (3*(TMath::Sin(qrn)-qrn*TMath::Cos(qrn))/TMath::Power(qrn,3))*TMath::Exp(-(qs*qs)/2.);  // diverso rispetto a plot di esclusione
    crsect = ((GF*GF*iQw*iQw*E_nu*E_nu*Fq*Fq*(1+cos_theta)*(2./d_theta))/(8*pi))*TMath::Power(iMeV_in_fm,2)*TMath::Power(10,-26);

    return crsect;
  }

  inline double num_eventi_per_energia(){

    int_const = (M_riv*imass_fract)/(iA_el*uma_in_ton*4*pi*dist_sn*dist_sn*kpc_in_cm*kpc_in_cm);
    intE_nu_e = int_const*((N_nu_e*TMath::Power(beta_nu_e,3))/2.)*E_nu*E_nu*crsect*TMath::Exp(-beta_nu_e*E_nu);
    intE_antinu_e = int_const*((N_antinu_e*TMath::Power(beta_antinu_e,3))/2.)*E_nu*E_nu*crsect*TMath::Exp(-beta_antinu_e*E_nu);
    intE_nu_x = int_const*((4*N_nu_x*TMath::Power(beta_nu_x,3))/2.)*E_nu*E_nu*crsect*TMath::Exp(-beta_nu_x*E_nu);
    dN_vs_dE = (intE_nu_e + intE_antinu_e + intE_nu_x);

    return dN_vs_dE;
  }
  
  inline double err_num_eventi_per_energia(){

    int_const = (M_riv*imass_fract)/(iA_el*uma_in_ton*4*pi*dist_sn*dist_sn*kpc_in_cm*kpc_in_cm);
    intE_nu_e = int_const*((N_nu_e*TMath::Power(beta_nu_e,3))/2.)*E_nu*E_nu*crsect*TMath::Exp(-beta_nu_e*E_nu)*0.2248*TMath::Sqrt(1+7.61*TMath::Power(E_nu/24.82 - 1,2));
    intE_antinu_e = int_const*((N_antinu_e*TMath::Power(beta_antinu_e,3))/2.)*E_nu*E_nu*crsect*TMath::Exp(-beta_antinu_e*E_nu)*0.2248*TMath::Sqrt(1+7.61*TMath::Power(E_nu/24.82 - 1,2));
    intE_nu_x = int_const*((4*N_nu_x*TMath::Power(beta_nu_x,3))/2.)*E_nu*E_nu*crsect*TMath::Exp(-beta_nu_x*E_nu)*0.2248*TMath::Sqrt(1+7.61*TMath::Power(E_nu/24.82 - 1,2));

    err_num_ev = (intE_nu_e + intE_antinu_e + intE_nu_x);
    //err_num_ev = dN_vs_dE*0.2248*TMath::Sqrt(1+7.61*TMath::Power(E_nu/24.82 - 1,2));

    return err_num_ev;
  }


  inline double num_eventi_bkg(){
    
    int_const_bkg = (M_riv*imass_fract)/(iA_el*uma_in_ton);
    dN_vs_dE_bkg = int_const_bkg*crsect*iflux_hep_vs_E;    
    return dN_vs_dE_bkg;
  }
  
  inline TCanvas* gr_crsect(){

    cmg = new TCanvas("cmg","cmg",600,600);

    for(int i=0;i<n_el;i++){
      gr_sigma[i]->SetTitle(el[i]);
      gr_sigma[i]->SetLineWidth(2);
      gr_sigma[i]->SetLineColor(i+1);
      //gr_sigma[i]->SetMarkerStyle(22);
      gr_sigma[i]->SetFillColor(0);
      mg->Add(gr_sigma[i]);
    }
    
    cmg->cd();
    cmg->SetLogy();
    mg->Draw("AL");
    mg->GetXaxis()->SetTitle("n");
    mg->GetYaxis()->SetTitle("#pi");
    cmg->BuildLegend();

    return cmg;
    
  }

  inline TCanvas* trk_len(){

  cc7 = new TCanvas("cc7","cc7",600,600);
  h7->Draw("");
  cc7->Update();
  TPaveStats *cc7statsa =(TPaveStats*)cc7->GetPrimitive("stats");
  cc7statsa->SetName(""); 
  cc7statsa->SetTextColor(1);
  hh3->Scale(num_ev/h7->Integral());
  hh3->Draw("sames");
  cc7->Update();
  TPaveStats *cc7statsaa =(TPaveStats*)cc7->GetPrimitive("stats");
  cc7statsaa->SetName(""); 
  cc7statsaa->SetTextColor(3);
  hh0->Scale(num_ev/h7->Integral());
  hh0->Draw("sames");
  cc7->Update();
  TPaveStats *cc7statsb =(TPaveStats*)cc7->GetPrimitive("stats");
  cc7statsb->SetName(""); 
  cc7statsb->SetTextColor(6);
  hh1->Scale(num_ev/h7->Integral());
  hh1->Draw("sames");
  cc7->Update();
  TPaveStats *cc7statsbb =(TPaveStats*)cc7->GetPrimitive("stats");
  cc7statsbb->SetName(""); 
  cc7statsbb->SetTextColor(8);
  hh2->Scale(num_ev/h7->Integral());
  hh2->Draw("sames");
  cc7->Update();
  TPaveStats *cc7statsc =(TPaveStats*)cc7->GetPrimitive("stats");
  cc7statsc->SetName(""); 
  cc7statsc->SetTextColor(42);
  hh4->Scale(num_ev/h7->Integral());
  hh4->Draw("sames");
  cc7->Update();
  TPaveStats *cc7statsd =(TPaveStats*)cc7->GetPrimitive("stats");
  cc7statsd->SetName(""); 
  cc7statsd->SetTextColor(2);
  hh5->Scale(num_ev/h7->Integral());
  hh5->Draw("sames");
  cc7->Update();
  TPaveStats *cc7statse =(TPaveStats*)cc7->GetPrimitive("stats");
  cc7statse->SetName(""); 
  cc7statse->SetTextColor(8);
  hh6->Scale(num_ev/h7->Integral());
  hh6->Draw("sames");
  cc7->Update();
  //h7->Scale(1./h7->Integral());
  TPaveStats *cc7statsf =(TPaveStats*)cc7->GetPrimitive("stats");
  cc7statsf->SetName(""); 
  cc7statsf->SetTextColor(4);
  h7->SetLineColor(1);
  hh3->SetLineColor(41);
  hh0->SetLineColor(2);
  hh1->SetLineColor(4);
  hh2->SetLineColor(3);
  hh4->SetLineColor(7);
  hh5->SetLineColor(8);
  hh6->SetLineColor(6);
  h7->SetLineWidth(2);
  hh3->SetLineWidth(2);
  hh0->SetLineWidth(2);
  hh1->SetLineWidth(2);
  hh2->SetLineWidth(2);
  hh4->SetLineWidth(2);
  hh5->SetLineWidth(2);
  hh6->SetLineWidth(2);
  TLegend * legP = new TLegend();
  legP = new TLegend(0.15,0.60,0.4,0.80);
  legP->SetTextSize(0.03);
  legP->SetBorderSize(0);
  legP->SetLineStyle(0);
  legP->SetTextSize(0.03);
  legP->SetFillStyle(1001);
  legP->SetFillColor(kWhite);
  legP->AddEntry(h7,"all targets","f");
  legP->AddEntry(hh0,"Ag","f");
  legP->AddEntry(hh1,"Br","f");
  legP->AddEntry(hh2,"I","f");
  legP->AddEntry(hh3,"C","f");
  legP->AddEntry(hh4,"O","f");
  legP->AddEntry(hh5,"N","f");
  legP->AddEntry(hh6,"S","f");
  legP->Draw();

  return cc7;
  
    }
  
 protected:

 private:
  
};



