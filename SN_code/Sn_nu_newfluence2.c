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

#include "Definitions_sn2.h"
#include "Functions.h"
#include "EnergyRangeCorrelation.C"
#include "TRandom3.h"
#include <TSystem.h>



using namespace std;

void Sn_nu_newfluence2() {

  Double_t threshold;
  cout << "Inserisci soglia [nm]" << endl;
  cin  >> threshold;

  ofstream log_file;
  log_file.open ("num_eventi.txt");
  log_file << "SN Neutrino induced recoil rate";
  log_file << endl; 

  TRandom3 *tt = new TRandom3();
  Int_t index=0;
  Functions* functions;

  functions = new Functions();
  
  gStyle->SetOptStat("eMR");
  //gStyle->SetStatW(0.2);
  //gStyle->SetStatH(0.15);

  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);                
  gStyle->SetStatW(0.3);                
  gStyle->SetStatH(0.15);                
  

  gr_spect[0]=new TGraph();
  gr_spect[1]=new TGraph();
  gr_spect[2]=new TGraph();

  h4->Sumw2();
  h44->Sumw2();
  h2->Sumw2();
  h22->Sumw2();

  TRandom3 *r = new TRandom3();
  Double_t t = 0;
  TRandom3 *rr = new TRandom3();
  Double_t a = 0;

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
    h11[i] = new TH1D(Form("phi_el_%d",i),"",100,-pi,pi);
    iel = el[i];
    iZ_el = Z_el[i];
    iA_el = A_el[i];
    ifract_el = atomic_fract_el[i]; // uma
    imass_fract = mass_fract[i];
    iM_el=iA_el*uma;
    iN_el = N_el[i];
    num_ev_el=0;
    err_num_ev_el=0;
    gr_sigma[i]=new TGraph();
    sum_E=0;
    index=0;
    
    for(int j=1; j<E_nu_max;j++){
	E_nu=j;
	sum_crsect=0;

	spect_nu_e = /*TMath::Power(10,15)**/(erg_in_MeV/(4*pi*dist_sn*dist_sn*kpc_in_fm*kpc_in_fm))*(tot_ene_e/mean_ene_e)*(TMath::Power(E_nu,alpha_e)/TMath::Gamma(alpha_e+1))*TMath::Power(((alpha_e+1)/mean_ene_e),alpha_e+1)*TMath::Exp(-E_nu*(alpha_e+1)/mean_ene_e);
	spect_antinu_e = /*TMath::Power(10,15)**/(erg_in_MeV/(4*pi*dist_sn*dist_sn*kpc_in_fm*kpc_in_fm))*(tot_ene_antie/mean_ene_antie)*(TMath::Power(E_nu,alpha_antie)/TMath::Gamma(alpha_antie+1))*TMath::Power(((alpha_antie+1)/mean_ene_antie),alpha_antie+1)*TMath::Exp(-E_nu*(alpha_antie+1)/mean_ene_antie);
	spect_nu_x = /*TMath::Power(10,15)**/(erg_in_MeV/(4*pi*dist_sn*dist_sn*kpc_in_fm*kpc_in_fm))*(tot_ene_x/mean_ene_x)*(TMath::Power(E_nu,alpha_x)/TMath::Gamma(alpha_x+1))*TMath::Power(((alpha_x+1)/mean_ene_x),alpha_x+1)*TMath::Exp(-E_nu*(alpha_x+1)/mean_ene_x);
		
      if(i==0){
	gr_spect[0]->SetPoint(j-1,E_nu,spect_nu_e);
	gr_spect[1]->SetPoint(j-1,E_nu,spect_antinu_e);
	gr_spect[2]->SetPoint(j-1,E_nu,spect_nu_x);
      }
      
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
	  
	  px_cm = E_nu_cm*TMath::Sqrt((1-TMath::Power(rnd_cos_theta_cm,2)))*TMath::Sin(rnd_phi_cm);  // punta al centro galattico
	  py_cm = E_nu_cm*TMath::Sqrt((1-TMath::Power(rnd_cos_theta_cm,2)))*TMath::Cos(rnd_phi_cm);  // punta alla direzione del Cigno
	  pz_cm = E_nu_cm*rnd_cos_theta_cm;                                                          // normale al piano galattico
	  p_mod_cm = TMath::Sqrt(px_cm*px_cm + py_cm*py_cm + pz_cm*pz_cm);
	  
	  beta_cm = E_nu/(E_nu + iM_el*1000);
	  gamma_cm = (E_nu + iM_el*1000)/(TMath::Sqrt(iM_el*iM_el*1000000 + 2*E_nu*iM_el*1000));
	  px_lab = beta_cm*gamma_cm*E_rec_cm + gamma_cm*px_cm;
	  E_tot_lab = gamma_cm*E_rec_cm + beta_cm*gamma_cm*py_cm;  //MeV
	  E_rec_lab = gamma_cm*E_rec_cm + beta_cm*gamma_cm*py_cm - iM_el*1000; //MeV
	  p_mod_lab = TMath::Sqrt(E_tot_lab*E_tot_lab - iM_el*iM_el*1000000);
	  the_rec = TMath::ATan2(TMath::Sqrt((1-TMath::Power(rnd_cos_theta_cm,2))),(gamma_cm*(beta_cm*(E_rec_cm/p_mod_cm)+rnd_cos_theta_cm)));  //buono
	  //phi_rec = TMath::ATan2(px_lab,px_cm);
	  phi_rec = TMath::ATan2(py_cm,px_lab);
	  
	  hpx_cm->Fill(px_cm);
	  hpy_cm->Fill(py_cm);
	  hpz_cm->Fill(pz_cm);
	  hpx_lab->Fill(px_lab);
	  ///////////////////////////////////////////////////////////////////////////////////////////////

	  h2->Fill(TMath::ACos(costheta_sc));
	  h22->Fill(the_rec);
	  h4->Fill(E_rec*1000);
	  h44->Fill(E_rec_lab*1000);
	  /*
	  if(i==0)henelenAg->Fill(E_rec,L_true/1000.);
	  if(i==1)henelenBr->Fill(E_rec,L_true/1000.);
	  if(i==3)henelenC->Fill(E_rec,L_true/1000.);
	  if(i==4)henelenO->Fill(E_rec,L_true/1000.);
	  if(i==5)henelenN->Fill(E_rec,L_true/1000.);*/

	  if(L_true>=threshold && L_true<=1000){
	    
	    num_ev_el += functions->num_eventi_per_energia_newfluence();
	    err_num_ev_el += functions->num_eventi_per_energia()*0.2248*TMath::Sqrt(1+7.61*TMath::Power(E_nu/24.82 - 1,2));
	    num_ev += functions->num_eventi_per_energia_newfluence();
	    err_num_ev += functions->num_eventi_per_energia()*0.2248*TMath::Sqrt(1+7.61*TMath::Power(E_nu/24.82 - 1,2));
	    
	    h1[i]->Fill(costheta_sc);
	    h11[i]->Fill(phi_rec);
	    hthe_thr->Fill(the_rec);
	    hcosthe_thr->Fill(TMath::Cos(the_rec));
	    hphi_thr->Fill(phi_rec);
	    hlen_thr->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);
	    hene_thr->Fill(E_rec*1000);

	    if(i==0)hh0->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);
	    if(i==1)hh1->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);
	    if(i==2)hh2->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);
	    if(i==3)hh3->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);
	    if(i==4)hh4->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);
	    if(i==5)hh5->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);
	    if(i==6)hh6->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);	    
	  }
	  ///////////////////////////////////	
	}
	gr_sigma[i]->SetPoint(index,E_nu,sum_crsect);
	index++;
    }

    log_file << "Threshold = " << threshold << " nm " << el[i] << "\t " << "Events= " << num_ev_el << " " << "Error= " << err_num_ev_el << endl;
    cout << i << " " << el[i] << " " <<  num_ev << " " << err_num_ev << " " << sum_crsect <<  " " << num_ev_miss <<  " " << " " << num_ev_el << endl;

for(int l=0;l<150;l++){
  //a = rr->Uniform(0,0.26);
  costheta_sc = h1[i]->GetRandom();
  a = TMath::ACos(costheta_sc);
    
    for(int k=0;k<150;k++){
      if(h1[i]->GetEntries()!=0){
	t = r->Uniform(-pi,pi);
	
	//costheta_sc = h1[i]->GetRandom();
        phi_sc = h11[i]->GetRandom();
	theta_sc = TMath::ACos(costheta_sc);
	pphi=-pi/2.;  // nu SN DAL CENTRO GALATTICO
	pthe=0.;  // nu SN NEL PIANO GALATTICO
	////////////////////////////////////////////////////
	//px_sum = TMath::Cos(phi_sc)*TMath::Cos(theta_sc) + TMath::Cos(pphi)*TMath::Cos(pthe);
	//py_sum = TMath::Sin(phi_sc)*TMath::Cos(theta_sc) + TMath::Sin(pphi)*TMath::Cos(pthe);
	//pz_sum = TMath::Sin(theta_sc) + TMath::Sin(pthe);
	//theta_sum  = TMath::ATan(pz_sum/TMath::Sqrt(px_sum*px_sum + py_sum*py_sum));
	//phi_sum = TMath::ATan2(py_sum,px_sum);

	/////////////////////////////////////////////////////
	theta_sum = theta_sc + pthe;
	//if(t<0)theta_sum=-theta_sum;
	phi_sum = phi_sc + pphi;
	///////////////////////////////////////////////////////////////////
	cos_theta_sum = TMath::Cos(theta_sum); 
	h3->Fill(theta_sum,atomic_fract_el[i]*A_el[i]*A_el[i]);
	h33->Fill(cos_theta_sum,atomic_fract_el[i]*A_el[i]*A_el[i]);
	h5->Fill(phi_sum,atomic_fract_el[i]*A_el[i]*A_el[i]);
	h55->Fill(phi_sc,atomic_fract_el[i]*A_el[i]*A_el[i]);

	/*yrec = -TMath::Cos(theta_sc);
	xrec = TMath::Sin(theta_sc)*TMath::Cos(t);
	zrec = TMath::Sin(theta_sc)*TMath::Sin(t);*/
	//cout << t << " " << a <<endl;
	yrec = -TMath::Cos(a);
	xrec = TMath::Sin(a)*TMath::Cos(t);
	zrec = TMath::Sin(a)*TMath::Sin(t);
	hxyzrec->Fill(xrec,yrec,zrec);
	therec = TMath::ATan(zrec/TMath::Sqrt(TMath::Power(xrec,2)+TMath::Power(yrec,2)));
	phirec = TMath::ATan2(yrec,xrec);
	hphirec->Fill(phirec,atomic_fract_el[i]*A_el[i]*A_el[i]);
	htherec->Fill(therec,atomic_fract_el[i]*A_el[i]*A_el[i]);
	hcostherec->Fill(TMath::Sin(therec));
      }
    }
 }
        
  }

  /*
  TCanvas *cc2 = new TCanvas("cc2","cc2",1200,600);
  //hthe->Draw();
  cc2->Divide(2,1);
  cc2->cd(1);
  //7hr->Draw("");
  //cc2->cd(2);
  hphiB->Draw("");
  cc2->cd(2);
  htheB->Draw("");
  */


  TCanvas *cenelen = new TCanvas("cenelen","cenelen",600,600);
  henelenAg->Draw("");
  henelenBr->Draw("sames");
  henelenC->Draw("sames");
  henelenO->Draw("sames");
  henelenN->Draw("sames");
  
  TCanvas *c0 = new TCanvas("c0","c0",1200,600);
  c0->Divide(2,1);
  c0->cd(1);
  h0->Draw("");
  c0->cd(2);
  h00->Draw("");
  
  
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

  
  TCanvas *c3 = new TCanvas("c3","c3",600,600);
  // c3->Divide(2,1);
  //c3->cd(1);
  h3->Draw("");
  //c3->cd(2);
  //h33->Draw("");

  TCanvas *c33 = new TCanvas("c33","c33",1200,600);
  // c3->Divide(2,1);
  //c3->cd(1);
  h33->Draw("");
  //c3->cd(2);
  //h33->Draw("");

  TCanvas *c5 = new TCanvas("c5","c5",600,600);
  h5->Draw("");

  TCanvas *cphi = new TCanvas("cphi","cphi",600,600);
  hphirec->Draw("");
  TCanvas *cthe = new TCanvas("cthe","cthe",600,600);
  htherec->Draw("");
  TCanvas *ccosthe = new TCanvas("ccosthe","ccosthe",600,600);
  hcostherec->Draw("");


  TCanvas *cxyz = new TCanvas("cxyz","cxyz",600,600);
  hxyzrec->Draw("ISO");
  hxyzrec->SetFillColor(kBlue+1);
  
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
  hlen_thr->Scale(num_ev/hlen_thr->Integral());
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
  
  // cc7 = (TCanvas*)functions->trk_len();

  cmg2 = new TCanvas("cmg2","cmg2",600,600);
  
  for(int i=0;i<3;i++){
    gr_spect[i]->SetTitle(flavour[i]);
    gr_spect[i]->SetLineWidth(2);
    gr_spect[i]->SetLineColor(i+1);
    //gr_spect[i]->SetMarkerStyle(22);
    gr_spect[i]->SetFillColor(0);
    mg2->Add(gr_spect[i]);
  }
  
  cmg2->cd();
  //cmg2->SetFillStyle(4000);
  cmg2->SetLogy();
  mg2->Draw("AL");
  mg2->GetXaxis()->SetTitle("E_{#nu} [MeV]");
  mg2->GetYaxis()->SetTitle("f(E) [ MeV^{-1} ]");
  cmg2->BuildLegend();

  cmg = (TCanvas*)functions->gr_crsect();

 
  
  delete functions;
}
