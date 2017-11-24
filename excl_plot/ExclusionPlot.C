//************************************//
//***    Author: A. Di Crescenzo   ***//
//***    Created: 15 Nov 2015      ***//
//***    Last modified: 3 Dec 2016 ***//
//************************************//
//
// Functions implemented:
//
// DrawRate(range_in_nm, SD=false)
// TotalRate(mass_in_gev, range_in_min, SD=false)
// DifferentialRate(mass_in_gev, lement_charge, Erecoil_in_kev, SD=false)
// Analysis_Exclusion()
// Analysis_Threshold()
// Analysis_Directional1(case)
// Analysis_Directional2(case)
// AntoVale()
// NeutrinoBkg() //Xe target
// NeutrinoFloor() //NEWS 
//
//************************************//

#include "TRandom3.h"
#include "TMath.h"
#include "Definitions.h"
#include "EnergyRangeCorrelation.C"
#include <TString.h>
#include <ctime>

void DrawRate(Float_t range_min, bool SD = false){

  cout.precision(15);
  t = time(0); now = localtime(&t);
  cout <<  now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec;
  cout << " *** STARTTTTTTTTTTT *** "<< endl;
  
  verbosity = false;
  Int_t cont = 0;

  for(Int_t iele=0; iele<Nelements; iele++){
    Emin[iele]  = EnergyFromCut(Charge[iele],range_min);
    Emax[iele]  = EnergyFromCut(Charge[iele],range_max);
    DiffRateE[iele] = 0;
  }
  
  char nameFile[1000];
  char nameFileExcl[255];
  char nameFileRate[255];
  char nameFileDM[255];
  if(element_selection==0){ // all but H
    sprintf(nameFile,"RESULTS/%0.0fnm_%dkg_%ddays_All.root",range_min,mass,time);
    sprintf(nameFileExcl,"RESULTS/Excl_%0.0fnm_%dkg_%ddays_All.pdf",range_min,mass,time);
    sprintf(nameFileRate,"RESULTS/Rate_%0.0fnm_%dkg_%ddays_All.pdf",range_min,mass,time);
    sprintf(nameFileDM,"RESULTS/Rate_%0.0fnm_%dkg_%ddays_All.txt",range_min,mass,time); }
  if(element_selection==1){ // all but H
    sprintf(nameFile,"RESULTS/%0.0fnm_%dkg_%ddays_All.root",range_min,mass,time);
    sprintf(nameFileExcl,"RESULTS/Excl_%0.0fnm_%dkg_%ddays_All.pdf",range_min,mass,time);
    sprintf(nameFileRate,"RESULTS/Rate_%0.0fnm_%dkg_%ddays_All.pdf",range_min,mass,time);
    sprintf(nameFileDM,"RESULTS/Rate_%0.0fnm_%dkg_%ddays_All.txt",range_min,mass,time); }
  if(element_selection==3){ // light nuclei
    sprintf(nameFile,"RESULTS/%0.0fnm_%dkg_%ddays_Light.root",range_min,mass,time);
    sprintf(nameFileExcl,"RESULTS/Excl_%0.0fnm_%dkg_%ddays_Light.pdf",range_min,mass,time);
    sprintf(nameFileRate,"RESULTS/Rate_%0.0fnm_%dkg_%ddays_Light.pdf",range_min,mass,time);
    sprintf(nameFileDM,"RESULTS/Rate_%0.0fnm_%dkg_%ddays_Light.txt",range_min,mass,time); }
  if(element_selection==2){ // heavy nuclei
    sprintf(nameFile,"RESULTS/%0.0fnm_%dkg_%ddays_Heavy.root",range_min,mass,time);
    sprintf(nameFileExcl,"RESULTS/Excl_%0.0fnm_%dkg_%ddays_Heavy.pdf",range_min,mass,time);
    sprintf(nameFileRate,"RESULTS/Rate_%0.0fnm_%dkg_%ddays_Heavy.pdf",range_min,mass,time);
    sprintf(nameFileDM,"RESULTS/Rate_%0.0fnm_%dkg_%ddays_Heavy.txt",range_min,mass,time); }
  if(element_selection==4){ // heavy nuclei
    sprintf(nameFile,"RESULTS/%0.0fnm_%dkg_%ddays_Solfur.root",range_min,mass,time);
    sprintf(nameFileExcl,"RESULTS/Excl_%0.0fnm_%dkg_%ddays_Sulfur.pdf",range_min,mass,time);
    sprintf(nameFileRate,"RESULTS/Rate_%0.0fnm_%dkg_%ddays_Sulfur.pdf",range_min,mass,time);
    sprintf(nameFileDM,"RESULTS/Rate_%0.0fnm_%dkg_%ddays_Sulfur.txt",range_min,mass,time); }
  
  TFile f(nameFile,"new");
  FILE *out = fopen (nameFileDM, "w");
	
  t = time(0); now = localtime(&t);
  cout <<  now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec;
  cout << " -> Investigating mass rage 1-10 GeV "<< endl;
  
  for(Int_t iMD=10; iMD<100; iMD++){
    MD = iMD*0.1;
    //cout << "Mass "<< MD<< endl;
    rate = TotalRate(MD, range_min,SD);
    if(rate<1.0E-15) continue;
    grRate->SetPoint(cont,MD,rate);
    //cout << "Rate "<< rate << endl;
    cross =  (rate_UL/rate) * sigmaP;
    //cout << "Cross Section "<< cross << endl;
    fprintf(out, " %i %e;" , MD, cross);
    grExcl->SetPoint(cont,MD,cross);
    cont ++;
  }

  t = time(0); now = localtime(&t);
  cout <<  now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec;
  cout << " -> Investigating mass rage 10-100 GeV "<< endl;
    
  for(Int_t iMD=10; iMD<100; iMD++){
    MD = iMD;
    //cout << "Mass "<< MD<< endl;
    rate = TotalRate(MD, range_min,SD);
    if(rate<1.0E-15) continue;
    grRate->SetPoint(cont,MD,rate);
    //cout << "Rate "<< rate << endl;
    cross =  (rate_UL/rate) * sigmaP;
    //cout << "Cross Section "<< cross << endl;
    fprintf(out, " %i %e;" , MD, cross);
    grExcl->SetPoint(cont,MD,cross);
    cont ++;
  }

  t = time(0); now = localtime(&t);
  cout <<  now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec; 
  cout << " -> Investigating mass rage 100-1000 GeV "<< endl;
  
  for(Int_t iMD=10; iMD<100; iMD++){
    MD = iMD*10;
    //cout << "Mass "<< MD<< endl;
    rate = TotalRate(MD, range_min,SD);
    if(rate<1.0E-15) continue;
    grRate->SetPoint(cont,MD,rate);
    //cout << "Rate "<< rate << endl;
    cross =  (rate_UL/rate) * sigmaP;
    fprintf(out, " %i %e;" , MD, cross);
    //cout << "Cross Section "<< cross << endl;
    grExcl->SetPoint(cont,MD,cross);
    cont ++;
  }

  t = time(0); now = localtime(&t);
  cout <<  now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec; 
  cout << " -> Investigating mass rage 1000-10000 GeV "<< endl;
  
  for(Int_t iMD=10; iMD<100; iMD++){
    MD = iMD*100;
    //cout << "Mass "<< MD<< endl;
    rate = TotalRate(MD, range_min,SD);
    if(rate<1.0E-15) continue;
    grRate->SetPoint(cont,MD,rate);
    //cout << "Rate "<< rate << endl;
    cross =  (rate_UL/rate) * sigmaP;
    fprintf(out, " %i %e;" , MD, cross);
    //cout << "Cross Section "<< cross << endl;
    grExcl->SetPoint(cont,MD,cross);
    cont ++;
  }
  
  TCanvas *cRate = new TCanvas("cRate","",600,600);
  TH1F *frRate = gPad->DrawFrame(0,1.0E-43,1000,1.0E-35);
  frRate->Draw();
  TSpline3 *sRate = new TSpline3("sRate",grRate->GetX(),grRate->GetY(),grRate->GetN());
  sRate->SetLineColor(kBlue);
  grRate->SetLineColor(kBlue);
  grRate->SetLineWidth(2);
  grRate->SetMarkerColor(kBlue);
  grRate->SetMarkerStyle(21);
  grRate->SetMarkerSize(0.1);
  grRate->GetYaxis()->SetTitleOffset(1.2);
  grRate->GetXaxis()->SetTitle("WIPM mass (GeV/c^{2})");
  grRate->GetYaxis()->SetTitle("Total rate");
  grRate->GetXaxis()->CenterTitle();
  grRate->GetYaxis()->CenterTitle();
  grRate->Draw("alp");
  sRate->Draw("l same");
  grRate->SetName("grRate");
  grRate->Write();
  sRate->Write();
  cRate->Write();
  cRate->SaveAs(nameFileRate);
    
  TCanvas *cExcl = new TCanvas("cExcl","",600,600);
  TH1F *frExcl = gPad->DrawFrame(0,0,1000,100);
  frExcl->Draw();
  TSpline3 *sExcl = new TSpline3("sExcl",grExcl->GetX(),grExcl->GetY(),grExcl->GetN());
  sExcl->SetLineColor(kBlue);
  grExcl->SetMarkerColor(kBlue);
  grExcl->SetLineColor(kBlue);
  grExcl->SetLineWidth(2);
  grExcl->SetMarkerStyle(21);
  grExcl->SetMarkerSize(0.1);
  grExcl->GetYaxis()->SetTitleOffset(1.7);
  grExcl->GetXaxis()->SetTitle("WIPM mass (GeV/c^{2})");
  grExcl->GetYaxis()->SetTitle("WIPM-nucleus cross section (cm^{2})");
  grExcl->GetXaxis()->CenterTitle();
  grExcl->GetYaxis()->CenterTitle();
  grExcl->Draw("alp");
  sExcl->Draw("l same");
  cExcl->SetLogx();
  cExcl->SetLogy();
  cExcl->SetGridx();
  cExcl->SetGridy();
  grExcl->SetName("grExcl");
  grExcl->Write();
  sExcl->Write();
  cExcl->Write();
  cExcl->SaveAs(nameFileExcl);

  t = time(0); now = localtime(&t);
  cout <<  now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec;
  cout << " *** END !!!!! *** "<< endl;

  fclose(out);
}

//-------------------------------------

Float_t TotalRate(Float_t MD, Float_t range_min, bool SD = false)
{

  rate = 0;

  if(verbosity==true)
    for(Int_t iele=0; iele<Nelements; iele++){
      Emin[iele]  = EnergyFromCut(Charge[iele],range_min);
      Emax[iele]  = EnergyFromCut(Charge[iele],range_max);
      DiffRateE[iele] = 0;
    }
  
  for(Int_t iele=0; iele<Nelements; iele++)
    DiffRateE[iele] = 0;
  
  for(Int_t iEnergy=1; iEnergy<1500; iEnergy++){ //from 1 keV to 1.5 MeV

    iEnergy = iEnergy;
    for(Int_t iele=0; iele<Nelements; iele++){
      if(iEnergy>Emin[iele] && iEnergy<Emax[iele]){
	diffRate = DifferentialRate(MD,Charge[iele],iEnergy,SD);
	if(diffRate>0){
	  DiffRateE[iele] += diffRate*MassFraction[iele];}
      }
    }// loop on elements
  }// loop on energies

   for(Int_t iele=0; iele<Nelements; iele++){
      RateE[iele]  = 1.*DiffRateE[iele] / MassFraction[iele];
      SigmaE[iele] = (rate_UL/DiffRateE[iele])*sigmaP;
   }// loop on elements
   
  for(Int_t iele=0; iele<Nelements; iele++){
    if(element_selection==0){ //all elements
      rate += DiffRateE[iele];
      crossInv += 1./SigmaE[iele];
    }
    if(element_selection==1 && Charge[iele]!=1){ //discard protons
      rate += DiffRateE[iele];
      if(DiffRateE[iele]>0) crossInv +=  1./ SigmaE[iele];
      //printf("CHECK %8.2e  %8.2e\n",SigmaE[iele],crossInv);
    }
    if(element_selection==2&&Charge[iele]!=1&&Charge[iele]!=16&&Charge[iele]!=6&&Charge[iele]!=7&&Charge[iele]!=8){ //discard light nucle
      rate += DiffRateE[iele];
      crossInv += 1./SigmaE[iele];
    }
    if(element_selection==3&&Charge[iele]!=1&&Charge[iele]!=16&&Charge[iele]!=35&&Charge[iele]!=47&&Charge[iele]!=53){ //discard heavy nucle
      rate += DiffRateE[iele];
      crossInv += 1./SigmaE[iele];
    }
    if(element_selection==4 && Charge[iele]==16){ //discard protons
      rate += DiffRateE[iele];
      if(DiffRateE[iele]>0) crossInv +=  1./ SigmaE[iele];
      //printf("CHECK %8.2e  %8.2e\n",SigmaE[iele],crossInv);
    }
  }
  
  rate  = rate * time * mass;
  cross = (rate_UL/rate) * sigmaP;
  crossBis = 1./crossInv;
  
  if(verbosity){
    printf("\n");
    printf("  Z   Emin(keV)  Emax(keV)  Rate\n");
    for(Int_t iele=0; iele<Nelements; iele++)
      printf(" %2d %8.1f   %8.1f   %8.2e   %8.2e  %8.2e \n",Charge[iele],Emin[iele],Emax[iele],DiffRateE[iele],RateE[iele],SigmaE[iele]);
    
    printf("\n");
    printf(" *****   PARAMETERS  *****\n");
    printf(" Threshold     = %1.0f nm\n",range_min);
    printf(" WIMP Mass     = %d GeV/c2\n",MD);
    printf(" Cross section = %0.1e cm2\n",sigmaP);
    printf(" Exposure time = %d days\n",time);
    printf(" Exposure mass = %d kg\n",mass);
    printf("\n");
    printf(" *****     RESULT    *****\n");
    printf(" Total Rate      = %0.2e\n",rate);
    printf(" CrossSec Limit  = %0.2e\n",cross);
    printf(" CrossSec Limit2 = %0.2e\n",crossBis);
    printf(" *************************\n");
    printf("\n");
  } //if verbosity
  
  return rate;
}

//-------------------------------------

Float_t DifferentialRate(Float_t MD, Int_t charge, Float_t Erecoil, bool SD = false){
  
  DiffRate = 0;

  Sigma = Sigma(MD,charge,SD);
  R0 = (2/sqrtPi)*(N0/A)*(rhoD/MD)*Sigma*(v0*c)*1.0E5*Nseconds; // differential rate/kg/day (tru)
  E0 = 0.5*MD*v0*v0*1.0E6; //in keV
  r  = 4*MD*MT/(MD+MT)**2; // dimensionless
  vmin = (v0*c)*sqrt(Erecoil/(E0*r)); // km/s
  vmin = vmin/c; //in units of c
  
  //diff rate, (vE,vEsc=inf)
  DiffRate0 = (R0/(E0*r)) * (sqrtPi/4) * (v0/vE) *
    (TMath::Erf((vmin+vE)/v0) - TMath::Erf((vmin-vE)/v0)); // differential rate/keV/kg/day

  //diff rate, (vE,vEsc)s
  DiffRate1 = DiffRate0 - (R0/(E0*r))*TMath::Exp(-z*z);
  DiffRate1 = (k0/k1)*DiffRate1;
  
  FormFactor = FormFactorSI(charge,Erecoil);
  if(SD==true) FormFactor = FormFactorSD(charge,Erecoil);

  DiffRate = DiffRate0*FormFactor;

  //printf("DM Mass %4.2f, charge %d, Energy %4.2f\n",MD,charge,Erecoil);
  //printf("CrossSection %f, R0 %4.5f, E0 %4.5f, r %4.5f, vmin %4.5f\n",Sigma,R0,E0,r,vmin);
  //printf("DiffRate %f %f\n",Sigma,FormFactor);
  
  return DiffRate;
}

//----------------------------------


Float_t Sigma(Float_t MD, Int_t charge, bool SD = false){

  Float_t Sigma = 0;

  A      = MassNumberFromCharge(charge);
  Lande  = LandeFromCharge(charge);
  MT     = A * Mnucleon;
  muT    = MT*MD/(MT+MD);
  muP    = Mproton*MD/(Mproton+MD);
  
  Sigma = sigmaP*(muT/muP)*(muT/muP)*A*A;
  if(SD==true)
    Sigma = sigmaP*(muT/muP)*(muT/muP)*(4./3.)*LandeFromCharge(charge);
  
  return Sigma;
}

//-------------------------------------

Float_t FormFactorSD(Int_t charge, Float_t Erecoil){
  //energy in keV
  
  Float_t FormFactor = 0;
  A = MassNumberFromCharge(charge);

  Float_t rn = 1.14*TMath::Power(A,1./3); //radius nuclei (fm)
  Float_t qr = 6.92E-3 * TMath::Sqrt(A*Erecoil) * rn; //dimensionless
  
  FormFactor = TMath::Sin(qr)/qr;
  FormFactor = TMath::Power(FormFactor,2);
  
  return FormFactor;
}

//----------------------------------

Float_t FormFactorSI(Int_t charge, Float_t Erecoil){
  //energy in keV
  
  Float_t FormFactor = 0;
  A = MassNumberFromCharge(charge);

  Float_t rn = 1.14*TMath::Power(A,1./3); //radius nuclei (fm)
  Float_t qr = 6.92E-3 * TMath::Sqrt(A*Erecoil) * rn; //dimensionless

  FormFactor = 3*(TMath::Sin(qr) - qr*TMath::Cos(qr))/TMath::Power(qr,3);
  FormFactor = TMath::Power(FormFactor,2);
  return FormFactor;
}

//----------------------------------

Float_t  EnergyFromCut(Int_t charge, Float_t cut){
  //cut in nm
  //energy in keV
  Float_t energy = 0;

  energy = GetEnergyFromRange(charge,cut);
  return energy;
}


//----------------------------------

Float_t LandeFromCharge(Int_t charge){
  Float_t _Lande = -1;
  
  if(charge==1) // H
    _Lande = 0.750;
  if(charge==6) // C
    _Lande = 0.0;
  if(charge==7) // N
    _Lande = 0.087;
  if(charge==8) // O
    _Lande = 0.0;
  if(charge==16) // S
    _Lande = 0.0;
  if(charge==53) // I
    _Lande = 0.007;
  if(charge==35) // Br
    _Lande = 0.027;
  if(charge==47) // Ag
    _Lande = 0.054;
  
  return _Lande;

}
//----------------------------------

Int_t MassNumberFromCharge(Int_t charge){
  Int_t MassNumber = -1;
  
  if(charge==1) // H
    MassNumber = 1;
  if(charge==6) // C
    MassNumber = 12;
  if(charge==7) // N
    MassNumber = 14;
  if(charge==8) // O
    MassNumber = 16;
  if(charge==16) // S
    MassNumber = 32;
  if(charge==53) // I
    MassNumber = 127;
  if(charge==35) // Br
    MassNumber = 80;
  if(charge==47) // Ag
    MassNumber = 108;
  
  return MassNumber;

}

//------------------------------------

void Analysis_Exclusion() {

  TMultiGraph *mgRate = new TMultiGraph();
  TMultiGraph *mgExcl = new TMultiGraph();
  
  TFile file100_1_365_all("RESULTS/100nm_10kg_3650days_All.root");
  TFile file100_1_365_heavy("RESULTS/100nm_10kg_3650days_Heavy.root");
  TFile file100_1_365_light("RESULTS/100nm_10kg_3650days_Light.root");
  TFile file100_1_365_sulfur("RESULTS/100nm_10kg_3650days_Sulfur.root");
  
  TGraph *gRate_100_1_365_all = (TGraph*)file100_1_365_all.Get("grRate");
  TGraph *gExcl_100_1_365_all = (TGraph*)file100_1_365_all.Get("grExcl");

  TGraph *gRate_100_1_365_heavy = (TGraph*)file100_1_365_heavy.Get("grRate");
  TGraph *gExcl_100_1_365_heavy = (TGraph*)file100_1_365_heavy.Get("grExcl");

  TGraph *gRate_100_1_365_light = (TGraph*)file100_1_365_light.Get("grRate");
  TGraph *gExcl_100_1_365_light = (TGraph*)file100_1_365_light.Get("grExcl");

  TGraph *gRate_100_1_365_sulfur = (TGraph*)file100_1_365_sulfur.Get("grRate");
  TGraph *gExcl_100_1_365_sulfur = (TGraph*)file100_1_365_sulfur.Get("grExcl");

  TFile file_Dama1("RESULTS/Dama1.root");
  TGraph *gDama1 = file_Dama1.Get("Graph");

  TFile file_Dama2("RESULTS/Dama2.root");
  TGraph *gDama2 = file_Dama2.Get("Graph");
  
  // DRAW RATE GRAPHS
  
  gRate_100_1_365_heavy->SetLineColor(kRed+1);
  gRate_100_1_365_heavy->SetLineWidth(2);
  gRate_100_1_365_light->SetLineColor(kGreen+1);
  gRate_100_1_365_light->SetLineWidth(2);
  gRate_100_1_365_all->SetLineColor(kBlack);
  gRate_100_1_365_all->SetLineWidth(2);
  gRate_100_1_365_all->SetLineStyle(2);
  gRate_100_1_365_sulfur->SetLineColor(kAzure+7);
  gRate_100_1_365_sulfur->SetLineWidth(2);

  TCanvas *cRate = new TCanvas("cRate","",700,700);
  TH1F *frRate = gPad->DrawFrame(0,1.0E-10,1000,100);
  frRate->Draw();
  //cRate->SetLogx();
  //cRate->SetLogy();
  cRate->SetGridx();
  cRate->SetGridy();
  
  mgRate->Add(gRate_100_1_365_heavy);
  mgRate->Add(gRate_100_1_365_light);
  mgRate->Add(gRate_100_1_365_sulfur);
  mgRate->Add(gRate_100_1_365_all);
  mgRate->Draw("al");

  cRate->Update();
  cRate->SetTitle("100nm - 10kg 10years #sigma=10^{-41} cm^{2}");
  mgRate->GetXaxis()->SetTitleOffset(1.1);
  mgRate->GetYaxis()->SetTitleOffset(1.5);
  mgRate->GetXaxis()->SetTitle("WIPM mass (GeV/c^{2})");
  //  mgRate->GetYaxis()->SetTitle("Total rate (kg^{-1} year^{-1})");
  mgRate->GetYaxis()->SetTitle("Total rate");
  mgRate->GetXaxis()->CenterTitle();
  mgRate->GetYaxis()->CenterTitle();
  cRate->Modified();

  TLegend * legE = new TLegend();
  legE = new TLegend(0.70,0.60,1.0,0.80);
  legE->SetTextSize(0.03);
  legE->SetBorderSize(1);
  legE->SetLineStyle(0);
  legE->SetTextSize(0.03);
  legE->SetFillStyle(1001);
  legE->SetFillColor(kWhite);
  legE->AddEntry(gRate_100_1_365_light,"C, N, O","l");
  legE->AddEntry(gRate_100_1_365_heavy,"Br, Ag, I","l");
  legE->AddEntry(gRate_100_1_365_sulfur,"S","l");
  legE->AddEntry(gRate_100_1_365_all,"ALL","l");
  legE->AddEntry(gDama1,"DAMA/LIBRA","f");
  legE->Draw();
  
  // DRAW EXCLUSION PLOTS
  
  gExcl_100_1_365_heavy->SetLineColor(kRed+1);
  gExcl_100_1_365_heavy->SetLineWidth(2);
  gExcl_100_1_365_light->SetLineColor(kGreen+1);
  gExcl_100_1_365_light->SetLineWidth(2);
  gExcl_100_1_365_sulfur->SetLineColor(kAzure+1);
  gExcl_100_1_365_sulfur->SetLineWidth(2);
  gExcl_100_1_365_all->SetLineColor(kBlack);
  gExcl_100_1_365_all->SetLineWidth(2);
  gExcl_100_1_365_all->SetLineStyle(2);
  gDama1->SetFillColor(kRed-10);
  gDama1->SetLineColor(kRed-10);
  gDama2->SetFillColor(kRed-10);
  gDama2->SetLineColor(kRed-10);
  
  TCanvas *cExcl = new TCanvas("cExcl","100nm - 10kg 10years",700,700);
  TH1F *frExcl = gPad->DrawFrame(0,1.0E-9,1000,100);
  frExcl->Draw();
  cExcl->SetLogx();
  cExcl->SetLogy();
  cExcl->SetGridx();
  cExcl->SetGridy();

  mgExcl->Add(gDama1,"F");
  mgExcl->Add(gDama2,"F");
  mgExcl->Add(gExcl_100_1_365_heavy);
  mgExcl->Add(gExcl_100_1_365_light);
  mgExcl->Add(gExcl_100_1_365_sulfur);
  mgExcl->Add(gExcl_100_1_365_all);
  mgExcl->Draw("al");
  
  cExcl->Update();
  mgExcl->SetTitle("100nm - 10kg 10years");
  mgExcl->SetName("100nm - 10kg 10years");
  mgExcl->GetYaxis()->SetTitleOffset(1.7);
  mgExcl->GetXaxis()->SetTitleOffset(1.1);
  mgExcl->GetYaxis()->SetTitle("WIPM-nucleus cross section (cm^{2})");
  mgExcl->GetXaxis()->SetTitle("WIPM mass (GeV/c^{2})");
  mgExcl->GetXaxis()->CenterTitle();
  mgExcl->GetYaxis()->CenterTitle();
  legE->Draw();
  cExcl->Modified();
  
}

//------------------------------------------------

void Analysis_Threshold() {

  TMultiGraph *mgExcl = new TMultiGraph();
  
  TFile file100_1_365_all("RESULTS/100nm_1kg_365days_All.root");
  TFile file200_1_365_all("RESULTS/200nm_1kg_365days_All.root");
  TFile file50_1_365_all("RESULTS/50nm_1kg_365days_All.root");
  
  TGraph *gExcl_100_1_365_all = (TGraph*)file100_1_365_all.Get("grExcl");
  TGraph *gExcl_200_1_365_all = (TGraph*)file200_1_365_all.Get("grExcl");
  TGraph *gExcl_50_1_365_all = (TGraph*)file50_1_365_all.Get("grExcl");

  TFile file_Dama1("RESULTS/Dama1.root");
  TGraph *gDama1 = file_Dama1.Get("Graph");

  TFile file_Dama2("RESULTS/Dama2.root");
  TGraph *gDama2 = file_Dama2.Get("Graph");
  
  TLegend * legE = new TLegend();
  legE = new TLegend(0.70,0.60,1.0,0.80);
  legE->SetTextSize(0.03);
  legE->SetBorderSize(1);
  legE->SetLineStyle(0);
  legE->SetTextSize(0.03);
  legE->SetFillStyle(1001);
  legE->SetFillColor(kWhite);
  legE->AddEntry(gExcl_200_1_365_all,"NEWS 1 kg - 1 year","");
  legE->AddEntry(gExcl_200_1_365_all,"200 nm","l");
  legE->AddEntry(gExcl_100_1_365_all,"100 nm","l");
  legE->AddEntry(gExcl_50_1_365_all,"50 nm","l");
  legE->AddEntry(gDama1,"DAMA/LIBRA","f");
  legE->Draw();
  
  // DRAW EXCLUSION PLOTS
  
  gExcl_100_1_365_all->SetLineColor(kRed+1);
  gExcl_100_1_365_all->SetLineWidth(3);
  gExcl_200_1_365_all->SetLineColor(kGreen+1);
  gExcl_200_1_365_all->SetLineWidth(3);
  gExcl_50_1_365_all->SetLineColor(kBlack);
  gExcl_50_1_365_all->SetLineWidth(3);
  gDama1->SetFillColor(kRed-10);
  gDama1->SetLineColor(kRed-10);
  gDama2->SetFillColor(kRed-10);
  gDama2->SetLineColor(kRed-10);
  
  TCanvas *cExcl = new TCanvas("cExcl","100nm - 1kg 1year",700,700);
  TH1F *frExcl = gPad->DrawFrame(0,1.0E-44,1000,100);
  frExcl->Draw();
  cExcl->SetLogx();
  cExcl->SetLogy();
  cExcl->SetGridx();
  cExcl->SetGridy();

  mgExcl->Add(gDama1,"F");
  mgExcl->Add(gDama2,"F");
  mgExcl->Add(gExcl_100_1_365_all);
  mgExcl->Add(gExcl_200_1_365_all);
  mgExcl->Add(gExcl_50_1_365_all);
  mgExcl->Draw("al");
  
  cExcl->Update();
  mgExcl->GetYaxis()->SetTitleOffset(1.7);
  mgExcl->GetXaxis()->SetTitleOffset(1.1);
  mgExcl->GetYaxis()->SetTitle("WIPM-nucleus cross section (cm^{2})");
  mgExcl->GetXaxis()->SetTitle("WIPM mass (GeV/c^{2})");
  mgExcl->GetXaxis()->CenterTitle();
  mgExcl->GetYaxis()->CenterTitle();
  legE->Draw();
  cExcl->Modified();
  
}

//------------------------------------------------

void AntoVale() {

  TMultiGraph *mgExcl = new TMultiGraph();
  
  TFile file_anto("RESULTS/100nm_1kg_365days_all.root");
  TFile file_vale("RESULTS/exclusion_valerio.root");

  TGraph *gExcl_anto = (TGraph*)file_anto.Get("grExcl");
  TCanvas *c = (TCanvas*) file_vale.Get("c1");
  TGraph *gExcl_vale;

  TList* l = c->GetListOfPrimitives();
  TIter next(l);
  TObject *found, *obj;
  while ((obj=next())) {
    if (obj->InheritsFrom(TGraph::Class())) {
      gExcl_vale = (TGraph*)obj;
    }
  }
  
  gExcl_anto->SetLineColor(kRed+1);
  gExcl_anto->SetLineWidth(3);
  gExcl_vale->SetLineColor(kAzure+7);
  gExcl_vale->SetLineWidth(3);
  gExcl_vale->SetLineStyle(10);
 
  TCanvas *cExcl = new TCanvas("cExcl","100nm - 1kg 1year",700,700);
  TH1F *frExcl = gPad->DrawFrame(0,1.0E-9,1000,100);
  frExcl->Draw();
  cExcl->SetLogx();
  cExcl->SetLogy();
  cExcl->SetGridx();
  cExcl->SetGridy();

  mgExcl->Add(gExcl_anto);
  mgExcl->Add(gExcl_vale);
  mgExcl->Draw("al");
  
  cExcl->Update();
  mgExcl->SetTitle("100nm - 1kg 1year");
  mgExcl->SetName("100nm - 1kg 1year");
  mgExcl->GetYaxis()->SetTitleOffset(1.7);
  mgExcl->GetXaxis()->SetTitleOffset(1.1);
  mgExcl->GetYaxis()->SetTitle("WIPM-nucleus cross section (cm^{2})");
  mgExcl->GetXaxis()->SetTitle("WIPM mass (GeV/c^{2})");
  mgExcl->GetXaxis()->CenterTitle();
  mgExcl->GetYaxis()->CenterTitle();
  cExcl->Modified(); 

  TLegend * legE = new TLegend();
  legE = new TLegend(0.70,0.60,1.0,0.80);
  legE->SetTextSize(0.03);
  legE->SetBorderSize(1);
  legE->SetLineStyle(0);
  legE->SetTextSize(0.03);
  legE->SetFillStyle(1001);
  legE->SetFillColor(kWhite);
  legE->AddEntry(gExcl_vale,"Valerio","l");
  legE->AddEntry(gExcl_anto,"Antonia","l");
  legE->Draw();
}

//--------------------------


void NeutrinoBkg() {

  TMultiGraph *mgExcl = new TMultiGraph();
  
  //TFile file_anto("RESULTS/100nm_1000kg_1825days_all.root");
  // PAP
  TFile file_anto("RESULTS/30nm_100kg_36500days_all.root");
  // PAP
  TFile file_anto1("RESULTS/50nm_1000kg_36500days_all.root");
  //TFile file_anto("RESULTS/100nm_1000kg_1825days_All.root");  
  //TFile file_anto1("RESULTS/50nm_1000kg_1825days_All.root");
  TFile file_neutrino("RESULTS/Neutrino.root");


  TGraph *gExcl_anto = (TGraph*)file_anto.Get("grExcl");
  TGraph *gExcl_anto1 = (TGraph*)file_anto1.Get("grExcl");
  TGraph *gExcl_neutrino = file_neutrino.Get("Graph");
  
  gExcl_anto->SetLineColor(kBlue+1);
  gExcl_anto->SetLineWidth(3);
  gExcl_anto1->SetLineColor(kRed);
  gExcl_anto1->SetLineStyle(9);
  gExcl_anto1->SetLineWidth(3);
  gExcl_neutrino->SetLineColor(kGray);
  gExcl_neutrino->SetFillColor(kWhite);
  gExcl_neutrino->SetFillStyle(3002);
  gExcl_neutrino->SetLineWidth(6);
  gExcl_neutrino->SetLineStyle(2);
 
  TCanvas *cExcl = new TCanvas("cExcl","100nm - 1kg 1year",700,700);
  TH1F *frExcl = gPad->DrawFrame(0,1.0E-9,1000,100);
  frExcl->Draw();
  cExcl->SetLogx();
  cExcl->SetLogy();
  cExcl->SetGridx();
  cExcl->SetGridy();
 
  mgExcl->Add(gExcl_neutrino,"F");
  mgExcl->Add(gExcl_neutrino);
  mgExcl->Add(gExcl_anto1);
  mgExcl->Add(gExcl_anto);
  mgExcl->Draw("al");
  
  cExcl->Update();
  mgExcl->SetTitle("100nm - 1kg 1year");
  mgExcl->SetName("100nm - 1kg 1year");
  mgExcl->GetYaxis()->SetTitleOffset(1.7);
  mgExcl->GetXaxis()->SetTitleOffset(1.1);
  mgExcl->GetYaxis()->SetTitle("WIPM-nucleus cross section (cm^{2})");
  mgExcl->GetXaxis()->SetTitle("WIPM mass (GeV/c^{2})");
  mgExcl->GetXaxis()->CenterTitle();
  mgExcl->GetYaxis()->CenterTitle();
  mgExcl->GetYaxis()->SetRangeUser(1.0E-50,1.0E-37);
  cExcl->Modified(); 

  TLegend * legE = new TLegend();
  legE = new TLegend(0.70,0.60,1.0,0.80);
  legE->SetTextSize(0.03);
  legE->SetBorderSize(1);
  legE->SetLineStyle(0);
  legE->SetTextSize(0.03);
  legE->SetFillStyle(1001);
  legE->SetFillColor(kWhite);
  legE->AddEntry(gExcl_neutrino,"Neutrino Background","l");
  // PAP
  legE->AddEntry(gExcl_anto,"30 nm - 10 ton x year","l");
  // PAP
  legE->AddEntry(gExcl_anto1,"50 nm - 100 ton x year","l");
  //legE->AddEntry(gExcl_anto,"100 nm - 5 ton x year","l");
  //legE->AddEntry(gExcl_anto1,"50 nm - 5 ton x year","l");
  legE->Draw();
}

//--------------------------

void NeutrinoFloor() {

  TMultiGraph *mgExcl = new TMultiGraph();
  
  //TFile file_anto("RESULTS/100nm_1000kg_1825days_all.root");
  // PAP TFile file_anto("RESULTS/30nm_100kg_36500days_all.root");
  // PAP TFile file_anto1("RESULTS/50nm_1000kg_36500days_all.root");
  TFile file_anto("RESULTS/100nm_1000kg_1825days_All.root");  
  TFile file_anto1("RESULTS/50nm_1000kg_1825days_All.root");
  TFile file_neutrino("RESULTS/Neutrino.root");


  TGraph *gExcl_anto = (TGraph*)file_anto.Get("grExcl");
  TGraph *gExcl_anto1 = (TGraph*)file_anto1.Get("grExcl");
  TGraph *gExcl_neutrino = file_neutrino.Get("Graph");
  
  gExcl_anto->SetLineColor(kBlue+1);
  gExcl_anto->SetLineWidth(3);
  gExcl_anto1->SetLineColor(kRed+1);
  gExcl_anto1->SetLineStyle(2);
  gExcl_anto1->SetLineWidth(3);
  gExcl_neutrino->SetLineColor(kOrange+2);
  gExcl_neutrino->SetFillColor(kOrange+2);
  gExcl_neutrino->SetFillStyle(3002);
  gExcl_neutrino->SetLineWidth(10);
  //gExcl_neutrino->SetLineStyle(10);
 
  TCanvas *cExcl = new TCanvas("cExcl","100nm - 1kg 1year",700,700);
  TH1F *frExcl = gPad->DrawFrame(0,1.0E-9,1000,100);
  frExcl->Draw();
  cExcl->SetLogx();
  cExcl->SetLogy();
  cExcl->SetGridx();
  cExcl->SetGridy();
 
  mgExcl->Add(gExcl_neutrino,"F");
  mgExcl->Add(gExcl_neutrino);
  mgExcl->Add(gExcl_anto1);
  mgExcl->Add(gExcl_anto);
  mgExcl->Draw("al");
  
  cExcl->Update();
  mgExcl->SetTitle("100nm - 1kg 1year");
  mgExcl->SetName("100nm - 1kg 1year");
  mgExcl->GetYaxis()->SetTitleOffset(1.7);
  mgExcl->GetXaxis()->SetTitleOffset(1.1);
  mgExcl->GetYaxis()->SetTitle("WIPM-nucleus cross section (cm^{2})");
  mgExcl->GetXaxis()->SetTitle("WIPM mass (GeV/c^{2})");
  mgExcl->GetXaxis()->CenterTitle();
  mgExcl->GetYaxis()->CenterTitle();
  mgExcl->GetYaxis()->SetRangeUser(1.0E-50,1.0E-37);
  cExcl->Modified(); 

  TLegend * legE = new TLegend();
  legE = new TLegend(0.70,0.60,1.0,0.80);
  legE->SetTextSize(0.03);
  legE->SetBorderSize(1);
  legE->SetLineStyle(0);
  legE->SetTextSize(0.03);
  legE->SetFillStyle(1001);
  legE->SetFillColor(kWhite);
  legE->AddEntry(gExcl_neutrino,"Neutrino Background","l");
  // PAP legE->AddEntry(gExcl_anto,"30 nm - 10 ton x year","l");
  // PAP legE->AddEntry(gExcl_anto1,"50 nm - 100 ton x year","l");
  legE->AddEntry(gExcl_anto,"100 nm - 5 ton x year","l");
  legE->AddEntry(gExcl_anto1,"50 nm - 5 ton x year","l");
  legE->Draw();
}




//-------------------------------

void Analysis_Directional1(Int_t case=1){

  // Upper limit on signal events

   if(case==10){ // Nbkg=10, Mass=10kg, Time=10year, scattering=YES // filled on 14/07/2016 -
    const Int_t n = 15;
    Double_t mass[n]   = {1000, 800, 600, 400, 300, 200, 150, 100,  80,  60,  40,  30,  20,  15,  10};
    //Double_t limit[n]= {5.53,5.68,5.68,5.57,5.37,5.33,4.87,5.10,5.17,5.25,5.07,5.01,4.73,4.37,3.54};
    Double_t limit[n]  = {6.05,6.05,5.99,6.02,5.87,5.82,5.68,5.66,5.70,5.53,5.43,5.43,5.36,5.29,5.06}; 
    Double_t s1h[n]    = {0};
    Double_t s1l[n]    = {0};
    Double_t s2h[n]    = {0};
    Double_t s2l[n]    = {0};
    Double_t limit3[n] = {0};
    Double_t s3l[n]    = {0};
    Double_t s3h[n]    = {0};
    for(Int_t i=0; i<n; i++){
      limit3[i]= 6.7;
      s3l[i]   = 2.2;
      s3h[i]   = 2.7;}}

   if(case==100){ // Nbkg=100, Mass=10kg, Time=10year, scattering=YES // filled on 14/07/2016 -
    const Int_t n = 15;
    Double_t mass[n]   = {1000,   800,  600,  400,  300,  200,  150,  100,   80,   60,   40,   30,   20,   15,  10};
    //Double_t limit[n]= {15.57,14.94,14.92,14.82,14.76,14.41,14.09,14.02,14.43,14.72,13.91,13.22,13.03,11.43,8.98};
    Double_t limit[n]  = {16.58,15.82,16.32,16.25,16.04,15.49,16.25,15.89,15.71,15.47,15.07,15.02,14.80,14.89,14.68}; 
    Double_t s1h[n]    = {0};
    Double_t s1l[n]    = {0};
    Double_t s2h[n]    = {0};
    Double_t s2l[n]    = {0};
    Double_t limit3[n] = {0};
    Double_t s3l[n]    = {0};
    Double_t s3h[n]    = {0};
    for(Int_t i=0; i<n; i++){
       limit3[i]= 18.02;
      s3l[i]   = 5.22;
      s3h[i]   = 7.48;}}
   
  if(case==1){ // Nbkg=10, Mass=10kg, Time=10year // filled on 10/03/2016
    const Int_t n = 19;
    Double_t mass[n]   = {1000, 800, 600, 400, 300, 200, 150, 100,  80,  60,  40,  30,  20,  15,  10,   9,   8,   7,   6};
    Double_t limit[n]  = {5.53,5.68,5.68,5.57,5.37,5.33,4.87,5.10,5.17,5.25,5.07,5.01,4.73,4.37,3.54,3.37,3.15,3.01,2.76}; 
    Double_t s1h[n]    = {2.63,2.33,2.32,2.50,2.29,2.48,2.87,2.39,2.61,2.53,2.47,2.34,2.19,2.01,1.69,1.60,1.44,1.29,1.70};
    Double_t s1l[n]    = {1.46,1.76,1.79,1.71,1.43,1.51,1.17,1.28,1.72,1.54,1.32,1.52,1.33,1.13,1.01,0.23,0.61,0.50,0.34};
    Double_t s2h[n]    = {6.22,6.01,6.13,5.94,5.83,4.23,5.95,5.63,5.87,6.00,5.60,5.91,5.40,4.79,4.05,3.90,3.40,3.20,3.00};
    Double_t s2l[n]    = {2.21,2.75,2.57,2.67,2.30,2.78,1.66,2.02,2.13,2.14,2.28,2.06,1.87,1.77,1.23,1.08,0.79,0.68,0.456};
    Double_t limit3[n] = {0};
    Double_t s3l[n]    = {0};
    Double_t s3h[n]    = {0};
    for(Int_t i=0; i<n; i++){
      limit3[i]= 6.7;
      s3l[i]   = 2.2;
      s3h[i]   = 2.7;}}
  
  if(case==2){ // Nbkg=100, Mass=10kg, Time=10year // filled on 14/03/2016
    const Int_t n = 19;
    Double_t mass[n]   = { 1000,  800,  600,  400,  300,  200,  150,  100,   80,   60,   40,   30,   20,   15,  10,   9,   8,   7,   6};
    Double_t limit[n]  = {15.57,14.94,14.92,14.82,14.76,14.41,14.09,14.02,14.43,14.72,13.91,13.22,13.03,11.43,8.98,8.08,7.47,6.93,5.94};
    Double_t s1h[n]    = { 6.70, 7.81, 7.42, 6.60, 6.60, 6.04, 6.05, 5.76, 6.41, 6.36, 6.37, 6.46, 5.02, 5.39,4.48,3.88,3.46,3.15,2.55};
    Double_t s1l[n]    = { 4.88, 3.59, 3.89, 4.10, 4.25, 3.09, 3.59, 4.40, 4.63, 4.49, 3.41, 3.26, 3.70, 3.68,2.75,2.44,2.22,2.13,1.81};
    Double_t s2h[n]    = {16.00,15.90,15.80,15.60,15.60,14.60,14.40,14.20,15.00,14.30,14.90,14.90,13.00,12.27,9.60,9.32,7.97,7.25,6.19};
    Double_t s2l[n]    = { 7.44, 6.14, 6.94, 7.20, 7.20, 6.70, 7.13, 7.35, 6.61, 7.54, 6.44, 5.22, 5.84, 5.30,4.78,3.04,3.15,3.75,2.78};
    Double_t limit3[n] = {0};
    Double_t s3l[n]    = {0};
    Double_t s3h[n]    = {0};
    for(Int_t i=0; i<n; i++){
      limit3[i]= 18.02;
      s3l[i]   = 5.22;
      s3h[i]   = 7.48;}}
  
  if(case==3){
    const Int_t n = 6;
    Double_t mass[n]   = {1000, 200, 100,  20,  10,   4};
    Double_t limit[n]  = { 3.4, 3.4, 3.4, 3.4, 3.4, 3.3};
    Double_t s1l[n]    = { 0.1, 0.3, 0.3, 0.5, 0.6, 0.2};
    Double_t s1h[n]    = { 1.3, 1.0, 0.9, 0.7, 0.7, 0.6};
    Double_t s2l[n]    = { 0.2, 0.4, 0.4, 0.5, 0.6, 0.2};
    Double_t s2h[n]    = { 3.6, 2.7, 2.7, 2.2, 2.1, 1.6};
    Double_t limit3[n] = {0};
    Double_t s3l[n]    = {0};
    Double_t s3h[n]    = {0};
    for(Int_t i=0; i<n; i++){
      limit3[i]= 3.44;
      s3l[i]   = 2.02;
      s3h[i]   = 1.81;}}

  if(case==4){
    const Int_t n = 7;
    Double_t mass[n]   = {1000, 200, 100,  20,  10,   6,  4};
    Double_t limit[n]  = {5.08,4.75,4.52,4.01,3.63,3.68,3.26};
    Double_t s1l[n]    = {1.00,0.85,0.62,0.52,0.44,0.44,0.40};
    Double_t s1h[n]    = {1.15,1.31,1.38,0.66,1.11,0.67,0.84};
    Double_t s2l[n]    = {1.10,2.30,0.72,0.61,0.49,0.52,0.49};
    Double_t s2h[n]    = {3.50,3.69,3.50,3.12,3.20,3.14,2.45};
    Double_t limit3[n] = {0};
    Double_t s3l[n]    = {0};
    Double_t s3h[n]    = {0};
    for(Int_t i=0; i<n; i++){
      limit3[i]= 5.58;
      s3l[i]   = 1.7;
      s3h[i]   = 1.5;}}

  if(case==5){
    const Int_t n = 7;
    Double_t mass[n]   = {1000, 200, 100,  20,  10,   6,   4};
    Double_t limit[n]  = {19.2,16.5,15.7,11.4,10.6,9.91,8.45};
    Double_t s1l[n]    = { 6.3, 5.1, 4.9, 2.6, 2.6, 2.8, 2.0};
    Double_t s1h[n]    = { 8.6, 6.4, 6.6, 5.0, 4.3, 3.8, 3.9};
    Double_t s2l[n]    = { 8.5, 7.4, 6.2, 3.6, 4.1, 3.3, 2.3};
    Double_t s2h[n]    = {20.2,15.5,14.9,12.5,10.6, 9.9, 9.3};
    Double_t limit3[n] = {0};
    Double_t s3l[n]    = {0};
    Double_t s3h[n]    = {0};
    for(Int_t i=0; i<n; i++){
      limit3[i]= 50.6;
      s3l[i]   = 3.8;
      s3h[i]   = 5.8;}}
  
  Double_t limit5[n] = {0};
  Double_t s5l[n]    = {0};
  Double_t s5h[n]    = {0};
  
  for(Int_t i=0; i<n; i++){
    limit5[i]= 2.44;}
  
  Double_t exl[n] = {0};
  Double_t exh[n] = {0};
  
  //----------------
  //   draw graphs
  //----------------
  
  TMultiGraph *mgExcl = new TMultiGraph();
  gr0 = new TGraph(n,mass,limit);
  gr1 = new TGraphAsymmErrors(n,mass,limit,exl,exh,s1l,s1h);
  gr2 = new TGraphAsymmErrors(n,mass,limit,exl,exh,s2l,s2h);
  gr3 = new TGraphAsymmErrors(n,mass,limit3,exl,exh,s3l,s3h);
  gr4 = new TGraph(n,mass,limit3);
  gr5 = new TGraph(n,mass,limit5);
  gr1->SetFillColor(kGreen);
  gr2->SetFillColor(kYellow);
  gr3->SetFillColor(kRed);
  gr0->SetMarkerStyle(20);
  gr0->SetMarkerColor(kRed+1);
  gr0->SetLineColor(kRed+1);
  gr0->SetLineStyle(2);
  gr0->SetMarkerSize(1.0);
  gr0->SetLineWidth(2.8);
  gr4->SetLineColor(kBlue);
  //gr4->SetLineStyle(1);
  //gr4->SetLineColor(4);
  gr4->SetLineWidth(2.8);
    
  TCanvas *c1 = new TCanvas("c1","",700,700);
  //TH1F *frExcl = gPad->DrawFrame(0,0,1000,500);
  //frExcl->Draw();
  c1->SetLogx();
  c1->SetGridx();
  c1->SetGridy();

  //mgExcl->Add(gr3);
  //mgExcl->Add(gr2);
  //mgExcl->Add(gr1);
  mgExcl->Add(gr0);
  mgExcl->Add(gr4);
  //mgExcl->Add(gr5);
  mgExcl->Draw("alp3");
 
  c1->Update();
  mgExcl->SetTitle("100nm - 10kg 10years");
  mgExcl->SetName("100nm - 10kg 10years");
  mgExcl->GetYaxis()->SetTitleOffset(1.1);
  mgExcl->GetXaxis()->SetTitleOffset(1.1);
  mgExcl->GetYaxis()->SetTitle("Number of signal events");
  mgExcl->GetXaxis()->SetTitle("WIPM mass (GeV/c^{2})");
  mgExcl->GetXaxis()->CenterTitle();
  mgExcl->GetYaxis()->CenterTitle();
  mgExcl->GetYaxis()->SetRangeUser(0,80);
  if(case==1) mgExcl->GetYaxis()->SetRangeUser(0,15);
  if(case==2) mgExcl->GetYaxis()->SetRangeUser(0,35);
  if(case==3) mgExcl->GetYaxis()->SetRangeUser(0,10);
  if(case==4) mgExcl->GetYaxis()->SetRangeUser(0,15);
  c1->Modified();

  TLegend * legE = new TLegend();
  legE = new TLegend(0.70,0.60,1.0,0.80);
  legE->SetTextSize(0.02);
  legE->SetBorderSize(1);
  legE->SetLineStyle(0);
  legE->SetTextSize(0.03);
  legE->SetFillStyle(1001);
  legE->SetFillColor(kWhite);
  legE->AddEntry(gr0,"Likelihood","lp");
  //legE->AddEntry(gr0,"Median 90% C.L.","lp");
  //legE->AddEntry(gr1,"1#sigma band ","f");
  //legE->AddEntry(gr2,"2#sigma band ","f");
  legE->AddEntry(gr4,"Counting only","l");
  //legE->AddEntry(gr4,"Median 90% C.L.","l");
  //legE->AddEntry(gr3,"1#sigma band ","f");
  //legE->AddEntry(gr5,"Zero Background","l");
  legE->Draw();
}

//------------------------------------------

void Analysis_Directional2(Int_t case=1){ 
  
  // Upper limit on signal events

  if(case==10){ // Nbkg=10, Mass=10kg, Time=10year, scattering=YES // filled on 14/07/2016 -
    const Int_t n = 15;
    Double_t mass[n]   = {1000, 800, 600, 400, 300, 200, 150, 100,  80,  60,  40,  30,  20,  15,  10};
    Double_t limitN[n] = {6.05,6.05,5.99,6.02,5.87,5.82,5.68,5.66,5.70,5.53,5.43,5.43,5.36,5.29,5.06};
    Double_t limit[n];
    Double_t s1h[n];
    Double_t s1l[n];
    Double_t s2h[n];
    Double_t s2l[n];
    Double_t N0[n] = { 115,134,158,186,196,190,175,159,155,147,114,77,28, 8,0.43};
    for(Int_t i=0; i<n; i++){
      limit[i]=limitN[i]*(sigmaP/N0[i]);
    }
    TFile fileCount("RESULTS/100nm_10kg_3650days_All_2.44.root");
    TFile fileCount0("RESULTS/100nm_10kg_3650days_All_6.7.root");
    TFile fileCountLow("RESULTS/100nm_10kg_3650days_All_9.4.root");
    TFile fileCountUp("RESULTS/100nm_10kg_3650days_All_4.5.root");}

   if(case==100){ // Nbkg=100, Mass=10kg, Time=10year, scattering=YES // filled on 14/07/2016 -
    const Int_t n = 15;
    Double_t mass[n]   = {1000, 800, 600, 400, 300, 200, 150, 100,  80,  60,  40,  30,  20,  15,  10};
    Double_t limitN[n]  = {16.58,15.82,16.32,16.25,16.04,15.49,16.25,15.89,15.71,15.47,15.07,15.02,14.80,14.89,14.68}; 
    Double_t limit[n];
    Double_t s1h[n];
    Double_t s1l[n];
    Double_t s2h[n];
    Double_t s2l[n];
    Double_t N0[n] = { 115,134,158,186,196,190,175,159,155,147,114,77,28, 8,0.43};
    for(Int_t i=0; i<n; i++){
      limit[i]=limitN[i]*(sigmaP/N0[i]);
    }
     TFile fileCount("RESULTS/100nm_10kg_3650days_All_2.44.root");
    TFile fileCount0("RESULTS/100nm_10kg_3650days_All_18.0.root");
    TFile fileCountLow("RESULTS/100nm_10kg_3650days_All_25.5.root");
    TFile fileCountUp("RESULTS/100nm_10kg_3650days_All_10.0.root");}

   if(case==1){  // Nbkg=10, Mass=10kg, Time=10year // filled on 10/03/2016
    const Int_t n = 19;
    Double_t mass[n]   = {    1000,     800,     600,     400,     300,     200,     150,     100,      80,      60,      40,      30,      20,      15,      10,       9,       8,       7,       6};
    Double_t limit[n]  = {4.81E-43,4.25E-43,3.60E-43,2.98E-43,2.73E-43,2.80E-43,2.77E-43,3.21E-43,3.34E-43,3.57E-43,4.46E-43,6.49E-43,1.70E-42,5.43E-42,8.25E-41,2.37E-40,9.63E-40,7.21E-39,1.40E-37};
    Double_t s1h[n]    = {2.29E-43,1.75E-43,1.47E-43,1.35E-43,1.16E-43,1.31E-43,1.64E-43,1.50E-43,1.68E-43,1.74E-43,2.16E-43,3.03E-43,0.79E-42,2.50E-42,3.93E-41,1.12E-40,4.40E-40,3.10E-39,0.54E-37};
    Double_t s1l[n]    = {1.31E-43,1.32E-43,1.14E-43,0.91E-43,0.73E-43,0.80E-43,0.67E-43,0.80E-43,1.11E-43,1.04E-43,1.16E-43,1.98E-43,0.48E-42,1.40E-42,2.36E-41,0.58E-40,1.90E-40,1.21E-39,0.17E-37};
    Double_t s2h[n]    = {5.39E-43,4.50E-43,3.87E-43,3.19E-43,2.97E-43,3.00E-43,3.39E-43,3.55E-43,3.78E-43,4.12E-43,4.92E-43,7.65E-43,1.94E-42,5.96E-42,9.43E-41,2.80E-40,10.4E-40,7.66E-39,1.53E-37};
    Double_t s2l[n]    = {1.91E-43,2.06E-43,1.63E-43,1.43E-43,1.18E-43,1.45E-43,0.94E-43,1.28E-43,1.37E-43,1.44E-43,2.00E-43,2.68E-43,0.67E-42,2.20E-42,2.87E-41,0.76E-40,2.40E-40,1.63E-39,0.23E-37};
    TFile fileCount("RESULTS/100nm_10kg_3650days_All_2.44.root");
    TFile fileCount0("RESULTS/100nm_10kg_3650days_All_6.7.root");
    TFile fileCountLow("RESULTS/100nm_10kg_3650days_All_9.4.root");
    TFile fileCountUp("RESULTS/100nm_10kg_3650days_All_4.5.root");}
   
  if(case==2){ //Nbkg=100, Mass=10kg, Time=10year // filled on 14/03/2016
    const Int_t n = 19;
    Double_t mass[n]   = {    1000,     800,     600,     400,     300,     200,     150,     100,      80,      60,      40,      30,      20,      15,      10,       9,       8,       7,       6};
    Double_t limit[n]  = {1.36E-42,1.12E-42,9.44E-43,7.96E-43,7.51E-43,7.57E-43,8.03E-43,8.81E-43,9.30E-43,10.0E-43,12.2E-43,17.1E-43,4.69E-42,14.2E-42,2.09E-40,5.67E-40,2.28E-39,1.66E-38,3.01E-37};
    Double_t s1h[n]    = {0.59E-42,0.58E-42,4.70E-43,3.53E-43,3.36E-43,3.18E-43,3.45E-43,3.63E-43,4.13E-43,4.33E-43,5.59E-43,8.37E-43,1.81E-42,6.70E-42,1.04E-40,2.72E-40,1.06E-39,0.75E-38,1.29E-37};
    Double_t s1l[n]    = {0.42E-42,0.27E-42,2.40E-43,2.10E-43,2.16E-43,2.05E-43,2.05E-43,2.76E-43,2.98E-43,3.07E-43,3.00E-43,4.22E-43,1.35E-42,4.58E-42,0.64E-40,1.71E-40,0.67E-39,0.51E-38,0.91E-37};
    Double_t s2h[n]    = {1.39E-42,1.19E-42,9.90E-43,8.37E-43,7.96E-43,7.67E-43,8.20E-43,8.95E-43,9.67E-43,9.78E-43,13.1E-43,19.4E-43,4.69E-42,15.2E-42,2.24E-40,6.53E-40,2.44E-39,1.74E-38,3.13E-37};
    Double_t s2l[n]    = {0.65E-42,0.46E-42,4.40E-43,3.81E-43,3.66E-43,3.52E-43,4.06E-43,4.62E-43,4.26E-43,5.14E-43,5.55E-43,6.77E-43,2.10E-42,6.60E-42,1.11E-40,2.13E-40,0.96E-39,0.86E-38,1.41E-37};
    TFile fileCount("RESULTS/100nm_10kg_3650days_All_2.44.root");
    TFile fileCount0("RESULTS/100nm_10kg_3650days_All_18.0.root");
    TFile fileCountLow("RESULTS/100nm_10kg_3650days_All_25.5.root");
    TFile fileCountUp("RESULTS/100nm_10kg_3650days_All_10.0.root");}

  if(case==3){ //Nbkg=10, Mass=1 kg, TIme=1 year // filled on 03/04/2016
    const Int_t n = 19;
    Double_t s1l[n] ={0};
    Double_t s1h[n] ={0};
    Double_t s2l[n] ={0};
    Double_t s2h[n] ={0};
    Double_t mass[n]    = {1000, 800, 600, 400, 300, 200, 150, 100,  80,  60,  40,  30,  20,  15,  10,   9,   8,   7,   6};
    Double_t limitEv[n] = {5.53,5.68,5.68,5.57,5.37,5.33,4.87,5.10,5.17,5.25,5.07,5.01,4.73,4.37,3.54,3.37,3.15,3.01,2.76};
    Double_t N0[n]      = { 115,134,158,186,196,190,175,159,155,147,114,77,28, 8,0.43,0.14,0.03,4.17E-3,1.98E-4}; 
    Double_t limit[n];
    for(Int_t i=0; i<n; i++){
      limit[i]=limitEv[i]*(sigmaP/(N0[i]/100)); //gr0
    }
    TFile fileCount("RESULTS/100nm_1kg_365days_All_2.44.root");  //gExCount
    TFile fileCount0("RESULTS/100nm_1kg_365days_All_6.7.root");  //gExCount0
    TFile fileCountLow("RESULTS/100nm_1kg_365days_All_3.22.root"); //gExCountLow
    TFile fileCountUp("RESULTS/100nm_1kg_365days_All_5.46.root");}

  if(case==4){ 
    const Int_t n = 7;
    Double_t mass[n]   = {    1000,     200,     100,      20,      10,       6,       4};
    Double_t limit[n]  = {4.42E-41,2.50E-41,2.84E-41,1.44E-40,8.44E-39,1.86E-35,1.21E-29};
    Double_t s1l[n]    = {0.88E-41,0.45E-41,0.39E-41,0.19E-40,1.02E-39,0.22E-35,0.15E-29};
    Double_t s1h[n]    = {1.00E-41,0.69E-41,0.87E-41,0.24E-40,2.56E-39,0.34E-35,0.31E-29};
    Double_t s2l[n]    = {0.94E-41,1.21E-41,0.45E-41,0.22E-40,1.14E-39,0.27E-35,0.19E-29};
    Double_t s2h[n]    = {3.04E-41,1.94E-41,2.19E-41,1.13E-40,7.36E-39,1.58E-35,0.89E-29};
    TFile fileCount("RESULTS/100nm_1kg_365days_All_2.44.root");
    TFile fileCount0("RESULTS/100nm_1kg_365days_All_50.6.root");
    TFile fileCountLow("RESULTS/100nm_1kg_365days_All_46.8.root");
    TFile fileCountUp("RESULTS/100nm_1kg_365days_All_56.4.root");}

  if(case==5){ 
    const Int_t n = 7;
    Double_t mass[n]   = {    1000,     200,     100,      20,      10,       6,       4};
    Double_t limit[n]  = {1.67E-40,0.87E-40,0.99E-40,4.10E-40,2.46E-38,5.00E-35,5.56E-29};
    Double_t s1l[n]    = {0.55E-40,0.27E-40,0.31E-40,0.94E-40,0.61E-38,1.61E-35,1.30E-29};
    Double_t s1h[n]    = {0.70E-40,0.33E-40,0.41E-40,1.80E-40,0.64E-38,2.70E-35,2.60E-29};
    Double_t s2l[n]    = {0.74E-40,0.39E-40,0.39E-40,1.29E-40,0.82E-38,2.04E-35,1.50E-29};
    Double_t s2h[n]    = {1.75E-40,0.93E-40,0.93E-40,4.50E-40,2.58E-38,6.25E-35,6.10E-29};
    TFile fileCount("RESULTS/100nm_1kg_365days_All_2.44.root");
    TFile fileCount0("RESULTS/100nm_1kg_365days_All_50.6.root");
    TFile fileCountLow("RESULTS/100nm_1kg_365days_All_46.8.root");
    TFile fileCountUp("RESULTS/100nm_1kg_365days_All_56.4.root");}

  //-----------------------------
  // Draw graphs
  //-----------------------------
  TFile file_Dama1("RESULTS/Dama1.root");
  TGraph *gDama1 = file_Dama1.Get("Graph");

  TFile file_Dama2("RESULTS/Dama2.root");
  TGraph *gDama2 = file_Dama2.Get("Graph");

  TGraph *gExCount    = (TGraph*)fileCount.Get("grExcl");
  TGraph *gExCount0   = (TGraph*)fileCount0.Get("grExcl");
  TGraph *gExCountLow = (TGraph*)fileCountLow.Get("grExcl");
  TGraph *gExCountUp  = (TGraph*)fileCountUp.Get("grExcl");
  
  Double_t exl[n] = {0};
  Double_t exh[n] = {0};
  
  TMultiGraph *mgExcl = new TMultiGraph();
  gr0 = new TGraph(n,mass,limit);
  gr1 = new TGraphAsymmErrors(n,mass,limit,exl,exh,s1l,s1h);
  gr2 = new TGraphAsymmErrors(n,mass,limit,exl,exh,s2l,s2h);
  gr1->SetFillColor(kGreen);
  gr2->SetFillColor(kYellow);
  //gr3->SetFillColor(kRed);
  if(case!=3)gr0->SetMarkerStyle(20);
  if(case!=3)gr0->SetMarkerColor(kRed+1);
  if(case!=3)gr0->SetMarkerSize(0.8);
  gr0->SetLineStyle(2);
  gr0->SetLineWidth(2.8);
  gr0->SetLineColor(kRed+1);
  if(case!=3) gExCount0->SetLineStyle(1);
  gExCount0->SetLineColor(kBlue);
  gExCount0->SetLineWidth(2.8);
  gExCountUp->SetLineWidth(5);
  gExCountUp->SetLineStyle(1);
  gExCountUp->SetLineColor(kRed);
  gExCountUp->SetMarkerStyle(0);
  gExCountLow->SetLineColor(kRed);
  gExCountLow->SetLineWidth(2);
  if(case!=3) gExCountLow->SetMarkerStyle(0);
  if(case!=3) gExCountLow->SetLineStyle(1);

  if(case==3){
    gr0->SetLineStyle(1);
    gr0->SetMarkerStyle(0);
    gExCount->SetMarkerStyle(0);
    gExCount->SetLineStyle(2);
    gExCount->SetLineWidth(2.8);
    gExCount->SetLineColor(kBlack);
    gExCount0->SetLineColor(kRed+1);
    gExCount0->SetLineWidth(2.8);
    gExCount0->SetMarkerStyle(0);
    gExCount0->SetLineStyle(2);
    gExCountLow->SetLineColor(kBlue+1);
    gExCountLow->SetLineWidth(2.8);
    gExCountLow->SetMarkerStyle(0);
  }
  
  gDama1->SetFillColor(kRed-10);
  gDama1->SetLineColor(kRed-10);
  gDama2->SetFillColor(kRed-10);
  gDama2->SetLineColor(kRed-10);
  
  TCanvas *c1 = new TCanvas("c1","",700,700);
  //TH1F *frExcl = gPad->DrawFrame(0,0,1000,500);
  //frExcl->Draw();
  c1->SetLogx();
  c1->SetLogy();
  c1->SetGridx();
  c1->SetGridy();
  
  //if(case!=2) mgExcl->Add(gDama1,"F");
  //if(case!=2) mgExcl->Add(gDama2,"F");
  //mgExcl->Add(gr2);
  //mgExcl->Add(gr1);
  if(case==3) mgExcl->Add(gExCount);
  if(case==3)mgExcl->Add(gExCountLow);
  mgExcl->Add(gr0);
  //if(case!=3) mgExcl->Add(gExCountLow);
  mgExcl->Add(gExCount0);
  mgExcl->Draw("alp3");
 
  c1->Update();
  mgExcl->SetTitle("100nm - 10kg 10years");
  mgExcl->SetName("100nm - 10kg 10years");
  mgExcl->GetYaxis()->SetTitleOffset(1.7);
  mgExcl->GetXaxis()->SetTitleOffset(1.1);
  mgExcl->GetYaxis()->SetTitle("WIPM-nucleus cross section (cm^{2})");
  mgExcl->GetXaxis()->SetTitle("WIPM mass (GeV/c^{2})");
  mgExcl->GetXaxis()->CenterTitle();
  mgExcl->GetYaxis()->CenterTitle();
  mgExcl->GetYaxis()->SetRangeUser(1.0E-43,1.0E-37);
  c1->Modified();

  TLegend * legE = new TLegend();
  legE = new TLegend(0.70,0.60,0.9,0.80);
  legE->SetTextSize(0.02);
  legE->SetBorderSize(1);
  legE->SetLineStyle(0);
  legE->SetTextSize(0.03);
  legE->SetFillStyle(1001);
  legE->SetFillColor(kWhite);
  //legE->AddEntry(gr0,"Profile Likelihood","lp");
  if(case==3) legE->AddEntry(gExCount,"#mu_{b}=0","l");
  //legE->AddEntry(gExCountLow,"#mu_{b}=1 - Likelihood Method","l");
  if(case==3) legE->AddEntry(gr0,"#mu_{b}=10 - Likelihood","l");
  legE->AddEntry(gr0,"#mu_{b}=100 - Likelihood ","lp");
  //legE->AddEntry(gr0,"Median 90% C.L.","lp");
  //legE->AddEntry(gr1,"1#sigma band ","f");
  //legE->AddEntry(gr2,"2#sigma band ","f");
  legE->AddEntry(gExCount0,"#mu_{b}=100 - Counting only","l");
  //if(case==10) legE->AddEntry(gExCount0,"#mu_{b}=10 - Poisson Method","l");
  //legE->AddEntry(gExCountUp,"1#sigma band ","l");
  //legE->AddEntry(gExCount,"Zero Background","l");
  //if(case!=2) legE->AddEntry(gDama1,"DAMA/LIBRA","f");
  legE->Draw();
}

//-----------------------------------------------------------

void Analysis_Significance1(){ 
  
  // Upper limit on signal events
  const Int_t n = 9;
  Double_t lambda[n] = {0.1,0.2,0.3,0.4,0.5,0.6, 0.7, 0.8, 0.9};
  Double_t sig100[n] = {1.4,2.7,4.5,5.8,7.4,9.3,11.6,14.4,18.4};
  Double_t sig50[n]  = {1.0,1.9,2.9,4.2,5.5,6.9, 8.2,10.3,13.0};
  Double_t sig10[n]  = {0.5,1.1,1.4,1.9,2.4,2.8, 3.6, 4.6, 6.3};

   TMultiGraph *mg = new TMultiGraph();
   gr100 = new TGraph(n,lambda,sig100);
   gr50 = new TGraph(n,lambda,sig50);
   gr10 = new TGraph(n,lambda,sig10);
   gr100->SetLineStyle(1);
   gr100->SetLineColor(kBlack);
   gr100->SetMarkerColor(kBlack);
   gr100->SetMarkerStyle(20);
   gr50->SetLineStyle(2);
   gr50->SetLineColor(kRed);
   gr50->SetMarkerColor(kRed);
   gr50->SetMarkerStyle(21);
   gr10->SetLineStyle(9);
   gr10->SetLineColor(kBlue);
   gr10->SetMarkerColor(kBlue);
   gr10->SetMarkerStyle(22);
   gr100->SetLineWidth(2);
   gr50->SetLineWidth(2);
   gr10->SetLineWidth(2);

   TCanvas *c1 = new TCanvas("c1","",700,700);
    
   mg->Add(gr100);
   mg->Add(gr50);
   mg->Add(gr10);
   mg->Draw("alp3");
   c1->Update();
   mg->GetYaxis()->SetTitleOffset(1.1);
   mg->GetXaxis()->SetTitleOffset(1.1);
   mg->GetXaxis()->SetTitle("Signal fraction #lambda");
   mg->GetYaxis()->SetTitle("Significance");
   mg->GetXaxis()->CenterTitle();
   mg->GetYaxis()->CenterTitle();
   mg->GetXaxis()->SetRangeUser(0.1,0.9);
   mg->GetYaxis()->SetRangeUser(0,24);
   c1->Modified();
   
   TLegend * legE = new TLegend();
   legE = new TLegend(0.70,0.60,0.85,0.75);
   legE->SetTextSize(0.02);
   legE->SetBorderSize(1);
   legE->SetLineStyle(0);
   legE->SetTextSize(0.03);
   legE->SetFillStyle(1001);
   legE->SetFillColor(kWhite);
   legE->AddEntry(gr100,"N_{TOT}=100","lp");
   legE->AddEntry(gr50,"N_{TOT}=50","lp");
   legE->AddEntry(gr10,"N_{TOT}=10","lp");
   legE->Draw();
}

//-------------------------------------------------

void Analysis_Significance2(){ 
  
  // Upper limit on signal events
  /*const Int_t n = 19;
  Double_t Mass[n]    = {1000,800,600,400,300,200,150,100, 80, 60, 40,30,20,15,  10,   9,   8,     7,      6};
  Double_t Nsig1[n]   = {   4,  4,  4,3.5,  3,  3,  3,  3,  3,  3,  3, 3, 3, 3,   3, 2.8, 2.5,     2,      2};
  Double_t Nsig10[n]  = {  10, 10, 10, 10, 10, 10,  9,  9,  9,  9,  9, 9, 8, 8,   6, 5.5,   5,     4,      4};
  Double_t Nsig100[n] = {  28, 28, 28, 28, 27, 25, 25, 25, 25, 25, 25,23,20,17,  15,  14,  13,    12,     11};
  // expected sigma0
  Double_t N0[n]      = { 115,134,158,186,196,190,175,159,155,147,114,77,28, 8,0.43,0.14,0.03,4.17E-3,1.98E-4}; 
  */

  const Int_t n = 15; // filled on 14/07/2016 scatterig=YES
  Double_t Mass[n]    = {1000,800,600,400,300,200,150,100,  80, 60, 40,30,20,15,  10};
  Double_t Nsig1[n]   = {   5,  5,  5,  5,  5,  5,  5,   4,  4,  4,  4, 4, 4, 4,   4};
  Double_t Nsig10[n]  = {  11, 10, 10, 10, 10, 10,  10, 10,  9,  9,  9, 9, 9, 9,   9};
  Double_t Nsig100[n] = {  29, 28, 28, 28, 27, 27,  26, 26, 26, 26, 26,26,26,25,  25};
  // expected sigma0
  Double_t N0[n]      = { 115,134,158,186,196,190,175,159,155,147,114,77,28, 8,0.43}; 
  
  Double_t Nsig1_Poi[n];
  Double_t Nsig10_Poi[n];
  Double_t Nsig100_Poi[n];

  Double_t Sigma1[n];
  Double_t Sigma10[n];
  Double_t Sigma100[n];
  Double_t Sigma1_Poi[n];
  Double_t Sigma10_Poi[n];
  Double_t Sigma100_Poi[n];
  
  for(Int_t i=0; i<n; i++){
    Nsig1_Poi[i]= 5;
    Nsig10_Poi[i]= 12;
    Nsig100_Poi[i]= 33;
    Sigma1[i]=Nsig1[i]*(sigmaP/N0[i]);
    Sigma10[i]=Nsig10[i]*(sigmaP/N0[i]);
    Sigma100[i]=Nsig100[i]*(sigmaP/N0[i]);
    Sigma1_Poi[i]=Nsig1_Poi[i]*(sigmaP/N0[i]);
    Sigma10_Poi[i]=Nsig10_Poi[i]*(sigmaP/N0[i]);
    Sigma100_Poi[i]=Nsig100_Poi[i]*(sigmaP/N0[i]);
  }

  
  //------------------------------------
  //   Graphs on Sensitivity curve
  //------------------------------------
  TFile file_Dama1("RESULTS/Dama1.root");
  TGraph *gDama1 = file_Dama1.Get("Graph");
  
  TFile file_Dama2("RESULTS/Dama2.root");
  TGraph *gDama2 = file_Dama2.Get("Graph");
  
  TMultiGraph *mgS = new TMultiGraph();
   gr1S = new TGraph(n,Mass,Sigma1);
   gr10S = new TGraph(n,Mass,Sigma10);
   gr100S = new TGraph(n,Mass,Sigma100);
   gr1_PoiS = new TGraph(n,Mass,Sigma1_Poi);
   gr10_PoiS = new TGraph(n,Mass,Sigma10_Poi);
   gr100_PoiS = new TGraph(n,Mass,Sigma100_Poi);
    
   gr1S->SetLineStyle(1);
   gr1S->SetLineColor(kBlack);
   gr1S->SetMarkerColor(kBlack);
   gr1S->SetMarkerStyle(20);
   gr1S->SetMarkerSize(0.5);
   gr1_PoiS->SetLineStyle(2);
   gr1_PoiS->SetLineColor(kBlack);
     
   gr10S->SetLineStyle(1);
   gr10S->SetLineColor(kRed);
   gr10S->SetMarkerColor(kRed);
   gr10S->SetMarkerStyle(21);
   gr10S->SetMarkerSize(0.5);
   gr10_PoiS->SetLineStyle(2);
   gr10_PoiS->SetLineColor(kRed);
   
   gr100S->SetLineStyle(1);
   gr100_PoiS->SetLineStyle(2);
   gr100S->SetLineColor(kBlue);
   gr100S->SetMarkerColor(kBlue);
   gr100S->SetMarkerStyle(22);
   gr100S->SetMarkerSize(0.5);
   gr100_PoiS->SetLineColor(kBlue);

   gDama1->SetFillColor(kRed-10);
   gDama1->SetLineColor(kRed-10);
   gDama2->SetFillColor(kRed-10);
   gDama2->SetLineColor(kRed-10);
   
   TCanvas *cS = new TCanvas("cS","",700,700);
   cS->SetLogx();
   cS->SetLogy();
   cS->SetGridx();
   cS->SetGridy();
   //mgS->Add(gDama1,"F");
   //mgS->Add(gDama2,"F");
   mgS->Add(gr1S);
   mgS->Add(gr10S);
   mgS->Add(gr100S);
   mgS->Add(gr1_PoiS);
   mgS->Add(gr10_PoiS);
   mgS->Add(gr100_PoiS);

   mgS->Draw("alp3");
   //cS->Update();
   mgS->GetYaxis()->SetTitleOffset(1.7);
   mgS->GetXaxis()->SetTitleOffset(1.1);
   mgS->GetXaxis()->SetTitle("WIMP mass (Gev/c^{2})");
   mgS->GetYaxis()->SetTitle("WIPM-nucleus cross section (cm^{2})");
   mgS->GetXaxis()->CenterTitle();
   mgS->GetYaxis()->CenterTitle();
   mgS->GetYaxis()->SetRangeUser(1.0E-43,1.0E-37);
   //mg->GetXaxis()->SetRangeUser(0.1,0.9);
   //mg->GetYaxis()->SetRangeUser(0,40);
   cS->Modified();
   
   TLegend * legE = new TLegend();
   legE = new TLegend(0.70,0.60,0.95,0.85);
   legE->SetLineStyle(0);
   legE->SetTextSize(0.03);
   legE->SetFillStyle(1001);
   legE->SetFillColor(kWhite);
   legE->AddEntry(gr100_PoiS,"#mu_{b}=100 - Counting only","lp");
   legE->AddEntry(gr100S,    "#mu_{b}=100 - Likelihood","lp");
   legE->AddEntry(gr10_PoiS, "#mu_{b}=10   - Counting only","lp");
   legE->AddEntry(gr10S,     "#mu_{b}=10   - Likelihood","lp");
   legE->AddEntry(gr1_PoiS,  "#mu_{b}=1     - Counting only","lp");
   legE->AddEntry(gr1S,      "#mu_{b}=1     - Likelihood","lp");
   //legE->AddEntry(gDama1,"DAMA/LIBRA","f");

   legE->Draw();
   
  //------------------------------------
  //   Graphs on number of events
  //------------------------------------
  
   TMultiGraph *mg = new TMultiGraph();
   gr1 = new TGraph(n,Mass,Nsig1);
   gr10 = new TGraph(n,Mass,Nsig10);
   gr100 = new TGraph(n,Mass,Nsig100);
   gr1_Poi = new TGraph(n,Mass,Nsig1_Poi);
   gr10_Poi = new TGraph(n,Mass,Nsig10_Poi);
   gr100_Poi = new TGraph(n,Mass,Nsig100_Poi);
    
   gr1->SetLineStyle(1);
   gr1_Poi->SetLineStyle(2);
   gr1->SetLineColor(kBlack);
   gr1->SetMarkerColor(kBlack);
   gr1->SetMarkerStyle(20);
   gr1->SetMarkerSize(0.5);
   gr1_Poi->SetLineColor(kBlack);
     
   gr10->SetLineStyle(1);
   gr10_Poi->SetLineStyle(2);
   gr10->SetLineColor(kRed);
   gr10->SetMarkerColor(kRed);
   gr10->SetMarkerStyle(21);
   gr10->SetMarkerSize(0.5);
   gr10_Poi->SetLineColor(kRed);
   
   gr100->SetLineStyle(1);
   gr100_Poi->SetLineStyle(2);
   gr100->SetLineColor(kBlue);
   gr100->SetMarkerColor(kBlue);
   gr100->SetMarkerStyle(22);
   gr100->SetMarkerSize(0.5);
   gr100_Poi->SetLineColor(kBlue);
  
   TSpline3 *sNsig1 = new TSpline3("sNsig1",gr1->GetX(),gr1->GetY(),gr1->GetN());
   sNsig1->SetLineColor(kRed);
   
   TCanvas *c1 = new TCanvas("c1","",700,700);
   c1->SetLogx();
   //c1->SetLogy();
   //c1->SetGridx();
    //c1->SetGridy();
   mg->Add(gr1);
   mg->Add(gr10);
   mg->Add(gr100);
   mg->Add(gr1_Poi);
   mg->Add(gr10_Poi);
   mg->Add(gr100_Poi);
   mg->Draw("alp");
   //sNsig1->Draw("l same");
   //c1->Update();
   mg->GetYaxis()->SetTitleOffset(1.1);
   mg->GetXaxis()->SetTitleOffset(1.1);
   mg->GetXaxis()->SetTitle("WIMP mass (Gev/c^{2})");
   mg->GetYaxis()->SetTitle("Number of events");
   mg->GetXaxis()->CenterTitle();
   mg->GetYaxis()->CenterTitle();
   //mg->GetXaxis()->SetRangeUser(0.1,0.9);
   mg->GetYaxis()->SetRangeUser(0,40);
   c1->Modified();
   
   TLegend * legE = new TLegend();
   legE = new TLegend(0.70,0.60,0.95,0.85);
   legE->SetLineStyle(0);
   legE->SetTextSize(0.03);
   legE->SetFillStyle(1001);
   legE->SetFillColor(kWhite);
   legE->AddEntry(gr100_Poi,"#mu_{b}=100 - Counting only","lp");
   legE->AddEntry(gr100,    "#mu_{b}=100 - Likelihood","lp");
   legE->AddEntry(gr10_Poi, "#mu_{b}=10   - Counting only","lp");
   legE->AddEntry(gr10,     "#mu_{b}=10   - Likelihood","lp");
   legE->AddEntry(gr1_Poi,  "#mu_{b}=1     - Counting only","lp");
   legE->AddEntry(gr1,      "#mu_{b}=1     - Likelihood","lp");

   legE->Draw();
}
