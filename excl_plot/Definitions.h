#include <TMath.h>
#include <iostream>

//*****   PARAMETERS TO SET   *****//

const Float_t time      = 1*365; //days
const Float_t mass      = 400; //kg
const Float_t range_max = 1000; //nm 
Double_t rate_UL     = 2.44;
Bool_t verbosity        = true;
Int_t element_selection = 1;
// 0=all, 1=all but H, 2=heavy, 3=light, 4=S

//********************************//
// Functions implemented:
//
// DrawSigmaVSMass()
// DrawSigma()
//
//********************************//

time_t t;
struct tm *now;

Float_t v0       = 230;  //km/s
Float_t vE       = 244;  //km/s (lab velocity)
Float_t vRot     = 29.8; //km/s
Float_t vEsc     = 544;  //km/s *** CORRECT ***
//Float_t vEsc     = 600;  //km/s *** WRONG ***
Float_t c        = 299792.458; // km/s
Float_t z        = vEsc/v0;
Float_t sigma_v  = v0/sqrt(2);
Float_t N0       = 6.02E26; // kg^(-1) (Avogadro Number)
Float_t rhoD     = 0.4; // GeV/c2/cm3
Float_t Mproton  = 0.938; // proton mass
Float_t Mnucleon = 0.932; // nucleon mass to build target mass as MT=A*Mnucleon
Float_t Nseconds = 86400.0; //seconds in a day
long double sigmaP   = 1.0E-41; //WIMP-proton cross section (cm^2)
//Float_t sigmaP = 1.0; //WIMP-proton cross section (cm^2)
Float_t hbar     = 6.582E-16; //plank constant (eV*sec)
Float_t Pi       = TMath::Pi();
Float_t sqrtPi   = TMath::Sqrt(Pi);
Double_t rate = 0.;
Double_t diffRate = 0.;
Double_t cross    = 0.;
Float_t crossBis = 0.;
Double_t crossInv= 0.;
Double_t fake    = 0.;
Float_t weight   = 0.;

// ALL ELEMENTS
const Int_t Nelements = 8;
//Char_t  Name[Nelements]        ={    H,    C,    N,    O,    S,   Br,   Ag,    I};
Int_t   AMass[Nelements]         ={    1,   12,   14,   16,   32,   80,  108,  127};
Int_t   Charge[Nelements]        ={    1,    6,    7,    8,   16,   35,   47,   53};
Float_t MassFraction[Nelements]  ={0.016,0.101,0.027,0.074,0.003,0.320,0.440,0.019};
Float_t AtmFraction[Nelements]   ={0.410,0.214,0.049,0.118,0.003,0.100,0.100,0.004};

TF1 *fC = new TF1("f","0.4888-0.000775*x+7.7E-7*x*x",0,500);

TSpline3 *Line_C;
TSpline3 *Line_N;
TSpline3 *Line_O;
TSpline3 *Line_S;
TSpline3 *Line_Ag;
TSpline3 *Line_Br;
TSpline3 *Line_I;

Float_t Emin[Nelements]      ={0};
Float_t Emax[Nelements]      ={0};
Float_t RateE[Nelements]     ={0};
Float_t SigmaE[Nelements]    ={0};
Float_t DiffRateE[Nelements] ={0};

Float_t k0   = TMath::Power(Pi*v0*v0,3./2);
Float_t k1   = k0 * TMath::Erf(z)-2*z*TMath::Exp(-z*z)/sqrtPi;

v0   = v0/c;    //in units of c
vE   = vE/c;    //in units of c
vRot = vRot/c;  //in units of c
vEsc = vEsc/c;  //in units of c

Float_t E, v, vi;
Float_t MD; // WIMP mass
Float_t MT; // target mass
Float_t r, R0, E0, Emin, vmin;
Float_t muT, muP;
Float_t sigma0, resl;
Int_t A, MassNumber; // mass number 
Float_t DiffRate0, DiffRate1, DiffRate;
Float_t Lande, FormFactor, Sigma;

Float_t angle, angleCUT;
Float_t angle0, angle1, angleDM, angleDM0, energy, angle0res, angle1res;
Float_t fvi, fv0, fv1, fv2, fv3;
Float_t theta, phi, theta0, phi0;
Float_t vDM_x, vDM_y, vDM_z, vDM;
Float_t vCM_x, vCM_y, vCM_z, vCM;
Float_t vT_CM_x, vT_CM_y, vT_CM_z, vT_CM;
Float_t vT_x, vT_y, vT_z, vT;
Float_t p_x, p_y, p_z, pT;
Double_t v_g;
Double_t vx_g, vy_g, vz_g;

TRandom3 *rn = new TRandom();
TRandom3 *rand = new TRandom();

TGraph *grRate = new TGraph();
TGraph *grExcl = new TGraph();

TH1F *hv       = new TH1F("hv","",500,0,1000);
TH1F *hv0      = new TH1F("hv0","",100,0,600);
TH1F *hv0_400  = new TH1F("hv0_400","",100,0,600);
TH1F *hv0_500  = new TH1F("hv0_500","",100,0,600);
TH1F *hv0_600  = new TH1F("hv0_600","",100,0,600);
TH1F *hvDM_x   = new TH1F("hvDM_x","",400,-800,800);
TH1F *hvDM_y   = new TH1F("hvDM_y","",400,-800,800);
TH1F *hvDM_z   = new TH1F("hvDM_z","",400,-800,800);
TH1F *hvDM     = new TH1F("hvDM","",400,0,800);
TH1F *htDM     = new TH1F("htDM","",400,-4,4);
TH1F *htDM0    = new TH1F("htDM0","",400,-4,4);
TH1F *hvT_x    = new TH1F("hvT_x","",400,-1000,1000);
TH1F *hvT_y    = new TH1F("hvT_y","",400,-1000,1000);
TH1F *hvT_z    = new TH1F("hvT_z","",400,-1000,1000);
TH1F *hvT      = new TH1F("hvT","",400,0,1000);

TH1F *hp_x    = new TH1F("hp_x","",100,-50,50);
TH1F *hp_y    = new TH1F("hp_y","",100,-50,50);
TH1F *hp_z    = new TH1F("hp_z","",100,-50,50);
TH1F *hpT      = new TH1F("hpT","",100,0,80);

TH1F *htT0     = new TH1F("htT0","",200,-4,4);
TH1F *htT1     = new TH1F("htT1","",200,-4,4);
TH1F *htT0CUT  = new TH1F("htT0CUT","",60,-4,4);
TH1F *htT1CUT  = new TH1F("htT1CUT","",60,-4,4);
TH1F *htT0RES  = new TH1F("htT0RES","",60,-4,4);
TH1F *htT1RES  = new TH1F("htT1RES","",60,-4,4);
TH1F *he       = new TH1F("he","",200,0,1000);
TH1F *heCUT    = new TH1F("heCUT","",200,0,1000);
TH1F *ht       = new TH1F("ht","",50,-4,4);
TH1F *htCUT    = new TH1F("htCUT","",50,-4,4);
TH1F *he_heavy = new TH1F("he_heavy","",200,0,1000);
TH1F *he_light = new TH1F("he_light","",200,0,1000);
TH1F *ht_heavy = new TH1F("ht_heavy","",50,-4,4);
TH1F *ht_light = new TH1F("ht_light","",50,-4,4);
TH1F *he_heavyCUT = new TH1F("he_heavyCUT","",200,0,1000);
TH1F *he_lightCUT = new TH1F("he_lightCUT","",200,0,1000);
TH1F *ht_heavyCUT = new TH1F("ht_heavyCUT","",50,-4,4);
TH1F *ht_lightCUT = new TH1F("ht_lightCUT","",50,-4,4);
TProfile *het  = new TProfile("het","",20,0,500,-4,4);

TH1F *heCUT_C  = new TH1F("heCUT_C","",200,0,1000);
TH1F *heCUT_N  = new TH1F("heCUT_N","",200,0,1000);
TH1F *heCUT_O  = new TH1F("heCUT_O","",200,0,1000);
TH1F *heCUT_S  = new TH1F("heCUT_S","",200,0,1000);
TH1F *heCUT_Br = new TH1F("heCUT_Br","",200,0,1000);
TH1F *heCUT_Ag = new TH1F("heCUT_Ag","",200,0,1000);
TH1F *heCUT_I = new TH1F("heCUT_I","",200,0,1000);

//-----------------------------------

void ScatteringEvaluation(){

   const Int_t nL = 8;
   Double_t energyL[nL] ={   30,   40,   50,   70,  100,  200,  300,   900};
   Double_t sigmaC[nL] = {0.473,0.481,0.479,0.460,0.437,0.366,0.325, 0.320};
   Double_t sigmaO[nL] = {0.404,0.445,0.465,0.474,0.455,0.395,0.356, 0.300};
   Double_t sigmaN[nL] = {0.430,0.466,0.471,0.461,0.435,0.367,0.332, 0.330};
   Double_t sigmaS[nL] = {0.430,0.466,0.471,0.461,0.435,0.367,0.332, 0.330};


  const Int_t nH = 10;
  Double_t energyH[nH] = {  100,  150,  200,  250,  300,  350,  400,  450,  500,   900};
  Double_t sigmaBr[nH] = {0.270,0.303,0.314,0.332,0.337,0.328,0.331,0.331,0.325, 0.320};
  Double_t sigmaAg[nH] = {0.205,0.244,0.272,0.280,0.288,0.291,0.298,0.290,0.293, 0.290};
  Double_t sigmaI[nH] =  {0.205,0.244,0.272,0.280,0.288,0.291,0.298,0.290,0.293, 0.290};
  
  gC = new TGraph(nL,energyL,sigmaC);
  gN = new TGraph(nL,energyL,sigmaN);
  gO = new TGraph(nL,energyL,sigmaO);
  gS = new TGraph(nL,energyL,sigmaS);
  gBr = new TGraph(nH,energyH,sigmaBr);
  gAg = new TGraph(nH,energyH,sigmaAg);
  gI = new TGraph(nH,energyH,sigmaI);
  
  Line_C  = new TSpline3("sC",gC->GetX(),gC->GetY(),gC->GetN());
  Line_N  = new TSpline3("sN",gN->GetX(),gN->GetY(),gN->GetN());
  Line_O  = new TSpline3("sO",gO->GetX(),gO->GetY(),gO->GetN());
  Line_S  = new TSpline3("sS",gS->GetX(),gS->GetY(),gS->GetN());
  Line_Br = new TSpline3("sBr",gBr->GetX(),gBr->GetY(),gBr->GetN());
  Line_Ag = new TSpline3("sAg",gAg->GetX(),gAg->GetY(),gAg->GetN());
  Line_I = new TSpline3("sI",gI->GetX(),gI->GetY(),gI->GetN());
}

//-----------------------------------

void DrawSigmaVSMass(){

   TGraph *g1 = new TGraph();
   TGraph *g2 = new TGraph();
   TGraph *g3 = new TGraph();
   TGraph *g4 = new TGraph();
      
  const Int_t n = 19;
  Double_t mass[n]   =  { 1000,  800,  600,  400,  300,  200,  150,  100,   80,   60,   40,   30,   20,   15,   10,    9,    8,    7,    6};
  // cut 50 nm
  Double_t sigma[n]   = {0.897,0.887,0.871,0.841,0.815,0.766,0.725,0.652,0.614,0.606,0.688,0.688,0.626,0.568,0.476,0.429,0.384,0.337,0.230};
  // cut 100 nm
  Double_t sigma1[n]  = {0.715,0.707,0.696,0.676,0.659,0.628,0.613,0.614,0.622,0.608,0.575,0.542,0.468,0.408,0.280,0.245,0.214,0.181,0.142};
  Double_t sigma2[n]  = {0.718,0.710,0.696,0.674,0.650,0.603,0.570,0.551,0.584,0.603,0.575,0.541,0.478,0.401,0.339,0.270,0.250,0.230,0.210};
  Double_t sigma3[n]  = {0.715,0.709,0.694,0.668,0.643,0.593,0.569,0.553,0.574,0.607,0.580,0.544,0.477,0.409,0.267,0.236,0.208,0.176,0.139};
  Double_t sigma4[n]  = {0.706,0.697,0.682,0.650,0.622,0.570,0.541,0.528,0.561,0.589,0.544,0.504,0.426,0.350,0.200,0.168,0.139,0.107,0.075};

  Double_t em[n]      = {0};
  Double_t err[n]     = {0.006,0.006,0.006,0.007,0.008,0.009,0.010,0.010,0.010,0.020,0.020,0.020,0.030,0.030,0.040,0.040,0.040,0.040,0.040};

  TMultiGraph *mg = new TMultiGraph();

  g  = new TGraphErrors(n,mass,sigma,em,err);
  g1 = new TGraphErrors(n,mass,sigma1,em,err);
  g2 = new TGraphErrors(n,mass,sigma2,em,err);
  g3 = new TGraphErrors(n,mass,sigma3,em,err);
  g4 = new TGraphErrors(n,mass,sigma4,em,err);

  g->SetMarkerStyle(20);
  g->SetLineColor(kBlack);
  g->SetMarkerColor(kBlack);
  g->SetMarkerSize(0.8);
  g1->SetMarkerStyle(20);
  g1->SetLineColor(kBlack);
  g1->SetMarkerColor(kBlack);
  g1->SetMarkerSize(0.8);
  g2->SetMarkerStyle(20);
  g2->SetLineColor(kBlue);
  g2->SetMarkerColor(kBlue);
  g2->SetMarkerSize(0.8);
  g3->SetMarkerStyle(20);
  g3->SetLineColor(kRed);
  g3->SetMarkerColor(kRed);
  g3->SetMarkerSize(0.8);
  g4->SetMarkerStyle(20);
  g4->SetLineColor(kGreen-2);
  g4->SetMarkerColor(kGreen-2);
  g4->SetMarkerSize(0.8);
  
  TCanvas *c1 = new TCanvas("c1","",700,700);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetLogx();
  mg->Add(g);
  //mg->Add(g1);
  //mg->Add(g2);
  //mg->Add(g3);
  mg->Add(g4);
  mg->Draw("alp");
   
  c1->Update();
  mg->GetYaxis()->SetTitle("#sigma_{#theta} (rad)");
  mg->GetXaxis()->SetTitle("WIPM mass (GeV/c^{2})");
  c1->Modified();

  TLegend * legE = new TLegend();
  legE = new TLegend(0.70,0.60,1.1,0.75);
  legE->SetTextSize(0.02);
  legE->SetBorderSize(1);
  legE->SetLineStyle(0);
  legE->SetTextSize(0.02);
  legE->SetFillStyle(1001);
  legE->SetFillColor(kWhite);
  //legE->AddEntry(g1,"w=MassFraction,    vEsc = 600","lp");
  //legE->AddEntry(g3,"w=A*A*AtmFraction, vEsc = 600","lp");
  //legE->AddEntry(g4,"w=A*A*AtmFraction, vEsc = 544","lp");
  legE->AddEntry(g,"cut = 50 nm","lp");
  legE->AddEntry(g4,"cut = 100 nm","lp");
  legE->Draw();
  
}

 //-------------------------

void DrawSigma(){
  TF1 *g10 = new TF1("g10","TMath::Gaus(x,0,0.200,1)",0,TMath::Pi()/2.);
  TF1 *g50 = new TF1("g50","TMath::Gaus(x,0,0.400,1)",0,TMath::Pi()/2.);
  TF1 *g100 = new TF1("g100","TMath::Gaus(x,0,0.600,1)",0,TMath::Pi()/2.);
  TF1 *g1000 = new TF1("g1000","TMath::Gaus(x,0,0.800,1)",0,TMath::Pi()/2.);

  g10->SetLineColor(kRed);
  g10->SetLineWidth(0.5);
  g10->SetMarkerColor(kRed);
  g10->SetMarkerStyle(20);
  g10->SetMarkerSize(0.9);
  g10->SetNpx(20);
  g50->SetLineColor(kGreen);
  g50->SetLineWidth(0.5);
  g50->SetMarkerColor(kGreen);
  g50->SetMarkerStyle(21);
  g50->SetMarkerSize(0.8);
  g50->SetNpx(20);
  g100->SetLineColor(kBlue);
  g100->SetLineWidth(0.5);
  g100->SetMarkerColor(kBlue);;
  g100->SetMarkerStyle(22);
  g100->SetMarkerSize(0.7);
  g100->SetNpx(20);
  g1000->SetLineColor(kMagenta);
  g1000->SetLineWidth(0.5);
  g1000->SetMarkerColor(kMagenta);
  g1000->SetMarkerStyle(23);
  g1000->SetMarkerSize(1);
  g1000->SetNpx(20);

  TGraph *f10 = new TGraph(g10);
  TGraph *f50 = new TGraph(g50);
  TGraph *f100 = new TGraph(g100);
  TGraph *f1000 = new TGraph(g1000);
  
  TMultiGraph *mgRate = new TMultiGraph();
  mgRate->Add(f10);
  mgRate->Add(f50);
  mgRate->Add(f100);
  mgRate->Add(f1000);
  
  TCanvas *c = new TCanvas("c","",400,400);
  mgRate->Draw("ALP");

  mgRate->GetXaxis()->SetTitle("2D recoil angle (rad)");
  mgRate->GetXaxis()->CenterTitle();
  
  TLegend * legE = new TLegend();
  legE = new TLegend(0.50,0.60,1.0,0.80);
  legE->SetTextSize(0.04);
  legE->SetBorderSize(1);
  legE->SetLineStyle(0);
  legE->SetFillStyle(1001);
  legE->SetFillColor(kWhite);
  legE->AddEntry(g10,"M_{W} = 10 GeV/c^{2}","lp");
  legE->AddEntry(g50,"M_{W} = 50 GeV/c^{2}","lp");
  legE->AddEntry(g100,"M_{W} = 100 GeV/c^{2}","lp");
  legE->AddEntry(g1000,"M_{W} = 1000 GeV/c^{2}","lp");
  legE->Draw();
}
