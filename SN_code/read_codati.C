#include "Riostream.h"
void read_codati() {

   TString dir = gSystem->UnixPathName(__FILE__);
   dir.ReplaceAll("read_codati.C","");
   dir.ReplaceAll("/./","/");
   ifstream in;
   in.open(Form("%scodati.dat",dir.Data()));

   Float_t thr,C_ev,N_ev,O_ev,S_ev,Br_ev,Ag_ev,I_ev;
   Int_t nlines = 0;
  
   TMultiGraph *mg = new TMultiGraph();
   TGraph *grC = new TGraph();
   TGraph *grN = new TGraph();
   TGraph *grO = new TGraph();
   TGraph *grS = new TGraph();
   TGraph *grBr = new TGraph();
   TGraph *grAg = new TGraph();
   TGraph *grI = new TGraph();

   grC->SetName("grC");
   grN->SetName("grN");
   grO->SetName("grO");
   grS->SetName("grS");
   grBr->SetName("grBr");
   grAg->SetName("grAg");
   grI->SetName("grI");
   

   while (1) {
     in >> thr >> C_ev >> N_ev >> O_ev >> S_ev >> Br_ev >> Ag_ev >> I_ev;
      if (!in.good()) break;
      if (nlines < 5) printf("x=%8f, y=%8f, z=%8f\n",thr,C_ev,N_ev);
      //h1->Fill(x);
      //ntuple->Fill(x,y,z);
      grC->SetPoint(nlines,thr,C_ev);
      grN->SetPoint(nlines,thr,N_ev);
      grO->SetPoint(nlines,thr,O_ev);
      grS->SetPoint(nlines,thr,S_ev);
      grBr->SetPoint(nlines,thr,Br_ev);
      grAg->SetPoint(nlines,thr,Ag_ev);
      grI->SetPoint(nlines,thr,I_ev);     
      nlines++;
   }
   printf(" found %d points\n",nlines);

   grC->SetLineColor(1);
   grN->SetLineColor(2);
   grO->SetLineColor(3);
   grS->SetLineColor(4);
   grBr->SetLineColor(5);
   grAg->SetLineColor(6);
   grI->SetLineColor(7);
   grC->SetFillColor(0);
   grN->SetFillColor(0);
   grO->SetFillColor(0);
   grS->SetFillColor(0);
   grBr->SetFillColor(0);
   grAg->SetFillColor(0);
   grI->SetFillColor(0);

   mg->Add(grC);
   mg->Add(grN);
   mg->Add(grO);
   mg->Add(grS);
   mg->Add(grBr);
   mg->Add(grAg);
   mg->Add(grI);

   mg->Draw("AL");

   in.close();

   //f->Write();
}
