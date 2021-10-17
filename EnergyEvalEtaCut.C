/*
- Energy plots
- Eta cuts
- Individual detector
- Standalone
- Output file - energy_verification_EtaCut_detector.root
- Plots made - [te/ge vs ge], [ce/ge vs ge], [ce-te vs ge], [(ce-ge)/ge vs ge], [counts for te], 
  [counts for ce], [counts for ge]. (debug not implemented)
*/

/*
  authors - Sagar Joshi      (ee190002054@iiti.ac.in)
            Siddhant Rathi   (me190003061@iiti.ac.in)
	    
  version - 1.0

*/


#include <iostream>
#include <eicqa_modules/EvalRootTTree.h>
#include <eicqa_modules/EvalHit.h>
#include "TMath.h"
#include "TStyle.h"

R__LOAD_LIBRARY(libeicqa_modules.so)

void EnergyEvalEtaCut(TString detector, int print = 1, int debug = 0)
{
  
  TFile *f1 = new TFile("merged_Eval_" + detector + ".root","READ"); //change the input evaluator root file here
 
  TTree* T = (TTree*)f1->Get("T");
  EvalRootTTree *evaltree = nullptr;
  
  TCanvas *c = new TCanvas();

  // Changing default plotting style 
  gStyle->SetTitleW(0.8);
  gStyle->SetTitleH(0.1);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0102);
  gStyle->SetTitleXOffset(1);
  gStyle->SetTitleYOffset(1);
  gStyle->SetLabelSize(0.05);  
  gStyle->SetTitleXSize(0.05);  
  gStyle->SetTitleYSize(0.05);

  c->SetTickx();
  c->SetTicky();

  // Initialsing Hisogram Objects  
  TH2D *te_by_ge_ge_EtaCut = new TH2D("te_by_ge_ge_EtaCut","te/ge vs ge",200,0,30,200,-0.5,1.5); 
  TH2D *ce_by_ge_ge_EtaCut = new TH2D("ce_by_ge_ge_EtaCut","ce/ge vs ge",200,0,30,200,-0.5,1.5);
  TH2D *ce_minus_te_ge_EtaCut = new TH2D("ce_minus_te_ge_EtaCut","ce-te vs ge",200,0,30,200,-2,0.5);
  TH2D *ce_minus_ge_by_ge_ge_EtaCut = new TH2D("ce_minus_ge_by_ge_ge_EtaCut","#frac{#Delta e}{truth e} vs truth e",200,0,30,200,-2,1);
  TH1D *counts_te_EtaCut = new TH1D("counts_te_EtaCut","tower_e",200,0,40);
  TH1D *counts_ce_EtaCut = new TH1D("counts_ce_EtaCut","cluster_e",200,0,40);
  TH1D *counts_ge_EtaCut = new TH1D("counts_ge_EtaCut","true_e",200,0,40);

  long double total_te = 0;
  long double total_ce = 0; 
  long double total_ge = 0;   

  Double_t te = 0;
  Double_t ce = 0;

  Double_t eta_min, eta_max;
  TString cut_text;

  // Applying detector-specific Eta Cuts

  if (detector == "FEMC"){
    eta_min = 1.3;
    eta_max = 3.3;
    cut_text = " {1.3 < geta < 3.3} ";
  }

  else if(detector == "EEMC"){
    eta_min = -3.5;
    eta_max = -1.7;
    cut_text = " {-3.5 < geta < -1.7} ";
  }

  else if(detector == "CEMC"){
    eta_min = -1.5;
    eta_max = 1.2;
    cut_text = " {-1.5 < geta < 1.2} ";
  }

  else if(detector == "FHCAL"){
    eta_min = 1.2;
    eta_max = 3.5;
    cut_text = " {1.2 < geta < 3.5} ";
  }

  else if(detector == "HCALIN"){
    eta_min = -1.1;
    eta_max = 1.1;
    cut_text = " {-1.1 < geta < 1.1} ";
  }

  else if(detector == "HCALOUT"){
    eta_min = -1.1;
    eta_max = 1.1;
    cut_text = " {-1.1 < geta < 1.1} ";
  }

  else{
    std::cout<<"Please try again.";
    return 1;
  }
 
  T->SetBranchAddress("DST#EvalTTree_" + detector,&evaltree);  

  for(int i=0; i<T->GetEntries(); i++){
    T->GetEntry(i);


    cout << "ge: " << evaltree->get_ge() << endl;
   
    Double_t ge = evaltree->get_ge();
    Double_t geta = evaltree->get_geta();
    te = 0; 
   
    for (int k=0; k<evaltree->get_ntowers(); k++)
      {
	EvalTower *twr = evaltree->get_tower(k);
	if (twr)
	  {;
	    te += twr->get_te();
	  }
      }
   
    ce = 0;
    for (int j=0; j<evaltree->get_nclusters(); j++) // changed iteration varaible from i to j 
      {
	EvalCluster *clus = evaltree->get_cluster(j);
	if (clus)
	  {
	    ce += clus->get_ce();           
	  }
      }
    if(geta>=eta_min && geta<=eta_max)
      {
	ce_minus_ge_by_ge_ge_EtaCut->Fill(ge, (ce-ge)/ge);

	te_by_ge_ge_EtaCut->Fill(ge, te/ge);

	ce_by_ge_ge_EtaCut->Fill(ge, ce/ge);

	ce_minus_te_ge_EtaCut->Fill(ge, ce-te);

	counts_te_EtaCut->Fill(te);
            
	counts_ce_EtaCut->Fill(ce);

	counts_ge_EtaCut->Fill(ge);
	   
	total_ce += ce;
	total_ge += ge;
	total_te += te;

      }          
  }


  te_by_ge_ge_EtaCut->GetXaxis()->SetTitle("Generated Energy (GeV)");
  te_by_ge_ge_EtaCut->GetXaxis()->SetLabelSize(0.05);  
  te_by_ge_ge_EtaCut->GetXaxis()->SetTitleSize(0.05);
  te_by_ge_ge_EtaCut->GetYaxis()->SetTitle("te_{agg}/ge");
  te_by_ge_ge_EtaCut->GetYaxis()->SetLabelSize(0.05);  
  te_by_ge_ge_EtaCut->GetYaxis()->SetTitleSize(0.05);
  te_by_ge_ge_EtaCut->SetTitle("te_{agg}/ge vs ge" + cut_text);
    
  ce_by_ge_ge_EtaCut->GetXaxis()->SetTitle("Generated Energy (GeV)");
  ce_by_ge_ge_EtaCut->GetXaxis()->SetLabelSize(0.05);  
  ce_by_ge_ge_EtaCut->GetXaxis()->SetTitleSize(0.05);
  ce_by_ge_ge_EtaCut->GetYaxis()->SetTitle("ce_{agg}/ge"); 
  ce_by_ge_ge_EtaCut->GetYaxis()->SetLabelSize(0.05);  
  ce_by_ge_ge_EtaCut->GetYaxis()->SetTitleSize(0.05);
  ce_by_ge_ge_EtaCut->SetTitle("ce_{agg}/ge vs ge" + cut_text);

  ce_minus_te_ge_EtaCut->GetXaxis()->SetTitle("ge (GeV)");
  ce_minus_te_ge_EtaCut->GetXaxis()->SetLabelSize(0.05);  
  ce_minus_te_ge_EtaCut->GetXaxis()->SetTitleSize(0.05);
  ce_minus_te_ge_EtaCut->GetYaxis()->SetTitle("ce-te");
  ce_minus_te_ge_EtaCut->GetYaxis()->SetLabelSize(0.05);  
  ce_minus_te_ge_EtaCut->GetYaxis()->SetTitleSize(0.05);
  ce_minus_te_ge_EtaCut->SetTitle("ce_{agg}-te_{agg} vs ge" + cut_text);

  counts_te_EtaCut->GetXaxis()->SetTitle("te 1d");
  counts_te_EtaCut->GetXaxis()->SetLabelSize(0.05);  
  counts_te_EtaCut->GetXaxis()->SetTitleSize(0.05);
  counts_te_EtaCut->GetYaxis()->SetTitle("counts");
  counts_te_EtaCut->GetYaxis()->SetLabelSize(0.05);  
  counts_te_EtaCut->GetYaxis()->SetTitleSize(0.05);
  counts_te_EtaCut->SetTitle("te_1d" + cut_text);

  counts_ce_EtaCut->GetXaxis()->SetTitle("ce 1d");
  counts_ce_EtaCut->GetXaxis()->SetLabelSize(0.05);  
  counts_ce_EtaCut->GetXaxis()->SetTitleSize(0.05);
  counts_ce_EtaCut->GetYaxis()->SetTitle("counts");
  counts_ce_EtaCut->GetYaxis()->SetLabelSize(0.05);  
  counts_ce_EtaCut->GetYaxis()->SetTitleSize(0.05);
  counts_ce_EtaCut->SetTitle("ce_1d" + cut_text);

  counts_ge_EtaCut->GetXaxis()->SetTitle("ge 1d");
  counts_ge_EtaCut->GetXaxis()->SetLabelSize(0.05);  
  counts_ge_EtaCut->GetXaxis()->SetTitleSize(0.05);
  counts_ge_EtaCut->GetYaxis()->SetTitle("counts");
  counts_ge_EtaCut->GetYaxis()->SetLabelSize(0.05);  
  counts_ge_EtaCut->GetYaxis()->SetTitleSize(0.05);
  counts_ge_EtaCut->SetTitle("ge_1d" + cut_text);

  ce_minus_ge_by_ge_ge_EtaCut->GetXaxis()->SetTitle("Generated Energy (GeV)");
  ce_minus_ge_by_ge_ge_EtaCut->GetXaxis()->SetLabelSize(0.05);  
  ce_minus_ge_by_ge_ge_EtaCut->GetXaxis()->SetTitleSize(0.05);
  ce_minus_ge_by_ge_ge_EtaCut->GetYaxis()->SetTitle("(ce_{agg}-ge)/ge");
  ce_minus_ge_by_ge_ge_EtaCut->GetYaxis()->SetLabelSize(0.05);  
  ce_minus_ge_by_ge_ge_EtaCut->GetYaxis()->SetTitleSize(0.05);
  ce_minus_ge_by_ge_ge_EtaCut->SetTitle("(#Delta e/ true_e) vs true_e" + cut_text);
    

  TFile *f = new TFile("energy_verification_EtaCut_" + detector + ".root","RECREATE"); //change the calorimeter name here
  f->GetList()->Add(ce_minus_ge_by_ge_ge_EtaCut);
  f->GetList()->Add(ce_by_ge_ge_EtaCut);
  f->GetList()->Add(te_by_ge_ge_EtaCut);
  f->GetList()->Add(ce_minus_te_ge_EtaCut);
  f->GetList()->Add(counts_te_EtaCut);
  f->GetList()->Add(counts_ce_EtaCut);
  f->GetList()->Add(counts_ge_EtaCut);     
  f->Write();

  if(print == 1){

    ce_minus_ge_by_ge_ge_EtaCut->Draw("colz");
    c->Print(detector + "ce_minus_ge_by_ge_ge_EtaCut.png");
    ce_by_ge_ge_EtaCut->Draw("colz");
    c->Print(detector + "ce_by_ge_ge_EtaCut.png");
    te_by_ge_ge_EtaCut->Draw("colz");
    c->Print(detector + "te_by_ge_ge_EtaCut.png");
    ce_minus_te_ge_EtaCut->Draw("colz");
    c->Print(detector + "ce_minus_te_ge_EtaCut.png");
    counts_te_EtaCut->Draw("colz");
    c->Print(detector + "counts_te_EtaCut.png");
    counts_ce_EtaCut->Draw("colz");
    c->Print(detector + "counts_ce_EtaCut.png");
    counts_ge_EtaCut->Draw("colz");
    c->Print(detector + "counts_ge_EtaCut.png");

  }

  std::cout << "The total ce is: " << total_ce << endl;
  std::cout << "The total te is: " << total_te << endl;
  std::cout << "The total ge is: " << total_ge << endl;

  c->Close();
}
