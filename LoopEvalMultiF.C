/*
- Energy plots including [sigma_e vs ge] plots
- Eta Cuts
- Combined FEMC+FHCAL
- Standalone
- Output root file - energy_verification_EtaCut_FHCAL_FEMC.root
- Plots made -  [te/ge vs ge], [ce/ge vs ge], [ce-te vs ge], [(ce-ge)/ge vs ge], [counts for te], [counts for ce], 
  [counts for ge], [sigma_e vs ge], [mean_e vs ge], [chi2_e vs ge], slices for the [(ce-ge)/ge vs ge] plot.
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

void LoopEvalMultiF(int print = 1, int debug = 0)
{

  int nSlices = 15; // Number of ge-axis slices fitted while making sigma_e vs e plot
  TString detector = "FHCAL_FEMC";
  TFile *f1 = new TFile("merged_Eval_FHCAL.root","READ"); //change the input evaluator root file here
  TFile *f2 = new TFile("merged_Eval_FEMC.root","READ"); //change the input evaluator root file here

  TTree* T1 = (TTree*)f1->Get("T");
  EvalRootTTree *evaltree1 = nullptr;
  
  TTree* T2 = (TTree*)f2->Get("T");
  EvalRootTTree *evaltree2 = nullptr;
  
  TCanvas *c = new TCanvas();
  c->SetTickx();
  c->SetTicky();
  
  // Modifying default plotting style
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(102);
  gStyle->SetTitleXOffset(1);
  gStyle->SetTitleYOffset(1);
  gStyle->SetLabelSize(0.05);  
  gStyle->SetTitleXSize(0.05);  
  gStyle->SetTitleYSize(0.05);


  // Initialising histogram objects
  TH2D *ce_minus_te_ge_EtaCut = new TH2D("ce_minus_te_ge_EtaCut","ce-te vs ge",200,0,30,200,-2,1);
  TH1D *counts_ge_EtaCut = new TH1D("counts_ge_EtaCut","true_e",200,0,40);

  TH2D *te_by_ge_ge_EtaCut = new TH2D("te_by_ge_ge_EtaCut","te/ge vs ge",200,0,30,200,-0.5,1.5); 
  TH2D *ce_by_ge_ge_EtaCut = new TH2D("ce_by_ge_ge_EtaCut","ce/ge vs ge",200,0,30,200,-0.5,1.5);
  TH2D *ce_minus_ge_by_ge_ge_EtaCut = new TH2D("ce_minus_ge_by_ge_ge_EtaCut","#frac{#Delta e}{truth e} vs truth e",200,0,30,200,-2,1);
  TH2D *ce_minus_ge_by_ge_ge_EtaCut_temp = new TH2D("ce_minus_ge_by_ge_ge_EtaCut_temp","#frac{#Delta e}{truth e} vs truth e",nSlices,0,30,2000,-1.5,0.5);

  TH1D *counts_te_EtaCut = new TH1D("counts_te_EtaCut","tower_e",200,0,40);
  TH1D *counts_ce_EtaCut = new TH1D("counts_ce_EtaCut","cluster_e",200,0,40);
  
  long double total_te = 0;
  long double total_ce = 0; 
  long double total_ge = 0; 

  Double_t eta_min = 1.3;
  Double_t eta_max = 3.3;
  TString cut_text = " {1.3 < geta < 3.3} ";
 
  T1->SetBranchAddress("DST#EvalTTree_FHCAL",&evaltree1);
  T2->SetBranchAddress("DST#EvalTTree_FEMC",&evaltree2);
 
  for(int i=0; i<T1->GetEntries(); i++){ // We assume same no. of entries, since no cuts are applied

      T1->GetEntry(i);
      T2->GetEntry(i);

      if(debug==1){
	std::cout<<"\n\n\n------------------------------------------\nParticle: "<<i<<"\n\n";
	std::cout<<"Initial Parameters "<<"\n";
      }

      Double_t geta1 = evaltree1->get_geta();

      if(debug==1){      
	std::cout<<"geta: "<<geta1<<"\n";
      }

      if(geta1>=eta_min && geta1<=eta_max){
    
	if(debug==1){
	  cout<<"\ngeta cut applied (1.3, 3.3)"<<"\n\n";
	}

	Double_t ge = evaltree1->get_ge();

	counts_ge_EtaCut->Fill(ge); 

	if(debug==1){
	  std::cout<<"ge: "<<ge<<"\n";
	}

	total_ge += ge;

	if(debug==1){
	  std::cout<<"total_ge till now = "<<total_ge<<"\n";
	}

	Double_t te_aggregate = 0;    
	for (int j=0; j<evaltree1->get_ntowers(); j++){
	  
	  if(debug==1){
	    std::cout<<"\nFHCAL Tower: "<<j<<"\n";
	  }

	  Double_t te = 0;

	  EvalTower *twr1 = evaltree1->get_tower(j);

	  if (twr1){

	    if(debug==1){
	      cout<<"non-empty FHCAL tower\n";
	    }

	    te = twr1->get_te();
	    te_aggregate += te;
	    if(debug==1){
	      cout<<"te_aggregate till now = "<<te_aggregate<<"\n";
	      std::cout<<"te += "<<twr1->get_te()<<"\n";
	    }
	  }
	}
       
	for (int j=0; j<evaltree2->get_ntowers(); j++){

	  if(debug==1){
	    std::cout<<"\nFEMC Tower: "<<j<<"\n";
	  }

	  Double_t te = 0;

	  EvalTower *twr2 = evaltree2->get_tower(j);

	  if (twr2){
	    
	    if(debug==1){
	      cout<<"non-empty FEMC tower\n";
	    }

	    te = twr2->get_te();
	    te_aggregate += te;

	    if(debug==1){
	      cout<<"te_aggregate till now = "<<te_aggregate<<"\n";
	      std::cout<<"te += "<<twr2->get_te()<<"\n";
	    }
	  }
	}  
      
	total_te += te_aggregate;
	if(debug==1){
	  cout<<"total_te till now = "<<total_te<<"\n";
	  cout<<"\n";
	}

	Double_t ce_aggregate = 0;
	
	if(debug==1){
	  cout<<"\nFHCAL clusters "<<evaltree1->get_nclusters()<<"\n";
	  cout<<"\nFEMC clusters "<<evaltree2->get_nclusters()<<"\n";
	}
	
	for (int j=0; j<evaltree1->get_nclusters(); j++){
	 
	  if(debug==1){
	    std::cout<<"\nFHCAL Cluster: "<<j<<"\n";
	  }
	  
	  Double_t ce = 0;

	  EvalCluster *clus1 = evaltree1->get_cluster(j);
	  
	  if (clus1){
	  
	    if(debug==1){
	      cout<<"non-empty FHCAL cluster\n";
	    }

	    ce = clus1->get_ce();
	    ce_aggregate += ce;

	    if(debug==1){
	      cout<<"ce_aggregate till now = "<<ce_aggregate<<"\n";
	      std::cout<<"ce += "<<clus1->get_ce()<<"\n";
	    }
	  }
	}

	for (int j=0; j<evaltree2->get_nclusters(); j++){
	 
	  if(debug==1){
	    std::cout<<"\nCluster: "<<j<<"\n";
	  }
	  Double_t ce = 0;
       
	  EvalCluster *clus2 = evaltree2->get_cluster(j);
	  if (clus2){

	    if(debug==1){
	      cout<<"non-empty FEMC cluster\n";
	    }

	    ce += clus2->get_ce();
	    ce_aggregate += ce;

	    if(debug==1){
	      cout<<"ce_aggregate till now = "<<ce_aggregate<<"\n";
	      std::cout<<"ce += "<<clus2->get_ce()<<"\n";
	    }
	  }
	}  

	total_ce += ce_aggregate;
	if(debug==1){
	  cout<<"total_ce till now = "<<total_ce<<"\n";
	}

	counts_te_EtaCut->Fill(te_aggregate);
	if(debug==1){
	  cout<<"te_aggregate: ("<<te_aggregate<<")\n";
	}

	te_by_ge_ge_EtaCut->Fill(ge, te_aggregate/ge);
	if(debug==1){
	  cout<<"(ge, te_aggregate/ge): ("<<ge<<", "<<te_aggregate/ge<<")\n";
	}

	counts_ce_EtaCut->Fill(ce_aggregate);
	if(debug==1){
	  cout<<"ce_aggregate: ("<<ce_aggregate<<")\n";
	}

	ce_by_ge_ge_EtaCut->Fill(ge, ce_aggregate/ge);
	if(debug==1){
	  cout<<"(ge, ce_aggregate/ge): ("<<ge<<", "<<ce_aggregate/ge<<")\n";
	}	

	ce_minus_ge_by_ge_ge_EtaCut_temp->Fill(ge, (ce_aggregate-ge)/ge);
	ce_minus_ge_by_ge_ge_EtaCut->Fill(ge, (ce_aggregate-ge)/ge);
	if(debug==1){
	  cout<<"(ge, (ce_aggregate-ge)/ge): ("<<ge<<", "<<(ce_aggregate-ge)/ge<<")\n";
	}

	ce_minus_te_ge_EtaCut->Fill(ge, ce_aggregate-te_aggregate);
	if(debug==1){
	  cout<<"(ge, (ce_aggregate-te_aggregate)): ("<<ge<<", "<<ce_aggregate-te_aggregate<<")\n";
	}      
      }
      
  }

  if(debug==1){
    cout<<"\neta cut if ends"<<"\n";
  }


  TString arr[nSlices]; // Used for naming slices
  for(int sno = 0; sno < nSlices; sno++){
    arr[sno] = TString::Itoa(sno + 1, 10);
  }

  // Generating sigma_e vs ge plot
  ce_minus_ge_by_ge_ge_EtaCut_temp->FitSlicesY(0, 0, -1, 0, "QN");
  TH2F *ce_minus_ge_by_ge_ge_EtaCut_temp_2 = (TH2F*)gDirectory->Get("ce_minus_ge_by_ge_ge_EtaCut_temp_2");
  TH2F *ce_minus_ge_by_ge_ge_EtaCut_temp_1 = (TH2F*)gDirectory->Get("ce_minus_ge_by_ge_ge_EtaCut_temp_1");
  TH2F *ce_minus_ge_by_ge_ge_EtaCut_temp_chi2 = (TH2F*)gDirectory->Get("ce_minus_ge_by_ge_ge_EtaCut_temp_chi2");

  TH1D* slices[nSlices];

  // Generating plots for individual slices
  for(int sno = 0; sno < nSlices; sno++){
    int plusOne = sno+1;
    TString sname = "slice " + arr[sno];
    slices[sno] = ce_minus_ge_by_ge_ge_EtaCut_temp->ProjectionY(sname, plusOne, plusOne);
  }

  // Initialising fit functions 
  TF1 *fit = new TF1("fit", "gaus");   
  TF1 *fExp = new TF1("fExp","0.1 + 0.5/sqrt(x)",0,30);
  TF1 *fTrue = new TF1("fTrue","[0] + [1]/sqrt(x)",0,30);   
     
  if(debug==1){
    std::cout<<"\nHistogram Formatting\n\n";
  }

  ce_minus_te_ge_EtaCut->GetXaxis()->SetTitle("Generated Energy (GeV)");
  ce_minus_te_ge_EtaCut->GetXaxis()->SetLabelSize(0.05);  
  ce_minus_te_ge_EtaCut->GetXaxis()->SetTitleSize(0.05);
  ce_minus_te_ge_EtaCut->GetYaxis()->SetTitle("ce_aggregate - te_aggregate (GeV)");
  ce_minus_te_ge_EtaCut->GetYaxis()->SetLabelSize(0.05);  
  ce_minus_te_ge_EtaCut->GetYaxis()->SetTitleSize(0.05);
  //ce_minus_te_ge_EtaCut->SetTitle("ce_agg - te_agg vs ge" + cut_text);

  counts_ge_EtaCut->GetXaxis()->SetTitle("ge 1d");
  counts_ge_EtaCut->GetXaxis()->SetLabelSize(0.05);  
  counts_ge_EtaCut->GetXaxis()->SetTitleSize(0.05);
  counts_ge_EtaCut->GetYaxis()->SetTitle("counts");
  counts_ge_EtaCut->GetYaxis()->SetLabelSize(0.05);  
  counts_ge_EtaCut->GetYaxis()->SetTitleSize(0.05);
  //counts_ge_EtaCut->SetTitle("ge_1d" + cut_text);

  te_by_ge_ge_EtaCut->GetXaxis()->SetTitle("Generated Energy (GeV)");
  te_by_ge_ge_EtaCut->GetXaxis()->SetLabelSize(0.05);  
  te_by_ge_ge_EtaCut->GetXaxis()->SetTitleSize(0.05);
  te_by_ge_ge_EtaCut->GetYaxis()->SetTitle("te_{agg}/ge");
  te_by_ge_ge_EtaCut->GetYaxis()->SetLabelSize(0.05);  
  te_by_ge_ge_EtaCut->GetYaxis()->SetTitleSize(0.05);
  //te_by_ge_ge_EtaCut->SetTitle("te_{agg}/ge vs ge" + cut_text);
    
  ce_by_ge_ge_EtaCut->GetXaxis()->SetTitle("Generated Energy (GeV)");
  ce_by_ge_ge_EtaCut->GetXaxis()->SetLabelSize(0.05);  
  ce_by_ge_ge_EtaCut->GetXaxis()->SetTitleSize(0.05);
  ce_by_ge_ge_EtaCut->GetYaxis()->SetTitle("ce_{agg}/ge"); 
  ce_by_ge_ge_EtaCut->GetYaxis()->SetLabelSize(0.05);  
  ce_by_ge_ge_EtaCut->GetYaxis()->SetTitleSize(0.05);
  //ce_by_ge_ge_EtaCut->SetTitle("ce_{agg}/ge vs ge" + cut_text);

  counts_te_EtaCut->GetXaxis()->SetTitle("te_{agg} 1d");
  counts_te_EtaCut->GetXaxis()->SetLabelSize(0.05);  
  counts_te_EtaCut->GetXaxis()->SetTitleSize(0.05);
  counts_te_EtaCut->GetYaxis()->SetTitle("counts");
  counts_te_EtaCut->GetYaxis()->SetLabelSize(0.05);  
  counts_te_EtaCut->GetYaxis()->SetTitleSize(0.05);
  //counts_te_EtaCut->SetTitle("te_{agg}_1d" + cut_text);

  counts_ce_EtaCut->GetXaxis()->SetTitle("ce_{agg} 1d");
  counts_ce_EtaCut->GetXaxis()->SetLabelSize(0.05);  
  counts_ce_EtaCut->GetXaxis()->SetTitleSize(0.05);
  counts_ce_EtaCut->GetYaxis()->SetTitle("counts");
  counts_ce_EtaCut->GetYaxis()->SetLabelSize(0.05);  
  counts_ce_EtaCut->GetYaxis()->SetTitleSize(0.05);
  //counts_ce_EtaCut->SetTitle("ce_{agg}_1d" + cut_text);

  ce_minus_ge_by_ge_ge_EtaCut->GetXaxis()->SetTitle("Generated Energy (GeV)");
  ce_minus_ge_by_ge_ge_EtaCut->GetXaxis()->SetLabelSize(0.05);  
  ce_minus_ge_by_ge_ge_EtaCut->GetXaxis()->SetTitleSize(0.05);
  ce_minus_ge_by_ge_ge_EtaCut->GetYaxis()->SetTitle("(ce_{agg}-ge)/ge");
  ce_minus_ge_by_ge_ge_EtaCut->GetYaxis()->SetLabelSize(0.05);  
  ce_minus_ge_by_ge_ge_EtaCut->GetYaxis()->SetTitleSize(0.05);
  //ce_minus_ge_by_ge_ge_EtaCut->SetTitle("(#Delta e_{agg}/ true_e) vs true_e" + cut_text);    

  ce_minus_ge_by_ge_ge_EtaCut_temp_2->GetXaxis()->SetTitle("Generated Energy (GeV)");
  ce_minus_ge_by_ge_ge_EtaCut_temp_2->GetXaxis()->SetLabelSize(0.05);  
  ce_minus_ge_by_ge_ge_EtaCut_temp_2->GetXaxis()->SetTitleSize(0.05);
  ce_minus_ge_by_ge_ge_EtaCut_temp_2->GetYaxis()->SetTitle("#sigma_{e_{agg}}");
  ce_minus_ge_by_ge_ge_EtaCut_temp_2->GetYaxis()->SetLabelSize(0.05);  
  ce_minus_ge_by_ge_ge_EtaCut_temp_2->GetYaxis()->SetTitleSize(0.05);
  //ce_minus_ge_by_ge_ge_EtaCut_2->SetTitle("#sigma_{e_{agg}} vs true_e" + cut_text);

  ce_minus_ge_by_ge_ge_EtaCut_temp_chi2->GetXaxis()->SetTitle("Generated Energy (GeV)");
  ce_minus_ge_by_ge_ge_EtaCut_temp_chi2->GetXaxis()->SetLabelSize(0.05);  
  ce_minus_ge_by_ge_ge_EtaCut_temp_chi2->GetXaxis()->SetTitleSize(0.05);
  ce_minus_ge_by_ge_ge_EtaCut_temp_chi2->GetYaxis()->SetTitle("reduced_#chi^{2}_{e_{agg}}");
  ce_minus_ge_by_ge_ge_EtaCut_temp_chi2->GetYaxis()->SetLabelSize(0.05);  
  ce_minus_ge_by_ge_ge_EtaCut_temp_chi2->GetYaxis()->SetTitleSize(0.04);
  //ce_minus_ge_by_ge_ge_EtaCut_chi2->SetTitle("#chi^{2}_{e_{agg}} vs true_e" + cut_text);

  ce_minus_ge_by_ge_ge_EtaCut_temp_1->GetXaxis()->SetTitle("Generated Energy (GeV)");
  ce_minus_ge_by_ge_ge_EtaCut_temp_1->GetXaxis()->SetLabelSize(0.05);  
  ce_minus_ge_by_ge_ge_EtaCut_temp_1->GetXaxis()->SetTitleSize(0.05);
  ce_minus_ge_by_ge_ge_EtaCut_temp_1->GetYaxis()->SetTitle("mean_{e_{agg}}");
  ce_minus_ge_by_ge_ge_EtaCut_temp_1->GetYaxis()->SetLabelSize(0.05);  
  ce_minus_ge_by_ge_ge_EtaCut_temp_1->GetYaxis()->SetTitleSize(0.05);
  //ce_minus_ge_by_ge_ge_EtaCut_1->SetTitle("mean_{e_{agg}} vs true_e" + cut_text);

  for(int sno = 0; sno < nSlices; sno++){
    slices[sno]->GetXaxis()->SetTitle("#Delta e^{agg}/ ge");
    slices[sno]->GetXaxis()->SetLabelSize(0.05);  
    slices[sno]->GetXaxis()->SetTitleSize(0.05);
    slices[sno]->GetYaxis()->SetTitle("Counts");
    slices[sno]->GetYaxis()->SetLabelSize(0.05);  
    slices[sno]->GetYaxis()->SetTitleSize(0.05);
  }

  if(debug==1){
    std::cout<<"\nWrite Histograms to File\n";
  }

  TFile *f = new TFile("energy_verification_EtaCut_" + detector + ".root","RECREATE"); //change the calorimeter name here
  f->GetList()->Add(ce_minus_te_ge_EtaCut);
  f->GetList()->Add(counts_ge_EtaCut);  
  f->GetList()->Add(ce_minus_ge_by_ge_ge_EtaCut);
  f->GetList()->Add(ce_minus_ge_by_ge_ge_EtaCut_temp_2);
  f->GetList()->Add(ce_minus_ge_by_ge_ge_EtaCut_temp_1);
  f->GetList()->Add(ce_minus_ge_by_ge_ge_EtaCut_temp_chi2);
  f->GetList()->Add(ce_by_ge_ge_EtaCut);
  f->GetList()->Add(te_by_ge_ge_EtaCut);
  f->GetList()->Add(counts_te_EtaCut);
  f->GetList()->Add(counts_ce_EtaCut);  
  
  for(int sno = 0; sno < nSlices; sno++){
    f->GetList()->Add(slices[sno]);
  }

  f->Write();

  if(print==1){
    if(debug==1){
      std::cout<<"\nSaving Histograms as .png\n";
    }
  
    gStyle -> SetOptStat(1);

    ce_minus_te_ge_EtaCut->Draw("colz");
    c->Print(detector + "ce_minus_te_ge_EtaCut.png");

    counts_ge_EtaCut->Draw("colz");
    c->Print(detector + "counts_ge_EtaCut.png");

    ce_minus_ge_by_ge_ge_EtaCut->Draw("colz");
    c->Print(detector + "ce_minus_ge_by_ge_ge_EtaCut.png");  

    gStyle -> SetOptStat(0);
    gStyle -> SetOptFit(0);
    fExp->SetLineColor(4); //38
    fTrue->SetLineColor(2); //46

    ce_minus_ge_by_ge_ge_EtaCut_temp_2->SetMarkerStyle(kFullCircle);
    ce_minus_ge_by_ge_ge_EtaCut_temp_2->SetMarkerColor(46); //30
    ce_minus_ge_by_ge_ge_EtaCut_temp_2->SetMarkerSize(0.75);
    ce_minus_ge_by_ge_ge_EtaCut_temp_2->SetAxisRange(0,1,"Y");
  
    ce_minus_ge_by_ge_ge_EtaCut_temp_2->Fit("fTrue", "M+");
    ce_minus_ge_by_ge_ge_EtaCut_temp_2->Draw("same");
    fExp->Draw("same");
  
    TLegend* legend = new TLegend(1.75,1.75);
    legend->SetHeader("Legend", "C");
    legend->AddEntry(ce_minus_ge_by_ge_ge_EtaCut_temp_2, "#sigma e_{agg} vs Generated Energy", "flep");
    legend->AddEntry((TObject*)0,"","");
    legend->AddEntry(fTrue, "p_{0} + p_{1}/#sqrt{ge} (Fitted)", "l");
    legend->AddEntry((TObject*)0,"","");
    legend->AddEntry(fExp, "0.1 + 0.5/#sqrt{ge} (Requirement)", "l");
    legend->SetTextSize(0.033);
    legend->Draw();

    cout<<"reduced_chi2 of fit: "<<fTrue->GetChisquare()<<"\n";
    cout<<"parameters: "<<fTrue->GetParameters()<<"\n";

    c->Print(detector + "sigmaE_ge_EtaCut.png");

    gStyle -> SetOptFit(102);

    ce_minus_ge_by_ge_ge_EtaCut_temp_1->Draw("colz");
    c->Print(detector + "meanE_ge_EtaCut.png");

    ce_minus_ge_by_ge_ge_EtaCut_temp_chi2->Draw("colz");
    c->Print(detector + "chi2E_ge_EtaCut.png");
  
    gStyle -> SetOptStat(1);

    ce_by_ge_ge_EtaCut->Draw("colz");
    c->Print(detector + "ce_by_ge_ge_EtaCut.png");

    te_by_ge_ge_EtaCut->Draw("colz");
    c->Print(detector + "te_by_ge_ge_EtaCut.png");

    counts_te_EtaCut->Draw("colz");
    c->Print(detector + "counts_te_EtaCut.png");
    counts_ce_EtaCut->Draw("colz");
    c->Print(detector + "counts_ce_EtaCut.png");

    gStyle -> SetOptStat(11);
  
    for(int sno = 0; sno < nSlices; sno++){
      TString plusOne = (TString)(sno + 1);
      TString nameF = detector + "_sigma_slice" + arr[sno] + ".png";
      slices[sno] -> Fit("fit", "M+");
      slices[sno] -> Draw("hist same");
      c->Print(nameF);
    }
  }

  std::cout << "The total ce is: " << total_ce << endl;
  std::cout << "The total te is: " << total_te << endl;
  std::cout << "The total ge is: " << total_ge << endl;
  c->Close();

  std::cout<<"\n\nDone\n----------------------------------------------------------------------\n\n";
}
