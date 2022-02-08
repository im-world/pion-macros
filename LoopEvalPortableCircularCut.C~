// Energy plots with manual clustering; Eta and Circular Cuts; Individual detectors; 
//Standalone; Plots made - [counts for twr_count], [te/ge vs ge], [tphi vs gphi],
// [theta vs gtheta], [dphi vs dtheta], [(te-ge)/ge vs ge], [sigma_e vs ge], 
//[mean_e vs ge], [chi2_e vs ge] both with and without circular cuts.

#include <iostream>
#include <eicqa_modules/EvalRootTTree.h>
#include <eicqa_modules/EvalHit.h>
#include "TMath.h"
#include "TStyle.h"

R__LOAD_LIBRARY(libeicqa_modules.so)

void LoopEvalPortableCircularCut(TString detector, int print = 1, int debug = 0){
 
  TFile *f1 = new TFile("merged_Eval_" + detector + ".root","READ"); 

  TTree* T1 = (TTree*)f1->Get("T");
  EvalRootTTree *evaltree1 = nullptr;
   
  TCanvas *c = new TCanvas();
  c->SetTickx();
  c->SetTicky();

  // Modifying the default plotting style  
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(102);
  gStyle->SetTitleXOffset(1);
  gStyle->SetTitleYOffset(1);
  gStyle->SetLabelSize(0.05);  
  gStyle->SetTitleXSize(0.05);  
  gStyle->SetTitleYSize(0.05);
 
  long double total_te = 0; 
  long double total_te_CircularCut = 0;
  long double total_ge = 0; 

  int nSlices = 15; // Number of ge-axis slices in sigma_e vs ge plot 

  //Initialising histogram variables

  TH2D *te_minus_ge_by_ge_ge_EtaCut = new TH2D("te_minus_ge_by_ge_ge_EtaCut","#frac{#Delta e_{agg}}{truth e} vs truth e",200,0,30,200,-2,1);
  TH2D *te_minus_ge_by_ge_ge_EtaCut_CircularCut = new TH2D("te_minus_ge_by_ge_ge_EtaCut_CircularCut","#frac{#Delta e_{agg}}{truth e} vs truth e",200,0,30,200,-2,1);
  TH2D *te_minus_ge_by_ge_ge_EtaCut_temp = new TH2D("te_minus_ge_by_ge_ge_EtaCut_temp","#frac{#Delta e_{agg}}{truth e} vs truth e",nSlices,0,30,2000,-0.99,1);

  TH2D *te_by_ge_ge_EtaCut = new TH2D("te_by_ge_ge_EtaCut","te_{agg}/ge vs ge",200,0,30,200,-0.5,1.5);
  TH2D *te_by_ge_ge_EtaCut_CircularCut = new TH2D("te_by_ge_ge_EtaCut_CircularCut","te_{agg}/ge vs ge",200,0,30,200,-0.5,1.5);

  TH2D *tphi_gphi_EtaCut = new TH2D("tphi_gphi_EtaCut","tphi vs gphi",200,-4,4,200,-4,4);
  TH2D *tphi_gphi_EtaCut_CircularCut = new TH2D("tphi_gphi_EtaCut_CircularCut","tphi vs gphi",200,-4,4,200,-4,4);

  TH2D *ttheta_gtheta_EtaCut = new TH2D("ttheta_gtheta_EtaCut","ttheta vs gtheta",200,0,0.6,200,0,0.6);
  TH2D *ttheta_gtheta_EtaCut_CircularCut = new TH2D("ttheta_gtheta_EtaCut_CircularCut","ttheta vs gtheta",200,0,0.6,200,0,0.6);

  TH2D *dphi_dtheta_EtaCut = new TH2D("dphi_dtheta_EtaCut","dphi vs dtheta",200,-1,1,200,-1,1);

  TH1D *counts_towerCounts_EtaCut = new TH1D("counts_towerCounts_EtaCut","n_towers",200,-1,400);
  TH1D *counts_towerCounts_EtaCut_CircularCut = new TH1D("counts_towerCounts_EtaCut_CircularCut","n_towers",200,-1,200);
  
  Double_t eta_min, eta_max, x_radius, y_radius;
  TString cut_text;

  if (detector == "FEMC"){
    x_radius = 0.15; //radius of circular cut
    y_radius = 0.4;
    eta_min = 1.3;
    eta_max = 3.3;
    cut_text = " {1.3 < true_eta < 3.3} ";
      
    if(debug==1){
      std::cout<<"FEMC eta cut & circular cut applied"<<"\n\n";
    }
  }

  else if(detector == "EEMC"){
    x_radius = 0.15; //radius of circular cut
    y_radius = 0.4;    
    eta_min = -3.5;
    eta_max = -1.7;
    cut_text = " {-3.5 < true_eta < -1.7} ";

    if(debug==1){
      std::cout<<"EEMC eta cut & circular cut applied"<<"\n\n";
    }
  }

  else if(detector == "CEMC"){
    x_radius = 0.15; //radius of circular cut
    y_radius = 0.4;
    eta_min = -1.5;
    eta_max = 1.2;
    cut_text = " {-1.5 < true_eta < 1.2} ";

    if(debug==1){
      std::cout<<"CEMC eta cut & circular cut applied"<<"\n\n";
    }
  }

  else if(detector == "FHCAL"){
    x_radius = 0.2; //radius of circular cut
    y_radius = 0.5;
    eta_min = 1.2;
    eta_max = 3.5;
    cut_text = " {1.2 < true_eta < 3.5} ";
     
    if(debug==1){
      std::cout<<"FHCAL eta cut & circular cut applied"<<"\n\n";
    }
  }

  else if(detector == "HCALIN"){
    x_radius = 0.15; //radius of circular cut
    y_radius = 0.4;
    eta_min = -1.1;
    eta_max = 1.1;
    cut_text = " {-1.1 < true_eta < 1.1} ";
    
    if(debug==1){
      std::cout<<"HCALIN eta cut & circular cut applied"<<"\n\n";
    }
  }

  else if(detector == "HCALOUT"){
    x_radius = 0.15; //radius of circular cut
    y_radius = 0.4;
    eta_min = -1.1;
    eta_max = 1.1;
    cut_text = " {-1.1 < true_eta < 1.1} ";
      
    if(debug==1){
      std::cout<<"HCALOUT eta cut & circular cut applied"<<"\n\n";
    }
  }

  else{
    std::cout<<"Please try again.";
    return 1;
  }

  T1->SetBranchAddress("DST#EvalTTree_" + detector,&evaltree1);
  
  for(int i=0; i<T1->GetEntries(); i++){

    T1->GetEntry(i);
      
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
      if(debug==1){
	std::cout<<"ge: "<<ge<<"\n";
      }	

      Double_t gphi = evaltree1->get_gphi();
      if(debug==1){
	std::cout<<"gphi: "<<gphi<<"\n";
      }
	
      Double_t gtheta = evaltree1->get_gtheta();
      if(debug==1){
	std::cout<<"gtheta: "<<gtheta<<"\n";
 }

      total_ge += ge;
      if(debug==1){
	std::cout<<"total_ge till now = "<<total_ge<<"\n";
      }
	
      int twr_count = 0;
      int twr_count_CircularCut = 0;
      Double_t te_aggregate = 0;
      Double_t te_aggregate_CircularCut = 0;    

      for (int j=0; j<evaltree1->get_ntowers(); j++){
	
	if(debug==1){
	  std::cout<<"\nFEMC Tower: "<<j<<"\n";
	}

	EvalTower *twr1 = evaltree1->get_tower(j);
	if (twr1){ 
	  twr_count += 1;
	  if(debug==1){
	    cout<<"non-empty "<<detector<<" tower\n";
	  }
	  Double_t tphi = twr1->get_tphi();
	  if(debug==1){
	    std::cout<<"tphi: "<<tphi<<"\n";
	  }
	  Double_t ttheta = twr1->get_ttheta();
	  if(debug==1){
	    std::cout<<"ttheta: "<<ttheta<<"\n";
	  }
	  Double_t dphi = tphi - gphi;
	  if(debug==1){
	    std::cout<<"tphi-gphi: "<<dphi<<"\n";
	  }
	  Double_t dtheta = ttheta - gtheta;
	  if(debug==1){
	    std::cout<<"ttheta-gtheta: "<<dtheta<<"\n";
	  }
	  Double_t te = twr1->get_te();
	  if(debug==1){
	    std::cout<<"te: "<<te<<"\n";
	  }
	  te_aggregate += te;
     
	  dphi_dtheta_EtaCut->Fill(dtheta, dphi);
	  tphi_gphi_EtaCut->Fill(gphi, tphi, te);
	  ttheta_gtheta_EtaCut->Fill(gtheta, ttheta, te); 

	  if(debug==1){
	    std::cout<<"pow(dphi/y_radius,2)+pow(dtheta/x_radius,2): "<<pow(dphi/y_radius,2)+pow(dtheta/x_radius,2)<<"\n";    
	  }

	  if (pow(dphi/y_radius,2)+pow(dtheta/x_radius,2)<=1){
	    if(debug==1){
	      cout<<"Tower included after circular cut\n";
	    }  
	    twr_count_CircularCut += 1;
	    te_aggregate_CircularCut += te;
	    tphi_gphi_EtaCut_CircularCut->Fill(gphi, tphi, te);
	    ttheta_gtheta_EtaCut_CircularCut->Fill(gtheta, ttheta, te); 
	  }

	  if(debug==1){
	    cout<<"te_aggregate till now = "<<te_aggregate<<"\n";
	    cout<<"te_aggregate_CircularCut till now = "<<te_aggregate_CircularCut<<"\n";
	    std::cout<<"te += "<<twr1->get_te()<<"\n";
	  }
	}
      }
       
      total_te += te_aggregate;
      total_te_CircularCut += te_aggregate_CircularCut;
      if(debug==1){
	cout<<"total_te till now = "<<total_te<<"\n";
	cout<<"total_te_CircularCut till now = "<<total_te_CircularCut<<"\n\n";
      }

      te_minus_ge_by_ge_ge_EtaCut->Fill(ge, (te_aggregate-ge)/ge);
      if(debug==1){     
	cout<<"(ge, (te_aggregate-ge)/ge): ("<<ge<<", "<<(te_aggregate-ge)/ge<<")\n";
      }
      te_minus_ge_by_ge_ge_EtaCut_CircularCut->Fill(ge, (te_aggregate_CircularCut-ge)/ge);
      if(debug==1){
	cout<<"(ge, (te_aggregate_CircularCut-ge)/ge): ("<<ge<<", "<<(te_aggregate_CircularCut-ge)/ge<<")\n";
      }
      te_minus_ge_by_ge_ge_EtaCut_temp->Fill(ge, (te_aggregate_CircularCut-ge)/ge);
   
      te_by_ge_ge_EtaCut->Fill(ge, te_aggregate/ge);
      if(debug==1){
	cout<<"(ge, te_aggregate/ge): ("<<ge<<", "<<te_aggregate/ge<<")\n";
      }
      te_by_ge_ge_EtaCut_CircularCut->Fill(ge, te_aggregate_CircularCut/ge);
      if(debug==1){
	cout<<"(ge, te_aggregate_CircularCut/ge): ("<<ge<<", "<<te_aggregate_CircularCut/ge<<")\n";   
      }

      counts_towerCounts_EtaCut->Fill(twr_count);
      if(debug==1){
	cout<<"(twr_count): ("<<twr_count<<")\n";
      }
     
      counts_towerCounts_EtaCut_CircularCut->Fill(twr_count_CircularCut);
      if(debug==1){
	cout<<"(twr_count_CircularCut): ("<<twr_count_CircularCut<<")\n";
      }      
    }
  }

  if(debug==1){  
    cout<<"\neta cut if ends"<<"\n";
  }

  
  TString arr[nSlices];
  for(int sno = 0; sno < nSlices; sno++){
    arr[sno] = TString::Itoa(sno + 1, 10); // Used for naming of slices
  }

  if(debug==1){
    std::cout<<"\nGenerating sigma_e vs ge plots\n\n";
  }

  te_minus_ge_by_ge_ge_EtaCut_temp->FitSlicesY(0, 0, -1, 0, "QN");
  TH2D *te_minus_ge_by_ge_ge_EtaCut_temp_2 = (TH2D*)gDirectory->Get("te_minus_ge_by_ge_ge_EtaCut_temp_2");
  TH2D *te_minus_ge_by_ge_ge_EtaCut_temp_1 = (TH2D*)gDirectory->Get("te_minus_ge_by_ge_ge_EtaCut_temp_1");
  TH2D *te_minus_ge_by_ge_ge_EtaCut_temp_chi2 = (TH2D*)gDirectory->Get("te_minus_ge_by_ge_ge_EtaCut_temp_chi2");
 

  // Generating individual slices

  TH1D* slices[nSlices];

  for(int sno = 0; sno < nSlices; sno++){
    int plusOne = sno+1;
    TString sname = "slice " + arr[sno];
    slices[sno] = te_minus_ge_by_ge_ge_EtaCut_temp->ProjectionY(sname, plusOne, plusOne);
  }

  // Initialising fit functions 
  TF1 *fit = new TF1("fit", "gaus");   
  TEllipse *el1 = new TEllipse(0,0,x_radius,y_radius);
  TF1 *fExp = new TF1("fExp","0.02 + 0.08/sqrt(x) + 0.02/x",0,30);
  TF1 *fTrue = new TF1("fTrue","[0] + [1]/sqrt(x) + [2]/x",0,30);  

  if(debug==1){
    std::cout<<"\nHistogram Formatting\n\n";
  }

  counts_towerCounts_EtaCut->GetXaxis()->SetTitle("Number of Towers");
  counts_towerCounts_EtaCut->GetXaxis()->SetLabelSize(0.045);  
  counts_towerCounts_EtaCut->GetXaxis()->SetTitleSize(0.045);
  counts_towerCounts_EtaCut->GetYaxis()->SetTitle("Counts");
  counts_towerCounts_EtaCut->GetYaxis()->SetLabelSize(0.045);  
  counts_towerCounts_EtaCut->GetYaxis()->SetTitleSize(0.045);

  counts_towerCounts_EtaCut_CircularCut->GetXaxis()->SetTitle("Number of Towers");
  counts_towerCounts_EtaCut_CircularCut->GetXaxis()->SetLabelSize(0.045);  
  counts_towerCounts_EtaCut_CircularCut->GetXaxis()->SetTitleSize(0.045);
  counts_towerCounts_EtaCut_CircularCut->GetYaxis()->SetTitle("Counts");
  counts_towerCounts_EtaCut_CircularCut->GetYaxis()->SetLabelSize(0.045);  
  counts_towerCounts_EtaCut_CircularCut->GetYaxis()->SetTitleSize(0.045);       

  te_minus_ge_by_ge_ge_EtaCut->GetXaxis()->SetTitle("Generated Energy (GeV)");
  te_minus_ge_by_ge_ge_EtaCut->GetXaxis()->SetLabelSize(0.05);  
  te_minus_ge_by_ge_ge_EtaCut->GetXaxis()->SetTitleSize(0.05);
  te_minus_ge_by_ge_ge_EtaCut->GetYaxis()->SetTitle("(te_{agg}-ge)/ge");
  te_minus_ge_by_ge_ge_EtaCut->GetYaxis()->SetLabelSize(0.05);  
  te_minus_ge_by_ge_ge_EtaCut->GetYaxis()->SetTitleSize(0.05);

  te_minus_ge_by_ge_ge_EtaCut_CircularCut->GetXaxis()->SetTitle("Generated Energy (GeV)");
  te_minus_ge_by_ge_ge_EtaCut_CircularCut->GetXaxis()->SetLabelSize(0.05);  
  te_minus_ge_by_ge_ge_EtaCut_CircularCut->GetXaxis()->SetTitleSize(0.05);
  te_minus_ge_by_ge_ge_EtaCut_CircularCut->GetYaxis()->SetTitle("(te_{agg}-ge)/ge");
  te_minus_ge_by_ge_ge_EtaCut_CircularCut->GetYaxis()->SetLabelSize(0.05);  
  te_minus_ge_by_ge_ge_EtaCut_CircularCut->GetYaxis()->SetTitleSize(0.05);

  te_by_ge_ge_EtaCut->GetXaxis()->SetTitle("Generated Energy (GeV)");
  te_by_ge_ge_EtaCut->GetXaxis()->SetLabelSize(0.05);  
  te_by_ge_ge_EtaCut->GetXaxis()->SetTitleSize(0.05);
  te_by_ge_ge_EtaCut->GetYaxis()->SetTitle("te_{agg}/ge");
  te_by_ge_ge_EtaCut->GetYaxis()->SetLabelSize(0.05);  
  te_by_ge_ge_EtaCut->GetYaxis()->SetTitleSize(0.05);

  te_by_ge_ge_EtaCut_CircularCut->GetXaxis()->SetTitle("Generated Energy (GeV)");
  te_by_ge_ge_EtaCut_CircularCut->GetXaxis()->SetLabelSize(0.05);  
  te_by_ge_ge_EtaCut_CircularCut->GetXaxis()->SetTitleSize(0.05);
  te_by_ge_ge_EtaCut_CircularCut->GetYaxis()->SetTitle("te_{agg}/ge");
  te_by_ge_ge_EtaCut_CircularCut->GetYaxis()->SetLabelSize(0.05);  
  te_by_ge_ge_EtaCut_CircularCut->GetYaxis()->SetTitleSize(0.05);

  tphi_gphi_EtaCut->GetXaxis()->SetTitle("Generated #Phi");
  tphi_gphi_EtaCut->GetXaxis()->SetLabelSize(0.05);  
  tphi_gphi_EtaCut->GetXaxis()->SetTitleSize(0.05);
  tphi_gphi_EtaCut->GetYaxis()->SetTitle("Tower #Phi");
  tphi_gphi_EtaCut->GetYaxis()->SetLabelSize(0.05);  
  tphi_gphi_EtaCut->GetYaxis()->SetTitleSize(0.05);

  tphi_gphi_EtaCut_CircularCut->GetXaxis()->SetTitle("Generated #Phi");
  tphi_gphi_EtaCut_CircularCut->GetXaxis()->SetLabelSize(0.05);  
  tphi_gphi_EtaCut_CircularCut->GetXaxis()->SetTitleSize(0.05);
  tphi_gphi_EtaCut_CircularCut->GetYaxis()->SetTitle("Tower #Phi");
  tphi_gphi_EtaCut_CircularCut->GetYaxis()->SetLabelSize(0.05);  
  tphi_gphi_EtaCut_CircularCut->GetYaxis()->SetTitleSize(0.05);

  ttheta_gtheta_EtaCut->GetXaxis()->SetTitle("Generated #theta");
  ttheta_gtheta_EtaCut->GetXaxis()->SetLabelSize(0.05);  
  ttheta_gtheta_EtaCut->GetXaxis()->SetTitleSize(0.05);
  ttheta_gtheta_EtaCut->GetYaxis()->SetTitle("Tower #theta");
  ttheta_gtheta_EtaCut->GetYaxis()->SetLabelSize(0.05);  
  ttheta_gtheta_EtaCut->GetYaxis()->SetTitleSize(0.05);
   
  ttheta_gtheta_EtaCut_CircularCut->GetXaxis()->SetTitle("Generated #theta");
  ttheta_gtheta_EtaCut_CircularCut->GetXaxis()->SetLabelSize(0.05);  
  ttheta_gtheta_EtaCut_CircularCut->GetXaxis()->SetTitleSize(0.05);
  ttheta_gtheta_EtaCut_CircularCut->GetYaxis()->SetTitle("Tower #theta");
  ttheta_gtheta_EtaCut_CircularCut->GetYaxis()->SetLabelSize(0.05);  
  ttheta_gtheta_EtaCut_CircularCut->GetYaxis()->SetTitleSize(0.05);

  dphi_dtheta_EtaCut->GetXaxis()->SetTitle("ttheta-gtheta");
  dphi_dtheta_EtaCut->GetXaxis()->SetLabelSize(0.05);  
  dphi_dtheta_EtaCut->GetXaxis()->SetTitleSize(0.05);
  dphi_dtheta_EtaCut->GetYaxis()->SetTitle("tphi-gphi");
  dphi_dtheta_EtaCut->GetYaxis()->SetLabelSize(0.05);  
  dphi_dtheta_EtaCut->GetYaxis()->SetTitleSize(0.05);

  te_minus_ge_by_ge_ge_EtaCut_temp_2->GetXaxis()->SetTitle("Generated Energy (GeV)");
  te_minus_ge_by_ge_ge_EtaCut_temp_2->GetXaxis()->SetLabelSize(0.05);  
  te_minus_ge_by_ge_ge_EtaCut_temp_2->GetXaxis()->SetTitleSize(0.05);
  te_minus_ge_by_ge_ge_EtaCut_temp_2->GetYaxis()->SetTitle("#sigma_{e_{agg}}");
  te_minus_ge_by_ge_ge_EtaCut_temp_2->GetYaxis()->SetLabelSize(0.05);  
  te_minus_ge_by_ge_ge_EtaCut_temp_2->GetYaxis()->SetTitleSize(0.05);
  //te_minus_ge_by_ge_ge_EtaCut_2->SetTitle("#sigma_{e_{agg}} vs true_e" + cut_text);

  te_minus_ge_by_ge_ge_EtaCut_temp_chi2->GetXaxis()->SetTitle("Generated Energy (GeV)");
  te_minus_ge_by_ge_ge_EtaCut_temp_chi2->GetXaxis()->SetLabelSize(0.05);  
  te_minus_ge_by_ge_ge_EtaCut_temp_chi2->GetXaxis()->SetTitleSize(0.05);
  te_minus_ge_by_ge_ge_EtaCut_temp_chi2->GetYaxis()->SetTitle("Reduced_#chi^{2}_{e_{agg}}");
  te_minus_ge_by_ge_ge_EtaCut_temp_chi2->GetYaxis()->SetLabelSize(0.05);  
  te_minus_ge_by_ge_ge_EtaCut_temp_chi2->GetYaxis()->SetTitleSize(0.04);
  //te_minus_ge_by_ge_ge_EtaCut_chi2->SetTitle("#chi^{2}_{e_{agg}} vs true_e" + cut_text);

  te_minus_ge_by_ge_ge_EtaCut_temp_1->GetXaxis()->SetTitle("Generated Energy (GeV)");
  te_minus_ge_by_ge_ge_EtaCut_temp_1->GetXaxis()->SetLabelSize(0.05);  
  te_minus_ge_by_ge_ge_EtaCut_temp_1->GetXaxis()->SetTitleSize(0.05);
  te_minus_ge_by_ge_ge_EtaCut_temp_1->GetYaxis()->SetTitle("Mean_{e_{agg}}");
  te_minus_ge_by_ge_ge_EtaCut_temp_1->GetYaxis()->SetLabelSize(0.05);  
  te_minus_ge_by_ge_ge_EtaCut_temp_1->GetYaxis()->SetTitleSize(0.05);
  //te_minus_ge_by_ge_ge_EtaCut_1->SetTitle("mean_{e_{agg}} vs true_e" + cut_text);
	
  for(int sno = 0; sno < nSlices; sno++){
    slices[sno]->GetXaxis()->SetTitle("#Delta e^{agg}/ ge");
    slices[sno]->GetXaxis()->SetLabelSize(0.05);  
    slices[sno]->GetXaxis()->SetTitleSize(0.05);
    slices[sno]->GetYaxis()->SetTitle("Counts");
    slices[sno]->GetYaxis()->SetLabelSize(0.05);  
    slices[sno]->GetYaxis()->SetTitleSize(0.05);
  }

  if(debug==1){
    std::cout<<"\nWriting Histograms to File\n";
  }

  TFile *f = new TFile("energy_verification_EtaCut_CircularCut" + detector + ".root","RECREATE"); 

  f->GetList()->Add(counts_towerCounts_EtaCut);
  f->GetList()->Add(counts_towerCounts_EtaCut_CircularCut);
  f->GetList()->Add(te_by_ge_ge_EtaCut);
  f->GetList()->Add(te_by_ge_ge_EtaCut_CircularCut);
  f->GetList()->Add(tphi_gphi_EtaCut);
  f->GetList()->Add(tphi_gphi_EtaCut_CircularCut);
  f->GetList()->Add(ttheta_gtheta_EtaCut);
  f->GetList()->Add(ttheta_gtheta_EtaCut_CircularCut);
  f->GetList()->Add(dphi_dtheta_EtaCut);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut_CircularCut);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut_temp_2);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut_temp_1);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut_temp_chi2);
  
  for(int sno = 0; sno < nSlices; sno++){
    f->GetList()->Add(slices[sno]);
  }

  f->Write();

  if(print==1){
    if(debug==1){
      std::cout<<"\nSaving Histograms as .png\n";
    }

    gStyle -> SetOptStat(0);

    te_minus_ge_by_ge_ge_EtaCut_temp_1->Draw("colz");
    c->Print(detector + "meanE_ge_EtaCut_CircularCut.png");
    te_minus_ge_by_ge_ge_EtaCut_temp_chi2->Draw("colz");
    c->Print(detector + "chi2E_ge_EtaCut_CircularCut.png");
 
    gStyle -> SetOptStat(0);
    gStyle -> SetOptFit(0);
    fExp->SetLineColor(4); //38
    fTrue->SetLineColor(2); //46
  
    te_minus_ge_by_ge_ge_EtaCut_temp_2->SetMarkerStyle(kFullCircle);
    te_minus_ge_by_ge_ge_EtaCut_temp_2->SetMarkerColor(46); //30
    te_minus_ge_by_ge_ge_EtaCut_temp_2->SetMarkerSize(0.75);
    te_minus_ge_by_ge_ge_EtaCut_temp_2->SetAxisRange(0,0.5,"Y");

    te_minus_ge_by_ge_ge_EtaCut_temp_2->Fit("fTrue", "M+");
    te_minus_ge_by_ge_ge_EtaCut_temp_2->Draw("same");
    fExp->Draw("same");
  
    TLegend* legend = new TLegend(1.75,1.75);
    legend->SetHeader("Legend", "C");
    legend->AddEntry(te_minus_ge_by_ge_ge_EtaCut_temp_2, "#sigma_{e_{agg}} vs Generated Energy", "flep");
    legend->AddEntry((TObject*)0,"","");
    legend->AddEntry(fTrue, "p_{0} + p_{1}/#sqrt{ge} + p_{2}/ge (Fitted)", "l");
    legend->AddEntry((TObject*)0,"","");
    legend->AddEntry(fExp, "0.02 + 0.08/#sqrt{ge} +0.02/ge", "l");
    legend->AddEntry((TObject*)0,"(Requirement)","" );
    legend->SetTextSize(0.033);
    legend->Draw();

    //cout<<"reduced_chi2 of fit: "<<fTrue->GetChisquare()<<"\n";
    //cout<<"parameters: "<<fTrue->GetParameters()<<"\n";

    c->Print(detector + "sigmeE_ge_EtaCut_CircularCut.png");

    gStyle -> SetOptStat(11);
    gStyle -> SetOptFit(102);
    for(int sno = 0; sno < nSlices; sno++){
      TString plusOne = (TString)(sno + 1);
      TString nameF = detector + "sigmaE_slice" + arr[sno] + "_EtaCut_CircularCut.png";
      slices[sno] -> Fit("fit", "M+");
      slices[sno] -> Draw("hist same");
      c->Print(nameF);
    }

    gStyle -> SetOptStat(1);

    counts_towerCounts_EtaCut->Draw("colz");
    c->Print(detector + "_counts_towerCounts_EtaCut.png");
    counts_towerCounts_EtaCut_CircularCut->Draw("colz");
    c->Print(detector + "_counts_towerCounts_EtaCut_CircularCut.png");

    te_by_ge_ge_EtaCut->Draw("colz");
    c->Print(detector + "_te_by_ge_ge_EtaCut.png");
    te_by_ge_ge_EtaCut_CircularCut->Draw("colz");
    c->Print(detector + "_te_by_ge_ge_EtaCut_CircularCut.png");

    tphi_gphi_EtaCut->Draw("colz");
    c->Print(detector + "_tphi_gphi_EtaCut.png");
    tphi_gphi_EtaCut_CircularCut->Draw("colz");
    c->Print(detector + "_tphi_gphi_EtaCut_CircularCut.png");
 
    ttheta_gtheta_EtaCut->Draw("colz");
    c->Print(detector + "_ttheta_gtheta_EtaCut.png");
    ttheta_gtheta_EtaCut_CircularCut->Draw("colz");
    c->Print(detector + "_ttheta_gtheta_EtaCut_CircularCut.png");
  
    te_minus_ge_by_ge_ge_EtaCut->Draw("colz");
    c->Print(detector + "_te_minus_ge_by_ge_ge_EtaCut.png");
    te_minus_ge_by_ge_ge_EtaCut_CircularCut->Draw("colz");
    c->Print(detector + "_te_minus_ge_by_ge_ge_EtaCut_CircularCut.png");

    c->SetCanvasSize(700,700);
    dphi_dtheta_EtaCut->Draw("colz");
    el1->SetFillColorAlpha(1, 0);
    el1->SetLineColorAlpha(2, 1);
    el1->SetLineWidth(3);
    el1->Draw("same");
    c->Print(detector + "_dphi_dtheta_EtaCut_CircularCut.png");

    c->SetCanvasSize(500,700);

    std::cout << "The total_te is: " << total_te << endl;
    std::cout << "The total_te_CircularCut is: " << total_te_CircularCut << endl;
    std::cout << "The total ge is: " << total_ge << endl; 
  }
  c->Close();
  std::cout<<"\n\nDone\n----------------------------------------------------------------------\n\n";
 
}
