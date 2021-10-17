/*
> LoopEvalMultiFCircularCut.C 
- Tower energy plots with manual clustering
- Eta and Circular Cuts
- Combined CEMC+HCALIN+HCALOUT
- Standalone
- Output file - energy_verification_EtaCut_CircularCut_HCALIN_HCALOUT_CEMC.root
- Plots made - [sigma_e vs ge], [mean_e vs ge], [chi2_e vs ge], slices for the
  [(te-ge)/ge vs ge], [te/ge vs ge], [(te-ge)/ge vs ge], both with and without circular cuts.
*/

/*
  authors - Sagar Joshi      (ee190002054@iiti.ac.in)
            Siddhant Rathi   (me190003061@iiti.ac.in)
	    
  version - 1.1
*/

#include <iostream>
#include <eicqa_modules/EvalRootTTree.h>
#include <eicqa_modules/EvalHit.h>
#include "TMath.h"
#include "TStyle.h"

R__LOAD_LIBRARY(libeicqa_modules.so)

void LoopEvalMultiHCircularCut(int print = 1, int debug = 0)
{
  
  TString detector = "HCALIN_HCALOUT_CEMC";
  TFile *f1 = new TFile("merged_Eval_HCALIN.root","READ"); 
  TFile *f2 = new TFile("merged_Eval_HCALOUT.root","READ"); 
  TFile *f3 = new TFile("merged_Eval_CEMC.root","READ"); 

  TTree* T1 = (TTree*)f1->Get("T");
  EvalRootTTree *evaltree1 = nullptr;
  
  TTree* T2 = (TTree*)f2->Get("T");
  EvalRootTTree *evaltree2 = nullptr;

  TTree* T3 = (TTree*)f3->Get("T");
  EvalRootTTree *evaltree3 = nullptr;
    
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

  long double total_te = 0; 
  long double total_te_CircularCut = 0;
  long double total_ge = 0; 

  int nSlicesx = 10; // Number of ge-axis slices taken for making sigma_e vs ge plot
  Double_t eta_min, eta_max; // Eta range of detectors
  Double_t x_radius_HCALIN = 0.15; // Length of semi-minor axis of circular (elliptical) in HCALIN
  Double_t y_radius_HCALIN = 0.25; // Length of semi-major axis of circular (elliptical) in HCALIN
  Double_t x_radius_HCALOUT = 0.2; // Length of semi-minor axis of circular (elliptical) in HCALOUT
  Double_t y_radius_HCALOUT = 0.3; // Length of semi-major axis of circular (elliptical) in HCALOUT
  Double_t x_radius_CEMC = 0.1; // Length of semi-minor axis of circular (elliptical) in CEMC
  Double_t y_radius_CEMC = 0.2; // Length of semi-minor axis of circular (elliptical) in CEMC
  Double_t fit_min, fit_max; // Fit range of the first slice of [(te-ge)/ge vs ge] plot
  Double_t sigma_min, sigma_max; // Range of Y-axis in sigma_e vs ge plot
  Double_t mean_min, mean_max; // Range of Y-axis in mean_e vs ge plot
  Double_t chi2_min, chi2_max; // Range of Y-axis in chi2_e vs ge plot
  Double_t recalibration_factor; // Number divided from entries of [(te-ge)/ge vs ge] plot for recalibration   
  
  fit_min = -0.8;
  fit_max = 1.0;
  eta_min = -1.1;
  eta_max = 1.1;
  sigma_min = 0;
  sigma_max = 1.2;
  mean_min = -0.6;
  mean_max = 0.2;
  chi2_min = 0;
  chi2_max = 1.60;
  recalibration_factor = 0.7461;

  TString cut_text = " {-1.1 < geta < 1.1} ";
  

  TH2D *te_minus_ge_by_ge_ge_EtaCut = new TH2D("te_minus_ge_by_ge_ge_EtaCut","#frac{#Delta e_{agg}}{truth e} vs truth e",200,0,30,200,-2,1);
  TH2D *te_minus_ge_by_ge_ge_EtaCut_CircularCut = new TH2D("te_minus_ge_by_ge_ge_EtaCut_CircularCut","#frac{#Delta e_{agg}}{truth e} vs truth e",200,0,30,200,-1.5,3);
  TH2D *te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated = new TH2D("te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated","#frac{#Delta e_{agg}}{truth e} vs truth e",200,0,30,200,-1.5,1.5);
  TH2D *te_minus_ge_by_ge_ge_EtaCut_temp = new TH2D("te_minus_ge_by_ge_ge_EtaCut_temp","#frac{#Delta e_{agg}}{truth e} vs truth e",nSlicesx,0,30,500,-0.99,1.3); // histogram from which mean vs ge, sigma vs ge, and reduced_chi2 vs ge plots are derived; Recalibration applied

  TH2D *te_by_ge_ge_EtaCut = new TH2D("te_by_ge_ge_EtaCut","te_{agg}/ge vs ge",200,0,30,200,-0.5,1.5);
  TH2D *te_by_ge_ge_EtaCut_CircularCut = new TH2D("te_by_ge_ge_EtaCut_CircularCut","te_{agg}/ge vs ge",200,0,30,200,-0.5,3.5);

  T1->SetBranchAddress("DST#EvalTTree_HCALIN",&evaltree1);
  T2->SetBranchAddress("DST#EvalTTree_HCALOUT",&evaltree2);
  T3->SetBranchAddress("DST#EvalTTree_CEMC",&evaltree3);

  for(int i=0; i<T1->GetEntries(); i++) // We assume same no. of entries, since no cuts are applied
    {
      T1->GetEntry(i);
      T2->GetEntry(i);
      T3->GetEntry(i);

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
	  cout<<"\ngeta cut applied (-1.1, 1.1)"<<"\n\n";
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

	Double_t te_aggregate = 0;
	Double_t te_aggregate_CircularCut = 0;    
        for (int j=0; j<evaltree1->get_ntowers(); j++){
	
	  if(debug==1){
	    std::cout<<"\nHCALIN Tower: "<<j<<"\n";
	  }

	  EvalTower *twr1 = evaltree1->get_tower(j);
	  if (twr1){ 
	  
	    if(debug==1){
	      cout<<"non-empty HCALIN tower\n";
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
     
	    if(debug==1){
	      std::cout<<"HCALIN: pow(dphi/y_radius,2)+pow(dtheta/x_radius,2) = "<<pow(dphi/y_radius_HCALIN,2)+pow(dtheta/x_radius_HCALIN,2)<<"\n";    
	    }

	    if (pow(dphi/y_radius_HCALIN,2)+pow(dtheta/x_radius_HCALIN,2)<=1){
	      if(debug==1){
		cout<<"HCALIN Tower included after circular cut\n";
	      }	      
	      te_aggregate_CircularCut += te;
	    }
	    if(debug==1){
	      cout<<"te_aggregate till now = "<<te_aggregate<<"\n";
	      std::cout<<"te += "<<twr1->get_te()<<"\n";
	    }
	  }
	}

	for (int j=0; j<evaltree2->get_ntowers(); j++){
	

	  if(debug==1){
	    std::cout<<"\nHCALOUT Tower: "<<j<<"\n";
	  }

	  EvalTower *twr2 = evaltree2->get_tower(j);
	  if (twr2){ 

	    if(debug==1){
	      cout<<"non-empty HCALOUT tower\n";
	    }

	    Double_t tphi = twr2->get_tphi();
	    if(debug==1){
	      std::cout<<"tphi: "<<tphi<<"\n";
	    }

	    Double_t ttheta = twr2->get_ttheta();
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

	    Double_t te = twr2->get_te();
	    if(debug==1){
	      std::cout<<"te: "<<te<<"\n";
	    }

	    te_aggregate += te;
     
	    if(debug==1){
	      std::cout<<"HCALOUT: pow(dphi/y_radius,2)+pow(dtheta/x_radius,2): "<<pow(dphi/y_radius_HCALOUT,2)+pow(dtheta/x_radius_HCALOUT,2)<<"\n";    
	    }

	    if (pow(dphi/y_radius_HCALOUT,2)+pow(dtheta/x_radius_HCALOUT,2)<=1){
	      if(debug==1){
		cout<<"HCALOUT Tower included after circular cut\n";
	      }
	      te_aggregate_CircularCut += te;
	    }

	    if(debug==1){
	      cout<<"te_aggregate till now = "<<te_aggregate<<"\n";
	      std::cout<<"te += "<<twr2->get_te()<<"\n";
	    }
	  }
	}

        for (int j=0; j<evaltree3->get_ntowers(); j++){
	

	  if(debug==1){
	    std::cout<<"\nCEMC Tower: "<<j<<"\n";
	  }

	  EvalTower *twr3 = evaltree3->get_tower(j);
	  if (twr3){ 

	    if(debug==1){
	      cout<<"non-empty CEMC tower\n";
	    }

	    Double_t tphi = twr3->get_tphi();
	    if(debug==1){
	      std::cout<<"tphi: "<<tphi<<"\n";
	    }

	    Double_t ttheta = twr3->get_ttheta();
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

	    Double_t te = twr3->get_te();
	    if(debug==1){
	      std::cout<<"te: "<<te<<"\n";
	    }

	    te_aggregate += te;
     
	    if(debug==1){
	      std::cout<<"CEMC: pow(dphi/y_radius,2)+pow(dtheta/x_radius,2): "<<pow(dphi/y_radius_CEMC,2)+pow(dtheta/x_radius_CEMC,2)<<"\n";    
	    }

	    if (pow(dphi/y_radius_CEMC,2)+pow(dtheta/x_radius_CEMC,2)<=1){
	      if(debug==1){
		cout<<"CEMC Tower included after circular cut\n";
	      }
	      te_aggregate_CircularCut += te;
	    }

	    if(debug==1){
	      cout<<"te_aggregate till now = "<<te_aggregate<<"\n";
	      std::cout<<"te += "<<twr3->get_te()<<"\n";
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
	te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated->Fill(ge, ((te_aggregate_CircularCut/recalibration_factor)-ge)/ge);

	if(debug==1){
	  cout<<"(ge, ((te_aggregate_CircularCut/recalibration_factor)-ge)/ge: ("<<ge<<", "<<((te_aggregate_CircularCut/recalibration_factor)-ge)/ge<<")\n";
	}
	te_minus_ge_by_ge_ge_EtaCut_temp->Fill(ge, ((te_aggregate_CircularCut/recalibration_factor)-ge)/ge);
   
	te_by_ge_ge_EtaCut->Fill(ge, te_aggregate/ge);

	if(debug==1){	
	  cout<<"(ge, te_aggregate/ge): ("<<ge<<", "<<te_aggregate/ge<<")\n";
	}

	te_by_ge_ge_EtaCut_CircularCut->Fill(ge, te_aggregate_CircularCut/ge);
	if(debug==1){
	  cout<<"(ge, te_aggregate_CircularCut/ge): ("<<ge<<", "<<te_aggregate_CircularCut/ge<<")\n";   
	}      
      }  
    }

  if(debug==1){
    cout<<"\neta cut if ends"<<"\n";
  }



  TString arr[nSlicesx]; // Used for naming fitted slices used in sigma_e vs ge 
  for(int sno = 0; sno < nSlicesx; sno++){
    arr[sno] = TString::Itoa(sno + 1, 10);
  }

  // Generating sigma_e vs ge plot
  te_minus_ge_by_ge_ge_EtaCut_temp->FitSlicesY(0, 2, -1, 0, "QN");
  TH2D *te_minus_ge_by_ge_ge_EtaCut_temp_2 = (TH2D*)gDirectory->Get("te_minus_ge_by_ge_ge_EtaCut_temp_2");
  TH2D *te_minus_ge_by_ge_ge_EtaCut_temp_1 = (TH2D*)gDirectory->Get("te_minus_ge_by_ge_ge_EtaCut_temp_1");
  TH2D *te_minus_ge_by_ge_ge_EtaCut_temp_chi2 = (TH2D*)gDirectory->Get("te_minus_ge_by_ge_ge_EtaCut_temp_chi2");

  TH1D* slices[nSlicesx];

  // Generating plots for individual slices
  for(int sno = 0; sno < nSlicesx; sno++){
    int plusOne = sno+1;
    TString sname = "slice " + arr[sno];
    slices[sno] = te_minus_ge_by_ge_ge_EtaCut_temp->ProjectionY(sname, plusOne, plusOne);
  }

  // Initialising fit functions 
  TF1 *fit = new TF1("fit", "gaus");   
  TF1 *fit1 = new TF1("fit1","gaus",fit_min,fit_max);
  TF1 *fExp = new TF1("fExp","0.1 + 1.0/sqrt(x)",0,30);
  TF1 *fTrue = new TF1("fTrue","[0] + [1]/sqrt(x)",0,30);


  if(debug==1){
    cout<<"\nHistogram Formatting\n\n";
  }

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

  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated->GetXaxis()->SetTitle("Generated Energy (GeV)");
  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated->GetXaxis()->SetLabelSize(0.05);  
  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated->GetXaxis()->SetTitleSize(0.05);
  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated->GetYaxis()->SetTitle("(te_{agg}-ge)/ge");
  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated->GetYaxis()->SetLabelSize(0.05);  
  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated->GetYaxis()->SetTitleSize(0.05);

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
	
  for(int sno = 0; sno < nSlicesx; sno++){
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

  TFile *f = new TFile("energy_verification_EtaCut_CircularCut" + detector + ".root","RECREATE"); 
    
  f->GetList()->Add(te_by_ge_ge_EtaCut);
  f->GetList()->Add(te_by_ge_ge_EtaCut_CircularCut);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut_CircularCut);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut_temp);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut_temp_2);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut_temp_1);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut_temp_chi2);
  
  for(int sno = 0; sno < nSlicesx; sno++){
    f->GetList()->Add(slices[sno]);
  }

    gStyle -> SetOptStat(11);
  gStyle -> SetOptFit(112);
  
  int sno = 0;
  TString plusOne = (TString)(sno + 1);
  TString nameF = detector + "_sigmaE_slice" + arr[sno] + "_EtaCut_CircularCut.png";
  slices[sno] -> Fit("fit1", "R+");
  slices[sno] -> Draw("hist same");
  c->Print(nameF);

  double_t mean = fit1->GetParameter(1);
  double_t mean_error = fit1->GetParError(1);
  double_t sigma = fit1->GetParameter(2);
  double_t sigma_error = fit1->GetParError(2); 
  double_t chi2 = (fit1->GetChisquare())/(fit1->GetNDF());
 
  TLine *Mean = new TLine(0,mean,30.0/nSlicesx,mean);
  TLine *Sigma = new TLine(0,sigma,30.0/nSlicesx,sigma);
  TLine *Chi2 = new TLine(0,chi2,30.0/nSlicesx,chi2);
  TLine *Periphery_Chi2 = new TLine(30.0/nSlicesx,0,30.0/nSlicesx,chi2);
  TLine *Sigma_error = new TLine(30.0/(2.0*nSlicesx),sigma-(sigma_error/2),30.0/(2.0*nSlicesx),sigma+(sigma_error/2));
  TLine *Mean_error = new TLine(30.0/(2.0*nSlicesx),mean-(mean_error/2),30.0/(2.0*nSlicesx),mean+(mean_error/2));
  TLine *Point = new TLine(30.0/(2.0*nSlicesx),sigma,30.0/(2.0*nSlicesx),sigma);

  if(print==1){
    if(debug==1){
      std::cout<<"\nSaving Histograms as .png\n";
    }

    gStyle -> SetOptStat(0);

    te_minus_ge_by_ge_ge_EtaCut_temp_1->SetAxisRange(mean_min,mean_max,"Y");
    te_minus_ge_by_ge_ge_EtaCut_temp_1->Draw("colz");
    Mean->SetLineColor(1);
    Mean->SetLineWidth(1); 
    Mean_error->SetLineColor(1);
    Mean_error->SetLineWidth(1);
    Mean->Draw("same");
    Mean_error->Draw("same");
    c->Print(detector + "_meanE_ge_EtaCut_CircularCut.png");

    te_minus_ge_by_ge_ge_EtaCut_temp_chi2->SetAxisRange(chi2_min,chi2_max,"Y");
    te_minus_ge_by_ge_ge_EtaCut_temp_chi2->Draw("colz");
    Chi2->SetLineColor(1);
    Chi2->SetLineWidth(1); 
    Periphery_Chi2->SetLineColor(1);
    Periphery_Chi2->SetLineWidth(1);
    Chi2->Draw("same");
    Periphery_Chi2->Draw("same");
    c->Print(detector + "_chi2E_ge_EtaCut_CircularCut.png");
 
    gStyle -> SetOptStat(0);
    gStyle -> SetOptFit(0);

    Sigma->SetLineColor(1);
    Sigma->SetLineWidth(1); 
    Sigma_error->SetLineColor(1);
    Sigma_error->SetLineWidth(1);
    Point->SetLineColor(2);
    Point->SetLineWidth(5);
    fExp->SetLineColor(4); //38
    fTrue->SetLineColor(2); //46
  
    te_minus_ge_by_ge_ge_EtaCut_temp_2->SetMarkerStyle(kFullCircle);
    te_minus_ge_by_ge_ge_EtaCut_temp_2->SetMarkerColor(46); //30
    te_minus_ge_by_ge_ge_EtaCut_temp_2->SetMarkerSize(0.75);
    te_minus_ge_by_ge_ge_EtaCut_temp_2->SetAxisRange(sigma_min,sigma_max,"Y");

    te_minus_ge_by_ge_ge_EtaCut_temp_2->Fit("fTrue", "M+");
    te_minus_ge_by_ge_ge_EtaCut_temp_2->Draw("same");
    fExp->Draw("same");
    Sigma->Draw("same");
    Sigma_error->Draw("same");
    Point->Draw("same");

    TLegend* legend = new TLegend(1.75,1.75);
    legend->SetHeader("Legend", "C");
    legend->AddEntry(te_minus_ge_by_ge_ge_EtaCut_temp_2, "#sigma_{e_{agg}} vs Generated Energy", "flep");
    legend->AddEntry((TObject*)0,"","");
    legend->AddEntry(fTrue, "p_{0} + p_{1}/#sqrt{ge} (Fitted)", "l");
    legend->AddEntry((TObject*)0,"","");
    legend->AddEntry(fExp, "0.1 + 1.0/#sqrt{ge} (Requirement)", "l");
    legend->SetTextSize(0.033);
    legend->Draw();

    cout<<"reduced_chi2 of fit: "<<fTrue->GetChisquare()/fTrue->GetNDF()<<"\n";

    c->Print(detector + "_sigmaE_ge_EtaCut_CircularCut.png");

    gStyle -> SetOptStat(11);
    gStyle -> SetOptFit(112);
    for(int sno = 1; sno < nSlicesx; sno++){
      TString plusOne = (TString)(sno + 1);
      TString nameF = detector + "_sigmaE_slice" + arr[sno] + "_EtaCut_CircularCut.png";
      slices[sno] -> Fit("fit", "M+");
      slices[sno] -> Draw("hist same");
      c->Print(nameF);
    }

    gStyle -> SetOptStat(1);

    te_by_ge_ge_EtaCut->Draw("colz");
    c->Print(detector + "_te_by_ge_ge_EtaCut.png");
    te_by_ge_ge_EtaCut_CircularCut->Draw("colz");
    c->Print(detector + "_te_by_ge_ge_EtaCut_CircularCut.png");

    te_minus_ge_by_ge_ge_EtaCut->Draw("colz");
    c->Print(detector + "_te_minus_ge_by_ge_ge_EtaCut.png");
    te_minus_ge_by_ge_ge_EtaCut_CircularCut->Draw("colz");
    c->Print(detector + "_te_minus_ge_by_ge_ge_EtaCut_CircularCut.png");
    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated->Draw("colz");
    c->Print(detector + "_te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated.png");
  }

  std::cout << "The total te is: " << total_te << endl;
  std::cout << "The total te_CircularCut is: " << total_te_CircularCut << endl;
  std::cout << "The total ge is: " << total_ge << endl; 
  c->Close();
  std::cout<<"\n\nDone\n----------------------------------------------------------------------\n\n";
 
}
