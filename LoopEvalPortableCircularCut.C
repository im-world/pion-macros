/*
  - Tower energy plots with manual clustering
  - Eta and Circular Cuts
  - Individual detectors
  - Standalone
  - Output file - energy_verification_EtaCut_CircularCut_detector.root
  - Plots made - [counts for twr_count], [te/ge vs ge], [tphi vs gphi], [theta vs gtheta], [dphi vs dtheta], 
  [(te-ge)/ge vs ge], [sigma_e vs ge], [mean_e vs ge], [chi2_e vs ge], slices for the [(te-ge)/ge vs ge] 
  plot, both with and without circular cuts.
*/

/*
  Authors:
 
  Sagar Joshi      (ee190002054@iiti.ac.in)
  Siddhant Rathi   (me190003061@iiti.ac.in)
	    
  version - 1.2
*/

#include <iostream>
#include <cmath>
#include <eicqa_modules/EvalRootTTree.h>
#include <eicqa_modules/EvalHit.h>
#include "TMath.h"
#include "TStyle.h"

R__LOAD_LIBRARY(libeicqa_modules.so)

void LoopEvalPortableCircularCut(TString detector, int print = 1, int debug = 0, Double_t energyCut = 0.0){
 
  TFile *f1 = new TFile("merged_Eval_" + detector + ".root","READ"); 

  TTree* T1 = (TTree*)f1->Get("T");
  EvalRootTTree *evaltree1 = nullptr;
   
  gStyle->SetCanvasPreferGL(kTRUE);
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

  int nSlicesx;                                                  // Number of ge-axis slices in sigma_e vs ge plot 
  int nSlicesy;                                                  // Number of sigma_e-axis slices in sigma_e vs ge plot
  Double_t eta_min, eta_max;                                     // Eta range of detectors
  Double_t x_radius, y_radius;                                   // Length of semi-minor and semi-major axes of circular (elliptical) cut
  Double_t fit_min, fit_max;                                     // Fit range of the first slice of [(te-ge)/ge vs ge] plot
  Double_t ttheta_gtheta_min, ttheta_gtheta_max;                 // Range of both axes in ttheta vs gtheta plot 
  Double_t towerCounts_max, towerCounts_CircularCut_max;         // Maximum value of x axis in 1d distribution showing number of towers
  Double_t sigma_min, sigma_max;                                 // Range of Y-axis in sigma_e vs ge plot
  Double_t mean_min, mean_max;                                   // Range of Y-axis in mean_e vs ge plot
  Double_t chi2_min, chi2_max;                                   // Range of Y-axis in chi2_e vs ge plot
  Double_t recalibration_firstSlice;                             // Recalibration factor of first slice (needed to be done manually because of low statistics)
  Double_t te_minus_ge_by_ge_ge_min, te_minus_ge_by_ge_ge_max;   // Range of y-axis in [(te-ge)/ge vs ge] plot
  TString fitting_function;                                      // Function fitted to slices of [(te-ge)/ge vs ge] plot
  TString cut_text, eRes, eRes1; 

  if (detector == "FEMC"){
    x_radius = 0.13; 
    y_radius = 0.40; 
    fit_min = -0.35;
    fit_max = 0.2;
    eta_min = 1.3;
    eta_max = 3.3;
    sigma_min = 0.;
    sigma_max = 0.5;
    mean_min = -0.1;
    mean_max = 0.1;
    chi2_min = 0;
    chi2_max = 1.82;
    ttheta_gtheta_min = 0;
    ttheta_gtheta_max = 0.6;
    towerCounts_max = 400;
    towerCounts_CircularCut_max = 200;
    te_minus_ge_by_ge_ge_min = -0.99;
    te_minus_ge_by_ge_ge_max = 1;
    recalibration_firstSlice = 0.7088;
    fitting_function = "gaus";
    cut_text = " {1.3 < true_eta < 3.3} ";
    eRes = "0.02 + 0.08/sqrt(x) + 0.02/x";
    eRes1 = "0.02 + 0.08/#sqrt{ge} + 0.02/ge";
    nSlicesx = 10;
    nSlicesy = 1000;

    if(debug==1){
      std::cout<<"FEMC eta cut & circular cut applied"<<"\n\n";
    }
  }

  else if(detector == "EEMC"){
    x_radius = 0.10; 
    y_radius = 0.40;   
    fit_min = 0.14;
    fit_max = 0.21; 
    eta_min = -3.5;
    eta_max = -1.7;
    sigma_min = 0;
    sigma_max = 0.3;
    mean_min = -0.1;
    mean_max = 0.1;
    chi2_min = 0;
    chi2_max = 3.1;
    ttheta_gtheta_min = 2.7;
    ttheta_gtheta_max = 3.1;
    towerCounts_max = 400;
    towerCounts_CircularCut_max = 200;
    te_minus_ge_by_ge_ge_min = -0.58; 
    te_minus_ge_by_ge_ge_max = 0.45; 
    recalibration_firstSlice = 0.75;
    fitting_function = "gaus";
    cut_text = " {-3.5 < true_eta < -1.7} ";
    eRes = "0.01 + 0.025/sqrt(x) + 0.01/x";
    eRes1 =  "0.01 + 0.025/#sqrt{ge} + 0.01/ge";
    nSlicesx = 10;
    nSlicesy = 350;

    if(debug==1){ 
      std::cout<<"EEMC eta cut & circular cut applied"<<"\n\n";
    }
  }

  else if(detector == "CEMC"){
    x_radius = 0.10; 
    y_radius = 0.20;
    fit_min = -0.4;
    fit_max = 0.4;
    eta_min = -1.5;
    eta_max = 1.2;
    sigma_min = 0;
    sigma_max = 0.5;
    mean_min = -0.1;
    mean_max = 0.1;
    chi2_min = 0;
    chi2_max = 1.42;
    ttheta_gtheta_min = 0.45;
    ttheta_gtheta_max = 2.6;
    towerCounts_max = 400;
    towerCounts_CircularCut_max = 200;
    te_minus_ge_by_ge_ge_min = -0.99;
    te_minus_ge_by_ge_ge_max = 1;
    recalibration_firstSlice = 1.414;
    fitting_function = "gaus";
    cut_text = " {-1.5 < true_eta < 1.2} ";
    eRes = "0.025 + 0.13/sqrt(x) + 0.02/x";
    eRes1 = "0.025 + 0.23/#sqrt{ge} + 0.02/ge";
    nSlicesx = 15;
    nSlicesy = 1000;

    if(debug==1){
      std::cout<<"CEMC eta cut & circular cut applied"<<"\n\n";
    }
  }

  else if(detector == "FHCAL"){
    x_radius = 0.2; 
    y_radius = 0.5;
    fit_min = -0.35;
    fit_max = -0.05;
    eta_min = 1.2;
    eta_max = 3.5;
    sigma_min = 0;
    sigma_max = 0.5;
    mean_min = 0.5;
    mean_max =0.95;
    chi2_min = 0;
    chi2_max = 1.82;
    ttheta_gtheta_min = 0;
    ttheta_gtheta_max = 0.6;
    towerCounts_max = 400;
    towerCounts_CircularCut_max = 200;
    te_minus_ge_by_ge_ge_min = -0.99;
    te_minus_ge_by_ge_ge_max = 1;
    recalibration_firstSlice = 1.414;
    fitting_function = "gaus";
    cut_text = " {1.2 < true_eta < 3.5} ";
    eRes = "1";
    eRes1 = "1";
    nSlicesx = 15;
    nSlicesy = 1000;

    if(debug==1){
      std::cout<<"FHCAL eta cut & circular cut applied"<<"\n\n";
    }
  }

  else if(detector == "HCALIN"){
    x_radius = 0.15; 
    y_radius = 0.4;
    fit_min = -0.35;
    fit_max = -0.05;
    eta_min = -1.1;
    eta_max = 1.1;
    sigma_min = 0;
    sigma_max = 0.5;
    mean_min = 0.5;
    mean_max = 0.95;
    chi2_min = 0;
    chi2_max = 1.82;
    ttheta_gtheta_min = 0;
    ttheta_gtheta_max = 0.6;
    towerCounts_max = 400;
    towerCounts_CircularCut_max = 200;
    te_minus_ge_by_ge_ge_min = -0.99;
    te_minus_ge_by_ge_ge_max = 1;
    recalibration_firstSlice = 1.414;
    fitting_function = "gaus";
    cut_text = " {-1.1 < true_eta < 1.1} ";
    eRes = "1";
    eRes1 = "1";
    nSlicesx = 15;
    nSlicesy = 1000;

    if(debug==1){
      std::cout<<"HCALIN eta cut & circular cut applied"<<"\n\n";
    }
  }

  else if(detector == "HCALOUT"){
    x_radius = 0.15; 
    y_radius = 0.4;
    fit_min = -0.35;
    fit_max = -0.05;
    eta_min = -1.1;
    eta_max = 1.1;
    sigma_min = 0;
    sigma_max = 0.5;
    mean_min = 0.5;
    mean_max = 0.95;
    chi2_min = 0;
    chi2_max = 1.82;
    ttheta_gtheta_min = 0;
    ttheta_gtheta_max = 0.6;
    towerCounts_max = 400;
    towerCounts_CircularCut_max = 200;
    te_minus_ge_by_ge_ge_min = -0.99;
    te_minus_ge_by_ge_ge_max = 1;
    recalibration_firstSlice = 1.414;
    fitting_function = "gaus";
    cut_text = " {-1.1 < true_eta < 1.1} ";
    eRes = "1";
    eRes1 = "1";
    nSlicesx = 15;
    nSlicesy = 1000;
      
    if(debug==1){
      std::cout<<"HCALOUT eta cut & circular cut applied"<<"\n\n";
    }
  }

  else{
    std::cout<<"Please try again.";
    return 1;
  }

  //Initialising histogram variables

  TH2D *te_minus_ge_by_ge_ge_EtaCut = new TH2D("te_minus_ge_by_ge_ge_EtaCut","#frac{#Delta e_{agg}}{truth e} vs truth e",200,0,30,200,-1.5,1);
  TH2D *te_minus_ge_by_ge_ge_EtaCut_CircularCut = new TH2D("te_minus_ge_by_ge_ge_EtaCut_CircularCut","#frac{#Delta e_{agg}}{truth e} vs truth e",200,0,30,200,-1.5,2);
  TH2D *te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated = new TH2D("te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated","#frac{#Delta e_{agg}}{truth e} vs truth e",200,0,30,200,-1.5,1);
  TH2D *te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp = new TH2D("te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp","#frac{#Delta e_{agg}}{truth e} vs truth e",nSlicesx,0,30,nSlicesy,te_minus_ge_by_ge_ge_min,te_minus_ge_by_ge_ge_max); // histogram from which mean vs ge, sigma vs ge, and reduced_chi2 vs ge plots are derived

  TH2D *te_by_ge_ge_EtaCut = new TH2D("te_by_ge_ge_EtaCut","te_{agg}/ge vs ge",200,0,30,200,-0.5,3.5);
  TH2D *te_by_ge_ge_EtaCut_CircularCut = new TH2D("te_by_ge_ge_EtaCut_CircularCut","te_{agg}/ge vs ge",200,0,30,200,-0.5,3.5);
  auto *mean_te_by_ge_ge_EtaCut_CircularCut = new TProfile("mean_te_by_ge_ge_EtaCut_CircularCut","Mean_{te/ge}",nSlicesx,0,30,-0.5,3.5);

  TH2D *tphi_gphi_EtaCut = new TH2D("tphi_gphi_EtaCut","tphi vs gphi",200,-4,4,200,-4,4);
  TH2D *tphi_gphi_EtaCut_CircularCut = new TH2D("tphi_gphi_EtaCut_CircularCut","tphi vs gphi",200,-4,4,200,-4,4);

  TH2D *ttheta_gtheta_EtaCut = new TH2D("ttheta_gtheta_EtaCut","ttheta vs gtheta",200,ttheta_gtheta_min,ttheta_gtheta_max,200,ttheta_gtheta_min,ttheta_gtheta_max);
  TH2D *ttheta_gtheta_EtaCut_CircularCut = new TH2D("ttheta_gtheta_EtaCut_CircularCut","ttheta vs gtheta",200,ttheta_gtheta_min,ttheta_gtheta_max,200,ttheta_gtheta_min,ttheta_gtheta_max);
 
  TH1D *counts_towerCounts_EtaCut = new TH1D("counts_towerCounts_EtaCut","n_towers",200,-1,towerCounts_max);
  TH1D *counts_towerCounts_EtaCut_CircularCut = new TH1D("counts_towerCounts_EtaCut_CircularCut","n_towers",200,-1,towerCounts_CircularCut_max);

  TH2D *dphi_dtheta_EtaCut = new TH2D("dphi_dtheta_EtaCut","dphi vs dtheta",200,-1,1,200,-1,1);

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
	  std::cout<<"\n"<<detector<<" Tower: "<<j<<"\n";
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

	  if(te > energyCut){
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
      //te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated->Fill(ge, ((te_aggregate_CircularCut/recalibration_factor)-ge)/ge);
      //if(debug==1){
      //cout<<"(ge, ((te_aggregate_CircularCut/recalibration_factor)-ge)/ge: ("<<ge<<", "<<((te_aggregate_CircularCut/recalibration_factor)-ge)/ge<<")\n";
      //}
      //te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp->Fill(ge, ((te_aggregate_CircularCut/recalibration_factor)-ge)/ge);
   
      te_by_ge_ge_EtaCut->Fill(ge, te_aggregate/ge);
      if(debug==1){
	cout<<"(ge, te_aggregate/ge): ("<<ge<<", "<<te_aggregate/ge<<")\n";
      }
     
      te_by_ge_ge_EtaCut_CircularCut->Fill(ge, te_aggregate_CircularCut/ge);
      mean_te_by_ge_ge_EtaCut_CircularCut->Fill(ge, te_aggregate_CircularCut/ge);
  
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

  
  TString arr[nSlicesx];
  for(int sno = 0; sno < nSlicesx; sno++){
    arr[sno] = TString::Itoa(sno + 1, 10); // Used for naming of slices
  }

  if(debug==1){
    std::cout<<"\nGenerating sigma_e vs ge plots\n\n";
  }


  Double_t recalibrationArr[nSlicesx];

  recalibrationArr[0] = recalibration_firstSlice;
  cout << "Recalibration factor for slice " << 1 << " is: " <<  recalibrationArr[0] << endl;

  for(int binIter = 2; binIter <= nSlicesx; binIter++){
    recalibrationArr[binIter-1] = mean_te_by_ge_ge_EtaCut_CircularCut->GetBinContent(binIter);
    cout << "Recalibration factor for slice " << binIter << " is: " <<  recalibrationArr[binIter-1] << endl;
  }

  for(int i=0; i<T1->GetEntries(); i++){

    T1->GetEntry(i);
          
    Double_t geta1 = evaltree1->get_geta();
    Double_t gphi = evaltree1->get_gphi();
    Double_t gtheta = evaltree1->get_gtheta();
    Double_t ge = evaltree1->get_ge();
    Double_t te_aggregate_CircularCut = 0;    

    if(geta1>=eta_min && geta1<=eta_max){
    	
      
      for (int j=0; j<evaltree1->get_ntowers(); j++){

	EvalTower *twr1 = evaltree1->get_tower(j);
	if (twr1){ 

	  Double_t tphi = twr1->get_tphi();
	  Double_t ttheta = twr1->get_ttheta();
	  Double_t dphi = tphi - gphi;
	  Double_t dtheta = ttheta - gtheta;	 
	  Double_t te = twr1->get_te();

          if(te > energyCut){
    
	    if(debug==1){
	      std::cout<<"pow(dphi/y_radius,2)+pow(dtheta/x_radius,2): "<<pow(dphi/y_radius,2)+pow(dtheta/x_radius,2)<<"\n";    
	    }

	    if (pow(dphi/y_radius,2)+pow(dtheta/x_radius,2)<=1){
	      te_aggregate_CircularCut += te;
	    }
	  }
   	}
      }
      int recalibration_factor = ceil((ge/30.0)*(Double_t)nSlicesx)- 1;
      te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated->Fill(ge, ((te_aggregate_CircularCut/recalibrationArr[recalibration_factor])-ge)/ge);
      te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp->Fill(ge, ((te_aggregate_CircularCut/recalibrationArr[recalibration_factor])-ge)/ge);
      if(debug==1){
        cout<<"(ge, ((te_aggregate_CircularCut/recalibration_factor)-ge)/ge: ("<<ge<<", "<<((te_aggregate_CircularCut/recalibrationArr[recalibration_factor])-ge)/ge<<")\n";
      }
    }    
  }


  // Initialising fit functions 
  TF1 *fit = new TF1("fit",fitting_function);   
  TF1 *fit1 = new TF1("fit1",fitting_function,fit_min,fit_max);
  TEllipse *el1 = new TEllipse(0,0,x_radius,y_radius);
  TF1 *fExp = new TF1("fExp",eRes,0,30);
  TF1 *fTrue = new TF1("fTrue","[0] + [1]/sqrt(x) + [2]/x",0,30); 

  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp->FitSlicesY(0, 2, -1, 0, "QN");
  TH2D *te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2 = (TH2D*)gDirectory->Get("te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2");
  TH2D *te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_1 = (TH2D*)gDirectory->Get("te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_1");
  TH2D *te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_chi2 = (TH2D*)gDirectory->Get("te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_chi2");


  // Generating individual slices

  TH1D* slices[nSlicesx];

  for(int sno = 0; sno < nSlicesx; sno++){
    int plusOne = sno+1;
    TString sname = "slice " + arr[sno];
    slices[sno] = te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp->ProjectionY(sname, plusOne, plusOne);
  }

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
 
  mean_te_by_ge_ge_EtaCut_CircularCut->GetXaxis()->SetTitle("Generated Energy (GeV)");
  mean_te_by_ge_ge_EtaCut_CircularCut->GetXaxis()->SetLabelSize(0.05);  
  mean_te_by_ge_ge_EtaCut_CircularCut->GetXaxis()->SetTitleSize(0.05);
  mean_te_by_ge_ge_EtaCut_CircularCut->GetYaxis()->SetTitle("Mean of te_{agg}/ge");
  mean_te_by_ge_ge_EtaCut_CircularCut->GetYaxis()->SetLabelSize(0.05);  
  mean_te_by_ge_ge_EtaCut_CircularCut->GetYaxis()->SetTitleSize(0.05);

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

  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2->GetXaxis()->SetTitle("Generated Energy (GeV)");
  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2->GetXaxis()->SetLabelSize(0.05);  
  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2->GetXaxis()->SetTitleSize(0.05);
  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2->GetYaxis()->SetTitle("#sigma_{e_{agg}}");
  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2->GetYaxis()->SetLabelSize(0.05);  
  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2->GetYaxis()->SetTitleSize(0.05);
  //te_minus_ge_by_ge_ge_EtaCut_2->SetTitle("#sigma_{e_{agg}} vs true_e" + cut_text);

  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_chi2->GetXaxis()->SetTitle("Generated Energy (GeV)");
  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_chi2->GetXaxis()->SetLabelSize(0.05);  
  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_chi2->GetXaxis()->SetTitleSize(0.05);
  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_chi2->GetYaxis()->SetTitle("Reduced_#chi^{2}_{e_{agg}}");
  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_chi2->GetYaxis()->SetLabelSize(0.05);  
  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_chi2->GetYaxis()->SetTitleSize(0.04);
  //te_minus_ge_by_ge_ge_EtaCut_chi2->SetTitle("#chi^{2}_{e_{agg}} vs true_e" + cut_text);

  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_1->GetXaxis()->SetTitle("Generated Energy (GeV)");
  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_1->GetXaxis()->SetLabelSize(0.05);  
  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_1->GetXaxis()->SetTitleSize(0.05);
  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_1->GetYaxis()->SetTitle("Mean_{e_{agg}}");
  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_1->GetYaxis()->SetLabelSize(0.05);  
  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_1->GetYaxis()->SetTitleSize(0.05);
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
    std::cout<<"\nWriting Histograms to File\n";
  }

  TFile *f = new TFile("Energy_verification_EtaCut_CircularCut" + detector + ".root","RECREATE"); 

  f->GetList()->Add(counts_towerCounts_EtaCut);
  f->GetList()->Add(counts_towerCounts_EtaCut_CircularCut);
  f->GetList()->Add(te_by_ge_ge_EtaCut);
  f->GetList()->Add(te_by_ge_ge_EtaCut_CircularCut);
  f->GetList()->Add(mean_te_by_ge_ge_EtaCut_CircularCut);
  f->GetList()->Add(tphi_gphi_EtaCut);
  f->GetList()->Add(tphi_gphi_EtaCut_CircularCut);
  f->GetList()->Add(ttheta_gtheta_EtaCut);
  f->GetList()->Add(ttheta_gtheta_EtaCut_CircularCut);
  f->GetList()->Add(dphi_dtheta_EtaCut);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut_CircularCut);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_1);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_chi2);
  
  for(int sno = 0; sno < nSlicesx; sno++){
    f->GetList()->Add(slices[sno]);
  }

  f->Write();

  gStyle -> SetOptStat(11);
  gStyle -> SetOptFit(112);
  
  
  int sno = 0;
  TString plusOne = (TString)(sno + 1);
  TString nameF = detector + "_sigmaE_slice" + arr[sno] + "_EtaCut_CircularCut.png";
  slices[sno] -> Fit("fit1", "R+");
  slices[sno] -> Draw("hist same");
  c->Print(nameF);

  /*
  ///////////////////////////////////////////////////////////////////////////////////////

  sno = 1;
  plusOne = (TString)(sno + 1);
  nameF = detector + "_ssj_sigmaE_slice" + arr[sno] + "_EtaCut_CircularCut.png";
  slices[sno] -> Fit("fit1", "M+");
  slices[sno] -> Draw("hist same");
  c->Print(nameF);

  sno = 2;
  plusOne = (TString)(sno + 1);
  nameF = detector + "_ssj_sigmaE_slice" + arr[sno] + "_EtaCut_CircularCut.png";
  slices[sno] -> Fit("fit1", "M+");
  slices[sno] -> Draw("hist same");
  c->Print(nameF);

  ///////////////////////////////////////////////////////////////////////////////////////
  */

  double_t mean = fit1->GetParameter(1);
  double_t mean_error = fit1->GetParError(1);
  double_t sigma = fit1->GetParameter(2);
  double_t sigma_error = fit1->GetParError(2); 
  double_t chi2 = (fit1->GetChisquare())/(fit1->GetNDF());
 
  TLine *Mean = new TLine(0,mean,30.0/nSlicesx,mean);
  TLine *Sigma = new TLine(0,sigma,30.0/nSlicesx,sigma);
  TLine *Chi2 = new TLine(0,chi2,30.0/nSlicesx,chi2);
  TLine *Periphery_Chi2 = new TLine(30.0/nSlicesx,0,30.0/nSlicesx,chi2);
  TLine *Sigma_error = new TLine(30.0/(2.0*nSlicesx),sigma-(sigma_error),30.0/(2.0*nSlicesx),sigma+(sigma_error));
  TLine *Mean_error = new TLine(30.0/(2.0*nSlicesx),mean-(mean_error),30.0/(2.0*nSlicesx),mean+(mean_error));
  TLine *Point = new TLine(30.0/(2.0*nSlicesx),sigma,30.0/(2.0*nSlicesx),sigma);

  if(print==1){
    if(debug==1){
      std::cout<<"\nSaving Histograms as .png\n";
    }

    gStyle -> SetOptStat(0);

    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_1->SetAxisRange(mean_min,mean_max,"Y");
    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_1->Draw();
    Mean->SetLineColor(1);
    Mean->SetLineWidth(1); 
    Mean_error->SetLineColor(1);
    Mean_error->SetLineWidth(1);
    Mean->Draw("same");
    Mean_error->Draw("same");
    c->Print(detector + "_meanE_ge_EtaCut_CircularCut.png");

    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_chi2->SetAxisRange(chi2_min,chi2_max,"Y");
    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_chi2->Draw();
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
    Point->SetLineColor(3);
    Point->SetLineWidth(5);
    fExp->SetLineColor(4); //38
    fTrue->SetLineColor(2); //46
  
    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2->SetMarkerStyle(kFullCircle);
    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2->SetMarkerColor(46); //30
    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2->SetMarkerSize(0.75);
    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2->SetAxisRange(sigma_min,sigma_max,"Y");

    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2->Fit("fTrue", "M+");
    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2->Draw("same");
    fExp->Draw("same");
    Sigma->Draw("same");
    Sigma_error->Draw("same");
    Point->Draw("same");
  
    TLegend* legend = new TLegend(1.75,1.75);
    legend->SetHeader("Legend", "C");
    legend->AddEntry(te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2, "#sigma_{e_{agg}} vs Generated Energy", "flep");
    legend->AddEntry((TObject*)0,"","");
    legend->AddEntry(fTrue, "p_{0} + p_{1}/#sqrt{ge} + p_{2}/ge (Fitted)", "l");
    legend->AddEntry((TObject*)0,"","");
    legend->AddEntry(fExp, eRes1, "l");
    legend->AddEntry((TObject*)0,"(Requirement)","" );
    legend->SetTextSize(0.033);
    legend->Draw();

    std::cout<<"reduced_chi2 of fit: "<<fTrue->GetChisquare()/fTrue->GetNDF()<<"\n";

    c->Print(detector + "_sigmaE_ge_EtaCut_CircularCut.png");

    gStyle -> SetOptStat(11);
    gStyle -> SetOptFit(112);
    // ofstream file1;
    // file1.open("gaus_1-4.txt"); 
    for(int sno = 1; sno < nSlicesx; sno++){  // sno = 1
      TString plusOne = (TString)(sno + 1);
      TString nameF = detector + "_sigmaE_slice" + arr[sno] + "_EtaCut_CircularCut.png";
      slices[sno] -> Fit("fit", "M+");
      slices[sno] -> Draw("hist same");
      c->Print(nameF);
      // double_t mean = fit->GetParameter(1);
      // double_t mean_error = fit->GetParError(1);
      // double_t sigma = fit->GetParameter(2);
      // double_t sigma_error = fit->GetParError(2); 
      // double_t chi2 = (fit->GetChisquare())/(fit->GetNDF());
      // file1 << mean << " " << mean_error << " " << sigma << " " << sigma_error << " " << chi2 << endl;
    }
    // file1.close();
    gStyle -> SetOptStat(1);

    counts_towerCounts_EtaCut->Draw("colz");
    c->Print(detector + "_counts_towerCounts_EtaCut.png");
    counts_towerCounts_EtaCut_CircularCut->Draw("colz");
    c->Print(detector + "_counts_towerCounts_EtaCut_CircularCut.png");

    te_by_ge_ge_EtaCut->Draw("colz");
    c->Print(detector + "_te_by_ge_ge_EtaCut.png");
    te_by_ge_ge_EtaCut_CircularCut->Draw("colz");
    c->Print(detector + "_te_by_ge_ge_EtaCut_CircularCut.png");

    gStyle -> SetOptStat(0);

    mean_te_by_ge_ge_EtaCut_CircularCut->Draw();
    c->Print(detector + "_mean_te_by_ge_ge_EtaCut_CircularCut.png");

    gStyle -> SetOptStat(1);

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
    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated->Draw("colz");
    c->Print(detector + "_te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated.png");

    c->SetCanvasSize(700,700);
    dphi_dtheta_EtaCut->Draw("colz");
    el1->SetFillColorAlpha(0, 0);
    el1->SetLineColorAlpha(2, 1);
    el1->SetLineWidth(3);
    el1->Draw("same");
    c->Print(detector + "_dphi_dtheta_EtaCut_CircularCut.png");

    c->SetCanvasSize(500,700);
    c->Close();
  }
  std::cout << "The total_te is: " << total_te << endl;
  std::cout << "The total_te_CircularCut is: " << total_te_CircularCut << endl;
  std::cout << "The total_ge is: " << total_ge << endl; 
  std::cout<<"\n\nDone\n----------------------------------------------------------------------\n\n";
 
}
