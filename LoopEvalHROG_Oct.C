/*
> LoopEvalMultiFCircularCut.C 
- Tower energy plots with manual clustering
- Eta and Circular Cuts
- Combined FEMC+FHCAL
- Standalone
- Output file - energy_verification_EtaCut_CircularCut_FHCAL_FEMC.root
- Plots made - [sigma_e vs ge], [mean_e vs ge], [chi2_e vs ge], slices for the
  [(te-ge)/ge vs ge], [te/ge vs ge], [(te-ge)/ge vs ge], both with and without circular cuts.
*/

/*
  authors - Sagar Joshi      (ee190002054@iiti.ac.in)
            Siddhant Rathi   (me190003061@iiti.ac.in)
	    
  version - 1.1

*/

#include <iostream>
#include <stdexcept>
#include <eicqa_modules/EvalRootTTree.h>
#include <eicqa_modules/EvalHit.h>
#include "TMath.h"
#include "TStyle.h"

R__LOAD_LIBRARY(libeicqa_modules.so)

void LoopEvalHROG_Oct(int print = 1, int debug = 0, Double_t energyCutAggregate = 0.0, Double_t energyCut = 0.0, int MIP_theta_parametrisation = 1){

  Double_t EMC_cut = 0.2;
  TF1 *mip_pmzn_energy_cut_ftheta = new TF1("mip_pmzn_energy_cut_ftheta", "(9.46093e-01) - 1.62771*x + 1.37776*(x^2) - (5.4996e-01)*(x^3) + (8.82673e-02)*(x^4)"); 
  

  if(MIP_theta_parametrisation == 1 && energyCut != 0){
    throw std::invalid_argument("We do not currently support theta-parametrized \nMIP cut on EMC simultaneously with individual tower cuts \non other detectors.:(;");
  } 


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

  gStyle->SetCanvasPreferGL(kTRUE); 
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
  int nSlicesy = 350;
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
  Double_t recalibration_firstSlice = 3.8147;
  Double_t recalibration_firstSlice_FHCAL = 0.0260;                             // Recalibration factor of first slice (needed to be done manually because of low statistics)
  Double_t recalibration_firstSlice_FEMC = 0.0300;
  Double_t te_minus_ge_by_ge_ge_min, te_minus_ge_by_ge_ge_max;   // Range of y-axis in [(te-ge)/ge vs ge] plot

  fit_min = -0.8;
  fit_max = 1.0;
  eta_min = -0.96;
  eta_max = 0.92;
  sigma_min = 0;
  sigma_max = 1.5;
  mean_min = -0.3;
  mean_max = 0.3;
  chi2_min = 0;
  chi2_max = 2.23;
  //  recalibration_factor = 0.7088;
  te_minus_ge_by_ge_ge_min = -0.99;
  te_minus_ge_by_ge_ge_max = 1.0;
   

  TString cut_text = " {-0.96 < geta < 0.92} ";
  

  TH2D *te_minus_ge_by_ge_ge_EtaCut = new TH2D("te_minus_ge_by_ge_ge_EtaCut","#frac{#Delta e_{agg}}{truth e} vs truth e",200,0,30,200,-2,1);
  TH2D *te_minus_ge_by_ge_ge_EtaCut_CircularCut = new TH2D("te_minus_ge_by_ge_ge_EtaCut_CircularCut","#frac{#Delta e_{agg}}{truth e} vs truth e",200,0,30,200,-1.5,2);;
  TH2D *te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated = new TH2D("te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated","#frac{#Delta e_{agg}}{truth e} vs truth e",200,0,30,200,-1.5,1.5);
  TH2D *te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp = new TH2D("te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp","#frac{#Delta e_{agg}}{truth e} vs truth e",nSlicesx,0,30,nSlicesy,te_minus_ge_by_ge_ge_min,te_minus_ge_by_ge_ge_max); // histogram from which mean vs ge, sigma vs ge, and reduced_chi2 vs ge plots are derived

  TH2D *te_minus_ge_by_ge_ge_EtaCut_temp = new TH2D("te_minus_ge_by_ge_ge_EtaCut_temp","#frac{#Delta e_{agg}}{truth e} vs truth e",nSlicesx,0,30,500,-0.99,1.3); // histogram from which mean vs ge, sigma vs ge, and reduced_chi2 vs ge plots are derived

  TH2D *te_by_ge_ge_EtaCut = new TH2D("te_by_ge_ge_EtaCut","te_{agg}/ge vs ge",200,0,30,200,-0.5,1.5);
  TH2D *te_by_ge_ge_EtaCut_CircularCut = new TH2D("te_by_ge_ge_EtaCut_CircularCut","te_{agg}/ge vs ge",200,0,30,200,-1,2);
  TH2D *te_by_ge_ge_EtaCut_CircularCut_HCALIN = new TH2D("te_by_ge_ge_EtaCut_CircularCut_HCALIN","te_{agg}/ge vs ge",200,0,30,200,-1,2);
  TH2D *te_by_ge_ge_EtaCut_CircularCut_HCALOUT = new TH2D("te_by_ge_ge_EtaCut_CircularCut_HCALOUT","te_{agg}/ge vs ge",200,0,30,200,-1,2);
  TH2D *te_by_ge_ge_EtaCut_CircularCut_CEMC = new TH2D("te_by_ge_ge_EtaCut_CircularCut_CEMC","te_{agg}/ge vs ge",200,0,30,200,-1,2);

  auto *mean_te_by_ge_ge_EtaCut_CircularCut = new TProfile("mean_te_by_ge_ge_EtaCut_CircularCut","Mean_{te/ge}",nSlicesx,0,30,-0.5,35);
  auto *mean_te_by_ge_ge_EtaCut_CircularCut_HCALIN = new TProfile("mean_te_by_ge_ge_EtaCut_CircularCut_HCALIN","Mean_{te/ge}",nSlicesx,0,30,-0.5,35);
  auto *mean_te_by_ge_ge_EtaCut_CircularCut_HCALOUT = new TProfile("mean_te_by_ge_ge_EtaCut_CircularCut_HCALOUT","Mean_{te/ge}",nSlicesx,0,30,-0.5,35);
  auto *mean_te_by_ge_ge_EtaCut_CircularCut_CEMC = new TProfile("mean_te_by_ge_ge_EtaCut_CircularCut_CEMC","Mean_{te/ge}",nSlicesx,0,30,-0.5,35);

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

	Double_t te_aggregate = 0;
	Double_t te_aggregate_CircularCut = 0;    
	Double_t te_aggregate_HCALIN_CircularCut = 0;
	Double_t te_aggregate_HCALOUT_CircularCut = 0;

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

	    if(te > energyCut){
	      te_aggregate += te;
     
	      if(debug==1){
		std::cout<<"HCALIN: pow(dphi/y_radius_HCALIN,2)+pow(dtheta/x_radius_HCALOUT,2) = "<<pow(dphi/y_radius_HCALIN,2)+pow(dtheta/x_radius_HCALIN,2)<<"\n";    
	      }

	      if (pow(dphi/y_radius_HCALIN,2)+pow(dtheta/x_radius_HCALIN,2)<=1){
		if(debug==1){
		  cout<<"HCALIN Tower included after circular cut\n";
		}	      

		te_aggregate_CircularCut += te;
		te_aggregate_HCALIN_CircularCut += te;

	      }

	      if(debug==1){
		cout<<"te_aggregate till now = "<<te_aggregate<<"\n";
		std::cout<<"te += "<<twr1->get_te()<<"\n";
		cout<<"te_aggregate_CircularCut till now = "<<te_aggregate_CircularCut<<"\n";
	      }
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

	    if(te > energyCut){
	      te_aggregate += te;
     
	      if(debug==1){
		std::cout<<"HCALOUT: pow(dphi/y_radius_HCALOUT,2)+pow(dtheta/x_radius_HCALOUT,2): "<<pow(dphi/y_radius_HCALOUT,2)+pow(dtheta/x_radius_HCALOUT,2)<<"\n";    
	      }

	      if (pow(dphi/y_radius_HCALOUT,2)+pow(dtheta/x_radius_HCALOUT,2)<=1){
		if(debug==1){
		  cout<<"HCALOUT Tower included after circular cut\n";
		}
		
		te_aggregate_HCALOUT_CircularCut += te;
		te_aggregate_CircularCut += te;
	      }

	      if(debug==1){
		cout<<"te_aggregate till now = "<<te_aggregate<<"\n";
		cout<<"te_aggregate_CircularCut till now = "<<te_aggregate_CircularCut<<"\n";
		std::cout<<"te += "<<twr2->get_te()<<"\n";
	      }
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
	    
	    if(MIP_theta_parametrisation == 1){
	      EMC_cut = mip_pmzn_energy_cut_ftheta->Eval(gtheta);
	    }

	    if(te > energyCut + EMC_cut){
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
		cout<<"te_aggregate_CircularCut till now = "<<te_aggregate_CircularCut<<"\n";
		std::cout<<"te += "<<twr3->get_te()<<"\n";
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

	if(te_aggregate_CircularCut > energyCutAggregate){

	  te_minus_ge_by_ge_ge_EtaCut->Fill(ge, (te_aggregate-ge)/ge);

	  if(debug==1){
	    cout<<"(ge, (te_aggregate-ge)/ge): ("<<ge<<", "<<(te_aggregate-ge)/ge<<")\n";
	  }

	  te_minus_ge_by_ge_ge_EtaCut_CircularCut->Fill(ge, (te_aggregate_CircularCut-ge)/ge);
	  if(debug==1){
	    cout<<"(ge, (te_aggregate_CircularCut-ge)/ge): ("<<ge<<", "<<(te_aggregate_CircularCut-ge)/ge<<")\n";
	  }
         

	  /*te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated->Fill(ge, ((te_aggregate_CircularCut/recalibration_factor)-ge)/ge);
	    r early on, it took me a while to realize what I need is best done by a weighted sum instead of trying to normalize the fit functions. if(debug==1){
	    cout<<"(ge, ((te_aggregate_CircularCut/recalibration_factor)-ge)/ge: ("<<ge<<", "<<((te_aggregate_CircularCut/recalibration_factor)-ge)/ge<<")\n";
	    }
	    te_minus_ge_by_ge_ge_EtaCut_Recalibrated_temp->Fill(ge, ((te_aggregate_CircularCut/recalibration_factor)-ge)/ge);*/
   
	  te_by_ge_ge_EtaCut->Fill(ge, te_aggregate/ge);

	  if(debug==1){	
	    cout<<"(ge, te_aggregate/ge): ("<<ge<<", "<<te_aggregate/ge<<")\n";
	  }

	  te_by_ge_ge_EtaCut_CircularCut->Fill(ge, te_aggregate_CircularCut/ge);
	  te_by_ge_ge_EtaCut_CircularCut_HCALIN->Fill(ge, te_aggregate_HCALIN_CircularCut/ge);
	  te_by_ge_ge_EtaCut_CircularCut_HCALOUT->Fill(ge, te_aggregate_HCALOUT_CircularCut/ge);
	  te_by_ge_ge_EtaCut_CircularCut_CEMC->Fill(ge, (te_aggregate_CircularCut-te_aggregate_HCALIN_CircularCut-te_aggregate_HCALOUT_CircularCut)/ge);

	  mean_te_by_ge_ge_EtaCut_CircularCut_HCALIN->Fill(ge, te_aggregate_HCALIN_CircularCut/ge);
	  mean_te_by_ge_ge_EtaCut_CircularCut_HCALOUT->Fill(ge, te_aggregate_HCALOUT_CircularCut/ge);
	  mean_te_by_ge_ge_EtaCut_CircularCut_CEMC->Fill(ge, (te_aggregate_CircularCut-te_aggregate_HCALIN_CircularCut-te_aggregate_HCALOUT_CircularCut)/ge);
     
     
	  if(debug==1){
	    cout<<"\n(ge, te_aggregate_CircularCut/ge): ("<<ge<<", "<<te_aggregate_CircularCut/ge<<")\n"; 
	    cout<<"(ge, te_aggregate_HCALIN_CircularCut/ge): ("<<ge<<", "<<te_aggregate_HCALIN_CircularCut/ge<<")\n"; 
	    cout<<"(ge, te_aggregate_HCALOUT_CircularCut/ge): ("<<ge<<", "<<te_aggregate_HCALOUT_CircularCut/ge<<")\n";
	    cout<<"(ge, (te_aggregate_CircularCut-te_aggregate_HCALIN_CircularCut - te_aggregate_HCALOUT_CircularCut)/ge): ("<<ge<<", "<<(te_aggregate_CircularCut-te_aggregate_HCALIN_CircularCut - te_aggregate_HCALOUT_CircularCut)/ge<<")\n\n"; 
	  }      
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

  if(debug==1){
    std::cout<<"\nGenerating sigma_e vs ge plots\n\n";
  }

  Double_t recalibrationArr1[nSlicesx];
  Double_t rf_integral_HCALIN = 0;
  Double_t weight_HCALIN = te_by_ge_ge_EtaCut_CircularCut_HCALIN->GetMean(2);

  //recalibrationArr1[0] = recalibration_firstSlice_FHCAL;
  //rf_integral_FHCAL += recalibrationArr1[0];
  //cout << "Recalibration factor for slice " << 1 << " of FHCAL is: " <<  recalibrationArr1[0] << endl;

  for(int binIter = 1; binIter <= nSlicesx; binIter++){
    recalibrationArr1[binIter-1] = mean_te_by_ge_ge_EtaCut_CircularCut_HCALIN->GetBinContent(binIter);    
    rf_integral_HCALIN += recalibrationArr1[binIter-1];

    if(debug == 1){
      cout<<"rf_integral_HCALIN += "<<recalibrationArr1[binIter-1]<<"\n";
    }

    cout << "Recalibration factor for slice " << binIter << " of HCALIN is: " <<  recalibrationArr1[binIter-1] << endl;
  }

  if(debug == 1){  
    cout<<"rf_integral_HCALIN = "<<rf_integral_HCALIN<<"\n\n";
  }

  Double_t recalibrationArr2[nSlicesx];
  Double_t rf_integral_HCALOUT = 0;
  Double_t weight_HCALOUT = te_by_ge_ge_EtaCut_CircularCut_HCALOUT->GetMean(2);

  //recalibrationArr2[0] = recalibration_firstSlice_HCALOUT;
  //rf_integral_HCALOUT += recalibrationArr2[0];
  //cout << "Recalibration factor for slice " << 1 << " of HCALOUT is: " <<  recalibrationArr2[0] << endl;

  for(int binIter = 1; binIter <= nSlicesx; binIter++){
    recalibrationArr2[binIter-1] = mean_te_by_ge_ge_EtaCut_CircularCut_HCALOUT->GetBinContent(binIter);
    rf_integral_HCALOUT += recalibrationArr2[binIter-1];

    if(debug == 1){   
      cout<<"rf_integral_HCALOUT += "<<recalibrationArr2[binIter-1]<<"\n";
    }

    cout << "Recalibration factor for slice " << binIter << " of HCALOUT is: " <<  recalibrationArr2[binIter-1] << endl;
  }

  if(debug == 1){
    cout<<"rf_integral_HCALOUT = "<<rf_integral_HCALOUT<<"\n\n";
  }


  Double_t recalibrationArr3[nSlicesx];
  Double_t rf_integral_CEMC = 0;
  Double_t weight_CEMC = te_by_ge_ge_EtaCut_CircularCut_CEMC->GetMean(2);

  //recalibrationArr3[0] = recalibration_firstSlice_CEMC;
  //rf_integral_CEMC += recalibrationArr3[0];
  //cout << "Recalibration factor for slice " << 1 << " of CEMC is: " <<  recalibrationArr3[0] << endl;

  for(int binIter = 1; binIter <= nSlicesx; binIter++){
    recalibrationArr3[binIter-1] = mean_te_by_ge_ge_EtaCut_CircularCut_CEMC->GetBinContent(binIter);
    rf_integral_CEMC += recalibrationArr3[binIter-1];

    if(debug == 1){   
      cout<<"rf_integral_CEMC += "<<recalibrationArr3[binIter-1]<<"\n";
    }

    cout << "Recalibration factor for slice " << binIter << " of CEMC is: " <<  recalibrationArr3[binIter-1] << endl;
  }

  if(debug == 1){
    cout<<"rf_integral_CEMC = "<<rf_integral_CEMC<<"\n\n";
  }

  for(int i=0; i<T1->GetEntries(); i++){

    T1->GetEntry(i);
    T2->GetEntry(i);
    T3->GetEntry(i);

    Double_t geta1 = evaltree1->get_geta();
    Double_t gphi = evaltree1->get_gphi();
    Double_t gtheta = evaltree1->get_gtheta();
    Double_t ge = evaltree1->get_ge();
    Double_t te_aggregate_CircularCut = 0;
    Double_t te_aggregate_CircularCut_normalised = 0;    

    int recalibration_factor = ceil((ge/30.0)*(Double_t)nSlicesx)- 1;

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
	      std::cout<<"pow(dphi/y_radius_HCALIN,2)+pow(dtheta/x_radius_HCALIN,2): "<<pow(dphi/y_radius_HCALIN,2)+pow(dtheta/x_radius_HCALIN,2)<<"\n";    
	    }

	    if (pow(dphi/y_radius_HCALIN,2)+pow(dtheta/x_radius_HCALIN,2)<=1){
              te_aggregate_CircularCut += te;
	      te_aggregate_CircularCut_normalised += te*weight_HCALIN/recalibrationArr1[recalibration_factor];
	    }
	  }
   	}
      }

      for (int j=0; j<evaltree2->get_ntowers(); j++){

	EvalTower *twr2 = evaltree2->get_tower(j);
	if (twr2){ 

	  Double_t tphi = twr2->get_tphi();
	  Double_t ttheta = twr2->get_ttheta();
	  Double_t dphi = tphi - gphi;
	  Double_t dtheta = ttheta - gtheta;	 
	  Double_t te = twr2->get_te();

          if(te > energyCut){
    
	    if(debug==1){
	      std::cout<<"pow(dphi/y_radius_HCALOUT,2)+pow(dtheta/x_radius_HCALOUT,2): "<<pow(dphi/y_radius_HCALOUT,2)+pow(dtheta/x_radius_HCALOUT,2)<<"\n";    
	    }

	    if (pow(dphi/y_radius_HCALOUT,2)+pow(dtheta/x_radius_HCALOUT,2)<=1){
	      te_aggregate_CircularCut += te;
	      te_aggregate_CircularCut_normalised += te*weight_HCALOUT/recalibrationArr2[recalibration_factor];
	    }
	  }
   	}
      }

      for (int j=0; j<evaltree3->get_ntowers(); j++){

	EvalTower *twr3 = evaltree3->get_tower(j);
	if (twr3){ 

	  Double_t tphi = twr3->get_tphi();
	  Double_t ttheta = twr3->get_ttheta();
	  Double_t dphi = tphi - gphi;
	  Double_t dtheta = ttheta - gtheta;	 
	  Double_t te = twr3->get_te();

	  if(MIP_theta_parametrisation == 1){
	    EMC_cut = mip_pmzn_energy_cut_ftheta->Eval(gtheta);
	  }

          if(te > energyCut + EMC_cut){
    
	    if(debug==1){
	      std::cout<<"pow(dphi/y_radius_CEMC,2)+pow(dtheta/x_radius_CEMC,2): "<<pow(dphi/y_radius_CEMC,2)+pow(dtheta/x_radius_CEMC,2)<<"\n";    
	    }

	    if (pow(dphi/y_radius_CEMC,2)+pow(dtheta/x_radius_CEMC,2)<=1){
	      te_aggregate_CircularCut += te;
	      te_aggregate_CircularCut_normalised += te*weight_CEMC/recalibrationArr3[recalibration_factor];
	    }
	  }
   	}
      }
 
      if(te_aggregate_CircularCut > energyCutAggregate){
	mean_te_by_ge_ge_EtaCut_CircularCut->Fill(ge, te_aggregate_CircularCut_normalised/ge);
      }

      if(debug == 1){
	cout<<"(ge, te_aggregate_CircularCut_normalised/ge): ("<<ge<<", "<<te_aggregate_CircularCut_normalised/ge<<")\n";
      }   
    }    
  }


  Double_t recalibrationArr[nSlicesx];

  // recalibrationArr[0] = recalibration_firstSlice;
  // cout << "Recalibration factor for slice " << 1 << " is: " <<  recalibrationArr[0] << endl;

  for(int binIter = 1; binIter <= nSlicesx; binIter++){
    recalibrationArr[binIter-1] = mean_te_by_ge_ge_EtaCut_CircularCut->GetBinContent(binIter);
    cout << "Recalibration factor for slice " << binIter << " is: " <<  recalibrationArr[binIter-1] << endl;
  }


  for(int i=0; i<T1->GetEntries(); i++){

    T1->GetEntry(i);
    T2->GetEntry(i);
    T3->GetEntry(i);

    Double_t geta1 = evaltree1->get_geta();
    Double_t gphi = evaltree1->get_gphi();
    Double_t gtheta = evaltree1->get_gtheta();
    Double_t ge = evaltree1->get_ge();
    Double_t te_aggregate_CircularCut = 0;
    Double_t te_aggregate_CircularCut_normalised = 0;    

    int recalibration_factor = ceil((ge/30.0)*(Double_t)nSlicesx)- 1;

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
	      std::cout<<"pow(dphi/y_radius_HCALIN,2)+pow(dtheta/x_radius_HCALIN,2): "<<pow(dphi/y_radius_HCALIN,2)+pow(dtheta/x_radius_HCALIN,2)<<"\n";    
	    }

	    if (pow(dphi/y_radius_HCALIN,2)+pow(dtheta/x_radius_HCALIN,2)<=1){
	      te_aggregate_CircularCut += te;
	      te_aggregate_CircularCut_normalised += te*weight_HCALIN/recalibrationArr1[recalibration_factor];
	    }
	  }
   	}
      }

      for (int j=0; j<evaltree2->get_ntowers(); j++){

	EvalTower *twr2 = evaltree2->get_tower(j);
	if (twr2){ 

	  Double_t tphi = twr2->get_tphi();
	  Double_t ttheta = twr2->get_ttheta();
	  Double_t dphi = tphi - gphi;
	  Double_t dtheta = ttheta - gtheta;	 
	  Double_t te = twr2->get_te();

          if(te > energyCut){
    
	    if(debug==1){
	      std::cout<<"pow(dphi/y_radius_HCALOUT,2)+pow(dtheta/x_radius_HCALOUT,2): "<<pow(dphi/y_radius_HCALOUT,2)+pow(dtheta/x_radius_HCALOUT,2)<<"\n";    
	    }

	    if (pow(dphi/y_radius_HCALOUT,2)+pow(dtheta/x_radius_HCALOUT,2)<=1){
	      te_aggregate_CircularCut += te;
	      te_aggregate_CircularCut_normalised += te*weight_HCALOUT/recalibrationArr2[recalibration_factor];
	    }
	  }
   	}
      }

      for (int j=0; j<evaltree3->get_ntowers(); j++){

	EvalTower *twr3 = evaltree3->get_tower(j);
	if (twr3){ 

	  Double_t tphi = twr3->get_tphi();
	  Double_t ttheta = twr3->get_ttheta();
	  Double_t dphi = tphi - gphi;
	  Double_t dtheta = ttheta - gtheta;	 
	  Double_t te = twr3->get_te();

	  if(MIP_theta_parametrisation == 1){
	    EMC_cut = mip_pmzn_energy_cut_ftheta->Eval(gtheta);
	  }

          if(te > energyCut + EMC_cut){
    
	    if(debug==1){
	      std::cout<<"pow(dphi/y_radius_CEMC,2)+pow(dtheta/x_radius_CEMC,2): "<<pow(dphi/y_radius_CEMC,2)+pow(dtheta/x_radius_CEMC,2)<<"\n";    
	    }

	    if (pow(dphi/y_radius_CEMC,2)+pow(dtheta/x_radius_CEMC,2)<=1){
	      te_aggregate_CircularCut += te;
	      te_aggregate_CircularCut_normalised += te*weight_CEMC/recalibrationArr3[recalibration_factor];
	    }
	  }
   	}
      }
      
      if(te_aggregate_CircularCut > energyCutAggregate){
	te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated->Fill(ge, ((te_aggregate_CircularCut_normalised/recalibrationArr[recalibration_factor])-ge)/ge);
	te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp->Fill(ge, ((te_aggregate_CircularCut_normalised/recalibrationArr[recalibration_factor])-ge)/ge);
      }

      if(debug==1){
        cout<<"(ge, ((te_aggregate_CircularCut_normalised/recalibration_factor)-ge)/ge: ("<<ge<<", "<<((te_aggregate_CircularCut_normalised/recalibrationArr[recalibration_factor])-ge)/ge<<")\n";
      }
    }    
  }


  // Initialising fit functions 
  TF1 *fit = new TF1("fit", "gaus");   
  TF1 *fit1 = new TF1("fit1","gaus",fit_min,fit_max);
  TF1 *fExp = new TF1("fExp","0.1 + 1.0/sqrt(x)",0,30);
  TF1 *fTrue = new TF1("fTrue","[0] + [1]/sqrt(x)",0,30);


  te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp->FitSlicesY(0, 1, -1, 0, "QN");
  TH2D *te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2 = (TH2D*)gDirectory->Get("te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2");
  TH2D *te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_1 = (TH2D*)gDirectory->Get("te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_1");
  TH2D *te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_chi2 = (TH2D*)gDirectory->Get("te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_chi2");



  TH1D* slices[nSlicesx];

  // Generating plots for individual slices
  for(int sno = 0; sno < nSlicesx; sno++){
    int plusOne = sno+1;
    TString sname = "slice " + arr[sno];
    slices[sno] = te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp->ProjectionY(sname, plusOne, plusOne);
  }


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

  mean_te_by_ge_ge_EtaCut_CircularCut->GetXaxis()->SetTitle("Generated Energy (GeV)");
  mean_te_by_ge_ge_EtaCut_CircularCut->GetXaxis()->SetLabelSize(0.05);  
  mean_te_by_ge_ge_EtaCut_CircularCut->GetXaxis()->SetTitleSize(0.05);
  mean_te_by_ge_ge_EtaCut_CircularCut->GetYaxis()->SetTitle("Mean of te_{agg}/ge");
  mean_te_by_ge_ge_EtaCut_CircularCut->GetYaxis()->SetLabelSize(0.05);  
  mean_te_by_ge_ge_EtaCut_CircularCut->GetYaxis()->SetTitleSize(0.05);

  mean_te_by_ge_ge_EtaCut_CircularCut_HCALIN->GetXaxis()->SetTitle("Generated Energy (GeV)");
  mean_te_by_ge_ge_EtaCut_CircularCut_HCALIN->GetXaxis()->SetLabelSize(0.05);  
  mean_te_by_ge_ge_EtaCut_CircularCut_HCALIN->GetXaxis()->SetTitleSize(0.05);
  mean_te_by_ge_ge_EtaCut_CircularCut_HCALIN->GetYaxis()->SetTitle("Mean of te_{agg}/ge (HCALIN)");
  mean_te_by_ge_ge_EtaCut_CircularCut_HCALIN->GetYaxis()->SetLabelSize(0.05);  
  mean_te_by_ge_ge_EtaCut_CircularCut_HCALIN->GetYaxis()->SetTitleSize(0.05);

  mean_te_by_ge_ge_EtaCut_CircularCut_HCALOUT->GetXaxis()->SetTitle("Generated Energy (GeV)");
  mean_te_by_ge_ge_EtaCut_CircularCut_HCALOUT->GetXaxis()->SetLabelSize(0.05);  
  mean_te_by_ge_ge_EtaCut_CircularCut_HCALOUT->GetXaxis()->SetTitleSize(0.05);
  mean_te_by_ge_ge_EtaCut_CircularCut_HCALOUT->GetYaxis()->SetTitle("Mean of te_{agg}/ge (HCALOUT)");
  mean_te_by_ge_ge_EtaCut_CircularCut_HCALOUT->GetYaxis()->SetLabelSize(0.05);  
  mean_te_by_ge_ge_EtaCut_CircularCut_HCALOUT->GetYaxis()->SetTitleSize(0.05);

  mean_te_by_ge_ge_EtaCut_CircularCut_CEMC->GetXaxis()->SetTitle("Generated Energy (GeV)");
  mean_te_by_ge_ge_EtaCut_CircularCut_CEMC->GetXaxis()->SetLabelSize(0.05);  
  mean_te_by_ge_ge_EtaCut_CircularCut_CEMC->GetXaxis()->SetTitleSize(0.05);
  mean_te_by_ge_ge_EtaCut_CircularCut_CEMC->GetYaxis()->SetTitle("Mean of te_{agg}/ge (CEMC)");
  mean_te_by_ge_ge_EtaCut_CircularCut_CEMC->GetYaxis()->SetLabelSize(0.05);  
  mean_te_by_ge_ge_EtaCut_CircularCut_CEMC->GetYaxis()->SetTitleSize(0.05);

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
    std::cout<<"\nWrite Histograms to File\n";
  }

  TFile *f = new TFile("energy_verification_EtaCut_CircularCut" + detector + ".root","RECREATE"); 
    
  f->GetList()->Add(te_by_ge_ge_EtaCut);
  f->GetList()->Add(te_by_ge_ge_EtaCut_CircularCut);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut_CircularCut);
  f->GetList()->Add(mean_te_by_ge_ge_EtaCut_CircularCut);
  f->GetList()->Add(mean_te_by_ge_ge_EtaCut_CircularCut_HCALIN);
  f->GetList()->Add(mean_te_by_ge_ge_EtaCut_CircularCut_HCALOUT);
  f->GetList()->Add(mean_te_by_ge_ge_EtaCut_CircularCut_CEMC);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_1);
  f->GetList()->Add(te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_chi2);
  
  for(int sno = 0; sno < nSlicesx; sno++){
    f->GetList()->Add(slices[sno]);
  }

  gStyle -> SetOptStat(11);
  gStyle -> SetOptFit(112);
  
  /*
  int sno = 0;
  TString plusOne = (TString)(sno + 1);
  TString nameF = detector + "_sigmaE_slice" + arr[sno] + "_EtaCut_CircularCut.png";
  slices[sno] -> Fit("fit1", "R+");
  slices[sno] -> Draw("hist same");
  c->Print(nameF);
  */

  /*  double_t mean = fit1->GetParameter(1);
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
  TLine *Point = new TLine(30.0/(2.0*nSlicesx),sigma,30.0/(2.0*nSlicesx),sigma);*/

  if(print==1){
    if(debug==1){
      std::cout<<"\nSaving Histograms as .png\n";
    }

    gStyle -> SetOptStat(0);

    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_1->SetAxisRange(mean_min,mean_max,"Y");
    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_1->Draw();
    /*Mean->SetLineColor(1);
    Mean->SetLineWidth(1); 
    Mean_error->SetLineColor(1);
    Mean_error->SetLineWidth(1);
    Mean->Draw("same");
    Mean_error->Draw("same");*/
    c->Print(detector + "_meanE_ge_EtaCut_CircularCut.png");

    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_chi2->SetAxisRange(chi2_min,chi2_max,"Y");
    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_chi2->Draw();
    /*Chi2->SetLineColor(1);
    Chi2->SetLineWidth(1); 
    Periphery_Chi2->SetLineColor(1);
    Periphery_Chi2->SetLineWidth(1);
    Chi2->Draw("same");
    Periphery_Chi2->Draw("same");*/
    c->Print(detector + "_chi2E_ge_EtaCut_CircularCut.png");
 
    gStyle -> SetOptStat(0);
    gStyle -> SetOptFit(0);

    /*Sigma->SetLineColor(1);
    Sigma->SetLineWidth(1); 
    Sigma_error->SetLineColor(1);
    Sigma_error->SetLineWidth(1);
    Point->SetLineColor(3);
    Point->SetLineWidth(5);*/
    fExp->SetLineColor(4); //38
    fTrue->SetLineColor(2); //46
  
    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2->SetMarkerStyle(kFullCircle);
    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2->SetMarkerColor(46); //30
    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2->SetMarkerSize(0.75);
    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2->SetAxisRange(sigma_min,sigma_max,"Y");

    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2->Fit("fTrue", "M+");
    te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2->Draw("same");
    fExp->Draw("same");
    /*Sigma->Draw("same");
    Sigma_error->Draw("same");
    Point->Draw("same");*/
  
    TLegend* legend = new TLegend(1.75,1.75);
    legend->SetHeader("Legend", "C");
    //legend->AddEntry(te_minus_ge_by_ge_ge_EtaCut_CircularCut_Recalibrated_temp_2, "#sigma_{e_{agg}} vs Generated Energy", "flep");
    //legend->AddEntry((TObject*)0,"","");
    legend->AddEntry(fTrue, "p_{0} + p_{1}/#sqrt{ge} (Fitted)", "l");
    legend->AddEntry((TObject*)0,"","");
    legend->AddEntry(fExp, "0.1 + 1.0/#sqrt{ge} (Requirement)", "l");
    legend->SetTextSize(0.033);
    legend->Draw();

    std::cout<<"reduced_chi2 of fit: "<<fTrue->GetChisquare()/fTrue->GetNDF()<<"\n";

    c->Print(detector + "_sigmaE_ge_EtaCut_CircularCut.png");

    gStyle -> SetOptStat(11);
    gStyle -> SetOptFit(112);

    for(int sno = 0; sno < nSlicesx; sno++){
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

    gStyle -> SetOptStat(0);

    mean_te_by_ge_ge_EtaCut_CircularCut->Draw();
    c->Print(detector + "_mean_te_by_ge_ge_EtaCut_CircularCut.png");

    mean_te_by_ge_ge_EtaCut_CircularCut_HCALIN->Draw();
    c->Print(detector + "_mean_te_by_ge_ge_EtaCut_CircularCut_HCALIN.png");

    mean_te_by_ge_ge_EtaCut_CircularCut_HCALOUT->Draw();
    c->Print(detector + "_mean_te_by_ge_ge_EtaCut_CircularCut_HCALOUT.png");

    mean_te_by_ge_ge_EtaCut_CircularCut_CEMC->Draw();
    c->Print(detector + "_mean_te_by_ge_ge_EtaCut_CircularCut_CEMC.png");

    gStyle -> SetOptStat(1);
    c->Close();

  }

  std::cout << "The total te is: " << total_te << endl;
  std::cout << "The total te_CircularCut is: " << total_te_CircularCut << endl;
  std::cout << "The total ge is: " << total_ge << endl; 
  std::cout<<"\n\nDone\n----------------------------------------------------------------------\n\n";
 
}
