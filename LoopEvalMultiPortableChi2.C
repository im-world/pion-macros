/*
- [sigma_e vs ge] plots 
- Eta Cuts 
- Individual detectors 
- Standalone 
- Output root file - sigma_e_EtaCut_detector.root
- Plots made - [sigma_e vs ge], [mean_e vs ge], [chi2_e vs ge], slices for the [(ce-ge)/ge vs ge] plot.
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

void LoopEvalMultiPortableChi2(TString detector, int print = 1, int debug = 0){
 
  TFile *f1 = new TFile("merged_Eval_" + detector + ".root","READ"); //change the  
  TTree* T1 = (TTree*)f1->Get("T");
  EvalRootTTree *evaltree1 = nullptr;

  T1->SetBranchAddress("DST#EvalTTree_" + detector,&evaltree1);

  TCanvas *c = new TCanvas();
  c->SetTickx();
  c->SetTicky();

  // Changing default plotting styles
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(102);
  gStyle->SetTitleXOffset(1);
  gStyle->SetTitleYOffset(1);
  gStyle->SetLabelSize(0.05);  
  gStyle->SetTitleXSize(0.05);  
  gStyle->SetTitleYSize(0.05);
  
  int nSlices = 21;  // number of slices for ge-axis
  
  TH2D *ce_minus_ge_by_ge_ge_EtaCut = new TH2D("ce_minus_ge_by_ge_ge_EtaCut","#frac{#Delta e_{agg}}{truth e} vs truth e",nSlices,0,30,2000,-2,1);
  
  Double_t eta_min, eta_max;
  TString cut_text, eRes, pRes;

  // Applying detector-specific eta cuts

  if (detector == "FEMC"){
    x_radius = 0.15; //radius of circular cut
    y_radius = 0.4;
    eta_min = 1.3;
    eta_max = 3.3;
    cut_text = " {1.3 < true_eta < 3.3} ";
    eRes = "0.02 + 0.08/sqrt(x) + 0.02/x";
    pRes = "(e-)";

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
    eRes = "0.01 + 0.025/sqrt(x) + 0.01/x";
    pRes = "(e-)";

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
    eRes = "0.025 + 0.13/sqrt(x) + 0.02/x";
    pRes = "(e-)";

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
    eRes = "1";
    pRes = "(e-)";

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
    eRes = "1";
    pRes = "(e-)";

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
    eRes = "1";
    pRes = "(e-)";
      
    if(debug==1){
      std::cout<<"HCALOUT eta cut & circular cut applied"<<"\n\n";
    }
  }

  else{
    std::cout<<"Please try again.";
    return 1;
  }

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
 
      Double_t te_aggregate = 0;    
      for (int j=0; j<evaltree1->get_ntowers(); j++){
	
	if(debug==1){
	  std::cout<<"\n Tower: "<<j<<"\n";
	}
       
	Double_t te = 0;

        EvalTower *twr1 = evaltree1->get_tower(j);
	if (twr1){
	  
	  if(debug==1){
	    cout<<"non-empty tower\n";
	  }

	  te = twr1->get_te();
	  te_aggregate += te;
	   
	  if(debug==1){ 
	    cout<<"te_aggregate till now = "<<te_aggregate<<"\n";
	    std::cout<<"te += "<<twr1->get_te()<<"\n";
	  }
	}
      }

      if(debug==1){
	cout<<"\n";
      }

      Double_t ce_aggregate = 0;
      for (int j=0; j<evaltree1->get_nclusters(); j++){

	if(debug==1){	
	  std::cout<<"\n Cluster: "<<j<<"\n";
	}

        Double_t ce = 0;

	EvalCluster *clus1 = evaltree1->get_cluster(j);
	if (clus1){
	  
	  if(debug==1){	  
	    cout<<"non-empty cluster\n";
	  }	  
	  
	  ce = clus1->get_ce();
	  ce_aggregate += ce;
	  
	  if(debug==1){
	    cout<<"ce_aggregate till now = "<<ce_aggregate<<"\n";
	    std::cout<<"ce += "<<clus1->get_ce()<<"\n";
	  }
	}
      }

      ce_minus_ge_by_ge_ge_EtaCut->Fill(ge, (ce_aggregate-ge)/ge);
      if(debug==1){
	cout<<"(ge, (ce_aggregate-ge)/ge): ("<<ge<<", "<<(ce_aggregate-ge)/ge<<")\n";
      }
      
    }
  }

  if(debug==1){
    cout<<"\neta cut if ends"<<"\n";
  }

  
  // Generating sigma_e plot from ce_minus_ge_by_ge_ge_EtaCut

  TString arr[nSlices]; // Used later for naming

  for(int sno = 0; sno < nSlices; sno++){
    arr[sno] = TString::Itoa(sno + 1, 10);
  }

  ce_minus_ge_by_ge_ge_EtaCut->FitSlicesY(0, 0, -1, 0, "QN");
  TH2D *ce_minus_ge_by_ge_ge_EtaCut_2 = (TH2D*)gDirectory->Get("ce_minus_ge_by_ge_ge_EtaCut_2");
  TH2D *ce_minus_ge_by_ge_ge_EtaCut_1 = (TH2D*)gDirectory->Get("ce_minus_ge_by_ge_ge_EtaCut_1");
  TH2D *ce_minus_ge_by_ge_ge_EtaCut_chi2 = (TH2D*)gDirectory->Get("ce_minus_ge_by_ge_ge_EtaCut_chi2");

  if(debug==1){
    std::cout<<"\nSigma, Mean, Chi2 Plots made\n";
  }

  //Generating individual slices for sigma_e plot

  TH1D* slices[nSlices];

  for(int sno = 0; sno < nSlices; sno++){
    int plusOne = sno+1;
    TString sname = "slice " + arr[sno];
    slices[sno] = ce_minus_ge_by_ge_ge_EtaCut->ProjectionY(sname, plusOne, plusOne);
  }
 
  if(debug==1){
    std::cout<<"Slices made\n";
  }

  TF1 *fit = new TF1("fit", "gaus");   
  TF1 *fExp = new TF1("fExp",eRes,0,30);
  TF1 *fTrue = new TF1("fTrue","[0] + [1]/sqrt(x) + [2]/x",0,30);  
 
  if(debug==1){
    std::cout<<"Fit functions defined\n\n";  
    std::cout<<"\nHistogram Formatting\n\n";
  }

  ce_minus_ge_by_ge_ge_EtaCut_2->GetXaxis()->SetTitle("Generated Energy (GeV)");
  ce_minus_ge_by_ge_ge_EtaCut_2->GetXaxis()->SetLabelSize(0.05);  
  ce_minus_ge_by_ge_ge_EtaCut_2->GetXaxis()->SetTitleSize(0.05);
  ce_minus_ge_by_ge_ge_EtaCut_2->GetYaxis()->SetTitle("#sigma_{e_{agg}}");
  ce_minus_ge_by_ge_ge_EtaCut_2->GetYaxis()->SetLabelSize(0.05);  
  ce_minus_ge_by_ge_ge_EtaCut_2->GetYaxis()->SetTitleSize(0.05);
  //hist_energy_agg_2->SetTitle("#sigma_{e_{agg}} vs true_e" + cut_text);

  ce_minus_ge_by_ge_ge_EtaCut_chi2->GetXaxis()->SetTitle("Generated Energy (GeV)");
  ce_minus_ge_by_ge_ge_EtaCut_chi2->GetXaxis()->SetLabelSize(0.05);  
  ce_minus_ge_by_ge_ge_EtaCut_chi2->GetXaxis()->SetTitleSize(0.05);
  ce_minus_ge_by_ge_ge_EtaCut_chi2->GetYaxis()->SetTitle("reduced_#chi^{2}_{e_{agg}}");
  ce_minus_ge_by_ge_ge_EtaCut_chi2->GetYaxis()->SetLabelSize(0.05);  
  ce_minus_ge_by_ge_ge_EtaCut_chi2->GetYaxis()->SetTitleSize(0.04);
  //hist_energy_agg_chi2->SetTitle("#chi^{2}_{e_{agg}} vs true_e" + cut_text);

  ce_minus_ge_by_ge_ge_EtaCut_1->GetXaxis()->SetTitle("Generated Energy (GeV)");
  ce_minus_ge_by_ge_ge_EtaCut_1->GetXaxis()->SetLabelSize(0.05);  
  ce_minus_ge_by_ge_ge_EtaCut_1->GetXaxis()->SetTitleSize(0.05);
  ce_minus_ge_by_ge_ge_EtaCut_1->GetYaxis()->SetTitle("mean_{e_{agg}}");
  ce_minus_ge_by_ge_ge_EtaCut_1->GetYaxis()->SetLabelSize(0.05);  
  ce_minus_ge_by_ge_ge_EtaCut_1->GetYaxis()->SetTitleSize(0.05);
  //hist_energy_agg_1->SetTitle("mean_{e_{agg}} vs true_e" + cut_text);
	
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
  
  TFile *f = new TFile("sigme_e_EtaCut_" + detector + ".root","RECREATE"); 
  f->GetList()->Add(ce_minus_ge_by_ge_ge_EtaCut);
  f->GetList()->Add(ce_minus_ge_by_ge_ge_EtaCut_2);
  f->GetList()->Add(ce_minus_ge_by_ge_ge_EtaCut_1);
  f->GetList()->Add(ce_minus_ge_by_ge_ge_EtaCut_chi2);
  
  for(int sno = 0; sno < nSlices; sno++){
    f->GetList()->Add(slices[sno]);
  }

  f->Write();

  if(print==1){
    if(debug==1){
      std::cout<<"\nSaving Histograms as .png\n";
    }

  
    gStyle -> SetOptStat(0);

    ce_minus_ge_by_ge_ge_EtaCut_1->Draw("colz");
    c->Print(detector + "meanE_ge_EtaCut.png");

    ce_minus_ge_by_ge_ge_EtaCut_chi2->Draw("colz");
    c->Print(detector + "chi2E_ge_EtaCut.png");
 
    gStyle -> SetOptStat(0);
    gStyle -> SetOptFit(0);
    fExp->SetLineColor(4); //38
    fTrue->SetLineColor(2); //46

    ce_minus_ge_by_ge_ge_EtaCut_2->SetMarkerStyle(kFullCircle);
    ce_minus_ge_by_ge_ge_EtaCut_2->SetMarkerColor(46); //30
    ce_minus_ge_by_ge_ge_EtaCut_2->SetMarkerSize(0.75);
    ce_minus_ge_by_ge_ge_EtaCut_2->SetAxisRange(0,0.5,"Y");

    ce_minus_ge_by_ge_ge_EtaCut_2->Fit("fTrue", "M+");
    ce_minus_ge_by_ge_ge_EtaCut_2->Draw("same");
    fExp->Draw("same");
  
    TLegend* legend = new TLegend(1.75,1.75);
    legend->SetHeader("Legend", "C");
    legend->AddEntry(ce_minus_ge_by_ge_ge_EtaCut_2, "#sigma e_{agg} vs Generated Energy", "flep");
    legend->AddEntry((TObject*)0,"","");
    legend->AddEntry(fTrue, "p_{0} + p_{1}/#sqrt{ge} + p_{2}/ge (Fitted)", "l");
    legend->AddEntry((TObject*)0,"","");
    legend->AddEntry(fExp, eRes + pRes, "l");
    legend->AddEntry((TObject*)0,"(Requirement)","");
    legend->SetTextSize(0.033);
    legend->Draw();

    cout<<"reduced_chi2 of fit: "<<fTrue->GetChisquare()<<"\n";
    cout<<"parameters: "<<fTrue->GetParameters()<<"\n";

    c->Print(detector + "sigmaE_ge_EtaCut.png");

    gStyle -> SetOptStat(11);
    gStyle -> SetOptFit(102);

    for(int sno = 0; sno < nSlices; sno++){
      TString plusOne = (TString)(sno + 1);
      TString nameF = detector + "_sigma_slice" + arr[sno] + ".png";
      slices[sno] -> Fit("fit", "M+");
      slices[sno] -> Draw("hist same");
      c->Print(nameF);
    }
  }
  c->Close();
  std::cout<<"\n\nDone\n----------------------------------------------------------------------\n\n";
 
}
