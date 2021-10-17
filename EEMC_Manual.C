#include <iostream>
#include <cmath>
#include <eicqa_modules/EvalRootTTree.h>
#include <eicqa_modules/EvalHit.h>
#include "TMath.h"
#include "TStyle.h"

R__LOAD_LIBRARY(libeicqa_modules.so)

void EEMC_Manual(){
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

  double_t mean[10];
  double_t mean_error[10];
  double_t sigma[10];
  double_t sigma_error[10];
  double_t chi2[10];
  double_t ge[10] = {2,5,8,11,14,17,20,23,26,29};
  double_t ge_error[10] = {1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5};
  double_t chi2_error[10] = {0,0,0,0,0,0,0,0,0,0};

  ifstream file1;
  file1.open("gaus_1-4.txt"); 
  for(int i = 0; i < 4; i++){
    file1 >> mean[i] >> mean_error[i] >> sigma[i] >> sigma_error[i] >> chi2[i];
  }
  file1.close();
 
  ifstream file2;
  file2.open("crystalball_5-10.txt");
  for(int i = 4; i < 10; i++){
    file2 >> mean[i] >> mean_error[i] >> sigma[i] >> sigma_error[i] >> chi2[i];
  }
  file2.close();

  TGraphErrors* Mean = new TGraphErrors(10,ge,mean,ge_error,mean_error);
  Mean->GetXaxis()->SetTitle("Generated Energy (GeV)");
  Mean->GetXaxis()->SetLabelSize(0.05);  
  Mean->GetXaxis()->SetTitleSize(0.05);
  Mean->GetYaxis()->SetTitle("Mean_{e_{agg}}");
  Mean->GetYaxis()->SetLabelSize(0.05);  
  Mean->GetYaxis()->SetTitleSize(0.05);

  Mean->GetXaxis()->SetLimits(0.,30.);
  Mean->GetHistogram()->SetMinimum(-0.1);
  Mean->GetHistogram()->SetMaximum(0.1);
  Mean->Draw("apz"); 
  c->Update();
  c->Modified();
  c->Print("EEMC_Mean.png");

  TGraphErrors* Chi2 = new TGraphErrors(10,ge,chi2,ge_error,chi2_error);
  Chi2->GetXaxis()->SetTitle("Generated Energy (GeV)");
  Chi2->GetXaxis()->SetLabelSize(0.05);  
  Chi2->GetXaxis()->SetTitleSize(0.05);
  Chi2->GetYaxis()->SetTitle("Reduced_#chi^{2}_{e_{agg}}");
  Chi2->GetYaxis()->SetLabelSize(0.05);  
  Chi2->GetYaxis()->SetTitleSize(0.04);

  Chi2->GetXaxis()->SetLimits(0.,30.);
  Chi2->Draw("apz");
  c->Print("EEMC_Chi2.png");

  TGraphErrors* Sigma = new TGraphErrors(10,ge,sigma,ge_error,sigma_error);
  Sigma->GetXaxis()->SetTitle("Generated Energy (GeV)");
  Sigma->GetXaxis()->SetLabelSize(0.05);  
  Sigma->GetXaxis()->SetTitleSize(0.05);
  Sigma->GetYaxis()->SetTitle("#sigma_{e_{agg}}");
  Sigma->GetYaxis()->SetLabelSize(0.05);  
  Sigma->GetYaxis()->SetTitleSize(0.05);

  TF1 *fExp = new TF1("fExp","0.01 + 0.025/sqrt(x) + 0.01/x",0,30);
  TF1 *fTrue = new TF1("fTrue","[0] + [1]/sqrt(x) + [2]/x",0,30);

  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(0);

  TLegend* legend = new TLegend(1.75,1.75);
  legend->SetHeader("Legend", "C");
  legend->AddEntry(Sigma, "#sigma_{e_{agg}} vs Generated Energy", "flep");
  legend->AddEntry((TObject*)0,"","");
  legend->AddEntry(fTrue, "p_{0} + p_{1}/#sqrt{ge} + p_{2}/ge (Fitted)", "l");
  legend->AddEntry((TObject*)0,"","");
  legend->AddEntry(fExp, "0.01 + 0.025/#sqrt{ge} + 0.01/ge", "l");
  legend->AddEntry((TObject*)0,"(Requirement)","" );
  legend->SetTextSize(0.033);

  Sigma->GetHistogram()->SetMinimum(0);
  Sigma->GetHistogram()->SetMaximum(0.3);
  Sigma->GetXaxis()->SetLimits(0.,30.);

  Sigma->SetMarkerSize(0.75);
  Sigma->SetMarkerStyle(kFullCircle);
  Sigma->SetMarkerColor(46);

  fExp->SetLineColor(4); //38
  fTrue->SetLineColor(2); //46

  Sigma->Fit("fTrue", "M+");
  Sigma->Draw("apz same");
  fExp->Draw("same");
  legend->Draw();
  c->Print("EEMC_Sigma.png");
}
