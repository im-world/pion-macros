/*
- Reference statistics to debug LoopEvalMultiH*.C
- Combined CEMC+HCALIN+HCALOUT
- Standalone
- Output file - None (LoopEvalMultiHVerify.txt when ran using script)
- Plots made - None
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

void LoopEvalMultiHVerify()
{
   

  TString detector = "HCALIN_HCALOUT_CEMC";
  TFile *f1 = new TFile("merged_Eval_HCALIN.root","READ"); //change the input evaluator root file here
  TFile *f2 = new TFile("merged_Eval_HCALOUT.root","READ");
  TFile *f3 = new TFile("merged_Eval_CEMC.root","READ"); //change the input evaluator root file here
  
  TCanvas *c = new TCanvas();
  gStyle->SetTitleW(0.8);
  gStyle->SetTitleH(0.1);

  TTree* T1 = (TTree*)f1->Get("T");
  EvalRootTTree *evaltree1 = nullptr;
  
  TTree* T2 = (TTree*)f2->Get("T");
  EvalRootTTree *evaltree2 = nullptr;

  TTree* T3 = (TTree*)f3->Get("T");
  EvalRootTTree *evaltree3 = nullptr;

  // ce/ge vs ge
  // ce - te vs ge
  // te/ge vs ge

  Double_t total_te = 0;
  Double_t total_ce = 0;
  Double_t total_ge = 0;
  Double_t eta_min, eta_max;
  TString cut_text;

     eta_min = -1.1;
     eta_max = 1.1;
     cut_text = " {-1.1 < true_eta < 1.1, meas_e > 0.2 GeV} ";
 
  
  T1->SetBranchAddress("DST#EvalTTree_HCALIN",&evaltree1);
  T2->SetBranchAddress("DST#EvalTTree_HCALOUT",&evaltree2);
  T3->SetBranchAddress("DST#EvalTTree_CEMC",&evaltree3);
  for(int i=0; i<T1->GetEntries(); i++)
  {
    T1->GetEntry(i);
    T2->GetEntry(i);
    T3->GetEntry(i);
   
    //   cout << "ge: " << evaltree->get_ge() << endl;
   
    std::cout<<"\n\n\n------------------------------------------\nParticle: "<<i<<"\n\n";
    std::cout<<"Initial Parameters "<<"\n";
    Double_t ge = evaltree1->get_ge();
    Double_t geta = evaltree1->get_geta();
    Double_t gphi = evaltree1->get_gphi();
    Double_t gtheta = evaltree1->get_gtheta();
    Double_t gpx = evaltree1->get_gpx(); 
    Double_t gpy = evaltree1->get_gpy(); 
    Double_t gpz = evaltree1->get_gpz();
    Double_t gp = sqrt(pow(gpx,2) + pow(gpy,2) + pow(gpz,2)); 

    std::cout<<"ge is: "<<ge<<"\n";
    std::cout<<"geta is: "<<geta<<"\n";
    std::cout<<"gphi is: "<<gphi<<"\n";
    std::cout<<"gtheta is: "<<gtheta<<"\n";
    std::cout<<"gpx is: "<<gpx<<"\n";
    std::cout<<"gpy is: "<<gpy<<"\n";
    std::cout<<"gpz is: "<<gpz<<"\n";
    std::cout<<"gp is: "<<gp<<"\n\n";
  
    // cut_text = "";
    if(geta>=eta_min && geta<=eta_max)
    {
        total_ge += ge;
        int clus_count_HCALIN = 0;
        int clus_count_HCALOUT = 0;
        int clus_count_CEMC = 0;
        Double_t ce_aggregate = 0;

        int twr_count_HCALIN = 0;
        int twr_count_HCALOUT = 0;
        int twr_count_CEMC = 0;
        Double_t te_aggregate = 0;
    
    //   cout << "Number of clusters: " << evaltree->get_nclusters() << endl;

    for (int k=0; k<evaltree1->get_ntowers(); k++) // changed iteration varaible from i to j 
    {
      std::cout<<"\nHCALIN Tower: "<<k<<"\n";
      Double_t te = 0;
      EvalTower *twr1 = evaltree1->get_tower(k);

      if (twr1)
      {
        cout<<"non-empty HCALIN Tower\n";
        te = twr1->get_te();
        std::cout<<"HCALIN te is: "<<te<<"\n";
        total_te += te;
        twr_count_HCALIN += 1;    
        te_aggregate += te;  
      }
      std::cout<<"te_aggregate is: "<<te_aggregate<<"\n\n";
    }

    for (int k=0; k<evaltree2->get_ntowers(); k++) // changed iteration varaible from i to j 
    {
      std::cout<<"\nHCALOUT Tower: "<<k<<"\n";
      Double_t te = 0;
      EvalTower *twr2 = evaltree2->get_tower(k);

      if (twr2)
      {
        cout<<"non-empty HCALOUT Tower\n";
        te = twr2->get_te();
        std::cout<<"HCALOUT te is: "<<te<<"\n";
        total_te += te;
        twr_count_HCALOUT += 1;    
        te_aggregate += te;  
      }
       std::cout<<"te_aggregate is: "<<te_aggregate<<"\n\n";
    }

    for (int k=0; k<evaltree3->get_ntowers(); k++) // changed iteration varaible from i to j 
    {
      std::cout<<"\nCEMC Tower: "<<k<<"\n";
      Double_t te = 0;
      EvalTower *twr3 = evaltree3->get_tower(k);
      if (twr3)
      {
        cout<<"non-empty CEMC Tower\n";
        te = twr3->get_te();
        std::cout<<"CEMC te is: "<<te<<"\n";
        total_te += te;
        twr_count_CEMC += 1;    
        te_aggregate += te;  
      }
       std::cout<<"te_aggregate is: "<<te_aggregate<<"\n\n";
    }
       std::cout<<"twr_count_HCALIN is: "<<twr_count_HCALIN<<"\n";
       std::cout<<"twr_count_HCALOUT is: "<<twr_count_HCALOUT<<"\n";
       std::cout<<"twr_count_CEMC is: "<<twr_count_CEMC<<"\n\n";

    for (int j=0; j<evaltree1->get_nclusters(); j++) // changed iteration varaible from i to j 
    {
      std::cout<<"\nHCALIN Cluster: "<<j<<"\n";

      EvalCluster *clus1 = evaltree1->get_cluster(j);
      if (clus1)
      {
        cout<<"non-empty HCALIN cluster\n";
	Double_t ce = clus1->get_ce();
      	Double_t ctowers = clus1->get_ctowers();
        std::cout<<"HCALIN ce is: "<<ce<<"\n";
        std::cout<<"HCALIN ctowers is: "<<ctowers<<"\n\n";
        total_ce += ce;
        clus_count_HCALIN += 1;    
        ce_aggregate += ce;  
      }
       std::cout<<"ce_aggregate is: "<<ce_aggregate<<"\n\n";
    }


    for (int j=0; j<evaltree2->get_nclusters(); j++) // changed iteration varaible from i to j 
    {
      std::cout<<"\nHCALOUT Cluster: "<<j<<"\n";
      EvalCluster *clus2 = evaltree2->get_cluster(j);
      if (clus2)
      {
        cout<<"non-empty HCALIN cluster\n";
	Double_t ce = clus2->get_ce();
      	Double_t ctowers = clus2->get_ctowers();
        std::cout<<"HCALOUT ce is: "<<ce<<"\n";
        std::cout<<"HCALOUT ctowers is: "<<ctowers<<"\n\n";
        total_ce += ce;
        clus_count_HCALOUT += 1;    
        ce_aggregate += ce;  
      }
       std::cout<<"ce_aggregate is: "<<ce_aggregate<<"\n\n";
    }

    for (int j=0; j<evaltree3->get_nclusters(); j++) // changed iteration varaible from i to j 
    {
      std::cout<<"\nCEMC Cluster: "<<j<<"\n";
      EvalCluster *clus3 = evaltree3->get_cluster(j);
      if (clus3)
      {
        cout<<"non-empty CEMC cluster\n";
	Double_t ce = clus3->get_ce();
      	Double_t ctowers = clus3->get_ctowers();
        std::cout<<"CEMC ce is: "<<ce<<"\n";
        std::cout<<"CEMC ctowers is: "<<ctowers<<"\n\n";
        total_ce += ce;
        clus_count_CEMC += 1;    
        ce_aggregate += ce;  
      }

       std::cout<<"ce_aggregate is: "<<ce_aggregate<<"\n\n";
    }
       std::cout<<"clus_count_HCALIN is: "<<clus_count_HCALIN<<"\n\n";
       std::cout<<"clus_count_HCALOUT is: "<<clus_count_HCALOUT<<"\n\n";
       std::cout<<"clus_count_CEMC is: "<<clus_count_CEMC<<"\n\n";            
   }
      std::cout<<"\n total ce = "<<total_ce<<"\n";
      std::cout<<"\n total te = "<<total_te<<"\n\n\n";
  }

  
  std::cout<<"total_ge is: "<<total_ge<<"\n";
  std::cout<<"total_ce is: "<<total_ce<<"\n";
  std::cout<<"total_te is: "<<total_te<<"\n";
  std::cout<<"\n\nDone\n----------------------------------------------------------------------\n\n";
}
