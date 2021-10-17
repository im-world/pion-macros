/*
- Reference statistics to debug LoopEvalMultiF*.C
- Combined FEMC+FHCAL
- Standalone
- Output file - None (LoopEvalMultiFVerify.txt when ran using script)
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

void LoopEvalMultiFVerify()
{
   

  TString detector = "FHCAL_FEMC";
  TFile *f1 = new TFile("merged_Eval_FHCAL.root","READ");
  TFile *f2 = new TFile("merged_Eval_FEMC.root","READ");
  
  TCanvas *c = new TCanvas();
  gStyle->SetTitleW(0.8);
  gStyle->SetTitleH(0.1);

  TTree* T1 = (TTree*)f1->Get("T");
  EvalRootTTree *evaltree1 = nullptr;
  
  TTree* T2 = (TTree*)f2->Get("T");
  EvalRootTTree *evaltree2 = nullptr;

  Double_t total_te = 0;
  Double_t total_ce = 0;
  Double_t total_ge = 0;
  Double_t eta_min, eta_max;
  TString cut_text;

  eta_min = 1.3;
  eta_max = 3.3;
  cut_text = " {1.3 < true_eta < 3.3, meas_e > 0.2 GeV} ";
 
  
  T1->SetBranchAddress("DST#EvalTTree_FHCAL",&evaltree1);
  T2->SetBranchAddress("DST#EvalTTree_FEMC",&evaltree2);
  for(int i=0; i<T1->GetEntries(); i++)
    {
      T1->GetEntry(i);
      T2->GetEntry(i);
   
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
	  int clus_count_FHCAL = 0;
	  int clus_count_FEMC = 0;
	  Double_t ce_aggregate = 0;

	  int twr_count_FHCAL = 0;
	  int twr_count_FEMC = 0;
	  Double_t te_aggregate = 0;
    
	  for (int k=0; k<evaltree1->get_ntowers(); k++) // changed iteration varaible from i to j 
	    {
	      std::cout<<"\nFHCAL Tower: "<<k<<"\n";
	      Double_t te = 0;
	      EvalTower *twr1 = evaltree1->get_tower(k);

	      if (twr1)
		{
		  cout<<"non-empty FHCAL Tower\n";
		  te = twr1->get_te();
		  std::cout<<"FHCAL te is: "<<te<<"\n";
		  total_te += te;
		  twr_count_FHCAL += 1;    
		  te_aggregate += te;  
		}
	      std::cout<<"te_aggregate is: "<<te_aggregate<<"\n\n";
	    }

	  for (int k=0; k<evaltree2->get_ntowers(); k++) // changed iteration varaible from i to j 
	    {
	      std::cout<<"\nFEMC Tower: "<<k<<"\n";
	      Double_t te = 0;
	      EvalTower *twr2 = evaltree2->get_tower(k);

	      if (twr2)
		{
		  cout<<"non-empty FEMC Tower\n";
		  te = twr2->get_te();
		  std::cout<<"FEMC te is: "<<te<<"\n";
		  total_te += te;
		  twr_count_FEMC += 1;    
		  te_aggregate += te;  
		}
	      std::cout<<"te_aggregate is: "<<te_aggregate<<"\n\n";
	    }
	  std::cout<<"twr_count_FHCAL is: "<<twr_count_FHCAL<<"\n";
	  std::cout<<"twr_count_FEMC is: "<<twr_count_FEMC<<"\n";

	  for (int j=0; j<evaltree1->get_nclusters(); j++) // changed iteration varaible from i to j 
	    {
	      std::cout<<"\nFHCAL Cluster: "<<j<<"\n";

	      EvalCluster *clus1 = evaltree1->get_cluster(j);
	      if (clus1)
		{
		  cout<<"non-empty FHCAL cluster\n";
		  Double_t ce = clus1->get_ce();
		  Double_t ctowers = clus1->get_ctowers();
		  std::cout<<"FHCAL ce is: "<<ce<<"\n";
		  std::cout<<"FHCAL ctowers is: "<<ctowers<<"\n\n";
		  total_ce += ce;
		  clus_count_FHCAL += 1;    
		  ce_aggregate += ce;  
		}
	      std::cout<<"ce_aggregate is: "<<ce_aggregate<<"\n\n";
	    }


	  for (int j=0; j<evaltree2->get_nclusters(); j++) // changed iteration varaible from i to j 
	    {
	      std::cout<<"\nFEMC Cluster: "<<j<<"\n";
	      EvalCluster *clus2 = evaltree2->get_cluster(j);
	      if (clus2)
		{
		  cout<<"non-empty FEMC cluster\n";
		  Double_t ce = clus2->get_ce();
		  Double_t ctowers = clus2->get_ctowers();
		  std::cout<<"FEMC ce is: "<<ce<<"\n";
		  std::cout<<"FEMC ctowers is: "<<ctowers<<"\n\n";
		  total_ce += ce;
		  clus_count_FEMC += 1;    
		  ce_aggregate += ce;  
		}
	      std::cout<<"ce_aggregate is: "<<ce_aggregate<<"\n\n";
	    }
	  std::cout<<"clus_count_FHCAL is: "<<clus_count_FHCAL<<"\n\n";
	  std::cout<<"clus_count_FEMC is: "<<clus_count_FEMC<<"\n\n";            
	}
      std::cout<<"\n total ce = "<<total_ce<<"\n";
      std::cout<<"\n total te = "<<total_te<<"\n\n\n";
    }

  
  std::cout<<"total_ge is: "<<total_ge<<"\n";
  std::cout<<"total_ce is: "<<total_ce<<"\n";
  std::cout<<"total_te is: "<<total_te<<"\n";
  std::cout<<"\n\nDone\n----------------------------------------------------------------------\n\n";
}
