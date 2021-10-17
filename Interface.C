// Interface for all other macros

/*
  authors - Sagar Joshi      (ee190002054@iiti.ac.in)
            Siddhant Rathi   (me190003061@iiti.ac.in)
	    
  version - 1.0

*/

#include <iostream>

#include <LoopEvalMultiPortableChi2.C>
#include <EnergyEvalEtaCut.C>
#include <LoopEvalPortableCircularCut.C>

#include <LoopEvalMultiF.C>
#include <LoopEvalMultiFCircularCut.C>
#include <LoopEvalMultiH.C>
#include <LoopEvalMultiHCircularCut.C>


void Interface(){

  int combine = 0;
  int detector = 0;
  int circularCut = 0;
  int print = 0;

  std::cout<<R"(
 ███████████                       █████ █████    █████████   ████  ████              ██████      █████████  
░░███░░░░░░█                      ░░███ ░░███    ███░░░░░███ ░░███ ░░███            ███░░░░███   ███░░░░░███ 
 ░███   █ ░  █████ ████ ████████   ░███  ░███ █ ░███    ░███  ░███  ░███           ███    ░░███ ░███    ░███ 
 ░███████   ░░███ ░███ ░░███░░███  ░███████████ ░███████████  ░███  ░███          ░███     ░███ ░███████████ 
 ░███░░░█    ░███ ░███  ░███ ░███  ░░░░░░░███░█ ░███░░░░░███  ░███  ░███          ░███   ██░███ ░███░░░░░███ 
 ░███  ░     ░███ ░███  ░███ ░███        ░███░  ░███    ░███  ░███  ░███          ░░███ ░░████  ░███    ░███ 
 █████       ░░████████ ████ █████       █████  █████   █████ █████ █████          ░░░██████░██ █████   █████
░░░░░         ░░░░░░░░ ░░░░ ░░░░░       ░░░░░  ░░░░░   ░░░░░ ░░░░░ ░░░░░             ░░░░░░ ░░ ░░░░░   ░░░░░ )";

  std::cout<<"\n\nPlease press the button corresponding to the option you\n want to select\n\n";

  std::cout<<"\n\nWould you like to use combined detectors or individual ones?\n\n1 - Combined \n2 - Individual\n\n";

  std::cout<<"Please select an option (1/2): ";
  std::cin>>combine;


  while(combine!=1 && combine !=2){
    std::cout<<"Invalid input";

    std::cout<<"\n\nWould you like to use combined detectors or individual ones?\n\n1 - Combined \n2 - Individual\n\n";

    std::cout<<"Please select an option (1/2): ";
    std::cin>>combine;
  }

  if(combine == 1){

    std::cout<<"\n\nWhich detector combination would you like to use?\n\n1 - FEMC + FHCAL\n2 - CEMC + HCALIN + HCALOUT\n\n";

    std::cout<<"Please select an option (1/2): ";
    std::cin>>detector;

    while(detector != 1 && detector != 2){
      std::cout<<"Invalid input";
      std::cout<<"\n\nWhich detector combination would you like to use?\n\n1 - FEMC + FHCAL\n2 - CEMC + HCALIN + HCALOUT\n\n";

      std::cout<<"Please select an option (1/2): ";
      std::cin>>detector;
    }
  }

  else if(combine == 2){

    std::cout<<"\n\nWhich detector would you like to use?\n\n1 - CEMC\n2 - EEMC\n3 - FEMC\n4 - FHCAL\n5 - HCALIN\n6 - HCALOUT\n\n";

    std::cout<<"Please select an option (1/2/3/4/5/6): ";
    std::cin>>detector;

    while(detector != 1 && detector != 2 && detector != 3 && detector != 4 && detector != 5 && detector != 6){
      std::cout<<"Invalid input";
      std::cout<<"\n\nWhich detector would you like to use?\n\n1 - CEMC\n2 - EEMC\n3 - FEMC\n4 - FHCAL\n5 - HCALIN\n6 - HCALOUT\n\n";

      std::cout<<"Please select an option (1/2/3/4/5/6): ";
      std::cin>>detector;
    }
  }

  std::cout<<"\n\nWould you like employ circular cuts (manual clustering)?\n\n1 - Yes\n2 - No\n\n";

  std::cout<<"Please select an option (1/2): ";
  std::cin>>circularCut;

  while(circularCut != 1 && circularCut != 2){
    std::cout<<"Invalid input";
    std::cout<<"\n\nWould you like employ circular cuts (manual clustering)?\n\n1 - Yes\n2 - No\n\n";

    std::cout<<"Please select an option (1/2): ";
    std::cin>>circularCut;
  }

  std::cout<<"\n\nWould you like to save the generated plots as .png files?\n\n1 - Yes\n2 - No\n\n";

  std::cout<<"Please select an option (1/2): ";
  std::cin>>print;

  while(print != 1 && print != 2){
    std::cout<<"Invalid input";
    std::cout<<"\n\nWould you like to save the generated plots as .png files?\n\n1 - Yes\n2 - No\n\n";

    std::cout<<"Please select an option (1/2): ";
    std::cin>>print;
  }

  

  ////////////////////
  //Processing section
  ////////////////////


  TString detectorName;

  if(detector == 1){
    detectorName = "CEMC";
  }
  
  else if(detector == 2){
    detectorName = "EEMC";
  }

  else if(detector == 3){
    detectorName = "FEMC";
  }

  else if(detector == 4){
    detectorName = "FHCAL";
  }

  else if(detector == 5){
    detectorName = "HCALIN";
  }

  else if(detector == 6){
    detectorName = "HCALOUT";
  }




  if(combine == 1 && detector == 1 && circularCut == 1){
    LoopEvalMultiFCircularCut(print);
  }

  else if(combine == 1 && detector == 1 && circularCut == 2){
    LoopEvalMultiF(print);
  }

  else if(combine == 1 && detector == 2 && circularCut == 1){
    LoopEvalMultiHCircularCut(print);
  }

  else if(combine == 1 && detector == 2 && circularCut == 2){
    LoopEvalMultiH(print);
  }

  else if(combine == 2 && circularCut == 1){
    LoopEvalPortableCircularCut(detectorName, print);
  }

  else if(combine == 2 && circularCut == 2){
    EnergyEvalEtaCut(detectorName, print);
    LoopEvalMultiPortableChi2(detectorName, print);
  }

}
