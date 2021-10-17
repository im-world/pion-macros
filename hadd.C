//- Used to combine the evaluator files while running jobs in batches

#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include <iostream> 

TList *FileList;
TFile *Target;
 
void MergeRootfile( TDirectory *target, TList *sourcelist );
 
 
void hadd(TString detector) {
  // Prepare the files to me merged

  //std::cout<<"Please enter the name of the detector (all caps): ";
  //std::cin>>detector;

  if(gSystem->AccessPathName("macros/Eval_" + detector + ".root")) {
    gSystem->CopyFile("hsimple.root", "macros1/Eval_" + detector + ".root");
    gSystem->CopyFile("hsimple.root", "macros2/Eval_" + detector + ".root");
    gSystem->CopyFile("hsimple.root", "macros3/Eval_" + detector + ".root");
    gSystem->CopyFile("hsimple.root", "macros4/Eval_" + detector + ".root");
    gSystem->CopyFile("hsimple.root", "macros5/Eval_" + detector + ".root");
    gSystem->CopyFile("hsimple.root", "macros6/Eval_" + detector + ".root");
    gSystem->CopyFile("hsimple.root", "macros7/Eval_" + detector + ".root");
    gSystem->CopyFile("hsimple.root", "macros8/Eval_" + detector + ".root");
    gSystem->CopyFile("hsimple.root", "macros9/Eval_" + detector + ".root");
    gSystem->CopyFile("hsimple.root", "macros10/Eval_" + detector + ".root");
    gSystem->CopyFile("hsimple.root", "macros11/Eval_" + detector + ".root");
    gSystem->CopyFile("hsimple.root", "macros12/Eval_" + detector + ".root");
    gSystem->CopyFile("hsimple.root", "macros13/Eval_" + detector + ".root");
    gSystem->CopyFile("hsimple.root", "macros14/Eval_" + detector + ".root");
  }
 
  // in an interactive ROOT session, edit the file names
  // Target and FileList, then
  // root > .L hadd.C
  // root > hadd()
 
  Target = TFile::Open("merged_Eval_" + detector + ".root", "RECREATE" );
 
  FileList = new TList();
  FileList->Add( TFile::Open("macros/Eval_" + detector + ".root") );
  FileList->Add( TFile::Open("macros1/Eval_" + detector + ".root") );
  FileList->Add( TFile::Open("macros2/Eval_" + detector + ".root") );
  FileList->Add( TFile::Open("macros3/Eval_" + detector + ".root") );
  FileList->Add( TFile::Open("macros4/Eval_" + detector + ".root") );
  FileList->Add( TFile::Open("macros5/Eval_" + detector + ".root") );
  FileList->Add( TFile::Open("macros6/Eval_" + detector + ".root") );
  FileList->Add( TFile::Open("macros7/Eval_" + detector + ".root") );
  FileList->Add( TFile::Open("macros8/Eval_" + detector + ".root") );
  FileList->Add( TFile::Open("macros9/Eval_" + detector + ".root") );
  FileList->Add( TFile::Open("macros10/Eval_" + detector + ".root") );
  FileList->Add( TFile::Open("macros11/Eval_" + detector + ".root") );
  FileList->Add( TFile::Open("macros12/Eval_" + detector + ".root") );
  FileList->Add( TFile::Open("macros13/Eval_" + detector + ".root") );
  FileList->Add( TFile::Open("macros14/Eval_" + detector + ".root") );


  MergeRootfile( Target, FileList );
 
}
 
void MergeRootfile( TDirectory *target, TList *sourcelist ) {
 
  //  cout << "Target path: " << target->GetPath() << endl;
  TString path( (char*)strstr( target->GetPath(), ":" ) );
  path.Remove( 0, 2 );
 
  TFile *first_source = (TFile*)sourcelist->First();
  first_source->cd( path );
  TDirectory *current_sourcedir = gDirectory;
  //gain time, do not add the objects in the list in memory
  Bool_t status = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
 
  // loop over all keys in this directory
  TChain *globChain = 0;
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key, *oldkey=0;
  while ( (key = (TKey*)nextkey())) {
 
    //keep only the highest cycle number for each key
    if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;
 
    // read object from first source file
    first_source->cd( path );
    TObject *obj = key->ReadObj();
 
    if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
      // descendant of TH1 -> merge it
 
      //      cout << "Merging histogram " << obj->GetName() << endl;
      TH1 *h1 = (TH1*)obj;
 
      // loop over all source files and add the content of the
      // correspondant histogram to the one pointed to by "h1"
      TFile *nextsource = (TFile*)sourcelist->After( first_source );
      while ( nextsource ) {
 
	// make sure we are at the correct directory level by cd'ing to path
	nextsource->cd( path );
	TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(h1->GetName());
	if (key2) {
	  TH1 *h2 = (TH1*)key2->ReadObj();
	  h1->Add( h2 );
	  delete h2;
	}
 
	nextsource = (TFile*)sourcelist->After( nextsource );
      }
    }
    else if ( obj->IsA()->InheritsFrom( TTree::Class() ) ) {
 
      // loop over all source files create a chain of Trees "globChain"
      const char* obj_name= obj->GetName();
 
      globChain = new TChain(obj_name);
      globChain->SetMaxTreeSize(200000000000LL); // set size to ~100 GB
      globChain->Add(first_source->GetName());
      TFile *nextsource = (TFile*)sourcelist->After( first_source );
      //      const char* file_name = nextsource->GetName();
      // cout << "file name  " << file_name << endl;
      while ( nextsource ) {
 
	globChain->Add(nextsource->GetName());
	nextsource = (TFile*)sourcelist->After( nextsource );
      }
 
    } else if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) ) {
      // it's a subdirectory
 
      cout << "Found subdirectory " << obj->GetName() << endl;
 
      // create a new subdir of same name and title in the target file
      target->cd();
      TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );
 
      // newdir is now the starting point of another round of merging
      // newdir still knows its depth within the target file via
      // GetPath(), so we can still figure out where we are in the recursion
      MergeRootfile( newdir, sourcelist );
 
    } else {
 
      // object is of no type that we know or can handle
      cout << "Unknown object type, name: "
	   << obj->GetName() << " title: " << obj->GetTitle() << endl;
    }
 
    // now write the merged histogram (which is "in" obj) to the target file
    // note that this will just store obj in the current directory level,
    // which is not persistent until the complete directory itself is stored
    // by "target->Write()" below
    if ( obj ) {
      target->cd();
 
      //!!if the object is a tree, it is stored in globChain...
      if(obj->IsA()->InheritsFrom( TTree::Class() ))
	globChain->Merge(target->GetFile(),0,"keep");
      else
	obj->Write( key->GetName() );
    }
 
  } // while ( ( TKey *key = (TKey*)nextkey() ) )
 
  // save modifications to target file
  target->SaveSelf(kTRUE);
  TH1::AddDirectory(status);
}
