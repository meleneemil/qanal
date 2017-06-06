//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 11 20:48:49 2017 by ROOT version 5.34/36
// from TTree vmm3_clu/vmm3_clu
// found on file: processed3_run_2670_2MDT24_03_2017_trigatpeak.root
//////////////////////////////////////////////////////////

#ifndef vmm3_clu_h
#define vmm3_clu_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class vmm3_clu {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   vector<vector<double> > *cl_stripID;
   vector<vector<double> > *cl_stripPdo;
   vector<vector<double> > *cl_stripTdo;
   Int_t           cl_event;
   vector<int>     *cl_ChmbId;
   vector<int>     *cl_BoardId;
   vector<int>     *cl_ChipId;
   vector<int>     *cl_nclu;
   vector<vector<double> > *cl_clsize;
   vector<vector<double> > *cl_clcharge;
   vector<vector<double> > *cl_clpos_strips;
   vector<vector<double> > *cl_cltime;

   // List of branches
   TBranch        *b_cl_stripID;   //!
   TBranch        *b_cl_stripPdo;   //!
   TBranch        *b_cl_stripTdo;   //!
   TBranch        *b_cl_event;   //!
   TBranch        *b_cl_ChmbId;   //!
   TBranch        *b_cl_BoardId;   //!
   TBranch        *b_cl_ChipId;   //!
   TBranch        *b_cl_nclu;   //!
   TBranch        *b_cl_clsize;   //!
   TBranch        *b_cl_clcharge;   //!
   TBranch        *b_cl_clpos_strips;   //!
   TBranch        *b_cl_cltime;   //!

   vmm3_clu(TTree *tree=0);
   virtual ~vmm3_clu();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef vmm3_clu_cxx
vmm3_clu::vmm3_clu(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("processed3_run_2670_2MDT24_03_2017_trigatpeak.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("processed3_run_2670_2MDT24_03_2017_trigatpeak.root");
      }
      f->GetObject("vmm3_clu",tree);

   }
   Init(tree);
}

vmm3_clu::~vmm3_clu()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t vmm3_clu::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t vmm3_clu::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void vmm3_clu::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   cl_stripID = 0;
   cl_stripPdo = 0;
   cl_stripTdo = 0;
   cl_ChmbId = 0;
   cl_BoardId = 0;
   cl_ChipId = 0;
   cl_nclu = 0;
   cl_clsize = 0;
   cl_clcharge = 0;
   cl_clpos_strips = 0;
   cl_cltime = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("cl_stripID", &cl_stripID, &b_cl_stripID);
   fChain->SetBranchAddress("cl_stripPdo", &cl_stripPdo, &b_cl_stripPdo);
   fChain->SetBranchAddress("cl_stripTdo", &cl_stripTdo, &b_cl_stripTdo);
   fChain->SetBranchAddress("cl_event", &cl_event, &b_cl_event);
   fChain->SetBranchAddress("cl_ChmbId", &cl_ChmbId, &b_cl_ChmbId);
   fChain->SetBranchAddress("cl_BoardId", &cl_BoardId, &b_cl_BoardId);
   fChain->SetBranchAddress("cl_ChipId", &cl_ChipId, &b_cl_ChipId);
   fChain->SetBranchAddress("cl_nclu", &cl_nclu, &b_cl_nclu);
   fChain->SetBranchAddress("cl_clsize", &cl_clsize, &b_cl_clsize);
   fChain->SetBranchAddress("cl_clcharge", &cl_clcharge, &b_cl_clcharge);
   fChain->SetBranchAddress("cl_clpos_strips", &cl_clpos_strips, &b_cl_clpos_strips);
   fChain->SetBranchAddress("cl_cltime", &cl_cltime, &b_cl_cltime);
   Notify();
}

Bool_t vmm3_clu::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void vmm3_clu::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t vmm3_clu::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef vmm3_clu_cxx
