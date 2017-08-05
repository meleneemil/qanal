//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jun 16 10:05:28 2017 by ROOT version 5.34/36
// from TTree vmm/vmm
// found on file: run_0035.root
//////////////////////////////////////////////////////////

#ifndef vmm_h
#define vmm_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
//#include "/Volumes/MacintoshHD2/build/electric_commander/build/workspace/root-release-5.34/BUILDTYPE/Release/COMPILER/native/LABEL/mac1011/sources/root_v5.34.36/root/cint/cint/lib/prec_stl/vector"

// Fixed size dimensions of array or collections stored in the TTree if any.

class vmm {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           eventFAFA;
   vector<int>     *triggerTimeStamp;
   vector<int>     *triggerCounter;
   vector<int>     *boardId;
   vector<int>     *chip;
   vector<int>     *eventSize;
   vector<vector<int> > *tdo;
   vector<vector<int> > *pdo;
   vector<vector<int> > *flag;
   vector<vector<int> > *threshold;
   vector<vector<int> > *bcid;
   vector<vector<int> > *relbcid;
   vector<vector<int> > *overflow;
   vector<vector<int> > *orbitCount;
   vector<vector<int> > *grayDecoded;
   vector<vector<int> > *channel;
   vector<vector<int> > *febChannel;
   vector<vector<int> > *mappedChannel;
   vector<int>     *art_valid;
   vector<int>     *art;

   // List of branches
   TBranch        *b_eventFAFA;   //!
   TBranch        *b_triggerTimeStamp;   //!
   TBranch        *b_triggerCounter;   //!
   TBranch        *b_boardId;   //!
   TBranch        *b_chip;   //!
   TBranch        *b_eventSize;   //!
   TBranch        *b_tdo;   //!
   TBranch        *b_pdo;   //!
   TBranch        *b_flag;   //!
   TBranch        *b_threshold;   //!
   TBranch        *b_bcid;   //!
   TBranch        *b_relbcid;   //!
   TBranch        *b_overflow;   //!
   TBranch        *b_orbitCount;   //!
   TBranch        *b_grayDecoded;   //!
   TBranch        *b_channel;   //!
   TBranch        *b_febChannel;   //!
   TBranch        *b_mappedChannel;   //!
   TBranch        *b_art_valid;   //!
   TBranch        *b_art;   //!

   vmm(TTree *tree=0);
   virtual ~vmm();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef vmm_cxx
vmm::vmm(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("run_0035.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("run_0035.root");
      }
      f->GetObject("vmm",tree);

   }
   Init(tree);
}

vmm::~vmm()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t vmm::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t vmm::LoadTree(Long64_t entry)
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

void vmm::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   triggerTimeStamp = 0;
   triggerCounter = 0;
   boardId = 0;
   chip = 0;
   eventSize = 0;
   tdo = 0;
   pdo = 0;
   flag = 0;
   threshold = 0;
   bcid = 0;
   relbcid = 0;
   overflow = 0;
   orbitCount = 0;
   grayDecoded = 0;
   channel = 0;
   febChannel = 0;
   mappedChannel = 0;
   art_valid = 0;
   art = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventFAFA", &eventFAFA, &b_eventFAFA);
   fChain->SetBranchAddress("triggerTimeStamp", &triggerTimeStamp, &b_triggerTimeStamp);
   fChain->SetBranchAddress("triggerCounter", &triggerCounter, &b_triggerCounter);
   fChain->SetBranchAddress("boardId", &boardId, &b_boardId);
   fChain->SetBranchAddress("chip", &chip, &b_chip);
   fChain->SetBranchAddress("eventSize", &eventSize, &b_eventSize);
   fChain->SetBranchAddress("tdo", &tdo, &b_tdo);
   fChain->SetBranchAddress("pdo", &pdo, &b_pdo);
   fChain->SetBranchAddress("flag", &flag, &b_flag);
   fChain->SetBranchAddress("threshold", &threshold, &b_threshold);
   fChain->SetBranchAddress("bcid", &bcid, &b_bcid);
   fChain->SetBranchAddress("relbcid", &relbcid, &b_relbcid);
   fChain->SetBranchAddress("overflow", &overflow, &b_overflow);
   fChain->SetBranchAddress("orbitCount", &orbitCount, &b_orbitCount);
   fChain->SetBranchAddress("grayDecoded", &grayDecoded, &b_grayDecoded);
   fChain->SetBranchAddress("channel", &channel, &b_channel);
   fChain->SetBranchAddress("febChannel", &febChannel, &b_febChannel);
   fChain->SetBranchAddress("mappedChannel", &mappedChannel, &b_mappedChannel);
   fChain->SetBranchAddress("art_valid", &art_valid, &b_art_valid);
   fChain->SetBranchAddress("art", &art, &b_art);
   Notify();
}

Bool_t vmm::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void vmm::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t vmm::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef vmm_cxx
