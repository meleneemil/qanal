//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 18 16:11:16 2017 by ROOT version 5.34/36
// from TTree vmm2/vmm2
// found on file: run_9028.root
//////////////////////////////////////////////////////////

#ifndef vmm2_h
#define vmm2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class vmm2 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           event;
   Int_t           eventCnt;
   Int_t           eventFAFA;
   vector<int>     *triggerTimestamp;
   vector<int>     *triggerCounter;
   vector<int>     *fec;
   vector<int>     *board;
   vector<int>     *chip;
   vector<int>     *eventSize;
   vector<vector<int> > *channel;
   vector<vector<int> > *flag;
   vector<vector<int> > *threshold;
   vector<vector<int> > *pdo;
   vector<vector<int> > *tdo;
   vector<vector<int> > *bcid;
   vector<vector<int> > *grayDecoded;
   vector<int>     *ARTHDMIChipID;
   vector<int>     *ARTTimeStamp;
   vector<vector<int> > *ARTChannel;
   vector<vector<int> > *ARTFlag;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_eventCnt;   //!
   TBranch        *b_eventFAFA;   //!
   TBranch        *b_triggerTimestamp;   //!
   TBranch        *b_triggerCounter;   //!
   TBranch        *b_fec;   //!
   TBranch        *b_board;   //!
   TBranch        *b_chip;   //!
   TBranch        *b_eventSize;   //!
   TBranch        *b_channel;   //!
   TBranch        *b_flag;   //!
   TBranch        *b_threshold;   //!
   TBranch        *b_pdo;   //!
   TBranch        *b_tdo;   //!
   TBranch        *b_bcid;   //!
   TBranch        *b_grayDecoded;   //!
   TBranch        *b_ARTHDMIChipID;   //!
   TBranch        *b_ARTTimeStamp;   //!
   TBranch        *b_ARTChannel;   //!
   TBranch        *b_ARTFlag;   //!

   vmm2(TTree *tree=0);
   virtual ~vmm2();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef vmm2_cxx
vmm2::vmm2(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
//   if (tree == 0) {
//      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("run_9028.root");
//      if (!f || !f->IsOpen()) {
//         f = new TFile("run_9028.root");
//      }
//      f->GetObject("vmm2",tree);

//   }
   Init(tree);
}

vmm2::~vmm2()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t vmm2::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t vmm2::LoadTree(Long64_t entry)
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

void vmm2::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   triggerTimestamp = 0;
   triggerCounter = 0;
   fec = 0;
   board = 0;
   chip = 0;
   eventSize = 0;
   channel = 0;
   flag = 0;
   threshold = 0;
   pdo = 0;
   tdo = 0;
   bcid = 0;
   grayDecoded = 0;
   ARTHDMIChipID = 0;
   ARTTimeStamp = 0;
   ARTChannel = 0;
   ARTFlag = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("eventCnt", &eventCnt, &b_eventCnt);
   fChain->SetBranchAddress("eventFAFA", &eventFAFA, &b_eventFAFA);
   fChain->SetBranchAddress("triggerTimestamp", &triggerTimestamp, &b_triggerTimestamp);
   fChain->SetBranchAddress("triggerCounter", &triggerCounter, &b_triggerCounter);
   fChain->SetBranchAddress("fec", &fec, &b_fec);
   fChain->SetBranchAddress("board", &board, &b_board);
   fChain->SetBranchAddress("chip", &chip, &b_chip);
   fChain->SetBranchAddress("eventSize", &eventSize, &b_eventSize);
   fChain->SetBranchAddress("channel", &channel, &b_channel);
   fChain->SetBranchAddress("flag", &flag, &b_flag);
   fChain->SetBranchAddress("threshold", &threshold, &b_threshold);
   fChain->SetBranchAddress("pdo", &pdo, &b_pdo);
   fChain->SetBranchAddress("tdo", &tdo, &b_tdo);
   fChain->SetBranchAddress("bcid", &bcid, &b_bcid);
   fChain->SetBranchAddress("grayDecoded", &grayDecoded, &b_grayDecoded);
   fChain->SetBranchAddress("ARTHDMIChipID", &ARTHDMIChipID, &b_ARTHDMIChipID);
   fChain->SetBranchAddress("ARTTimeStamp", &ARTTimeStamp, &b_ARTTimeStamp);
   fChain->SetBranchAddress("ARTChannel", &ARTChannel, &b_ARTChannel);
   fChain->SetBranchAddress("ARTFlag", &ARTFlag, &b_ARTFlag);
   Notify();
}

Bool_t vmm2::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void vmm2::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t vmm2::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef vmm2_cxx
