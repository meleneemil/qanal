//basic C++ commands
#include <vector>
#include <cstring>
#include <algorithm>
//#include <vector<vector>>
#include <iostream>
using namespace std;
//#include "vmm1_new.C"
// ROOT header
#include "Riostream.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMarker.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLine.h"
#include "TBranch.h"
#include "TTree.h"
#include "TPad.h"
#include "TChain.h"
#include "TSystem.h"
#include "vmm.C"

#include "emiltools.cpp"
//TFile *root_file;
//vmm2 *cl_tree;

//MMegas Clusterization script for the July 2016 configuration for 2 MMegas chambers with 8 VMM2 chips per MMFE8
//reading only the xstrips from T1, T2, T3, T4, T5, T6, T7, T8


//variables: 
Int_t strips;
Int_t iEn, nEntries;
Int_t cluster;
Int_t oversizeCluster=0;
const int clthr=30;
//int thresh = 0;
Int_t tmp1,tmp2,flag;
Int_t clusterSize[clthr];
Int_t strOverThr=0;
Float_t clustCharge;
Float_t clustTime;
Float_t clustPos;


int sumClust;

//vector <double> clusteringChamber;
//vector <double> clusteringChargeChamber;
//vector <double> clusteringTimeChamber;
//vector <double> clusteringPosChamber;

vector <double> clusteringMMFE8;
vector <double> clusteringChargeMMFE8;
vector <double> clusteringTimeMMFE8;
vector <double> clusteringPosMMFE8;

//vector <double> clusteringChipVMM2;
//vector <double> clusteringChargeChipVMM2;
//vector <double> clusteringTimeChipVMM2;
//vector <double> clusteringPosChipVMM2;

//vectors for separating the different chambers
int nT1, nT2, nT3, nT4, nT5, nT6, nT7, nT8;
int nStrips;
//vector for separating the different MMFE8:
int MMFE8_1, MMFE8_2, MMFE8_3, MMFE8_4, MMFE8_5, MMFE8_6, MMFE8_7, MMFE8_8;
//vector for separating the different VMM2 of each MMFE8:
int VMM2_1, VMM2_2, VMM2_3, VMM2_4, VMM2_5, VMM2_6, VMM2_7, VMM2_8;  

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||";
#define PBWIDTH 60

void printProgress (double percentage)
{
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    //printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
}

bool saturatedEvent(vector<float> *charge)
{
    bool isSaturated = false;
    for(uint i=0; i<charge->size(); i++)
    {
        cout << "saturatedEvent::charge size: " << charge->size() << endl;
        if(charge->at(i)>=1){
            cout << "saturatedEvent::charge->at("<< i << "): " << charge->at(i) << endl;
            isSaturated =true;
            return 	isSaturated;
        }
    }

}

//clearing the variables and vectors used for the temporary information
void initialize()	{	

    cluster=0;
    flag=0;
    tmp1=2000;
    sumClust=0;
    strOverThr=0;
    clustCharge=0;

    //----CHAMBER:
    //	   clusteringChamber.clear();
    // clusteringChargeChamber.clear();
    //    clusteringPosChamber.clear();
    //   clusteringTimeChamber.clear();
    //----BOARD MMFE8:
    clusteringMMFE8.clear();
    clusteringChargeMMFE8.clear();
    clusteringPosMMFE8.clear();
    clusteringTimeMMFE8.clear();
    //----VMM2 CHIP:
    //	  clusteringChipVMM2.clear();
    //clusteringChargeChipVMM2.clear();
    //   clusteringPosChipVMM2.clear();
    //  clusteringTimeChipVMM2.clear();
    //cout<<"Initializing ok"<<endl;
}	

//clusterization algorithm: Weighted cluster position & cluster charge as the sum of strip charges
void clusterize(int counter, vector<double> charge_ptr, vector<double> time_ptr, vector<double> hitLoc_ptr, int nStrip, float thresh)	{	

    thresh=0.00000000000000000001;
    //cout << "charge_ptr.at(" << counter << ") > thresh = " << charge_ptr.at(counter) << " > " << thresh << endl;

    if(charge_ptr.at(counter)>thresh) {
        strOverThr++;
        tmp2=hitLoc_ptr.at(counter);
        //cout << "tmp1 = " << tmp1 << endl;
        //cout << "tmp2 = hitLoc_ptr.at(" << counter << ") = " << hitLoc_ptr.at(counter) << endl;
        if(cluster==0) {cluster=1;flag=1;}
        if ( (tmp2-tmp1)==1 ) cluster++;
        if (flag==1 && counter==(nStrip-1) ) {
            //std::cout<<"stop cluster end of array"<<std::endl;
            clustCharge=0;clustPos=0;clustTime=0;
            for(Int_t clc=0;clc<cluster;clc++){
                //cout << "cluster = " << cluster << endl;
                //cout << "counter+1-cluster+clc = " << counter+1-cluster+clc << endl;
                clustCharge+=charge_ptr.at(counter+1-cluster+clc);
                //cout << "clustCharge = clustCharge + charge_ptr.at(counter+1-cluster+clc) = " << clustCharge << endl;
                clustPos+=charge_ptr.at(counter+1-cluster+clc)*hitLoc_ptr.at(counter+1-cluster+clc);
                clustTime+=charge_ptr.at(counter+1-cluster+clc)*time_ptr.at(counter+1-cluster+clc);
            }

            if(clustCharge>thresh){
                clustPos/=clustCharge;
                clustTime/=clustCharge;
                //cout << "clustPos = clustPos / clustCharge = " << clustPos << endl;
                //cout << "clustTime = clustTime / clustCharge = " << clustTime << endl;
            } //cout << "\n" << clustPos;

            //clusteringChargeChamber.push_back(clustCharge);
            //      clusteringChamber.push_back(cluster);
            //   clusteringPosChamber.push_back(clustPos);
            //  clusteringTimeChamber.push_back(clustTime);
            //				std::cout<<"In cluster 0"<<std::endl;

            clusteringChargeMMFE8.push_back(clustCharge);
            clusteringMMFE8.push_back(cluster);
            clusteringPosMMFE8.push_back(clustPos);
            clusteringTimeMMFE8.push_back(clustTime);

            //clusteringChargeChipVMM2.push_back(clustCharge);
            //      clusteringChipVMM2.push_back(cluster);
            //   clusteringPosChipVMM2.push_back(clustPos);
            //  clusteringTimeChipVMM2.push_back(clustTime);

            flag=0;
            sumClust++;

            if (clthr > cluster) {clusterSize[cluster]++;}
            else if (clthr < cluster) oversizeCluster++;
            cluster=0;
        }

        else if ( (tmp2-tmp1>1) && flag==1 && counter<(nStrip-1) ) {

            //std::cout<<"stop cluster"<<std::endl;
            clustCharge=0;clustPos=0;clustTime=0;
            for(Int_t clc=0;clc<cluster;clc++){
                //std::cout<<counter-cluster+clc<<std::endl;
                clustCharge+=charge_ptr.at(counter-cluster+clc);
                clustPos+=charge_ptr.at(counter-cluster+clc)*hitLoc_ptr.at(counter-cluster+clc);
                clustTime+=charge_ptr.at(counter-cluster+clc)*time_ptr.at(counter-cluster+clc);
            }
            if(clustCharge>0.000000001) {clustPos/=clustCharge; clustTime/=clustCharge;} //cout << "\n" << clustPos;

            //clusteringChargeChamber.push_back(clustCharge);
            //clusteringChamber.push_back(cluster);
            //clusteringPosChamber.push_back(clustPos);
            //clusteringTimeChamber.push_back(clustTime);

            clusteringChargeMMFE8.push_back(clustCharge);
            clusteringMMFE8.push_back(cluster);
            clusteringPosMMFE8.push_back(clustPos);
            clusteringTimeMMFE8.push_back(clustTime);

            flag=1;
            sumClust++;

            if (clthr > cluster) {clusterSize[cluster]++;}
            else if (clthr < cluster) oversizeCluster++;
            cluster=0;
        }
        else if ( (tmp2-tmp1>1) && flag==1 && counter==(nStrip-1) ) {
            //				std::cout<<"stop cluster for non consecutive strips"<<std::endl;
            clustCharge=0;clustPos=0;clustTime=0;
            for(Int_t clc=0;clc<cluster;clc++){
                clustCharge+=charge_ptr.at(counter-cluster+clc);
                clustPos+=charge_ptr.at(counter-cluster+clc)*hitLoc_ptr.at(counter-cluster+clc);
                clustTime+=charge_ptr.at(counter-cluster+clc)*time_ptr.at(counter-cluster+clc);
            }
            if(clustCharge>0.000000001) {clustPos/=clustCharge;clustTime/=clustCharge;} //cout << "\n" << clustPos;

            //clusteringChargeChamber.push_back(clustCharge);
            //clusteringChamber.push_back(cluster);
            //clusteringPosChamber.push_back(clustPos);
            //clusteringTimeChamber.push_back(clustTime);

            clusteringChargeMMFE8.push_back(clustCharge);
            clusteringMMFE8.push_back(cluster);
            clusteringPosMMFE8.push_back(clustPos);
            clusteringTimeMMFE8.push_back(clustTime);

            flag=1;
            sumClust++;
            if (clthr > cluster) {clusterSize[cluster]++;}
            else if (clthr < cluster) oversizeCluster++;
            cluster=1;
            clustCharge=charge_ptr.at(counter);
            clustPos=hitLoc_ptr.at(counter);
            //clustTime=hitLoc_ptr.at(counter);
            clustTime=time_ptr.at(counter);

            //clusteringChargeChamber.push_back(clustCharge);
            //clusteringChamber.push_back(cluster);
            //clusteringPosChamber.push_back(clustPos);
            //clusteringTimeChamber.push_back(clustTime);

            clusteringChargeMMFE8.push_back(clustCharge);
            clusteringMMFE8.push_back(cluster);
            clusteringPosMMFE8.push_back(clustPos);
            clusteringTimeMMFE8.push_back(clustTime);

            cluster=0;
            flag=0;
        }
        tmp1=hitLoc_ptr.at(counter);
    }

    else if(charge_ptr.at(counter)<=thresh && flag==1){
        //std::cout<<"stop cluster for low charge strips"<<std::endl;
        clustCharge=0;clustPos=0;clustTime=0;
        for(Int_t clc=0;clc<cluster;clc++){
            clustCharge+=charge_ptr.at(counter-cluster+clc);
            clustPos+=charge_ptr.at(counter-cluster+clc)*hitLoc_ptr.at(counter-cluster+clc);
            clustTime+=charge_ptr.at(counter-cluster+clc)*time_ptr.at(counter-cluster+clc);
        }
        if(clustCharge!=0) {clustPos/=clustCharge; clustTime/=clustCharge;} //cout << "\n" << clustPos;

        //clusteringChargeChamber.push_back(clustCharge);
        //clusteringChamber.push_back(cluster);
        //clusteringPosChamber.push_back(clustPos);
        //clusteringTimeChamber.push_back(clustTime);

        clusteringChargeMMFE8.push_back(clustCharge);
        clusteringMMFE8.push_back(cluster);
        clusteringPosMMFE8.push_back(clustPos);
        clusteringTimeMMFE8.push_back(clustTime);

        flag=0;
        sumClust++;
        if (clthr > cluster) {clusterSize[cluster]++;}
        else if (clthr < cluster) oversizeCluster++; cluster=0;
    }

    else if(charge_ptr.at(counter) > thresh && flag==1 && counter==nStrip-1){
        //		    std::cout<<"stop cluster"<<std::endl;
        clustCharge=0;clustPos=0;clustTime=0;
        for(Int_t clc=0;clc<cluster;clc++){
            clustCharge+=charge_ptr.at(counter-cluster+1+clc);
            clustPos+=charge_ptr.at(counter-cluster+1+clc)*hitLoc_ptr.at(counter-cluster+1+clc);
            clustTime+=charge_ptr.at(counter-cluster+clc)*time_ptr.at(counter-cluster+clc);
        }
        if(clustCharge!=0) {clustPos/=clustCharge; clustTime/=clustCharge;} //cout << "\n" << clustPos;

        //clusteringChargeChamber.push_back(clustCharge);
        //clusteringChamber.push_back(cluster);
        //clusteringPosChamber.push_back(clustPos);
        //clusteringTimeChamber.push_back(clustTime);

        clusteringChargeMMFE8.push_back(clustCharge);
        clusteringMMFE8.push_back(cluster);
        clusteringPosMMFE8.push_back(clustPos);
        clusteringTimeMMFE8.push_back(clustTime);

        sumClust++;if (clthr > cluster) {clusterSize[cluster]++;}
        else if (clthr < cluster) oversizeCluster++;cluster=0;

    }
}

//loading the file and the trees: RAW for strip#, mm_id, apv_id and DATA for apv_qqmax
void process(TString file, float thresh){	

    //cout << "Line = " << __LINE__ << endl;

    //vectors for the clustering algorithm
    int clEvent=0;
    //vector< TString > *clChmbId=0;
    vector< Int_t > *clChmbId=0;
    vector< Int_t > *clBoardId=0;
    vector< Int_t > *clChipId=0;
    vector< Int_t > *clNclu=0;
    vector< vector<double> > *clSize=0;
    vector< vector<double> > *clCharge=0;
    vector< vector<double> > *clTime=0;
    vector< vector<double> > *clPos=0;
    vector< vector<double> > *clStripID=0;
    vector< vector<double> > *clStripPdo=0;
    vector< vector<double> > *clStripTdo=0;
    //vector< vector<double> > *boardID=0;
    //vector< vector<double> > *chipID=0;

    //cout << "Line = " << __LINE__ << endl;

    //Also open the file as update-able to write in it the new cluster variables (efficient???)

    //TString fWriteName = file+"_clu";

    Int_t           eventFAFA;
    //vector<TString> *chmbID=0;
    vector<Int_t> *boardId=0;
    vector<Int_t> *chip=0;
    //vector<vector<double> > *stripID=0;
    vector<vector<double> > *channel=0;
    vector<vector<double> > *tdo=0;
    vector<vector<double> > *pdo=0;

    // List of branches
    TBranch        *b_eventFAFA;   //!
    //TBranch        *b_chmbID;   //!
    TBranch        *b_boardId;   //!
    TBranch		   *b_chip; //!
    //TBranch        *b_stripID;   //!
    TBranch        *b_channel;   //!
    TBranch        *b_tdo;   //!
    TBranch        *b_pdo;   //!

    //cout << "Line = " << __LINE__ << endl;
    //TFile *root_file = new TFile(file.Data(),"READ");
    TFile *fin = new TFile(file.Data(),"READ");

    //cout << "Line = " << __LINE__ << endl;
    //TTree *tr = (TTree*) fin->Get("vmm1_new");
    TTree *tr = (TTree*) fin->Get("vmm");
    //TTree *tr = (TTree*) root_file->Get("vmm2");

    //cout << "Line = " << __LINE__ << endl;
    //    tr->Print();
    //cout << "Line = " << __LINE__ << endl;
    tr->SetBranchStatus("*",1);
    //cout << "Line = " << __LINE__ << endl;
    //############################################################################## comenting out here
    tr->SetBranchAddress("eventFAFA", &eventFAFA, &b_eventFAFA);
    //cout << "Line = " << __LINE__ << endl;
    //tr->SetBranchAddress("chmbID", &chmbID, &b_chmbID);
    //##############################################################################
    tr->SetBranchAddress("boardId", &boardId, &b_boardId);
    //cout << "Line = " << __LINE__ << endl;
    tr->SetBranchAddress("chip", &chip, &b_chip);
    //cout << "Line = " << __LINE__ << endl;
    tr->SetBranchAddress("channel", &channel, &b_channel);
    //cout << "Line = " << __LINE__ << endl;
    //tr->SetBranchAddress("stripID", &stripID, &b_stripID);
    tr->SetBranchAddress("tdo", &tdo, &b_tdo);
    //cout << "Line = " << __LINE__ << endl;
    tr->SetBranchAddress("pdo", &pdo, &b_pdo);

    //cout << "Line = " << __LINE__ << endl;
    //// List of branches
    //TBranch *b_event;   //!
    //TBranch *b_eventCnt;   //!
    //TBranch *b_eventFAFA;   //!
    //TBranch *b_triggerTimestamp;   //!
    //TBranch *b_triggerCounter;   //!
    //TBranch *b_fec;   //!
    //TBranch *b_board;   //!
    //TBranch *b_chip;   //!
    //TBranch *b_eventSize;   //!
    //TBranch *b_channel;   //!
    //TBranch *b_flag;   //!
    //TBranch *b_threshold;   //!
    //TBranch *b_pdo;   //!
    //TBranch *b_tdo;   //!
    //TBranch *b_bcid;   //!
    //TBranch *b_grayDecoded;   //!
    //TBranch *b_ARTHDMIChipID;   //!
    //TBranch *b_ARTTimeStamp;   //!
    //TBranch *b_ARTChannel;   //!
    //TBranch *b_ARTFlag;   //!

    //cout << "Line = " << __LINE__ << endl;
    TString foutname = "processed3_"+file;
    TFile *fout = new TFile(foutname.Data(),"RECREATE");

    //TTree *tclone = tr->CloneTree();
    TTree *tout = new TTree("vmm3_clu","vmm3_clu");

    tout->Branch("cl_stripID",&clStripID);
    tout->Branch("cl_stripPdo",&clStripPdo);
    tout->Branch("cl_stripTdo",&clStripTdo);
    tout->Branch("cl_event", &clEvent);
    tout->Branch("cl_ChmbId",&clChmbId);
    tout->Branch("cl_BoardId",&clBoardId);
    tout->Branch("cl_ChipId",&clChipId);
    tout->Branch("cl_nclu", &clNclu);
    tout->Branch("cl_clsize", &clSize);
    tout->Branch("cl_clcharge", &clCharge);
    tout->Branch("cl_clpos_strips", &clPos);
    tout->Branch("cl_cltime", &clTime);

    //tr->LoadTree(0);
    nEntries = tr->GetEntries();
    cout << nEntries << endl;


    //cout << "Line = " << __LINE__ << endl;
    for (Long64_t iEntry = 0; iEntry < nEntries; iEntry++) {	//loop on tree entries

        //cout << "Line = " << __LINE__ <<" entry = "<<iEntry<< endl;
        //change these for the new trees
        tr->GetEntry(iEntry);
        //std::cout<<"Entry : "<<iEntry<<std::endl;
        clEvent = eventFAFA;

        printPercentage(iEntry,nEntries);//in 10s of %

        //for(uint iChamber = 0; iChamber < chmbID->size(); iChamber++)	{
        for(uint iBoard = 0; iBoard < boardId->size(); iBoard++)	{

            //cout << "Line = " << __LINE__ << endl;
            initialize();	//intializing
            //cout << "Line = " << __LINE__ << endl;
            //clChmbId->push_back(chmbID->at(iChamber));
            //cout << "Line = " << __LINE__ << endl;
            clBoardId->push_back(boardId->at(iBoard));
            //clStripID->push_back(stripID->at(iChamber));
            clChipId->push_back(chip->at(iBoard));
            //for(uint iChan=0; iChan<channel->at(iBoard).size(); iChan++){
            //	clStripID->push_back(channel->at(iBoard).at(iChan));
            //}

            clStripID->push_back(channel->at(iBoard));
            clStripPdo->push_back(pdo->at(iBoard));
            clStripTdo->push_back(tdo->at(iBoard));
            //nStrips = (stripID->at(iChamber)).size();

            nStrips = (channel->at(iBoard)).size();
            //cout << "nStrips = (channel->at(iBoard)).size() = " << nStrips << endl;
            if(nStrips==0){ continue; }

            for(Int_t i=0;i<nStrips;i++){ //loop on strips
                //clusterize(i,pdo->at(iChamber), tdo->at(iChamber), stripID->at(iChamber), nStrips, thresh);
                //cout << "* nStrips = " << nStrips << endl;
                clusterize(i,pdo->at(iBoard), tdo->at(iBoard), channel->at(iBoard), nStrips, thresh);

            }//end of loop on strips

            //clNclu->push_back(sumClust);
            //clSize->push_back(clusteringChamber);
            //clPos->push_back(clusteringPosChamber);
            //clCharge->push_back(clusteringChargeChamber);
            //clTime->push_back(clusteringTimeChamber);
            clNclu->push_back(sumClust);
            clSize->push_back(clusteringMMFE8);
            clPos->push_back(clusteringPosMMFE8);
            clCharge->push_back(clusteringChargeMMFE8);
            clTime->push_back(clusteringTimeMMFE8);

        }//end of loop on Boards

        tout->Fill();
        clChmbId->clear();
        clBoardId->clear();
        clChipId->clear();
        clStripID->clear();
        clStripPdo->clear();
        clStripTdo->clear();
        clNclu->clear();
        clSize->clear();
        clPos->clear();
        clCharge->clear();
        clTime->clear();
    }

    fout->cd();
    //tclone->Write("", TObject::kOverwrite);
    tout->Write("", TObject::kOverwrite);
    fout->Write();

    delete tr;
    delete tout;

    fout->Close();
    fin->Close();
    //root_file->Close();

}


