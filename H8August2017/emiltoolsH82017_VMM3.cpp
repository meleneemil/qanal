#include <vector>
#include <iostream>
#include <string>
#include <map>
using namespace std;
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TH1D.h"
#include <TStyle.h>
#include <TCanvas.h>
#include "TLatex.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TMarker.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TBranch.h"
#include "TTree.h"
#include "TPad.h"
#include "TChain.h"
#include "TBox.h"
#include "TMinuit.h"
#include "TSpectrum.h"
#include "TSystem.h"
#include "Riostream.h"
#include "TROOT.h"
#include "TProfile.h"
#include "TLine.h"
#include "TFitResult.h"
#include "TAxis.h"

//+++++++++++include vmm trees here
#include "vmm.C"
vmm *vmm_tree;
//++++++++++++++++++++++++++++++++++++++

//clusters *cl_tree;//this is the object from which we get the tree
TString s_runnumber;//here we will store the run number from the filename.

void load_tree_objects(TString file)
{
    //Get file
    TFile *root_file = new TFile(file.Data(),"READ");
    //Load all cluster trees from file
    vmm_tree = new vmm((TTree*)root_file->Get("vmm"));
}

TString get_run_number_from_file(TString file) {

    int index=file.Index("run");
    TString s_runnumber = file(index,5);
    return s_runnumber;
}
int findMax(double *a, int size)
{
    int max=0;

    for(int i=1;i<size;i++)
    {
        if(a[i]>a[max])
            max=i;
    }
    return max;
}
int findMin(double *a, int size)
{
    int min=0;

    for(int i=1;i<size;i++)
    {
        if(a[i]<a[min])
            min=i;
    }
    return min;
}
void printPercentage(int index, int total)
{
    if(index%(total/10)==0)
        cout << "\r"<<(int)(((double)index)/((double)total)*100.0+1) << "% of "<<total<<flush;
}
void reset(vector<double> &vector)
{
    for(uint i=0;i<vector.size();i++)
    {
        vector.at(i)=0.0;
    }
}
void reset(vector<bool> &vector)
{
    for(uint i=0;i<vector.size();i++)
    {
        vector.at(i)=false;
    }
}
int indexOf(vector<TString> vector,TString ch)
{

    for(uint i=0;i<vector.size();i++)
    {
        if(vector.at(i)==ch)
            return i;
    }
}
bool is2Dchamber(TString chamber)
{
    if(chamber.Contains("Tmm"))
        return true;
    else if(chamber.Contains("NTUA_MM"))
        return true;
    if(chamber.Contains("Tmm"))
        return true;

    return false;
}

//void histZoomFit(vector<vector<TH1*>*> histo)
//{
//    for(int i=0;i<histo->size();i++)
//        for(int j=0;j<histo->at(i)->size();j++)
//        {
//            histo->
//                    at(i)->
//                    at(j)->
//                    GetXaxis()->
//                    SetRange(
//                        histo->
//                        at(i)->
//                        at(j)->FindFirstBinAbove(),
//                        histo->
//                        at(i)->
//                        at(j)->FindLastBinAbove());
//        }
//}

void histZoomFit(TH1D* histo)
{
//    for(int i=0;i<histo->size();i++)
//        for(int j=0;j<histo->at(i)->size();j++)
        {
            histo->
//                    at(i)->
//                    at(j)->
                    GetXaxis()->
                    SetRange(
                        histo->
//                        at(i)->
//                        at(j)->
                        FindFirstBinAbove(),
                        histo->
//                        at(i)->
//                        at(j)->
                        FindLastBinAbove());
        }
}
