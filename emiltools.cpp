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
#include "TFitResult.h"
#include "TAxis.h"
#include "TPaveStats.h"
#include "TGaxis.h"

//+++++++++++include cluster trees here
//#include "vmm.C"
//++++++++++++++++++++++++++++++++++++++

TString get_run_number_from_file(TString file) {

    //remove ".root"
    int rootIndex=file.Index(".root");
    file=file(0,rootIndex);

    //remove "run" (as many as they are)
    while(file.Index("run") != -1)
        file = file(file.Index("run")+3,999);

    //check for "_"
    if(file(0,1)=="_")
        file=file(1,999);

    return file;
}
TString get_run_number_from_APVfile(TString file) {

    int index=file.Index("run");
    TString s_runnumber = file(index,8);
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

int findMax(int *a, int size)
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

int findMin(int *a, int size)
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
void histZoomFit(TH1* histo)
{

    histo->
            GetXaxis()->
            SetRange(
                histo->FindFirstBinAbove(),
                histo->FindLastBinAbove());
}

void histZoomFit(vector<vector<TH1*>*> *&histo)
{
    for(int i=0;i<histo->size();i++)
        for(int j=0;j<histo->at(i)->size();j++)
        {
            histo->
                    at(i)->
                    at(j)->
                    GetXaxis()->
                    SetRange(
                        histo->
                        at(i)->
                        at(j)->FindFirstBinAbove(),
                        histo->
                        at(i)->
                        at(j)->FindLastBinAbove());
        }
}


void saveCanvasAs(TCanvas *&c, TString extension)
{

    //add the dot if not there
    if(!extension.Contains("."))
        extension="."+extension;

    //set thinner lines for the plot
    gStyle->SetLineScalePS(1.0);

    if(extension==".png")
        c->SetWindowSize(1920,1080);

    c->SaveAs(c->GetName()+extension);
}
void saveCanvasAsAll(TCanvas *&c)
{
    saveCanvasAs(c,"png");
    saveCanvasAs(c,"pdf");
    saveCanvasAs(c,"eps");
}

void transpad(TH1D *&h1,TH1D *&h2, TString cTitle="canvas") {

   TCanvas *c1 = new TCanvas(cTitle,cTitle,200,10,700,500);
   TPad *pad1 = new TPad("pad1","",0,0,1,1);
   TPad *pad2 = new TPad("pad2","",0,0,1,1);

   pad2->SetFillStyle(4000); //will be transparent
   pad1->Draw();
   pad1->cd();

//   TH1F *h1 = new TH1F("h1","h1",100,-3,3);
//   TH1F *h2 = new TH1F("h2","h2",100,-3,3);


//   TRandom r;
//   for (Int_t i=0;i<100000;i++) {
//      Double_t x1 = r.Gaus(-1,0.5);
//      Double_t x2 = r.Gaus(1,1.5);
//      if (i <1000) h1->Fill(x1);
//      h2->Fill(x2);
//   }

   h1->Draw();

   pad1->Update(); //this will force the generation of the "stats" box
   TPaveStats *ps1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
   ps1->SetX1NDC(0.4); ps1->SetX2NDC(0.6);

   pad1->Modified();
   c1->cd();

   //compute the pad range with suitable margins

   Double_t ymin = 0;
   Double_t ymax = (  h1->GetMaximum()>h2->GetMaximum() ? h1->GetMaximum() : h2->GetMaximum()  );
   Double_t dy = (ymax-ymin)/0.8; //10 per cent margins top and bottom
   Double_t xmin = 0;
   Double_t xmax = 650;
   Double_t dx = (xmax-xmin)/0.8; //10 per cent margins left and right
   pad2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);


   pad2->Draw();
   pad2->cd();
   h2->SetLineColor(kRed);
   h2->Draw("][sames");
   pad2->Update();

   TPaveStats *ps2 = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
   ps2->SetX1NDC(0.65); ps2->SetX2NDC(0.85);
   ps2->SetTextColor(kRed);
   // draw axis on the right side of the pad

   TGaxis *axis = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,50510,"+L");
   axis->SetLabelColor(kRed);
   axis->Draw();
}
