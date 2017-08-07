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


double fit_dgaus_sameMean(TH1D* histogram)
{
    //here we do a quick silent fit to determine the offsets and sigmas for the
    TF1 *f_gausD= new TF1("f_gausD","gaus",-50,50);
    histogram->Fit("f_gausD","QN0");
    double peak_hist = f_gausD->GetParameter(0);
    double offset_hist = f_gausD->GetParameter(1);
    double sigma_hist  = f_gausD->GetParameter(2);
    //end of quick fit

    //create the fit function
    TF1 *total = new TF1("total","[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[1])/[4])**2)",-1*sigma_hist+offset_hist,1*sigma_hist+offset_hist);

    total->SetParNames("p_{tails}","#mu_{tails}","#sigma_{tails}","p_{core}","#sigma_{core}");

    total->SetParameters(peak_hist/3,offset_hist,sigma_hist*3,peak_hist,sigma_hist);

    total->SetLineColor(kBlue);
    //limit the mean values
    total->SetParLimits(1, offset_hist-.5, offset_hist+.5);
    //limit the sigmas
    total->SetParLimits(2, 0.01, 10.);
    total->SetParLimits(4, 0.01, 10.);
    //limit the peak values
    total->SetParLimits(0, 0., peak_hist*2);
    total->SetParLimits(3, 0., peak_hist*2);

    total->SetLineWidth(2);

    gStyle->SetOptFit(1);
    histogram->Fit(total,"QR");

    //    TF1 *total = new TF1("total","[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[1])/[4])**2)",-7*sigma_hist+offset_hist,7*sigma_hist+offset_hist);
    TF1  *tf1 = new TF1("tf1","gaus",-7*sigma_hist+offset_hist,7*sigma_hist+offset_hist);
    TF1  *tf2 = new TF1("tf2","gaus",-7*sigma_hist+offset_hist,7*sigma_hist+offset_hist);

    tf1->SetParameters(total->GetParameter(0),total->GetParameter(1),total->GetParameter(2));
    tf2->SetParameters(total->GetParameter(3),total->GetParameter(1),total->GetParameter(4));
    tf1->Draw("sames");
    tf2->SetLineColor(kBlack);
    tf2->Draw("sames");

    //    TF1* core = new TF1("gaus");
    //    TF1* tails = new TF1("gaus");

    //    core->SetParameter(0,total->GetParameter(0));
    //    core->SetParameter(1,total->GetParameter(3));
    //    core->SetParameter(2,total->GetParameter(3));

    float res_fin=0.0;
    float res_fin_weight;
    TLatex t_resolution_weight;
    TLatex t_resolution;

    t_resolution.SetNDC();
    t_resolution_weight.SetNDC();

    if(total->GetParameter(0)>total->GetParameter(3))
        res_fin=total->GetParameter(2);
    else
        res_fin=total->GetParameter(4);

    double sigma_core;
    if(total->GetParameter(2)<total->GetParameter(4))
    {
        sigma_core = total->GetParameter(2);
    }
    else
    {
        sigma_core = total->GetParameter(5);
    }
//    TF1 *total = new TF1("total","[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[1])/[4])**2)"

    res_fin_weight=sqrt(
                (total->GetParameter(0) * TMath::Power(total->GetParameter(2),3) + total->GetParameter(3) * TMath::Power(total->GetParameter(4),3)  )
                /
                (   total->GetParameter(0) * total->GetParameter(2)   + total->GetParameter(3) * total->GetParameter(4) )
                );

    cout << "res_fin = "<< res_fin << endl ;
    cout << "res_fin_weight = "<< res_fin_weight << endl;

    t_resolution.DrawLatex(0.2,0.5,Form("#scale[0.8]{#color[1]{#sigma_{core}=%3.0f #mum}}",1,TMath::Abs(1000.0*res_fin)));
    t_resolution_weight.DrawLatex(0.2,0.45,Form("#scale[0.8]{#color[1]{#sigma_{weight}=%3.0f #mum}}",TMath::Abs(1000.0*res_fin_weight)));

    t_resolution.DrawLatex(0.2,0.4,Form("#scale[0.7]{#color[1]{res_{core}=%3.0f #mum}}",1,TMath::Abs(1000.0*res_fin)/sqrt(2)));
    t_resolution_weight.DrawLatex(0.2,0.35,Form("#scale[0.8]{#color[1]{res_{weight}=%3.0f #mum}}",TMath::Abs(1000.0*res_fin_weight)/sqrt(2)));


    return res_fin;
}

double fit_dgaus(TH1D* histogram)
{
    TLegend *leg = new TLegend(.1285,.529,.43,.88);

    //here we do a quick silent fit to determine the offsets and sigmas for the
    TF1 *f_gausD= new TF1("f_gausD","gaus",-50,50);
    histogram->Fit("f_gausD","QN0");
    double peak_hist = f_gausD->GetParameter(0);
    double offset_hist = f_gausD->GetParameter(1);
    double sigma_hist  = f_gausD->GetParameter(2);
    //end of quick fit

    //create the fit functions.
    TF1 *g1= new TF1("g1","gaus",-3*sigma_hist+offset_hist,3*sigma_hist+offset_hist);
    TF1 *g2= new TF1("g2","gaus",-3*sigma_hist+offset_hist,3*sigma_hist+offset_hist);

    g1->SetParameters(peak_hist*0.4, offset_hist, sigma_hist*2);
    g2->SetParameters(peak_hist, offset_hist, sigma_hist);

    TF1 *total = new TF1("total","g1+g2",-3*sigma_hist+offset_hist,3*sigma_hist+offset_hist);

    total->SetParNames("p_{tails}","#mu_{tails}","#sigma_{tails}","p_{core}","#mu_{core}","#sigma_{core}");

    total->SetLineColor(kBlue);
    //limit the mean values
    total->SetParLimits(1, offset_hist-.5, offset_hist+.5);
    total->SetParLimits(4, offset_hist-.5, offset_hist+.5);
    //limit the sigmas
    total->SetParLimits(2, 0.1, 12.);
    total->SetParLimits(5, 0.01, 7.);
    //limit the peak values
    total->SetParLimits(0, 0., peak_hist*0.6);
    total->SetParLimits(3, 0., peak_hist*1.5);

    g1->SetLineWidth(3);
    g2->SetLineWidth(3);
    total->SetLineWidth(3);

    gStyle->SetOptFit(1);
    histogram->Fit(total,"QR");

    g1->SetParameters(total->GetParameter(0), total->GetParameter(1), total->GetParameter(2));
    g2->SetParameters(total->GetParameter(3), total->GetParameter(4), total->GetParameter(5));

    float res_fin=0.0;
    float res_fin_weight;
    TLatex t_resolution_weight;
    TLatex t_resolution;
    //    int color_core;

    t_resolution.SetNDC();

    t_resolution_weight.SetNDC();
    if(total->GetParameter(0)>total->GetParameter(3))
    {
        res_fin=total->GetParameter(2);
        //        color_core=2;
        g1->SetLineColor(kOrange);
        g2->SetLineColor(kRed);
    }
    else
    {
        res_fin=total->GetParameter(5);
        //        color_core=3;
        g1->SetLineColor(kRed);
        g2->SetLineColor(kOrange);
    }
    double sigma_core;
    if(total->GetParameter(2)<total->GetParameter(5))
    {
        sigma_core = total->GetParameter(2);
    }
    else
    {
        sigma_core = total->GetParameter(5);
    }
    g1->Draw("sames");
    g2->Draw("sames");

    res_fin_weight=sqrt(
                (total->GetParameter(0) * TMath::Power(total->GetParameter(2),3) + total->GetParameter(3) * TMath::Power(total->GetParameter(5),3)  )
                /
                (   total->GetParameter(0) * total->GetParameter(2)   + total->GetParameter(3) * total->GetParameter(5) )
                );

    cout << "res_fin = "<< res_fin << endl ;
    cout << "res_fin_weight = "<< res_fin_weight << endl;

    t_resolution.DrawLatex(0.2,0.5,Form("#scale[0.8]{#color[1]{#sigma_{core}=%3.0f #mum}}",1,TMath::Abs(1000.0*res_fin)));
    t_resolution_weight.DrawLatex(0.2,0.45,Form("#scale[0.8]{#color[1]{#sigma_{weight}=%3.0f #mum}}",TMath::Abs(1000.0*res_fin_weight)));

    t_resolution.DrawLatex(0.2,0.4,Form("#scale[0.7]{#color[1]{res_{core}=%3.0f #mum}}",1,TMath::Abs(1000.0*res_fin)/sqrt(2)));
    t_resolution_weight.DrawLatex(0.2,0.35,Form("#scale[0.8]{#color[1]{res_{weight}=%3.0f #mum}}",TMath::Abs(1000.0*res_fin_weight)/sqrt(2)));


    return res_fin;
}
