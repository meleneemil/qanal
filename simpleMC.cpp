#include "emiltools.cpp"
#include "TRandom3.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "TGraph2DErrors.h"
#include "TLegend.h"

TH1D *singleHist;
TRandom3 *myRand = new TRandom3();
float xMin;
float xMax;

float DxValue;
float SigValue;
float SigError;

TGraph *tg;
TGraph *tg_overlord;

TCanvas *c_multi;
TCanvas *c_mcs;

vector<float> *v_all_DxOverSigma;
vector<float> *v_all_Iter;
bool runAll=false;

TF1 * f1 = new TF1("f1","0.412062*TMath::Log(x)+4.01841",100,10000000);
TF1 * f2 = new TF1("f2","0.412062*TMath::Log10(x)+4.01841",100,10000000);


void runSingle(double noOfIterations, float MCSigma, long histoBins=200)
{
    DxValue=0.0;
    SigValue=0.0;
    SigError=0.0;
    xMin=999.0;
    xMax=-999.0;

    float lowLimit = -10*MCSigma;
    float upLimit = 10*MCSigma;

    singleHist = new TH1D("single","single",histoBins,lowLimit,upLimit);

    cout << "running loop for " << noOfIterations << " iterations" <<endl;

    for(int i=0;i<noOfIterations;i++)
    {
        float x = myRand->Gaus(0.0,MCSigma);
        if(x<xMin)
            xMin=x;
        else if(x>xMax)
            xMax=x;

        singleHist->Fill(x);
    }

    DxValue = xMax-xMin;
    TF1 *f1 = new TF1("myfit","gaus", lowLimit, upLimit);
    singleHist->Fit(f1,"Q0");
    SigValue = singleHist->GetRMS();
//    SigValue = f1->GetParameter(2);
//    SigError = f1->GetParError(2);
    cout << "noIters="<<noOfIterations<<" | FIT SigValue=" <<SigValue<<" | SET Sigma = "<<MCSigma << " | DX = "<<DxValue;
    cout <<endl;
}

void runAllN(float MCSigma, long histBins = 200, bool plot=false)
{
    int size=11;
    float *iterArray = new float[size];
    float *DxOverSigma = new float[size];
    float *DxOverSigmaERROR = new float[size];

    iterArray[0]=100;
    iterArray[1]=500;
    iterArray[2]=1000;
    iterArray[3]=5000;
    iterArray[4]=10000;
    iterArray[5]=50000;
    iterArray[6]=100000;
    iterArray[7]=500000;
    iterArray[8]=1000000;
    iterArray[9]=5000000;
    iterArray[10]=10000000;

    if(plot)
    {
        c_mcs = new TCanvas("mcs","mcs",400,100,700,700);
        c_mcs->Divide(4,3);

    }

    for(int i=0;i<size;i++)
    {
        runSingle(iterArray[i],MCSigma,histBins);
        DxOverSigma[i] = DxValue/SigValue;

        if(runAll)
        {
            if(iterArray[i] > histBins)
            {
                v_all_DxOverSigma->push_back(DxValue/SigValue);
                v_all_Iter->push_back(iterArray[i]);
            }
        }

        DxOverSigmaERROR[i] = 1/(SigValue*SigValue)*SigError*DxValue;

        if(plot)
        {
            c_mcs->cd(i+1);
            singleHist->Draw();
            c_mcs->Modified();
            c_mcs->Update();
        }
    }

    tg = new TGraphErrors(size,iterArray,DxOverSigma,0,/*DxOverSigmaERROR*/0);

    if(plot)
    {
        TCanvas *t = new TCanvas("ds","ds",100,100,600,600);
        t->SetLogx();
        tg->SetTitle("(xMax-xMin)/sigma VS noOfBins");
        tg->Draw("AL*");

    }
}

void runManySigmas(int TimesToRepeatEachSigma=1, long histoBins=500,bool saveCanvas=false)
{

    >>> kane average se kathe "column" stoixeiwn gia na einai pio katharo

    TMultiGraph *mg = new TMultiGraph();

    TLegend *legend = new TLegend(0.5,0.15,0.85,0.45);

    vector<float> *v_mcsigmas = new vector<float>();
    v_mcsigmas->push_back(1);
    //v_mcsigmas->push_back(1.5);
    //v_mcsigmas->push_back(2);
    //v_mcsigmas->push_back(5);
    //v_mcsigmas->push_back(10);
    //v_mcsigmas->push_back(20);

    for(int i=0;i<v_mcsigmas->size();i++)
    {
        float mcsigma = v_mcsigmas->at(i);
        for(int sigmaRepeats=0;sigmaRepeats<TimesToRepeatEachSigma;sigmaRepeats++)
        {
            runAllN(mcsigma,histoBins);
            tg->SetLineColor(1+i);
            tg->SetLineWidth(3);
            tg->SetTitle(Form("mc_sigma = %f",mcsigma));
            mg->Add(tg);

        }
        legend->AddEntry(tg,Form("Sigma=%0.f",mcsigma));

        //these 3 delete:
        //        mg->Draw("AL*");
        //        c_multi->Modified();
        //        c_multi->Update();
    }

    c_multi = new TCanvas("dsa2","dsa2",200,100,1000,900);
    c_multi->SetLogx();
    mg->SetTitle(Form("MC(gaus) of (xMax-xMin)/sigma VS noOfIterations - Run %i times foreach Sigma - bins=%i",TimesToRepeatEachSigma,histoBins));
    mg->Draw("AL*");
    c_multi->Modified();
    c_multi->Update();

    //-------------------
    f1->SetLineColor(kBlue);
    f1->SetLineWidth(5);
    f1->SetLineStyle(9);
    f1->Draw("sames");
    legend->AddEntry(f1,Form("0.41*Ln(x)+4"));
    //-------------------
    f2->SetLineColor(kRed);
    f2->SetLineWidth(5);
    f2->SetLineStyle(9);
    f2->Draw("sames");
    legend->AddEntry(f2,Form("0.41*Log10(x)+4"));
//0.412062*TMath::Log(x)+4.01841
    legend->Draw();

    if(saveCanvas)
        c_multi->SaveAs(Form("MC(gaus) of (xMax-xMin) over sigma VS noOfIterations - Run %i times foreach Sigma - bins=%i.png",TimesToRepeatEachSigma,histoBins));

}

void runManyBins()
{
    v_all_DxOverSigma = new vector<float>();
    v_all_Iter = new vector<float>();
    //runAll=true;

    vector<long> *v_bins = new vector<long>();
    v_bins->push_back(100);
    v_bins->push_back(1000);
    v_bins->push_back(10000);
    v_bins->push_back(100000);

    for(int i=0;i<v_bins->size();i++)
    {
        runManySigmas(10,v_bins->at(i),true);
    }

    if(runAll)
    {

        long size = v_all_DxOverSigma->size();
        float *iterArray = new float[size];
        float *dxOSigmaArray = new float[size];

        TCanvas *c_allFit = new TCanvas("allfit","allfit",100,100,500,500);

        for(int i=0;i<size;i++)
        {
            iterArray[i] = v_all_Iter->at(i);
            dxOSigmaArray[i] = v_all_DxOverSigma->at(i);
        }

        tg_overlord = new TGraph(size,iterArray,dxOSigmaArray);

        tg_overlord->Draw();
    }


}





