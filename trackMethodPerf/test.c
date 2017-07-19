//#include "emiltoolsH82017.cpp"
#include "TRandom3.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"

TRandom3 *myRand = new TRandom3();

double xR;
double xR1;
double xRT;
double xR2;

double xMC_1;
double xMC_T;
double xMC_2;

double res= 0.060;//[mm]
double sigma= res*sqrt(2.0);//[mm]

//offset is how much higher in x the system is located (/displaced)
//positive offset means lower values read
double offset1 = 0.0;
double offsetT = 0.0;
double offset2 = 0.0;
double z1=0.0;//leave this 0.0 by definition, for simplicity
double zT=450;
double z2=500;

//the track equation: x=az+b
//to reconstruct: we have the z of the coord.system, and the a,b.

//z is horizontal
//x is vertical

//first we assume straight tracks
double track_a=0.000;//mm
double track_b=1.1;//mm
//the reconstructed track's info
double rec_track_a;
double rec_track_b;

//---------
int hbp_bins = 100;
double hbp_min = -.5;
double hbp_max = .5;

int    hres_bins = 100;
double hres_min = -2;
double hres_max = 2;

TH1D * h_real1   ;//= new TH1D("h_real1","h_real1",hbp_bins,track_b-offset1+hbp_min,track_b-offset1+hbp_max);
TH1D * h_realT   ;//= new TH1D("h_realT","h_realT",hbp_bins,track_b-offsetT+hbp_min,track_b-offsetT+hbp_max);
TH1D * h_real2   ;//= new TH1D("h_real2","h_real2",hbp_bins,track_b-offset2+hbp_min,track_b-offset2+hbp_max);

TH1D * h_MC1     ;//= new TH1D("h_MC1","h_MC1",hbp_bins,track_b-offset1+hbp_min,track_b-offset1+hbp_max);
TH1D * h_MC2     ;//= new TH1D("h_MC2","h_MC2",hbp_bins,track_b-offset2+hbp_min,track_b-offset2+hbp_max);

TH1D * h_res_MC1 ;//= new TH1D("h_res_MC1","h_res_MC1",hres_bins,hres_min,hres_max);
TH1D * h_res_MC2 ;//= new TH1D("h_res_MC2","h_res_MC2",hres_bins,hres_min,hres_max);

TH1D * h_res_tr1 ;//= new TH1D("h_res_tr1","Diffs xR_T-xReco_T - reco:1",hres_bins,hres_min,hres_max);
TH1D * h_res_tr2 ;//= new TH1D("h_res_tr2","Diffs xR_T-xReco_T - reco:simple diff",hres_bins,hres_min,hres_max);
TH1D * h_res_tr3 ;//= new TH1D("h_res_tr3","Diffs xR_T-xReco_T - reco:prop tracking",hres_bins,hres_min,hres_max);

void setups()
{
    h_MC1->SetLineColor(kRed);
    h_MC2->SetLineColor(kRed);

    gStyle->SetOptFit(1);
}
void initHistos()
{
    h_real1   = new TH1D("h_real1","h_real1",hbp_bins,track_b-offset1+hbp_min,track_b-offset1+hbp_max);
    h_realT   = new TH1D("h_realT","h_realT",hbp_bins,track_b-offsetT+hbp_min,track_b-offsetT+hbp_max);
    h_real2   = new TH1D("h_real2","h_real2",hbp_bins,track_b-offset2+hbp_min,track_b-offset2+hbp_max);

    h_MC1     = new TH1D("h_MC1","h_MC1",hbp_bins,track_b-offset1+hbp_min,track_b-offset1+hbp_max);
    h_MC2     = new TH1D("h_MC2","h_MC2",hbp_bins,track_b-offset2+hbp_min,track_b-offset2+hbp_max);

    h_res_MC1 = new TH1D("h_res_MC1","h_res_MC1",hres_bins,hres_min,hres_max);
    h_res_MC2 = new TH1D("h_res_MC2","h_res_MC2",hres_bins,hres_min,hres_max);

    h_res_tr1 = new TH1D("h_res_tr1","Diffs xR_T-xReco_T - reco:1",hres_bins,hres_min,hres_max);
    h_res_tr2 = new TH1D("h_res_tr2","Diffs xR_T-xReco_T - reco:simple diff",hres_bins,hres_min,hres_max);
    h_res_tr3 = new TH1D("h_res_tr3","Diffs xR_T-xReco_T - reco:prop tracking",hres_bins,hres_min,hres_max);
}

void run(long noOfEvents=1000)
{
    initHistos();
    setups();

    for(int i=0;i<noOfEvents;i++)
    {

        /*
             * foreach event, in this loop the track is the same
             *
             * 0) [unnecessary] we assume that system1 is on 0,0, so that track_b=xR1
             * 1) we use the track info to generate the real X
             * 2) from the real X we create (MC) the xMC 1,2
             * 3) from the xMC1,2 we create the xMC_T (the reconstructed one)
             * 4) then we fill all histos, and save the data somehow...
             * 5) then we can run again with different track and z parameters.
             *
             **/


        //make real X
        xR1 = z1*track_a+track_b-offset1;//z1=0 in general
        xRT = zT*track_a+track_b-offsetT;
        xR2 = z2*track_a+track_b-offset2;

        h_real1->Fill(xR1);
        h_realT->Fill(xRT);
        h_real2->Fill(xR2);

        //create (MC)
        xMC_1 = xR1 + myRand->Gaus(0,sigma);
        xMC_2 = xR2 + myRand->Gaus(0,sigma);

        h_MC1->Fill(xMC_1);
        h_MC2->Fill(xMC_2);

        h_res_MC1->Fill(xMC_1-xR1);
        h_res_MC2->Fill(xMC_2-xR2);



        //***************************************************


        //NOW THE RECONSTRUCTION METHODS to create the x of Test chamber
        //1111 simply the first chamber-----------------
        xMC_T = xMC_1;
        //need to adjust offset
        h_res_tr1->Fill(xMC_T-xRT+1.5);

        //2222 average of 2 chambers--------------------
        xMC_T = (xMC_1+xMC_2)/2;
        h_res_tr2->Fill(xMC_T-xRT+1.5);

        //3333 proper tracking -------------------------

        // 1)calculate the reconstructed track's a,b
        //the system to solve is this:
        //xMC_1 = rec_track_a*z1 + rec_track_b;
        //xMC_2 = rec_track_a*z2 + rec_track_b;
        rec_track_a = (xMC_1-xMC_2) / (z1-z2);
        rec_track_b = xMC_1-rec_track_a*z1;

        //so the reconstructed xT is:
        xMC_T = rec_track_a*zT+rec_track_b;

        h_res_tr3->Fill(xMC_T-xRT+1.5);

    }//event loop

    //Beam profiles (real, with MC overlaid)
    TCanvas * c_bp = new TCanvas("c_bp","c_bp",200,200,500,500);
    c_bp->Divide(1,3);
    c_bp->cd(1);
    h_MC1->Draw();
    h_real1->Draw("sames");
    c_bp->cd(2);
    h_realT->Draw();
    c_bp->cd(3);
    h_MC2->Draw();
    h_real2->Draw("sames");

    //MC1,2 resolutions (to confirm the one set in the setup)
    TCanvas * c_resMC = new TCanvas("c_resMC","c_resMC",200,200,500,500);
    c_resMC->Divide(2,1);
    c_resMC->cd(1);
    h_res_MC1->Fit("gaus","Q");
    c_resMC->cd(2);
    h_res_MC2->Fit("gaus","Q");


    //the method resolutions
    TCanvas * c_resTr = new TCanvas("c_resTr","c_resTr",200,200,600,600);
    c_resTr->Divide(3,1);
    c_resTr->cd(1);
    h_res_tr1->Fit("gaus","Q");
    c_resTr->cd(2);
    h_res_tr2->Fit("gaus","Q");
    c_resTr->cd(3);
    h_res_tr3->Fit("gaus","Q");

}

