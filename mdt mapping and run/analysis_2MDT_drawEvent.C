#define analysis_cxx
#include "analysis.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
using namespace std;
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TH1D.h"
#include <TStyle.h>
#include <TCanvas.h>
#include "TROOT.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraph.h"
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
#include "TRandom.h"
#include "TBox.h"
#include "TLatex.h"
#include "TArc.h"

int run(const char *  file, int nEntries, bool debug=false) {
   cout << "====>>>>>> Start of analysis <<<<<<==== \n";
   
   TChain *chain = new TChain("vmm");
   chain->Add(file);
   analysis sample(chain);
   // Run Analysis
   sample.Loop(nEntries,debug);
   
   cout << "====>>>>>> End of analysis <<<<<<====" << endl;
   
   return 0;
   
} //end of main loop

// plots for cardMDT0
TH1D *times0=new TH1D("times0","corrected_times0",100,-1000.,1000.);
TH1D *allChannels0=new TH1D("allChannels0","allChannels0",64,0,64.);
TH1D *allCharge0=new TH1D("allCharge0","allCharge0",1024,0,1024.);
TH1D *allHits0=new TH1D("allHits0","allHits0",30,0,30.);
TH1D *allTime0=new TH1D("allTime0","allTime0",500,0,500.);
TH1D *artTime0=new TH1D("artTime0","artTime0",500,0,500.);
TH2D *artTimeVSallTimes0=new TH2D("artTimeVSallTimes0","artTimeVSallTimes0",200,0,200.,200,0.,200.);
TH1D *h_diffGray0=new TH1D("h_diffGray0","h_diffGray0",10000,-5000.,+5000.);
TH1D *h_allHits0_aftercut =new TH1D("h_allHits0_aftercut","h_allHits0_aftercut",30,0,30.);
TH1D *h_allHits0_outOfRange_aftercut =new TH1D("h_allHits0_outOfRange_aftercut","h_allHits0_outOfRange_aftercut",30,0,30.);


// plots for cardMDT1
TH1D *times1=new TH1D("times1","corrected_times1",100,-1000.,1000.);
TH1D *allChannels1=new TH1D("allChannels1","allChannels1",64,0,64.);
TH1D *allCharge1=new TH1D("allCharge1","allCharge1",1024,0,1024.);
TH1D *allHits1=new TH1D("allHits1","allHits1",30,0,30.);
TH1D *allTime1=new TH1D("allTime1","allTime1",500,0,500.);
TH1D *artTime1=new TH1D("artTime1","artTime1",500,0,500.);
TH2D *artTimeVSallTimes1=new TH2D("artTimeVSallTimes1","artTimeVSallTimes1",200,0,200.,200,0.,200.);
TH1D *h_diffGray1=new TH1D("h_diffGray1","h_diffGray1",10000,-5000.,+5000.);
TH1D *h_allHits1_aftercut =new TH1D("h_allHits1_aftercut","h_allHits1_aftercut",30,0,30.);
TH1D *h_allHits1_outOfRange_aftercut =new TH1D("h_allHits1_outOfRange_aftercut","h_allHits1_outOfRange_aftercut",30,0,30.);

TH2D *allHitsCorrelation=new TH2D("hitCorrelation","hitCorrelation", 20, 0,20,20,0,20);
TH1D *allHits=new TH1D("allHits","allHits",64,0,64.);


//TH1D *allHits=new TH1D("allHits","allHits",30,0,30.);



int mapChannelTube(int channel)
{
   int tube=1000;
   
   if(channel == 0) tube=4;
   else if(channel == 1) tube=1;
   else if(channel == 2) tube=5;
   else if(channel == 3) tube=2;
   else if(channel == 4) tube=6;
   else if(channel == 5) tube=3;
   else if(channel == 6) tube=10;
   else if(channel == 7) tube=7;
   else if(channel == 8) tube=11;
   else if(channel == 9) tube=8;
   else if(channel == 10) tube=12;
   else if(channel == 11) tube=9;
   else if(channel == 12) tube=16;
   else if(channel == 13) tube=13;
   else if(channel == 14) tube=17;
   else if(channel == 15) tube=14;
   else if(channel == 16) tube=18;
   else if(channel == 17) tube=15;
   else if(channel == 18) tube=22;
   else if(channel == 19) tube=19;
   else if(channel == 20) tube=23;
   else if(channel == 21) tube=20;
   else if(channel == 22) tube=24;
   else if(channel == 23) tube=21;
   

   if(tube == 1) tube=6;
   else if(tube == 2) tube=5;
   else if(tube == 3) tube=4;
   else if(tube == 4) tube=3;
   else if(tube == 5) tube=2;
   else if(tube == 6) tube=1;
   else if(tube == 7) tube=12;
   else if(tube == 8) tube=11;
   else if(tube == 9) tube=10;
   else if(tube == 10) tube=9;
   else if(tube == 11) tube=8;
   else if(tube == 12) tube=7;
   else if(tube == 13) tube=18;
   else if(tube == 14) tube=17;
   else if(tube == 15) tube=16;
   else if(tube == 16) tube=15;
   else if(tube == 17) tube=14;
   else if(tube == 18) tube=13;
   else if(tube == 19) tube=24;
   else if(tube == 20) tube=23;
   else if(tube == 21) tube=22;
   else if(tube == 22) tube=21;
   else if(tube == 23) tube=20;
   else if(tube == 24) tube=19;

   
   return tube-1;
    
}


void drawEvent_bothBoards(int evt,
                          vector<int> channelId_beforeCut_boardId0,
                          vector<int> channelId_beforeCut_boardId1,
                          vector<int> channelId_afterCut_boardId0,
                          vector<int> channelId_afterCut_boardId1)
{
   TCanvas c_event_display("c_event_display","c_event_display",10,10,500,675);

   c_event_display.Clear();
   
   std::vector<TArc> drift_circles;
   
   //specs for BIS chambers
   double spacer_thickness = 6.5;
   double pitch = 30.0350;//distance between consecutive wires in the same layer
   double radius = 15.;//tube radius 14.6mm + 0.4mm tube wall thickness
   
   double x_start_l1 = 15.0;//x positions staggering
   double x_start_l2 = 30.0175;
   double x_start_l3 = 15.0;
   double x_start_l4 = 30.0175;
   
   double y_start_l1 = 45.7675;//y positions
   double y_start_l2 = 71.7786;
   double y_start_l3 = 97.7896;
   double y_start_l4 = 123.8007;
   
   double x=0.;
   double y=0.;
   double offset=0.;
   
   std::vector <double> x_wires;
   std::vector <double> y_wires;
   
   for(int i=1; i<=2; i++)
   {
      for(int j=1; j<=4; j++)
      {
         for(int k=0; k<6; k++)
         {
            if(i==2)
            {
               offset = spacer_thickness+y_start_l4;
            }
            if(i==1){
               if(j==1) { x = x_start_l1; y = y_start_l1+offset;}
               if(j==2) { x = x_start_l2; y = y_start_l2+offset;}
               if(j==3) { x = x_start_l3; y = y_start_l3+offset;}
               if(j==4) { x = x_start_l4; y = y_start_l4+offset;}
            }
            else if(i==2)
            {
               if(j==1) { x = x_start_l4; y = y_start_l1+offset;}
               if(j==2) { x = x_start_l3; y = y_start_l2+offset;}
               if(j==3) { x = x_start_l2; y = y_start_l3+offset;}
               if(j==4) { x = x_start_l1; y = y_start_l4+offset;}
            }
            x = x+k*pitch;
            
            x_wires.push_back(x);
            y_wires.push_back(y);
            
            drift_circles.push_back(TArc(x,y,radius));
            //get_chamb_MdtHits().at(it)->get_hit_pos_arc()->SetFillColor(38);
            drift_circles.back().SetLineColor(38);
            //            if((it)->isOnSegment)
            //               drift_circles.back()->SetFillColor(38);
         }
      }
   }
   
   
   TGraph chamb_display_wires_graph(x_wires.size(),&x_wires[0],&y_wires[0]);
   
   chamb_display_wires_graph.SetTitle(Form("Event %i",evt));
   chamb_display_wires_graph.GetXaxis()->SetTitle("[mm]");
   chamb_display_wires_graph.GetYaxis()->SetTitle("[mm]");
   
   c_event_display.cd();
   
   chamb_display_wires_graph.Draw("AP");
   

   for(int i=0; i<channelId_beforeCut_boardId0.size(); i++) {
          cout<< "in drawEvent: cannelID_beforeCut_boardId0= " << channelId_beforeCut_boardId0.at(i)   <<endl;
          drift_circles.at(24*0+mapChannelTube(channelId_beforeCut_boardId0.at(i))).SetFillColor(38);
   }
    for(int i=0; i<channelId_beforeCut_boardId1.size(); i++) {
            if(channelId_beforeCut_boardId1.at(i)!=24) {
                cout<< "in drawEvent: cannelID_beforeCut_boardId1= " << channelId_beforeCut_boardId1.at(i)   <<endl;
                drift_circles.at(24*1+mapChannelTube(channelId_beforeCut_boardId1.at(i))).SetFillColor(38);
            }
    }
    cout << "-------------------------------" << endl;

   for(int i=0; i<channelId_afterCut_boardId0.size(); i++) {
         drift_circles.at(24*0+mapChannelTube(channelId_afterCut_boardId0.at(i))).SetFillColor(kRed);
   }
    for(int i=0; i<channelId_afterCut_boardId1.size(); i++) {
                drift_circles.at(24*1+mapChannelTube(channelId_afterCut_boardId1.at(i))).SetFillColor(kRed);
    }
    
    for(int i=0; i<drift_circles.size(); i++)
   {
      drift_circles.at(i).Draw("Psame");
   }

   c_event_display.Update();
   c_event_display.cd();
   c_event_display.WaitPrimitive();
   
   
}  // end of drawEvent_bothBoards()


void drawEvent(int evt, vector<vector<int>> channelId_beforeCut, vector<vector<int>> channelId_afterCut)
{
    TCanvas c_event_display("c_event_display","c_event_display",10,10,500,675);
    
    c_event_display.Clear();
    
    std::vector<TArc> drift_circles;
    
    //specs for BIS chambers
    double spacer_thickness = 6.5;
    double pitch = 30.0350;//distance between consecutive wires in the same layer
    double radius = 15.;//tube radius 14.6mm + 0.4mm tube wall thickness
    
    double x_start_l1 = 15.0;//x positions staggering
    double x_start_l2 = 30.0175;
    double x_start_l3 = 15.0;
    double x_start_l4 = 30.0175;
    
    double y_start_l1 = 45.7675;//y positions
    double y_start_l2 = 71.7786;
    double y_start_l3 = 97.7896;
    double y_start_l4 = 123.8007;
    
    double x=0.;
    double y=0.;
    double offset=0.;
    
    std::vector <double> x_wires;
    std::vector <double> y_wires;
    
    for(int i=1; i<=2; i++)
    {
        for(int j=1; j<=4; j++)
        {
            for(int k=0; k<6; k++)
            {
                if(i==2)
                {
                    offset = spacer_thickness+y_start_l4;
                }
                if(i==1){
                    if(j==1) { x = x_start_l1; y = y_start_l1+offset;}
                    if(j==2) { x = x_start_l2; y = y_start_l2+offset;}
                    if(j==3) { x = x_start_l3; y = y_start_l3+offset;}
                    if(j==4) { x = x_start_l4; y = y_start_l4+offset;}
                }
                else if(i==2)
                {
                    if(j==1) { x = x_start_l4; y = y_start_l1+offset;}
                    if(j==2) { x = x_start_l3; y = y_start_l2+offset;}
                    if(j==3) { x = x_start_l2; y = y_start_l3+offset;}
                    if(j==4) { x = x_start_l1; y = y_start_l4+offset;}
                }
                x = x+k*pitch;
                
                x_wires.push_back(x);
                y_wires.push_back(y);
                
                drift_circles.push_back(TArc(x,y,radius));
                //get_chamb_MdtHits().at(it)->get_hit_pos_arc()->SetFillColor(38);
                drift_circles.back().SetLineColor(38);
                //            if((it)->isOnSegment)
                //               drift_circles.back()->SetFillColor(38);
            }
        }
    }
    
    
    TGraph chamb_display_wires_graph(x_wires.size(),&x_wires[0],&y_wires[0]);
    
    chamb_display_wires_graph.SetTitle(Form("Event %i",evt));
    chamb_display_wires_graph.GetXaxis()->SetTitle("[mm]");
    chamb_display_wires_graph.GetYaxis()->SetTitle("[mm]");
    
    c_event_display.cd();
    
    chamb_display_wires_graph.Draw("AP");
    
    
    for(int i=0; i<channelId_beforeCut.size(); i++)
    {
        for(int j=0; j<channelId_beforeCut.at(i).size(); j++)
        {
            if(channelId_beforeCut.at(i).at(j)!=24) {
                cout<< "in drawEvent: cannelID_beforeCut= " << channelId_beforeCut.at(i).at(j)   <<endl;
            drift_circles.at(23*i+mapChannelTube(channelId_beforeCut.at(i).at(j))).SetFillColor(38);
            }
        }
    }
    cout << "-------------------------------" << endl;
    
    for(int i=0; i<channelId_afterCut.size(); i++)
    {
        for(int j=0; j<channelId_afterCut.at(i).size(); j++)
        {
            if(channelId_afterCut.at(i).at(j)!=24) {
                //            continue;
                //std::cout<<channelId_afterCut.at(i).at(j)<<" -> "<<23*i+mapChannelTube(channelId_afterCut.at(i).at(j))<<std::endl;
                drift_circles.at(23*i+mapChannelTube(channelId_afterCut.at(i).at(j))).SetFillColor(kRed);
            }
        }
    }
    
    for(int i=0; i<drift_circles.size(); i++)
    {
        drift_circles.at(i).Draw("Psame");
    }
    
    c_event_display.Update();
    c_event_display.cd();
    c_event_display.WaitPrimitive();
    
    
    
}  // end of drawEvent()


void analysis::Loop(int givenEntries, bool debug)
{
    
    //   bool debug=true;
    
    gROOT->Reset();
    
   Long64_t nentries;
   if (fChain == 0) return;
   if (givenEntries==-1){
      nentries = fChain->GetEntries();
   }else nentries = givenEntries;
   
   
   
   int numEventsBothBoards=0;
   
   for (Long64_t iEntry=0; iEntry <nentries; iEntry++) {
      
      double eventT0gray0=0., eventT0tdo0=0., artSignal0=0.;
      int numHits0=0;
      double eventT0gray1=0, eventT0tdo1=0, artSignal1=0.;
      int numHits1=0, numHits=0 ;
      double conversionFactor = 3.78;  // in adc counts/ns
      double timeWindow=70.;
      double grayDecodedChannel24_0=-5000., grayDecodedChannel24_1=-5000;
      
      if (iEntry % 1000 == 0){
         cout << "====>>>>>> " << iEntry << " processed so far <<<<<<==== \n";
      }
      fChain->GetEntry(iEntry);
      
      vector<int> &boardNum = *boardId;
      vector<vector<int> > &time = *tdo;
      vector<vector<int> > &amp = *pdo;
      vector<vector<int> > &channelId = *channel;
      vector<vector<int> > &gray = *grayDecoded;
      
      
      vector<vector<int> > channelId_afterCut;
      vector<int> channelId_afterCut_boardId0;
      vector<int> channelId_afterCut_boardId1;
       vector<vector<int> > channelId_beforeCut;
       vector<int> channelId_beforeCut_boardId0;
       vector<int> channelId_beforeCut_boardId1;
       
       channelId_afterCut.clear(); channelId_afterCut_boardId0.clear(); channelId_afterCut_boardId1.clear();
       channelId_beforeCut.clear(); channelId_beforeCut_boardId0.clear(); channelId_beforeCut_boardId1.clear();

       
      
      // assume two MDT boards (0 and 1) or (2 and 3)
      int boardID0=0, boardID1=1;
      
      
      // get the art signal from channel==24
      bool saw_board0 = false;
      bool saw_board1 = false;
      numHits0=0, numHits1=0, numHits=0, grayDecodedChannel24_0=-5000., grayDecodedChannel24_1=-5000.;
      for(int num1=0;num1<amp.size();num1++) {
         for(int num2=0;num2<amp[num1].size();num2++) {
            
            
            if (channelId[num1][num2]==24) {
               if (boardNum[num1]==boardID0) {
                  artSignal0=time[num1][num2]; artTime0 -> Fill(artSignal0);
                  eventT0gray0=gray[num1][num2]*25.;
                  eventT0tdo0=eventT0gray0+timeWindow-time[num1][num2]/conversionFactor;
                  grayDecodedChannel24_0 = gray[num1][num2];
               } else if (boardNum[num1]==boardID1) {
                  artSignal1=time[num1][num2]; artTime1 -> Fill(artSignal1);
                  eventT0gray1=gray[num1][num2]*25.;
                  eventT0tdo1=eventT0gray1+timeWindow-time[num1][num2]/conversionFactor;
                  grayDecodedChannel24_1 = gray[num1][num2];
               }
            } // plot art times

             if (boardNum[num1]==boardID0) {
                 allChannels0 -> Fill(channelId[num1][num2]);
             } else if (boardNum[num1]==boardID1) {
                 allChannels1 -> Fill(channelId[num1][num2]);
             }
             
            if (channelId[num1][num2]!=24 && boardNum[num1]==boardID0) { saw_board0 = true; numHits0++;
//                cout << "Saw boards0, numHit0: " << numHits0  << "   channelNum= " << channelId.at(num1).at(num2)  << endl;
                channelId_beforeCut_boardId0.push_back(channelId.at(num1).at(num2));
            } // number of hits per event
            if (channelId[num1][num2]!=24 && boardNum[num1]==boardID1) { saw_board1 = true; numHits1++;
//                cout << "Saw boards1, numHit1: " << numHits1  << "   channelNum= " << channelId.at(num1).at(num2)  << endl;
                channelId_beforeCut_boardId1.push_back(channelId.at(num1).at(num2));
            } // number of hits per event
         }
      }
      if(saw_board0 && saw_board1) {
         numEventsBothBoards++;
         cout << "Saw both boards in same event, numHit0: " << numHits0 << " numHit1: " << numHits1 << " numEventsBothBoards= " << numEventsBothBoards   << endl;
         numHits = numHits0 + numHits1;
         allHitsCorrelation->Fill(numHits1,numHits0);
         
         allHits->Fill(numHits); // plot it
      }
      allHits0->Fill(numHits0); // plot num of hits in boardID0
      allHits1->Fill(numHits1); // plot num of hits in boardID1
      
      
      double diffGray0=-5000, diffGray1=-5000.;
      int numHits0_aftercut=0, numHits1_aftercut=0, numHits0_outOfRange_aftercut=0, numHits1_outOfRange_aftercut=0;
      //        double lowLimit=-6., highLimit=10.;
      double lowLimit=-15., highLimit=10.;
      for(int num1=0;num1<amp.size();num1++){
         for(int num2=0;num2<amp[num1].size();num2++){
            double stripTime0=0, stripTime1=0;
//             if (boardNum[num1]==boardID0 && grayDecodedChannel24_0!=-5000.) {
                 if (boardNum[num1]==boardID0) {
//               allChannels0 -> Fill(channelId[num1][num2]); // plot beam profile w/ also channel 24==artSignal
               if (channelId[num1][num2]!=24) {
                  allCharge0 -> Fill(amp[num1][num2]); // plot pulse heights, pdo
                  allTime0 -> Fill(time[num1][num2]);
                  artTimeVSallTimes0 -> Fill(time[num1][num2],artSignal0);
                  stripTime0=gray[num1][num2]*25.+(timeWindow-(time[num1][num2]/conversionFactor))-eventT0tdo0;
                  double latency = -70.;
                  if(stripTime0<latency)
                     stripTime0=latency-stripTime0+latency;
                  times0->Fill(stripTime0);
                  diffGray0 = gray[num1][num2]-grayDecodedChannel24_0; h_diffGray0 ->Fill(diffGray0);
 ////                 if (diffGray0>lowLimit && diffGray0<highLimit && grayDecodedChannel24_0!=-5000.)
                      if (diffGray0>lowLimit && diffGray0<highLimit)
                      {
                     numHits0_aftercut++;
                     channelId_afterCut_boardId0.push_back(channelId.at(num1).at(num2));
                  }
                  if ((diffGray0<lowLimit || diffGray0<highLimit)) numHits0_outOfRange_aftercut++;
               }
            } else if (boardNum[num1]==boardID1) {
//               allChannels1 -> Fill(channelId[num1][num2]); // plot beam profile w/ also channel 24==artSignal
               if (channelId[num1][num2]!=24) {
                  allCharge1 -> Fill(amp[num1][num2]); // plot pulse heights, pdo
                  allTime1 -> Fill(time[num1][num2]);
                  artTimeVSallTimes1 -> Fill(time[num1][num2],artSignal1);
                  stripTime1=gray[num1][num2]*25.+(timeWindow-(time[num1][num2]/conversionFactor))-eventT0tdo1;
                  double latency = -70.;
                  if(stripTime0<latency)
                     stripTime1=latency-stripTime1+latency;
                  times1->Fill(stripTime1);
                  diffGray1 = gray[num1][num2]-grayDecodedChannel24_1; h_diffGray1 ->Fill(diffGray1);
                //if (diffGray1>lowLimit && diffGray1<highLimit && grayDecodedChannel24_1!=-5000.)
                       if (diffGray1>lowLimit && diffGray1<highLimit)
                  {
                     numHits1_aftercut++;
                     channelId_afterCut_boardId1.push_back(channelId.at(num1).at(num2));
                  }
//            if ((diffGray1<lowLimit || diffGray1<highLimit) && grayDecodedChannel24_1!=-5000.)
                   if ((diffGray1<lowLimit || diffGray1<highLimit) && grayDecodedChannel24_1!=-5000.) numHits1_outOfRange_aftercut++;
               }
               
            }
            
         }
      }
      
      channelId_afterCut.push_back(channelId_afterCut_boardId0);channelId_afterCut.push_back(channelId_afterCut_boardId1);
      channelId_beforeCut.push_back(channelId_beforeCut_boardId0);channelId_beforeCut.push_back(channelId_beforeCut_boardId1);
       
      if(debug && saw_board0 && saw_board1 && numHits0>1 && numHits1>1) {
//          drawEvent(iEntry, channelId_beforeCut, channelId_afterCut);
        drawEvent_bothBoards(numEventsBothBoards, channelId_beforeCut_boardId0, channelId_beforeCut_boardId1, channelId_afterCut_boardId0, channelId_afterCut_boardId1);
 //         drawEvent(iEntry, channelId, channelId_afterCut);
      }

       
      h_allHits0_aftercut->Fill(numHits0_aftercut); // plot it
      h_allHits1_aftercut->Fill(numHits1_aftercut); // plot it
      
      h_allHits0_outOfRange_aftercut->Fill(numHits0_outOfRange_aftercut); // plot it
      h_allHits1_outOfRange_aftercut->Fill(numHits1_outOfRange_aftercut); // plot it
      
      
      
   }  // end of iEntry
   
   // Canvas for boardMDT0
   TCanvas *boardMDT0 = new TCanvas("boardMDT0","boardMDT0",1500,1000);
   gStyle->SetOptStat(kTRUE);
   gStyle->SetOptStat(1);
   boardMDT0 -> Clear();
   boardMDT0 -> Divide(3,3);
   
   boardMDT0 -> cd(1); times0->Draw();
//   TF1 *test2 = new TF1("gg2","landau(0)",-1.,+1);
//   test2->SetParameters(25.,-400.,50.);
//   times0->Fit(test2);
   boardMDT0 -> cd(2); allChannels0 -> Draw();
   boardMDT0 -> cd(3); allHits0 -> Draw();
   boardMDT0 -> cd(4); allCharge0 -> Draw();
   boardMDT0 -> cd(5); allTime0 -> Draw();
   boardMDT0 -> cd(6); artTime0 -> Draw();
   boardMDT0 -> cd(7); artTimeVSallTimes0 -> Draw();
   boardMDT0 -> cd(8); h_diffGray0 -> Draw();
   boardMDT0 -> cd(9); h_allHits0_aftercut -> Draw();
   boardMDT0 -> Update();
   
   // Canvas for boardMDT1
   TCanvas *boardMDT1 = new TCanvas("boardMDT1","boardMDT1",1500,1000);
   gStyle->SetOptStat(kTRUE);
   gStyle->SetOptStat(1);
   boardMDT1 -> Clear();
   boardMDT1 -> Divide(3,3);
   
   boardMDT1 -> cd(1); times1->Draw();
//   TF1 *test3 = new TF1("gg2","landau(0)",-1.,+1);
//   test3->SetParameters(25.,-400.,50.);
//   times1->Fit(test3);
   boardMDT1 -> cd(2); allChannels1 -> Draw();
   boardMDT1 -> cd(3); allHits1 -> Draw();
   boardMDT1 -> cd(4); allCharge1 -> Draw();
   boardMDT1 -> cd(5); allTime1 -> Draw();
   boardMDT1 -> cd(6); artTime1 -> Draw();
   boardMDT1 -> cd(7); artTimeVSallTimes1 -> Draw();
   boardMDT1 -> cd(8); h_diffGray1 -> Draw();
   boardMDT1 -> cd(9); h_allHits1_aftercut -> Draw();
   
   // Canvas for All boardMD1
   TCanvas *boardMDT = new TCanvas("boardMDT","boardMDT",1500,1000);
   gStyle->SetOptStat(kTRUE);
   gStyle->SetOptStat(1);
   boardMDT -> Clear();
   boardMDT -> Divide(3,3);
   
   boardMDT -> cd(1); allHits -> Draw();
   boardMDT ->cd(2); allHitsCorrelation->Draw("colz");
   
   boardMDT -> cd(4); h_allHits0_outOfRange_aftercut -> Draw();
   boardMDT -> cd(5); h_allHits1_outOfRange_aftercut -> Draw();
   
   boardMDT -> Update();
   
   
   
} // end of analysis::Loop

