#include "emiltoolsH82017_VMM3.cpp"
//-------SETUP
int g_boardId0=0;
int g_boardId1=1;
int g_boardDistance_mm=37;
int g_PdoPedestal = 32;
double g_correctionSlope = -0.0061;
double g_MultiClusterSearchOffset01_mm=-0.6194;//this is taken from our single cluster event run, to know the misalignment of our chambers
double g_MultiClusterSearchEpsilon_mm=1.;//this is the error allowed when looking for the pair-cluster 1.5mm away e.g..
bool eventByEvent=false;
bool doClusterisation=true;
bool doSingleClusterEventAnalysis=true;
bool doMultiClusterEventAnalysis=true;
bool doEtaCorrection=false;
bool isEtaCorrectionHistoFilled=false;//leave it false
TLine *l1;

//--------------------------------------GLOBAL strip cuts config
int cut_strip_minPerChipEvent=2;

//PDO TDO cuts
int cut_strip_minPdo=32;//don't put this lower than g_PdoPedestal
int cut_strip_minTdo=10;

//edge cuts
int cut_strip_maxStrip= 60;
int cut_strip_minStrip= 2;//was 2

vector<int> *cut_strip_noisySet = new vector<int>();//remove these strips from the ones to be clusterised

//--------------------------------------CLUSTER CUTS config
int cut_cluster_maxNEmptyStrips=1;//inbetween, not total
int cut_cluster_maxClusterCharge=12000;
int cut_cluster_minClusterCharge=50;
int cut_cluster_maxClusterTimeWidth=50;
int cut_cluster_minClusterTimeWidth=-1;
int cut_cluster_minClusterSize = 2;
int cut_cluster_maxClusterSize = 8;
double cut_cluster_etaRatio_highCut=1.;
double cut_cluster_etaRatio_lowCut=0.0;

//vector<int> cut_cluster_stripContainBlacklist;//kill cluster if contains any of these strips | Let's leave this for now...
vector<int> *cut_cluster_stripNeighbourBlacklist = new vector<int>();//kill cluster if is 1st neighbour with any of these strips
vector<int> *cut_cluster_stripNeighbourBlacklist0 = new vector<int>();//kill cluster if is 1st neighbour with any of these strips
vector<int> *cut_cluster_stripNeighbourBlacklist1 = new vector<int>();//kill cluster if is 1st neighbour with any of these strips

//--------------------------------------EVENT CUTS config
int cut_event_maxNClustersPerBoard=3;
//-----------------------------------  Histograms

TH1D* h_bp_board0 = new TH1D("h_bp_board0","h_bp_board0",64,0,64);
TH1D* h_pdo_board0 = new TH1D("h_pdo_board0","h_pdo_board0",1024,0,1024);
TH1D* h_tdo_board0 = new TH1D("h_tdo_board0","h_tdo_board0",256,0,256);
TH1D* h_bcid_board0 = new TH1D("h_bcid_board0","h_bcid_board0",4096,0,4096);
TH1D* h_stripPerEv_board0 = new TH1D("h_stripPerEv_board0","h_stripPerEv_board0",64,0,64);
int i_stripPerEvent0=0;

TH2D* h_pdoVSchannel_board0 = new TH2D("h_pdoVSchannel_board0","h_pdoVSchannel_board0",64,0,64,1024,0,1024);
TH2D* h_tdoVSchannel_board0 = new TH2D("h_tdoVSchannel_board0","h_tdoVSchannel_board0",64,0,64,256,0,256);
TH2D* h_bcidVSchannel_board0 = new TH2D("h_bcidVSchannel_board0","h_bcidVSchannel_board0",64,0,64,256,0,256);
TH2D* h_pdoVStdo_board0 = new TH2D("h_pdoVStdo_board0","h_pdoVStdo_board0",256,0,256,1024,0,1024);
TH2D* h_pdoVSbcid_board0 = new TH2D("h_pdoVSbcid_board0","h_pdoVSbcid_board0",4096,0,4096,1024,0,1024);

TH1D* h_bp_board1 = new TH1D("h_bp_board1","h_bp_board1",64,0,64);
TH1D* h_pdo_board1 = new TH1D("h_pdo_board1","h_pdo_board1",1024,0,1024);
TH1D* h_tdo_board1 = new TH1D("h_tdo_board1","h_tdo_board1",256,0,256);
TH1D* h_bcid_board1 = new TH1D("h_bcid_board1","h_bcid_board1",4096,0,4096);
TH1D* h_stripPerEv_board1 = new TH1D("h_stripPerEv_board1","h_stripPerEv_board1",64,0,64);
int i_stripPerEvent1=0;

TH2D* h_pdoVSchannel_board1 = new TH2D("h_pdoVSchannel_board1","h_pdoVSchannel_board1",64,0,64,1024,0,1024);
TH2D* h_tdoVSchannel_board1 = new TH2D("h_tdoVSchannel_board1","h_tdoVSchannel_board1",64,0,64,256,0,256);
TH2D* h_bcidVSchannel_board1 = new TH2D("h_bcidVSchannel_board1","h_bcidVSchannel_board1",64,0,64,256,0,256);
TH2D* h_pdoVStdo_board1 = new TH2D("h_pdoVStdo_board1","h_pdoVStdo_board1",256,0,256,1024,0,1024);
TH2D* h_pdoVSbcid_board1 = new TH2D("h_pdoVSbcid_board1","h_pdoVSbcid_board1",4096,0,4096,1024,0,1024);

//------- Reverse Beam profile
//single cluster on one
//no cluster on the other
TH1D* h_reverseBeamProfile_board0 = new TH1D("h_reverseBeamProfile_board0","h_reverseBeamProfile_board0",240,0,64*.4);
TH1D* h_reverseBeamProfile_board1 = new TH1D("h_reverseBeamProfile_board1","h_reverseBeamProfile_board1",240,0,64*.4);

//------- SINGLE CLUSTER

TH1D* h_nClu_board0 = new TH1D("h_nClu_board0","h_nClu_board0",10,0,10);
TH1D* h_nClu_board1 = new TH1D("h_nClu_board1","h_nClu_board1",10,0,10);

TH1D* h_clPos_board0 = new TH1D("h_clPos_board0","h_clPos_board0",120,0,64*.4);
TH1D* h_clPos_board1 = new TH1D("h_clPos_board1","h_clPos_board1",120,0,64*.4);

TH1D* h_clCharge_board0 = new TH1D("h_clCharge_board0","h_clCharge_board0",800,0,8000);
TH1D* h_clCharge_board1 = new TH1D("h_clCharge_board1","h_clCharge_board1",800,0,8000);

TH1D* h_clSize_board0 = new TH1D("h_clSize_board0","h_clSize_board0",20,0,20);
TH1D* h_clSize_board1 = new TH1D("h_clSize_board1","h_clSize_board1",20,0,20);

TH2D* h_nCluVSnClu_boards = new TH2D("h_nCluVSnClu_boards","h_nCluVSnClu_boards",10,0,10,10,0,10);

TH1D* h_res = new TH1D("h_res","h_res",200,-2,2);

//------- MULTI CLUSTER

//TH1D* h_multiCluster_nClu_board0       = new TH1D("h_multiCluster_nClu_board0"      ,"h_multiCluster_nClu_board0",10,0,10);
//TH1D* h_multiCluster_nClu_board1       = new TH1D("h_multiCluster_nClu_board1"      ,"h_multiCluster_nClu_board1",10,0,10);
TH1D* h_multiCluster_clPos_board0      = new TH1D("h_multiCluster_clPos_board0"     ,"h_multiCluster_clPos_board0",120,0,64*.4);
TH1D* h_multiCluster_clPos_board1      = new TH1D("h_multiCluster_clPos_board1"     ,"h_multiCluster_clPos_board1",120,0,64*.4);
TH1D* h_multiCluster_clCharge_board0   = new TH1D("h_multiCluster_clCharge_board0"  ,"h_multiCluster_clCharge_board0",800,0,8000);
TH1D* h_multiCluster_clCharge_board1   = new TH1D("h_multiCluster_clCharge_board1"  ,"h_multiCluster_clCharge_board1",800,0,8000);
TH1D* h_multiCluster_clSize_board0     = new TH1D("h_multiCluster_clSize_board0"    ,"h_multiCluster_clSize_board0",20,0,20);
TH1D* h_multiCluster_clSize_board1     = new TH1D("h_multiCluster_clSize_board1"    ,"h_multiCluster_clSize_board1",20,0,20);
//TH2D* h_multiCluster_nCluVSnClu_boards = new TH2D("h_multiCluster_nCluVSnClu_boards","h_multiCluster_nCluVSnClu_boards",10,0,10,10,0,10);
TH1D* h_multiCluster_res               = new TH1D("h_multiCluster_multiCluster_res" ,"h_multiCluster_multiCluster_res",200,-2,2);

TH2D* h_multiCluster_diff01_VS_pos0 = new TH2D("h_multiCluster_diff01_VS_pos0","h_multiCluster_diff01_VS_pos0",120,0,64*.4,200,-2,2);
TH2D* h_multiCluster_diff01_VS_pos1 = new TH2D("h_multiCluster_diff01_VS_pos1","h_multiCluster_diff01_VS_pos1",120,0,64*.4,200,-2,2);

TH1D* h_multiCluster_diff_over_d01 = new TH1D("h_multiCluster_diff_over_d01","h_multiCluster_diff_over_d01",200,-2/g_boardDistance_mm,2/g_boardDistance_mm);

//--------------Event by event histos
TCanvas* c_ev;

TH1D* h_ev_pdo_board0 = new TH1D("h_ev_pdo_board0","h_ev_pdo_board0",64,0,64);
TH1D* h_ev_tdo_board0 = new TH1D("h_ev_tdo_board0","h_ev_tdo_board0",64,0,64);
TGraph *tg_ev_tdo_board0 = new TGraph(64);
TH1D* h_ev_bcid_board0 = new TH1D("h_ev_bcid_board0","h_ev_bcid_board0",64,0,64);

TH1D* h_ev_pdo_board1 = new TH1D("h_ev_pdo_board1","h_ev_pdo_board1",64,0,64);
TH1D* h_ev_tdo_board1 = new TH1D("h_ev_tdo_board1","h_ev_tdo_board1",64,0,64);
TGraph *tg_ev_tdo_board1 = new TGraph(64);
TH1D* h_ev_bcid_board1 = new TH1D("h_ev_bcid_board1","h_ev_bcid_board1",64,0,64);

//------------------------------ ETA CORRECTION (when cl_size>=2strip, check the ratios)
TH1D* h_etaRatio_board0 = new TH1D("h_etaRatio_board0","h_etaRatio_board0",200,0,1);
TH1D* h_etaRatio_board1 = new TH1D("h_etaRatio_board1","h_etaRatio_board1",200,0,1);
TH2D* h_eta0_vs_diff = new TH2D("h_eta0_vs_diff","h_eta0_vs_diff",200,-2/g_boardDistance_mm,2/g_boardDistance_mm,200,0,1);
TH2D* h_eta1_vs_diff = new TH2D("h_eta1_vs_diff","h_eta1_vs_diff",200,-2/g_boardDistance_mm,2/g_boardDistance_mm,200,0,1);

TH1D* h_qL = new TH1D("h_qL","h_L",256,0,1024);
TH1D* h_qR = new TH1D("h_qR","h_R",256,0,1024);
TH2D* h_qL_vs_qR = new TH2D("h_qL_vs_qR","h_qL_vs_qR",256,0,1024,256,0,1024);
TH2D* h_qSmall_vs_qLarge = new TH2D("h_qSmall_vs_qLarge","h_qSmall_vs_qLarge",256,0,1024,256,0,1024);
TH2D* h_eta_vs_smallQ = new TH2D("h_eta_vs_smallQ","h_eta_vs_smallQ",256,0,1024,200,0,1.0);
//------------------- Clusterization vars


//this set of vectors: change every Event.
//one event initially has subevents, which are split chip-wise
//after doing the strip cuts
//we go for clusterization.
//Each cluster belongs to a chip, and has a set of channel,pdo,tdo,bcid values.
//So one cluster has one chip value, and other vectors.
//So each event, has a vector of clusters.
//At the end of the event we go through the clusters to analyze.
//so these vectors need to be emptied at the beginning of each event.

vector<int>            *v_cl_boardId     = new vector < int >();
vector< vector <int> > *vv_cl_channel = new vector < vector <int> >();
vector< vector <int> > *vv_cl_pdo     = new vector < vector <int> >();
vector< vector <int> > *vv_cl_tdo     = new vector < vector <int> >();
vector< vector <int> > *vv_cl_bcid    = new vector < vector <int> >();

//while applying cuts, we gather the strips to be clusterized
//here we do it per board
//vector<int> *v_toCl_boardId;
vector<int> *v_toCl_channel_board0;
vector<int> *v_toCl_pdo_board0;
vector<int> *v_toCl_tdo_board0;
vector<int> *v_toCl_bcid_board0;
vector<int> *v_toCl_channel_board1;
vector<int> *v_toCl_pdo_board1;
vector<int> *v_toCl_tdo_board1;
vector<int> *v_toCl_bcid_board1;
vector<double>* cl_sizes0   ;
vector<double>* cl_sizes1   ;
vector<double>* cl_charges0 ;
vector<double>* cl_charges1 ;
vector<double>* cl_poss0    ;
vector<double>* cl_poss1    ;
vector<double>* cl_etaRatio_values0;
vector<double>* cl_etaRatio_values1;


//############################################################
string cmdInput="";
TString hname;
TString sRunNumber;

void initialize()
{
    l1 = new TLine();
    //    l1
    cut_strip_noisySet->push_back(63);

    //    cut_cluster_stripNeighbourBlacklist->push_back(8);
    //    cut_cluster_stripNeighbourBlacklist->push_back(17);
    //    cut_cluster_stripNeighbourBlacklist->push_back(21);
    //    cut_cluster_stripNeighbourBlacklist->push_back(48);

    //    cut_cluster_stripNeighbourBlacklist0->push_back(21);

    cut_cluster_stripNeighbourBlacklist1->push_back(8);
    //    cut_cluster_stripNeighbourBlacklist1->push_back(17);
    cut_cluster_stripNeighbourBlacklist1->push_back(47);

    gStyle->SetLabelSize(0.1,"Y");
    gStyle->SetOptStat(111111);
    //    gStyle->SetStatH(0.5);
    gStyle->SetTitleFontSize(0.1);

    h_ev_pdo_board0->SetLineColor(kRed);
    h_ev_pdo_board0->SetLineWidth(2);
    h_ev_tdo_board0->SetLineColor(kRed);
    h_ev_tdo_board0->SetLineWidth(2);
    h_ev_bcid_board0->SetLineColor(kRed);
    h_ev_bcid_board0->SetLineWidth(2);
    h_ev_bcid_board0->SetMaximum(4096);


    h_ev_pdo_board1->SetLineColor(kBlue);
    h_ev_pdo_board1->SetLineWidth(2);
    h_ev_tdo_board1->SetLineColor(kBlue);
    h_ev_tdo_board1->SetLineWidth(2);
    h_ev_bcid_board1->SetLineColor(kBlue);
    h_ev_bcid_board1->SetLineWidth(2);
    h_ev_bcid_board1->SetMaximum(4096);
}
void drawHistos()
{
    TCanvas *c_stat_boards = new TCanvas("c_stat_boards","c_stat_boards",200,200,500,500);
    c_stat_boards->Divide(2,5);
    c_stat_boards->cd(1);
    h_bp_board0->Draw();
    c_stat_boards->cd(2);
    h_bp_board1->Draw();
    c_stat_boards->cd(3);
    h_pdo_board0->Draw();
    c_stat_boards->cd(4);
    h_pdo_board1->Draw();
    c_stat_boards->cd(5);
    h_tdo_board0->Draw();
    c_stat_boards->cd(6);
    h_tdo_board1->Draw();
    c_stat_boards->cd(7);
    h_bcid_board0->Draw();
    c_stat_boards->cd(8);
    h_bcid_board1->Draw();
    c_stat_boards->cd(9);
    h_stripPerEv_board0->Draw();
    c_stat_boards->cd(10);
    h_stripPerEv_board1->Draw();


    //-----------------------------GENERAL VS STATS

    TCanvas *c_VSstat_boards = new TCanvas("c_VSstat_boards","c_VSstat_boards",300,200,500,500);
    c_VSstat_boards->Divide(2,5);
    c_VSstat_boards->cd(1);
    h_pdoVSchannel_board0->Draw("colz");
    c_VSstat_boards->cd(2);
    h_pdoVSchannel_board1->Draw("colz");

    c_VSstat_boards->cd(3);
    h_tdoVSchannel_board0->Draw("colz");
    c_VSstat_boards->cd(4);
    h_tdoVSchannel_board1->Draw("colz");

    c_VSstat_boards->cd(5);
    h_bcidVSchannel_board0->Draw("colz");
    c_VSstat_boards->cd(6);
    h_bcidVSchannel_board1->Draw("colz");

    c_VSstat_boards->cd(7);
    h_pdoVStdo_board0->Draw("colz");
    c_VSstat_boards->cd(8);
    h_pdoVStdo_board1->Draw("colz");

    c_VSstat_boards->cd(9);
    h_pdoVSbcid_board0->Draw("colz");
    c_VSstat_boards->cd(10);
    h_pdoVSbcid_board1->Draw("colz");


    //---------------------------- SINGLE CLUSTER ANALYSIS

    TCanvas *c_singleClusterStatistics_boards = new TCanvas("c_singleClusterStatistics_boards","c_singleClusterStatistics_boards",100,100,500,500);
    c_singleClusterStatistics_boards->Divide(2,5);

    c_singleClusterStatistics_boards->cd(1);
    h_nClu_board0->Draw();
    c_singleClusterStatistics_boards->cd(2);
    h_nClu_board1->Draw();

    c_singleClusterStatistics_boards->cd(3);
    h_clPos_board0->Draw();
    for(int iStrip=0;iStrip<cut_cluster_stripNeighbourBlacklist0->size();iStrip++)
    {
        l1 = new TLine(cut_cluster_stripNeighbourBlacklist0->at(iStrip)*0.4+.4,-0,cut_cluster_stripNeighbourBlacklist0->at(iStrip)*.4+.4,h_clPos_board0->GetMaximum());
        l1->SetLineWidth(5);
        l1->SetLineStyle(9);
        l1->Draw();
    }
    c_singleClusterStatistics_boards->cd(4);
    h_clPos_board1->Draw();
    for(int iStrip=0;iStrip<cut_cluster_stripNeighbourBlacklist1->size();iStrip++)
    {
        l1 = new TLine(cut_cluster_stripNeighbourBlacklist1->at(iStrip)*.4+.4,-0,cut_cluster_stripNeighbourBlacklist1->at(iStrip)*.4+.4,h_clPos_board1->GetMaximum());
        l1->SetLineWidth(5);
        l1->SetLineStyle(9);
        l1->Draw();
    }

    c_singleClusterStatistics_boards->cd(5);
    h_clCharge_board0->Draw();
    c_singleClusterStatistics_boards->cd(6);
    h_clCharge_board1->Draw();

    c_singleClusterStatistics_boards->cd(7);
    h_clSize_board0->Draw();
    c_singleClusterStatistics_boards->cd(8);
    h_clSize_board1->Draw();


    c_singleClusterStatistics_boards->cd(9);
    h_nCluVSnClu_boards->Draw("colz");
    c_singleClusterStatistics_boards->cd(10);
    fit_dgaus_sameMean(h_res);
    //h_res->Fit("gaus");

    TCanvas *tc=new TCanvas("resH8July2017","d",50,50,700,700);
    fit_dgaus_sameMean(h_res);


    //---------------------------- MULTI CLUSTER ANALYSIS

    TCanvas *c_multiClusterStatistics_boards = new TCanvas("c_multiClusterStatistics_boards","c_multiClusterStatistics_boards",100,150,500,500);
    c_multiClusterStatistics_boards->Divide(2,5);

    c_multiClusterStatistics_boards->cd(1);
    h_nClu_board0->Draw();
    c_multiClusterStatistics_boards->cd(2);
    h_nClu_board1->Draw();

    c_multiClusterStatistics_boards->cd(3);
    h_multiCluster_clPos_board0->Draw();
    for(int iStrip=0;iStrip<cut_cluster_stripNeighbourBlacklist0->size();iStrip++)
    {
        l1 = new TLine(cut_cluster_stripNeighbourBlacklist0->at(iStrip)*0.4+.4,-0,cut_cluster_stripNeighbourBlacklist0->at(iStrip)*.4+.4,h_clPos_board0->GetMaximum());
        l1->SetLineWidth(5);
        l1->SetLineStyle(9);
        l1->Draw();
    }
    c_multiClusterStatistics_boards->cd(4);
    h_multiCluster_clPos_board1->Draw();
    for(int iStrip=0;iStrip<cut_cluster_stripNeighbourBlacklist1->size();iStrip++)
    {
        l1 = new TLine(cut_cluster_stripNeighbourBlacklist1->at(iStrip)*.4+.4,-0,cut_cluster_stripNeighbourBlacklist1->at(iStrip)*.4+.4,h_clPos_board1->GetMaximum());
        l1->SetLineWidth(5);
        l1->SetLineStyle(9);
        l1->Draw();
    }

    c_multiClusterStatistics_boards->cd(5);
    h_multiCluster_clCharge_board0->Draw();
    c_multiClusterStatistics_boards->cd(6);
    h_multiCluster_clCharge_board1->Draw();

    c_multiClusterStatistics_boards->cd(7);
    h_multiCluster_clSize_board0->Draw();
    c_multiClusterStatistics_boards->cd(8);
    h_multiCluster_clSize_board1->Draw();


    c_multiClusterStatistics_boards->cd(9);
    h_nCluVSnClu_boards->Draw("colz");
    c_multiClusterStatistics_boards->cd(10);
    fit_dgaus(h_multiCluster_res);
    //    fit_dgaus_sameMean(h_multiCluster_res);

    TCanvas *tcMulti=new TCanvas("resH8July2017_Multiple","resH8July2017_Multiple",150,50,700,700);
    //    fit_dgaus(h_multiCluster_res);
    fit_dgaus_sameMean(h_multiCluster_res);

    TCanvas *c_Multi_diff_vs_pos = new TCanvas("c_Multi_diff_vs_pos","c_Multi_diff_vs_pos",250,10,600,600);
    c_Multi_diff_vs_pos->Divide(2,1);
    c_Multi_diff_vs_pos->cd(1);
    h_multiCluster_diff01_VS_pos0->Draw("colz");
    h_multiCluster_diff01_VS_pos0->Fit("pol1");
    c_Multi_diff_vs_pos->cd(2);
    h_multiCluster_diff01_VS_pos1->Draw("colz");

    TCanvas *c_angle = new TCanvas("c_angle","c_angle",40,40,400,400);
    h_multiCluster_diff_over_d01->Draw();

    TCanvas *c_eta = new TCanvas("c_eta","c_eta",60,70,600,600);
    c_eta->Divide(2,2);
    c_eta->cd(1);
    h_etaRatio_board0->Draw();
    c_eta->cd(2);
    h_etaRatio_board1->Draw();
    c_eta->cd(3);
    h_eta0_vs_diff->Draw("colz");
    c_eta->cd(4);
    h_eta1_vs_diff->Draw("colz");
    TCanvas *c_eta2 = new TCanvas("h_eta_vs_smallQ","h_eta_vs_smallQ",40,40,400,500);
    c_eta2->Divide(2,2);
    c_eta2->cd(1);
    h_eta_vs_smallQ->Draw("colz");
    c_eta2->cd(2);
    h_qL_vs_qR->Draw("colz");
    c_eta2->cd(3);
    h_qSmall_vs_qLarge->Draw("colz");
    //    c_eta->cd(6);
    //    h_qR->Draw();

    TCanvas *c_reverseBeamProfile = new TCanvas("c_reverseBeamProfile","c_reverseBeamProfile",60,60,600,600);
    c_reverseBeamProfile->Divide(2,1);
    c_reverseBeamProfile->cd(1);
    h_reverseBeamProfile_board0->Draw();
    c_reverseBeamProfile->cd(2);
    h_reverseBeamProfile_board1->Draw();

}

//the offset will be hardcoded, taken from the single cluster event run
//the 0 in the name, means that the clusterIndices returned are referring to the board0 clusters
void getClusterCombinationIndices_wSearchOffset(vector<double>* cl_poss_board0,
                                                vector<double>* cl_poss_board1,
                                                vector<int>& r_cl_indices_0,
                                                vector<int>& r_cl_indices_1
                                                )
{
    //we will loop over the cluster poss of board0
    //we check all other clusters for compatibility (diffs)

    //loop over board0 clusters
    for(int iCl0=0;iCl0<cl_poss_board0->size();iCl0++)
    {
        for(int iCl1=0;iCl1<cl_poss_board1->size();iCl1++)
        {
            double rawDiff = cl_poss_board0->at(iCl0)-cl_poss_board1->at(iCl1);
            //cout << "pos0 = "<<cl_poss_board0->at(iCl0)<<endl;
            //cout << "pos1 = "<<cl_poss_board1->at(iCl1)<<endl;
            //cout << "posDiff = "<<rawDiff<<endl;
            //cout << "|raw-off| = "<<rawDiff-g_MultiClusterSearchOffset01_mm<<endl;
            if(TMath::Abs(rawDiff-g_MultiClusterSearchOffset01_mm) > g_MultiClusterSearchEpsilon_mm)
            {
                //throw it away (do nothing)
            }
            else
            {
                //                cout << "#################################### KEEPING this pair\n";
                //they match within offset and epsilon, so we add both indices to the vectors
                r_cl_indices_0.push_back(iCl0);
                r_cl_indices_1.push_back(iCl1);
            }
        }//1 loop
    }//0 loop
}
int dsa=0;
vector<vector<int>*>* getClusterIndices(vector<int> *toCl_channel,
                                        vector<int> *toCl_pdo,
                                        vector<int> *toCl_tdo,
                                        vector<int> *toCl_bcid,
                                        int currentBoardId
                                        )
{
    vector <double> *test_pdo = new vector<double>();
    for(int i=0;i<toCl_pdo->size();i++)
    {
        test_pdo->push_back((double) toCl_pdo->at(i));
    }

    //try eta here:
    if(toCl_channel->size()==2 && false)
    {
        if(toCl_channel->at(0)==toCl_channel->at(1)-1)
        {
            double qL = toCl_pdo->at(0)*1.0;
            double qR = toCl_pdo->at(1)*1.0;
            double eta = qR/(qL+qR);
            if(currentBoardId==0)
                h_etaRatio_board0->Fill(eta);
            if(currentBoardId==1)
                h_etaRatio_board1->Fill(eta);
//            h_qLvsqR->Fill(qR,qL);
        }
        else if(toCl_channel->at(0)-1==toCl_channel->at(1))
            cout << "## Found non sorted\n";//never goes here
        else if(TMath::Abs(toCl_channel->at(0)-toCl_channel->at(1))==2)
        {
            cout <<"## clsize=2, but with 1 space #"<<dsa<<endl;
            dsa++;
        }
    }

    //try eta here:
    if(toCl_channel->size()==3 && false)
    {
        int i1=toCl_channel->at(0);
        int i2=toCl_channel->at(1);
        int i3=toCl_channel->at(2);

        if(i1==i2-1
                && i2==i3-1)
        {
            double qL = toCl_pdo->at(0)*1.0;
            double qR = toCl_pdo->at(2)*1.0;
            double qM = toCl_pdo->at(1)*1.0;

            if(qM>qL && qM>qR)
            {
                double eta;
                if(qL>qR)
                    eta = qM / (qM+qL);
                else
                    eta = qR/(qM+qR);


                if(currentBoardId==0)
                    h_etaRatio_board0->Fill(eta);
                if(currentBoardId==1)
                    h_etaRatio_board1->Fill(eta);
//                h_qLvsqR->Fill(qR,qL);
            }
        }
    }

    if(toCl_channel->size()==2 && false)
    {
        double eta;
        double shit;

        GiveEtaCorrection(toCl_channel,test_pdo,eta,shit);


        if(eta!=0)
        {
            if(currentBoardId==0)
                h_etaRatio_board0->Fill(eta);
            if(currentBoardId==1)
                h_etaRatio_board1->Fill(eta);
        }
        //        h_qLvsqR->Fill(qR,qL);
    }

    //this function returns a vector foreach cluster that is made.
    //each vector contains the corresponding indices from the input vector

    vector<int> *cutStripNeighbours;
    if(currentBoardId==0)
        cutStripNeighbours = cut_cluster_stripNeighbourBlacklist0;
    else if(currentBoardId==1)
        cutStripNeighbours = cut_cluster_stripNeighbourBlacklist1;

    if(eventByEvent)
        cout << "--- Clusters Board"<<currentBoardId<<endl;

    vector< vector < int >* > *clusters=new vector<vector<int>*>();
    vector<int> *newClusterIndices;
    int clusterCharge;
    int minBcid;
    int maxBcid;

    int newClusterTries=0;


    //sort the input vectors with ascending channel no
    for(size_t iStrip=0; iStrip<toCl_channel->size(); iStrip++)
        for(size_t jStrip=toCl_channel->size()-1; jStrip>iStrip; jStrip--)
            if(toCl_channel->at(jStrip-1)>toCl_channel->at(jStrip))
            {
                swap(toCl_channel->at(jStrip-1),toCl_channel->at(jStrip));
                swap(toCl_pdo->at(jStrip-1)    ,toCl_pdo->at(jStrip));
                swap(toCl_tdo->at(jStrip-1)    ,toCl_tdo->at(jStrip));
                swap(toCl_bcid->at(jStrip-1)   ,toCl_bcid->at(jStrip));
            }

    //now we are sorted in channel#

    bool startNewCluster=true;
    bool closeCluster=false;

    for(int i=0;i<toCl_channel->size();i++)
    {
        //        cout <<"--ON: "<<toCl_channel->at(i) <<endl;

        if(startNewCluster)
        {
            //just logging
            newClusterTries++;

            //            cout << "NEW Cluster\n";
            //start new cluster
            newClusterIndices = new vector<int>();
            clusterCharge=0;
            minBcid=4097;
            maxBcid=-1;
            startNewCluster=false;
        }


        //add the first available  strip
        //        cout <<"Adding: "<<toCl_channel->at(i)<<endl;
        newClusterIndices->push_back(i);
        //add the stripCharge to cl_charge
        clusterCharge+=toCl_pdo->at(i);
        if(toCl_bcid->at(i)>maxBcid)
            maxBcid=toCl_bcid->at(i);
        if(toCl_bcid->at(i)<minBcid)
            minBcid=toCl_bcid->at(i);

        //-----------------------IF LAST STRIP
        if(i==toCl_channel->size()-1)
            closeCluster=true;

        //-----------------------IF NEXT STRIP is too far
        else if(toCl_channel->at(i+1)-toCl_channel->at(i) - 1> cut_cluster_maxNEmptyStrips)
        {
            //            cout << "*** next strip too far\n";
            closeCluster = true;
        }

        if(closeCluster)
        {


            if(eventByEvent)
                cout << "cluster #"<< newClusterTries <<" | ";
            //------------------------------------- Cluster-wide Cuts
            bool clusterPassedCuts=true;
            int clusterTimeWidth = maxBcid-minBcid;
            int clusterSize = newClusterIndices->size();


//                        if(clusterSize==2)
//                        {
//                            double qL = toCl_pdo->at(newClusterIndices->at(0))*1.0;
//                            double qR = toCl_pdo->at(newClusterIndices->at(1))*1.0;

//                            //cout << "#### "<<newClusterIndices->at(0)<< " has q="<<qL<<" and "<<newClusterIndices->at(1)<<" has q="<<qR<<endl;

//                            double eta = qR/(qL+qR);

//                            if(currentBoardId==0)
//                                h_etaRatio_board0->Fill(eta);
//                            else if(currentBoardId==1)
//                                h_etaRatio_board1->Fill(eta);

//                            h_qL->Fill(qL);
//                            h_qR->Fill(qR);
//                            h_qLvsqR->Fill(qR,qL);
//                        }



            //check cluster charge and Size and TimeWidth(?)

            if(
                    clusterCharge > cut_cluster_maxClusterCharge
                    || clusterCharge < cut_cluster_minClusterCharge)
            {
                if(eventByEvent)
                    cout << "\t\t FAIL: ClusterCharge="<<clusterCharge;
                clusterPassedCuts=false;
            }

            if(clusterPassedCuts)
                if(
                        clusterSize   < cut_cluster_minClusterSize
                        || clusterSize   > cut_cluster_maxClusterSize)
                {
                    if(eventByEvent)
                        cout << "\t\t FAIL: ClusterSize="<<clusterSize;
                    clusterPassedCuts=false;
                }

            //             if(
            //                    clusterTimeWidth < cut_cluster_minClusterTimeWidth ||
            //                    clusterTimeWidth > cut_cluster_maxClusterTimeWidth
            //                    )
            //            {
            //                cout << "Cluster FAIL: ClusterTimeWidth\n";
            //                clusterPassedCuts=false;
            //            }

            //scan for bad strips
            if(clusterPassedCuts)
                for(int j=0;j<cutStripNeighbours->size();j++)
                {
                    int iStrip = cutStripNeighbours->at(j);
                    if(
                            toCl_channel->at(newClusterIndices->at(0))-1<=iStrip //a bad strip is 1st low neighbour or in the cluster
                            && toCl_channel->at(newClusterIndices->back())+1>=iStrip//a bad strip is 1st hi neghbour   or -- // --
                            )
                    {

                        if(eventByEvent)
                            cout << "\t\t FAIL: BadStrip: " <<
                                    (iStrip) <<
                                    " is between "<<
                                    (toCl_channel->at(newClusterIndices->at(0))-1)<<
                                    " and "<<
                                    (toCl_channel->at(newClusterIndices->back())+1)
                                    ;

                        clusterPassedCuts=false;
                        break;
                    }
                }


            //-----------------------ETA CORRECTION:

            if(clusterPassedCuts)
            {
                if(clusterSize==2)
                {
                    double eta = toCl_pdo->at(newClusterIndices->at(0))*1.0 / (toCl_pdo->at(newClusterIndices->at(0))*1.0 + toCl_pdo->at(newClusterIndices->at(1))*1.0 );

                    if(eta>cut_cluster_etaRatio_highCut)
                        clusterPassedCuts=false;
                    if(eta<cut_cluster_etaRatio_lowCut)
                        clusterPassedCuts=false;
                }
            }

            //-----------------------ETA CORRECTION: HISTO FILL
            if(clusterPassedCuts && true)
            {
                //fill the histo only in the first run (after which the flag 'filled' will goto 1)
                                if(clusterSize==2 && isEtaCorrectionHistoFilled==0)
                                {
                                    double qL = toCl_pdo->at(newClusterIndices->at(0))*1.0;
                                    double qR = toCl_pdo->at(newClusterIndices->at(1))*1.0;

                //                    cout << "#### "<<newClusterIndices->at(0)<< " has q="<<qL<<" and "<<newClusterIndices->at(1)<<" has q="<<qR<<endl;

                                    double eta = qR/(qL+qR);

                                    if(currentBoardId==0)
                                        h_etaRatio_board0->Fill(eta);
                                    else if(currentBoardId==1)
                                        h_etaRatio_board1->Fill(eta);

                                    h_qL_vs_qR->Fill(qR,qL);

                                    if(qL<qR)
                                    {
                                        h_eta_vs_smallQ->Fill(qL,eta);
                                        h_qSmall_vs_qLarge->Fill(qR,qL);
                                    }
                                    else
                                    {
                                        h_eta_vs_smallQ->Fill(qR,eta);
                                        h_qSmall_vs_qLarge->Fill(qL,qR);
                                    }
                                }
            }

            //-----IF ALL TESTS PASSED, ADD TO vector
            if(clusterPassedCuts)
            {//add the cluster
                clusters->push_back(newClusterIndices);

                //Printout
                if(eventByEvent)
                {

                    for(int k=0;k<clusters->back()->size();k++)
                    {
                        cout << toCl_channel->at(clusters->back()->at(k))<<" ";
                    }

                }//if ev by ev
            }

            if(eventByEvent)
                cout << endl;

            startNewCluster=true;
            closeCluster=false;
        }


    }//for

    return clusters;
}

void subEventLoop(int iEntry)
{
    i_stripPerEvent0=0;
    i_stripPerEvent1=0;

    v_toCl_channel_board0 = new vector<int>();
    v_toCl_pdo_board0     = new vector<int>();
    v_toCl_tdo_board0     = new vector<int>();
    v_toCl_bcid_board0    = new vector<int>();
    v_toCl_channel_board1 = new vector<int>();
    v_toCl_pdo_board1     = new vector<int>();
    v_toCl_tdo_board1     = new vector<int>();
    v_toCl_bcid_board1    = new vector<int>();

    tg_ev_tdo_board0 = new TGraph(64);
    tg_ev_tdo_board0->Set(0);
    tg_ev_tdo_board1 = new TGraph(64);
    tg_ev_tdo_board1->Set(0);

    //loop over the chips of the event
    for(int iSubEvent=0;iSubEvent < vmm_tree->chip->size();iSubEvent++)
    {

        int chipEventSize= vmm_tree->channel->at(iSubEvent).size();

        //Skip chip subevent if few strips
        if(chipEventSize<cut_strip_minPerChipEvent)
            continue;

        for(int iHit=0;iHit<chipEventSize;iHit++)
        {
            int channel = vmm_tree->channel->at(iSubEvent).at(iHit);
            int pdo = vmm_tree->pdo->at(iSubEvent).at(iHit);
            int tdo = vmm_tree->tdo->at(iSubEvent).at(iHit);
            int bcid = vmm_tree->bcid->at(iSubEvent).at(iHit);

            //Apply GLOBAL strip cuts
            if(
                    pdo < cut_strip_minPdo    ||
                    tdo < cut_strip_minTdo    ||
                    channel < cut_strip_minStrip ||
                    channel > cut_strip_maxStrip ||
                    find(cut_strip_noisySet->begin(), cut_strip_noisySet->end(), channel) != cut_strip_noisySet->end()
                    )
            {
                continue;
            }

            //Fill

            if(vmm_tree->boardId->at(iSubEvent)==g_boardId0)
            {

                h_bp_board0->Fill(channel);
                h_pdo_board0->Fill(pdo);
                h_tdo_board0->Fill(tdo);
                h_bcid_board0->Fill(bcid);
                i_stripPerEvent0++;

                h_pdoVSchannel_board0->Fill(channel,pdo);
                h_tdoVSchannel_board0->Fill(channel,tdo);
                h_bcidVSchannel_board0->Fill(channel,bcid);
                h_pdoVStdo_board0->Fill(tdo,pdo);
                h_pdoVSbcid_board0->Fill(bcid,pdo);

                if(eventByEvent)
                {
                    h_ev_pdo_board0->SetBinContent(channel+1,pdo);
                    tg_ev_tdo_board0->SetPoint(tg_ev_tdo_board0->GetN(),channel+1,tdo);
                    h_ev_tdo_board0->SetBinContent(channel+1,tdo);
                    h_ev_bcid_board0->SetBinContent(channel+1,bcid);
                }

                if(doClusterisation)
                {
                    //                    cout << "### channel = "<<channel<<endl;
                    v_toCl_channel_board0->push_back(channel);
                    v_toCl_pdo_board0->push_back(pdo);
                    v_toCl_tdo_board0->push_back(tdo);
                    v_toCl_bcid_board0->push_back(bcid);
                }
            }
            else if(vmm_tree->boardId->at(iSubEvent)==g_boardId1)
            {
                h_bp_board1->Fill(channel);
                h_pdo_board1->Fill(pdo);
                h_tdo_board1->Fill(tdo);
                h_bcid_board1->Fill(bcid);
                i_stripPerEvent1++;

                h_pdoVSchannel_board1->Fill(channel,pdo);
                h_tdoVSchannel_board1->Fill(channel,tdo);
                h_bcidVSchannel_board1->Fill(channel,bcid);
                h_pdoVStdo_board1->Fill(tdo,pdo);
                h_pdoVSbcid_board1->Fill(bcid,pdo);

                if(eventByEvent)
                {
                    h_ev_pdo_board1->SetBinContent(channel+1,pdo);
                    h_ev_tdo_board1->SetBinContent(channel+1,tdo);
                    tg_ev_tdo_board1->SetPoint(tg_ev_tdo_board1->GetN(),channel+1,tdo);
                    h_ev_bcid_board1->SetBinContent(channel+1,bcid);
                }
                if(doClusterisation)
                {
                    v_toCl_channel_board1->push_back(channel);
                    v_toCl_pdo_board1->push_back(pdo);
                    v_toCl_tdo_board1->push_back(tdo);
                    v_toCl_bcid_board1->push_back(bcid);
                }
            }


        }//hitLoop

    }//subEvent loop

}
void eventLoop_singleClusterFills()
{

    //--------try: Reverse Beam Profiles, to check for the pillars perhaps (normally needs tracking, but what the hey...)

    if(doSingleClusterEventAnalysis
            && cl_poss0->size()==1
            && cl_poss1->size()==0)
        h_reverseBeamProfile_board0->Fill(cl_poss0->at(0));
    if(doSingleClusterEventAnalysis
            && cl_poss0->size()==0
            && cl_poss1->size()==1)
        h_reverseBeamProfile_board1->Fill(cl_poss1->at(0));


    //---------NORMAL SINGLE CLUSTER EVENTS

    if(doSingleClusterEventAnalysis && cl_poss0->size()==1 && cl_poss1->size()==1)
    {
        h_clPos_board0->Fill(cl_poss0->at(0));
        h_clPos_board1->Fill(cl_poss1->at(0));

        h_clCharge_board0->Fill(cl_charges0->at(0));
        h_clCharge_board1->Fill(cl_charges1->at(0));

        h_clSize_board0->Fill(cl_sizes0->at(0));
        h_clSize_board1->Fill(cl_sizes1->at(0));

        h_res->Fill(cl_poss0->at(0)-cl_poss1->at(0));
    }
}
void eventLoop_multiClusterFills()
{

    if(doMultiClusterEventAnalysis &&
            cl_poss0->size()>=1 &&
            cl_poss1->size()>=1 &&
            cl_poss0->size()<= cut_event_maxNClustersPerBoard &&
            cl_poss1->size()<= cut_event_maxNClustersPerBoard
            )
    {
        //we need to know the pairs of indices to take the diffs of
        vector<int> goodClusters0;// = new vector<int>();
        vector<int> goodClusters1;// = new vector<int>();

        //            cout << "--- checking for\n";
        //            cout << "cl pos0: ";
        //            for(int k=0;k<cl_poss0->size();k++)
        //                cout << cl_poss0->at(k)<<" ";
        //            cout <<endl;
        //            cout << "cl pos1: ";
        //            for(int k=0;k<cl_poss1->size();k++)
        //                cout << cl_poss1->at(k)<<" ";
        //            cout <<endl;

        getClusterCombinationIndices_wSearchOffset(cl_poss0,cl_poss1,goodClusters0,goodClusters1);
        //            cout << "goodClusters0.size()="<<goodClusters0.size()<<endl;
        //now we have the good pairs = 2 vectors of same size

        for(int i=0;i<goodClusters0.size();i++)
        {
            //NOT plotting the clusters that are not paired...!!!

            h_multiCluster_clPos_board0   ->Fill(cl_poss0   ->at(goodClusters0.at(i)));
            h_multiCluster_clPos_board1   ->Fill(cl_poss1   ->at(goodClusters1.at(i)));
            h_multiCluster_clCharge_board0->Fill(cl_charges0->at(goodClusters0.at(i)));
            h_multiCluster_clCharge_board1->Fill(cl_charges1->at(goodClusters1.at(i)));
            h_multiCluster_clSize_board0  ->Fill(cl_sizes0  ->at(goodClusters0.at(i)));
            h_multiCluster_clSize_board1  ->Fill(cl_sizes1  ->at(goodClusters1.at(i)));

            double pos0=cl_poss0->at(goodClusters0.at(i));
            double pos1=cl_poss1->at(goodClusters1.at(i));
            double diff = pos0-pos1;
            double newDiff = diff - g_correctionSlope*pos0;
            h_multiCluster_res->Fill(newDiff-g_MultiClusterSearchOffset01_mm);


            h_multiCluster_diff01_VS_pos0->Fill(pos0,
                                                newDiff);
            h_multiCluster_diff01_VS_pos1->Fill(pos1,
                                                newDiff);


            h_multiCluster_diff_over_d01->Fill(newDiff/g_boardDistance_mm);

            if(doEtaCorrection)
            {
                h_eta0_vs_diff->Fill(newDiff,cl_etaRatio_values0->at(goodClusters0.at(i)));
                h_eta1_vs_diff->Fill(newDiff,cl_etaRatio_values1->at(goodClusters1.at(i)));
            }
        }

    }


}
void eventLoop_eventDisplay(int iEntry)
{

    if(eventByEvent)
    {
        bool isEmptyEvent = (h_ev_pdo_board0->GetEntries()==0 &&
                             h_ev_pdo_board1->GetEntries()==0 //&&
                             //                             h_ev_tdo_board0->GetEntries()==0 &&
                             //                             h_ev_tdo_board1->GetEntries()==0
                             );

        if(!isEmptyEvent)
        {
            c_ev->Clear();
            c_ev->Divide(3,3);
            c_ev->cd(1);
            h_ev_pdo_board0->Draw();
            for(int iStrip=0;iStrip<cut_cluster_stripNeighbourBlacklist0->size();iStrip++)
            {
                l1 = new TLine(cut_cluster_stripNeighbourBlacklist0->at(iStrip)+.4,-0,cut_cluster_stripNeighbourBlacklist0->at(iStrip)+.4,h_ev_pdo_board0->GetMaximum());
                l1->SetLineWidth(5);
                l1->SetLineStyle(9);
                l1->Draw();
            }
            c_ev->cd(2);
//            h_ev_tdo_board0->Draw();
            tg_ev_tdo_board0->Draw("AL*");
//            tg_ev_tdo_board0->GetXaxis()->SetRange(0,64);
            c_ev->cd(3);
            gPad->SetLogy();
            h_ev_bcid_board0->Draw();

            c_ev->cd(4);
            h_ev_pdo_board1->Draw();
            for(int iStrip=0;iStrip<cut_cluster_stripNeighbourBlacklist1->size();iStrip++)
            {
                l1 = new TLine(cut_cluster_stripNeighbourBlacklist1->at(iStrip)+.4,-0,cut_cluster_stripNeighbourBlacklist1->at(iStrip)+.4,h_ev_pdo_board0->GetMaximum());
                l1->SetLineWidth(5);
                l1->SetLineStyle(9);
                l1->Draw();
            }
            c_ev->cd(5);
//            h_ev_tdo_board1->Draw();
            tg_ev_tdo_board1->Draw("AL*");
            c_ev->cd(6);
            gPad->SetLogy();
            h_ev_bcid_board1->Draw();

            c_ev->cd(7);
            if(h_ev_pdo_board0->GetMaximum()<h_ev_pdo_board1->GetMaximum())
            {
                //                h_ev_pdo_board0->SetMaximum(h_ev_pdo_board1->GetMaximum()*1.1);
                //                h_ev_pdo_board1->SetMaximum(h_ev_pdo_board1->GetMaximum()*1.1);
            }
            h_ev_pdo_board0->Draw();
            h_ev_pdo_board1->Draw("sames");
            //            for(int iStrip=0;iStrip<cut_cluster_stripNeighbourBlacklist->size();iStrip++)
            //            {
            //                l1 = new TLine(cut_cluster_stripNeighbourBlacklist->at(iStrip)+.4,-0,cut_cluster_stripNeighbourBlacklist->at(iStrip)+.4,h_ev_pdo_board0->GetMaximum());
            //                l1->SetLineWidth(5);
            //                l1->SetLineStyle(9);
            //                l1->Draw();
            //            }
            c_ev->cd(8);
            h_ev_tdo_board0->Draw();
            h_ev_tdo_board1->Draw("sames");
            c_ev->cd(9);
            gPad->SetLogy();
            h_ev_bcid_board0->Draw();
            h_ev_bcid_board1->Draw("sames");

            c_ev->SetTitle(Form("trigCount = %i",iEntry));
            //        c_ev->Resize();
            //        c_ev->Modified();
            c_ev->Update();
            c_ev->WaitPrimitive();
            //        getline(cin, cmdInput);

        }//if ev by ev

        //reset eventByEvent histos anyway
        h_ev_pdo_board0->Reset();
        h_ev_pdo_board1->Reset();
        h_ev_tdo_board0->Reset();
        h_ev_tdo_board1->Reset();
        h_ev_bcid_board0->Reset();
        h_ev_bcid_board1->Reset();
    }
}

void eventLoop(int iEntry)
{
    h_stripPerEv_board0->Fill(i_stripPerEvent0);
    h_stripPerEv_board1->Fill(i_stripPerEvent1);

    //now we went through the strips
    //time to clusterize

    //board by board:

    //DEBUGGING!!!
    /*
        cout << "\n-------Debug Event: "<<iEntry<<endl;

        cout << "                                 Channels0: ";
        for(int i=0;i<v_toCl_channel_board0->size();i++)
            cout << v_toCl_channel_board0->at(i) <<" ";
        cout <<endl;
        cout << "                                 Channels1: ";
        for(int i=0;i<v_toCl_channel_board1->size();i++)
            cout << v_toCl_channel_board1->at(i) <<" ";
        cout <<endl;
*/

    vector<vector<int>*> *clusterIndices_board0 = getClusterIndices(
                v_toCl_channel_board0,
                v_toCl_pdo_board0,
                v_toCl_tdo_board0,
                v_toCl_bcid_board0,
                0
                );

    vector<vector<int>*> *clusterIndices_board1 = getClusterIndices(
                v_toCl_channel_board1,
                v_toCl_pdo_board1,
                v_toCl_tdo_board1,
                v_toCl_bcid_board1,
                1
                );

    h_nClu_board0->Fill(clusterIndices_board0->size());
    h_nClu_board1->Fill(clusterIndices_board1->size());
    h_nCluVSnClu_boards->Fill(clusterIndices_board0->size(),clusterIndices_board1->size());

    //####################################################
    //############## CLUSTER-WIDE CUTS
    //####################################################

    //first we fill these vectors, for any kind of event with clusters
    cl_sizes0   = new vector<double>();
    cl_sizes1   = new vector<double>();
    cl_charges0 = new vector<double>();
    cl_charges1 = new vector<double>();
    cl_poss0    = new vector<double>();
    cl_poss1    = new vector<double>();
    cl_etaRatio_values0    = new vector<double>();
    cl_etaRatio_values1    = new vector<double>();

    //---------CALCULATE CLUSTER POSITIONS
    for(int iCluster=0;iCluster<clusterIndices_board0->size();iCluster++)
    {
        double size=clusterIndices_board0->at(iCluster)->size();
        double charge=0;
        double over=0;

        for(int i=0;i<clusterIndices_board0->at(iCluster)->size();i++)
        {
            int index = clusterIndices_board0->at(iCluster)->at(i);
            charge+= v_toCl_pdo_board0->at(index)-g_PdoPedestal;
            over+= v_toCl_channel_board0->at(index) * (v_toCl_pdo_board0->at(index)-g_PdoPedestal);
        }
        double pos = 0.4 * over / charge;

        cl_sizes0->push_back(size);
        cl_charges0->push_back(charge);
        cl_poss0->push_back(pos);
    }
    for(int iCluster=0;iCluster<clusterIndices_board1->size();iCluster++)
    {
        double size=clusterIndices_board1->at(iCluster)->size();
        double charge=0;
        double over=0;

        for(int i=0;i<clusterIndices_board1->at(iCluster)->size();i++)
        {
            int index = clusterIndices_board1->at(iCluster)->at(i);
            charge+= v_toCl_pdo_board1->at(index)-g_PdoPedestal;
            over+= v_toCl_channel_board1->at(index) * (v_toCl_pdo_board1->at(index)-g_PdoPedestal);
        }
        double pos = 0.4 * over / charge;

        cl_sizes1->push_back(size);
        cl_charges1->push_back(charge);
        cl_poss1->push_back(pos);
    }

    //------ check clsize==2, to recalculate the pos
    //------------ETA CORRECTION (POS CORRECTING!!!)

    if(doEtaCorrection && isEtaCorrectionHistoFilled==1)
    {
        for(int iCluster=0;iCluster<clusterIndices_board0->size();iCluster++)
        {
            double size=clusterIndices_board0->at(iCluster)->size();
            if(size>=2)
            {
                double chLeft = v_toCl_pdo_board0->at(clusterIndices_board0->at(iCluster)->at(0));
                double chRight = v_toCl_pdo_board0->at(clusterIndices_board0->at(iCluster)->at(1));
                double P=chLeft/(chLeft+chRight);
                double N0 = h_etaRatio_board0->GetEntries();

                double etaIntegral;
                {
                    double xmin = 0;
                    double xmax=P;
                    TAxis *axis = h_etaRatio_board0->GetXaxis();
                    int bmin = axis->FindBin(xmin); //in your case xmin=-1.5
                    int bmax = axis->FindBin(xmax); //in your case xmax=0.8
                    etaIntegral = h_etaRatio_board0->Integral(bmin,bmax);
                    etaIntegral-= h_etaRatio_board0->GetBinContent(bmin) *  ( xmin-axis->GetBinLowEdge(bmin))    /axis->GetBinWidth(bmin);
                    etaIntegral-= h_etaRatio_board0->GetBinContent(bmax) *  ( axis->GetBinUpEdge(bmax)-xmax )    /axis->GetBinWidth(bmax);
                }
                double eta_cor_cluster_pos = P / N0 * etaIntegral;
                cl_etaRatio_values0->push_back(P);
                cl_poss0->at(iCluster)=eta_cor_cluster_pos;

            }
        }//0
        for(int iCluster=0;iCluster<clusterIndices_board1->size();iCluster++)
        {
            double size=clusterIndices_board1->at(iCluster)->size();
            if(size>=2)
            {
                double chLeft = v_toCl_pdo_board1->at(clusterIndices_board1->at(iCluster)->at(0));
                double chRight = v_toCl_pdo_board1->at(clusterIndices_board1->at(iCluster)->at(1));
                double P=chLeft/(chLeft+chRight);
                double N0 = h_etaRatio_board1->GetEntries();

                double etaIntegral;
                {
                    double xmin = 0;
                    double xmax=P;
                    TAxis *axis = h_etaRatio_board1->GetXaxis();
                    int bmin = axis->FindBin(xmin); //in your case xmin=-1.5
                    int bmax = axis->FindBin(xmax); //in your case xmax=0.8
                    etaIntegral = h_etaRatio_board1->Integral(bmin,bmax);
                    etaIntegral-= h_etaRatio_board1->GetBinContent(bmin) *  ( xmin-axis->GetBinLowEdge(bmin))    /axis->GetBinWidth(bmin);
                    etaIntegral-= h_etaRatio_board1->GetBinContent(bmax) *  ( axis->GetBinUpEdge(bmax)-xmax )    /axis->GetBinWidth(bmax);
                }
                double eta_cor_cluster_pos = P / N0 * etaIntegral;
                cl_etaRatio_values1->push_back(P);
                cl_poss1->at(iCluster)=eta_cor_cluster_pos;
            }
        }//1
    }//    if(doEtaCorrection && isEtaCorrectionHistoFilled==1)

    if(
            (doEtaCorrection && isEtaCorrectionHistoFilled)
            ||
            (!doEtaCorrection && !isEtaCorrectionHistoFilled)
            )
    {
        //now the vectors are filled.
        //so now we can easily check for single or multi cluster events.

        eventLoop_singleClusterFills();
        eventLoop_multiClusterFills();
        eventLoop_eventDisplay(iEntry);
    }
}
void run(TString file, bool doEventByEvent=false)
{
    if(g_PdoPedestal>cut_strip_minPdo)
    {
        cout << "cut_strip_minPdo="<<cut_strip_minPdo<<" is lower than g_PdoPedestal="<<g_PdoPedestal<<" ...exiting\n";
        exit(0);
    }
    //load file
    load_tree_objects(file);

    initialize();
    eventByEvent = doEventByEvent;

    size_t nEntries = vmm_tree->fChain->GetEntries();
    cout << "run#"<<sRunNumber<<" | Entries = " << nEntries << endl;

    if(eventByEvent)
        c_ev = new TCanvas("c_event","c_event",100,100,700,900);

    double timesToRun=1;
    if(doEtaCorrection)
        timesToRun=2;

    //just to stupidly force it to run again, to first fill the eta histos
    //and then to rerun, in order to correct the 2-strip cluster positions

    for(int i=0;i<timesToRun;i++)
    {
        for(int iEntry=0;iEntry<nEntries;iEntry++)
        {
            vmm_tree->GetEntry(iEntry);
            printPercentage(iEntry,nEntries);

            if(eventByEvent)
                cout << "\n--------------EVENT "<<iEntry<<" ------\n";
            //foreach entry (event) we loop over the subEvents
            subEventLoop(iEntry);

            eventLoop(iEntry);
        }

        if(doEtaCorrection)
            isEtaCorrectionHistoFilled=true;
    }

    cout <<endl;

    drawHistos();

    //    cout << "DONE (singleCluster): ADD -32 FOR PDSTAL\n";
    //    cout << "DONE (multiCluster): ADD -32 FOR PDSTAL\n";
    //    cout << "TODO: dgaus without mean\n";
    //    cout << "DONE: bad neighobur strips, split for both boards";
    //    cout << "DONE: multiCluster events";
    //    cout << "DONE: add angle plot (Dx/d) \n";
    cout << "TODO: reverse beam profiles for TLP etc\n";
}
/*
 *
        //print clusters
        if(eventByEvent && false)
        {
            //            cout << "== board"<<g_boardId0<<" has "<<clusterIndices_board0->size()<<" clusters"<<endl;
            //            cout << "== board"<<g_boardId1<<" has "<<clusterIndices_board1->size()<<" clusters"<<endl;


            cout << "== clusters - board"<<g_boardId0<<endl;
            for(int iCl=0;iCl<clusterIndices_board0->size();iCl++)
            {
                cout << "cluster #"<<iCl<<"  |  ";

                for(int iChan=0;iChan<clusterIndices_board0->at(iCl)->size();iChan++)
                {
                    cout << v_toCl_channel_board0->at(clusterIndices_board0->at(iCl)->at(iChan)) << " ";
                }
                cout <<endl;
            }
            cout << "== clusters - board"<<g_boardId1<<endl;

            for(int iCl=0;iCl<clusterIndices_board1->size();iCl++)
            {
                cout << "cluster #"<<iCl<<"  |  ";

                for(int iChan=0;iChan<clusterIndices_board1->at(iCl)->size();iChan++)
                {
                    cout << v_toCl_channel_board1->at(clusterIndices_board1->at(iCl)->at(iChan)) << " ";
                }
                cout <<endl;
            }
        }

*/
/*
        if(
                clusterIndices_board0->size()==1
                && clusterIndices_board1->size()==1
                && doSingleClusterEventAnalysis
                )
        {
            double size0=clusterIndices_board0->at(0)->size();
            double size1=clusterIndices_board1->at(0)->size();
            double charge0=0;
            double charge1=0;

            double over0=0;
            double over1=0;

            for(int i=0;i<clusterIndices_board0->at(0)->size();i++)
            {
                int index = clusterIndices_board0->at(0)->at(i);
                charge0+= v_toCl_pdo_board0->at(index)-g_PdoPedestal;
                over0+= v_toCl_channel_board0->at(index) * (v_toCl_pdo_board0->at(index)-g_PdoPedestal);
                //                if(eventByEvent)
                //                    cout << "PDO_0 = "<<v_toCl_pdo_board0->at(index)<<endl;
            }
            for(int i=0;i<clusterIndices_board1->at(0)->size();i++)
            {
                int index = clusterIndices_board1->at(0)->at(i);
                charge1+= v_toCl_pdo_board1->at(index)-g_PdoPedestal;
                over1+= v_toCl_channel_board1->at(index) * (v_toCl_pdo_board1->at(index)-g_PdoPedestal);
                //                if(eventByEvent)
                //                    cout << "PDO_1 = "<<v_toCl_pdo_board1->at(index)<<endl;
            }

            double pos0 = 0.4 * over0 / charge0;
            //            if(eventByEvent)
            //                cout << "POS0    --   "<<pos0<< "  |  over="<<over0<<"  charge0="<<charge0<<endl;
            double pos1 = 0.4 * over1 / charge1;
            //            if(eventByEvent)
            //                cout << "POS1    --   "<<pos1<< "  |  over="<<over1<<"  charge1="<<charge1<<endl;

            h_clPos_board0->Fill(pos0);
            h_clPos_board1->Fill(pos1);

            h_clCharge_board0->Fill(charge0);
            h_clCharge_board1->Fill(charge1);

            h_clSize_board0->Fill(size0);
            h_clSize_board1->Fill(size1);

            h_res->Fill(pos0-pos1);
        }//if SINGLE CLUSTER EVENTS oneach board

*/
/*
        //---------------MULTI-CLUSTER EVENTS (with nStrips>1)
//        if(false && doMultiClusterEventAnalysis)
//        {
//            //we will fill the resolution histo(s), but carefully combining clusters from both boards

//            //for this it helps to just make vectors with the strip numbers, and not indices
//            //board0
//            vector<vector<int>*>* clusterStrips_board0 = new vector<vector<int>*>();
//            vector<vector<int>*>* clusterPdo_board0 = new vector<vector<int>*>();
//            for(int iCl=0;iCl<clusterIndices_board0->size();iCl++)
//            {
//                vector<int>* thisClStrips = new vector<int>();
//                vector<int>* thisClPdo = new vector<int>();
//                for(int jStrip=0;jStrip<clusterIndices_board0->at(iCl)->size();jStrip++)
//                {
//                    thisClStrips->push_back(v_toCl_channel_board0->at(clusterIndices_board0->at(iCl)->at(jStrip)));
//                    thisClPdo->push_back(v_toCl_pdo_board0->at(clusterIndices_board0->at(iCl)->at(jStrip)));
//                }

//                clusterStrips_board0->push_back(thisClStrips);
//                clusterPdo_board0->push_back(thisClPdo);
//            }//board1
//            vector<vector<int>*>* clusterStrips_board1 = new vector<vector<int>*>();
//            vector<vector<int>*>* clusterPdo_board1 = new vector<vector<int>*>();
//            for(int iCl=0;iCl<clusterIndices_board1->size();iCl++)
//            {
//                vector<int>* thisClStrips = new vector<int>();
//                vector<int>* thisClPdo = new vector<int>();
//                for(int jStrip=0;jStrip<clusterIndices_board1->at(iCl)->size();jStrip++)
//                {
//                    thisClStrips->push_back(v_toCl_channel_board1->at(clusterIndices_board1->at(iCl)->at(jStrip)));
//                    thisClPdo->push_back(v_toCl_pdo_board1->at(clusterIndices_board1->at(iCl)->at(jStrip)));
//                }

//                clusterStrips_board1->push_back(thisClStrips);
//                clusterPdo_board1->push_back(thisClPdo);
//            }

//            //find cluster positions, and put them in a vector
//            vector<int> v_i_multiClPos0;
//            vector<int> v_i_multiClPos1;
//            for(int indexCluster=0;indexCluster<clusterStrips_board0.size();indexCluster++)
//            {
//                double size0=clusterStrips_board0->at(indexCluster)->size();
//                double size1=clusterStrips_board1->at(indexCluster)->size();
//                double charge0=0;
//                double charge1=0;
//                double over0=0;
//                double over1=0;

//                for(int iStrip=0;iStrip<clusterStrips_board0->at(indexCluster)->size();iStrip++)
//                {
//                    //                    int index = clusterIndices_board0->at(0)->at(i);
//                    charge0+= clusterPdo_board0->at(indexCluster)->at(iStrip)-g_PdoPedestal;
//                    over0+= clusterStrips_board0->at(indexCluster)->at(iStrip) * (clusterPdo_board0->at(indexCluster)->at(iStrip)-g_PdoPedestal);
//                }
//                double pos0 = 0.4 * over0 / charge0;

//                for(int iStrip=0;iStrip<clusterStrips_board1->at(indexCluster)->size();iStrip++)
//                {
//                    //                    int index = clusterIndices_board1->at(0)->at(i);
//                    charge1+= clusterPdo_board1->at(indexCluster)->at(iStrip)-g_PdoPedestal;
//                    over1+= clusterStrips_board1->at(indexCluster)->at(iStrip) * (clusterPdo_board1->at(indexCluster)->at(iStrip)-g_PdoPedestal);
//                }
//                double pos1 = 0.4 * over1 / charge1;
//            }
//            //check the 0+1 pos vectors, to match them





//            //OMG, this is implemented so badly!!!
//            //-------------------Find good CLUSTER PAIRS -----------------------
//            vector<int> comb_clIndices_0;
//            vector<int> comb_clIndices_1;
//            getClusterCombinationIndices_wSearchOffset(clusterStrips_board0,
//                                                       clusterStrips_board1,
//                                                       clusterPdo_board0,
//                                                       clusterPdo_board1,
//                                                       &comb_clIndices_0,
//                                                       &comb_clIndices_1
//                                                       );



//            //now we know which clusters to use, so we run again the pos-finder etc
//            //comb_clIndices_0 and 1 have same size, by definition. They contain the indices of which clusters are pairs inside them
//            //in reference to clusterStrips_board[0,1]

//            //se  we run to find pos diffs, and charges foreach pair of clusters we found
//            for(int indexCluster=0;indexCluster<comb_clIndices_0.size();indexCluster++)
//            {
//                double size0=clusterStrips_board0->at(indexCluster)->size();
//                double size1=clusterStrips_board1->at(indexCluster)->size();
//                double charge0=0;
//                double charge1=0;
//                double over0=0;
//                double over1=0;

//                for(int iStrip=0;iStrip<clusterStrips_board0->at(indexCluster)->size();iStrip++)
//                {
//                    //                    int index = clusterIndices_board0->at(0)->at(i);
//                    charge0+= clusterPdo_board0->at(indexCluster)->at(iStrip)-g_PdoPedestal;
//                    over0+= clusterStrips_board0->at(indexCluster)->at(iStrip) * (clusterPdo_board0->at(indexCluster)->at(iStrip)-g_PdoPedestal);
//                }
//                double pos0 = 0.4 * over0 / charge0;

//                for(int iStrip=0;iStrip<clusterStrips_board1->at(indexCluster)->size();iStrip++)
//                {
//                    //                    int index = clusterIndices_board1->at(0)->at(i);
//                    charge1+= clusterPdo_board1->at(indexCluster)->at(iStrip)-g_PdoPedestal;
//                    over1+= clusterStrips_board1->at(indexCluster)->at(iStrip) * (clusterPdo_board1->at(indexCluster)->at(iStrip)-g_PdoPedestal);
//                }
//                double pos1 = 0.4 * over1 / charge1;

//                //                //here we have both positions of the current cluster pair

//                //                h_multiCluster_clPos_board0->Fill(pos0);
//                //                h_multiCluster_clPos_board1->Fill(pos1);

//                //                h_multiCluster_clCharge_board0->Fill(charge0);
//                //                h_multiCluster_clCharge_board1->Fill(charge1);
//                //                h_multiCluster_clSize_board0->Fill(size0);
//                //                h_multiCluster_clSize_board1->Fill(size1);

//                //                h_multiCluster_res->Fill(pos0-pos1);


//            }



//        }//doMultiCluster*/
