//#include "vmm.C"
#include "emiltoolsH82017_VMM3.cpp"
//######## SETUP #########
bool drawRawChips = false;

//vector<int> *v_boardIds;
int g_boardId0=0;
int g_boardId1=1;

//########################
//--------------------------------------GLOBAL strip cuts config
int cut_strip_minPerChipEvent=2;

//PDO TDO cuts
int cut_strip_minPdo=60;
int cut_strip_minTdo=20;

//edge cuts
int cut_strip_hiLim= 63;
int cut_strip_loLim= 0;

vector<int> *cut_strip_noisySet = new vector<int>();

//--------------------------------------CLUSTER CUTS config
int cut_cluster_maxNEmptyStrips=2;
int cut_cluster_maxClusterCharge=12000;
int cut_cluster_minClusterCharge=50;
int cut_cluster_maxClusterTimeWidth=50;
int cut_cluster_minClusterTimeWidth=-1;

//########################

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



//--------------Event by event histos
TCanvas* c_ev;

TH1D* h_ev_pdo_board0 = new TH1D("h_ev_pdo_board0","h_ev_pdo_board0",64,0,64);
TH1D* h_ev_tdo_board0 = new TH1D("h_ev_tdo_board0","h_ev_tdo_board0",64,0,64);

TH1D* h_ev_pdo_board1 = new TH1D("h_ev_pdo_board1","h_ev_pdo_board1",64,0,64);
TH1D* h_ev_tdo_board1 = new TH1D("h_ev_tdo_board1","h_ev_tdo_board1",64,0,64);


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

//############################################################
string cmdInput="";
TString hname;
TString sRunNumber;

void initialize()
{

    //    v_boardIds->push_back(0);
    //    v_boardIds->push_back(1);

    cut_strip_noisySet->push_back(63);

    gStyle->SetStatH(0.5);
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


}

vector<vector<int>*>* getClusterIndices(vector<int> *toCl_channel,
                                        vector<int> *toCl_pdo,
                                        vector<int> *toCl_tdo,
                                        vector<int> *toCl_bcid
                                        )
{
    //this function returns a vector foreach cluster that is made.
    //each vector contains the corresponding indices from the input vector

    vector< vector < int >* > *clusters=new vector<vector<int>*>();
    vector<int> *newCluster;
    int clusterCharge;
    int minBcid;
    int maxBcid;


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
        cout <<"--ON: "<<toCl_channel->at(i) <<endl;

        if(startNewCluster)
        {
            cout << "NEW Cluster\n";
            //start new cluster
            newCluster = new vector<int>();
            clusterCharge=0;
            minBcid=4097;
            maxBcid=-1;
            startNewCluster=false;
        }


        //add the first available  strip
        cout <<"Adding: "<<toCl_channel->at(i)<<endl;
        newCluster->push_back(i);
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
            cout << "*** next strip too far\n";
            closeCluster = true;
        }

        if(closeCluster)
        {
            //Check Cluster-wide Cuts
            int clusterTimeWidth = maxBcid-minBcid;
            if(
                    clusterCharge > cut_cluster_maxClusterCharge ||
                    clusterCharge < cut_cluster_minClusterCharge
//                    clusterTimeWidth < cut_cluster_minClusterTimeWidth ||
//                    clusterTimeWidth > cut_cluster_maxClusterTimeWidth
                    )
            {
                cout << "*** CLUSTER REJECTED from CUTS #############################\n";
            }
            else
            {//add the cluster
                cout << "$$$ closing cluster\n";
                clusters->push_back(newCluster);
            }
            startNewCluster=true;
            closeCluster=false;
        }


    }//for

    return clusters;
}

void subEventLoop(int iEntry, bool eventByEvent, bool doClusterisation)
{
    i_stripPerEvent0=0;
    i_stripPerEvent1=0;


    //    v_cl_boardId->clear();
    if(doClusterisation)
    {
        v_toCl_channel_board0 = new vector<int>();
        v_toCl_pdo_board0     = new vector<int>();
        v_toCl_tdo_board0     = new vector<int>();
        v_toCl_bcid_board0    = new vector<int>();
        v_toCl_channel_board1 = new vector<int>();
        v_toCl_pdo_board1     = new vector<int>();
        v_toCl_tdo_board1     = new vector<int>();
        v_toCl_bcid_board1    = new vector<int>();
    }

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
                    channel < cut_strip_loLim ||
                    channel > cut_strip_hiLim ||
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
                    h_ev_tdo_board0->SetBinContent(channel+1,tdo);
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

    h_stripPerEv_board0->Fill(i_stripPerEvent0);
    h_stripPerEv_board1->Fill(i_stripPerEvent1);


    if(doClusterisation)
    {
        //now we went through the strips
        //time to clusterize

        //board by board:

        //DEBUGGING!!!

        cout << "\n-------Debug Event: "<<iEntry<<endl;

        cout << "                                 Channels0: ";
        for(int i=0;i<v_toCl_channel_board0->size();i++)
            cout << v_toCl_channel_board0->at(i) <<" ";
        cout <<endl;
        cout << "                                 Channels1: ";
        for(int i=0;i<v_toCl_channel_board1->size();i++)
            cout << v_toCl_channel_board1->at(i) <<" ";
        cout <<endl;



        //        cout << "\ngofor Clusterization \n";
        vector<vector<int>*> *clusterIndices_board0 = getClusterIndices(
                    v_toCl_channel_board0,
                    v_toCl_pdo_board0,
                    v_toCl_tdo_board0,
                    v_toCl_bcid_board0
                    );

        //        cout << "DONEDONEDONE \n";
        vector<vector<int>*> *clusterIndices_board1 = getClusterIndices(
                    v_toCl_channel_board1,
                    v_toCl_pdo_board1,
                    v_toCl_tdo_board1,
                    v_toCl_bcid_board1
                    );



        cout << "== board"<<g_boardId0<<" has "<<clusterIndices_board0->size()<<" clusters"<<endl;
        cout << "== board"<<g_boardId1<<" has "<<clusterIndices_board1->size()<<" clusters"<<endl;

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

    }//doClusterization


    if(eventByEvent)
    {
        bool isEmptyEvent = (h_ev_pdo_board0->GetEntries()==0 &&
                             h_ev_pdo_board1->GetEntries()==0 &&
                             h_ev_tdo_board0->GetEntries()==0 &&
                             h_ev_tdo_board1->GetEntries()==0
                             );

        if(!isEmptyEvent)
        {
            c_ev->Clear();
            c_ev->Divide(2,2);
            c_ev->cd(1);
            h_ev_pdo_board0->Draw();
            //        TLine *l1 = new TLine(cut_strip_loLim,-0,cut_strip_loLim,h_ev_pdo_board0->GetMaximum());
            //        l1->Draw();
            //        TLine *l2 = new TLine(cut_strip_hiLim,-0,cut_strip_hiLim,h_ev_pdo_board0->GetMaximum());
            //        l2->Draw();
            c_ev->cd(2);
            h_ev_tdo_board0->Draw();
            c_ev->cd(3);
            h_ev_pdo_board1->Draw();
            c_ev->cd(4);
            h_ev_tdo_board1->Draw();

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
    }



    //after clusterization, we can analyze them, inside the event.

}

void run(TString file, bool eventByEvent=false, bool doClusterisation=false)
{
    //load file
    load_tree_objects(file);

    initialize();

    size_t nEntries = vmm_tree->fChain->GetEntries();
    cout << "run#"<<sRunNumber<<" | Entries = " << nEntries << endl;

    if(eventByEvent)
        c_ev = new TCanvas("c_event","c_event",500,200,700,900);

    for(int iEntry=0;iEntry<nEntries;iEntry++)
    {
        vmm_tree->GetEntry(iEntry);
        printPercentage(iEntry,nEntries);

        //foreach entry (event) we loop over the subEvents
        subEventLoop(iEntry,eventByEvent,doClusterisation);

        //        cout <<endl;

    }//entry loop
    cout <<endl;

    drawHistos();

}

