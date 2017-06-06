#include "vmm3_clu.h"
#include "emiltools.cpp"
//######## SETUP ################
int nBoards=2;
int nChipsPerBoard=1;

bool drawRawChips = false;
bool drawRawBoards = true;

int stripsPerChip = 64;
int stripPitch = 0.4;

//########################
vector<vector<TH1D*>*> *h_bp_chips;
vector<vector<TH1D*>*> *h_charge_chips;
vector<vector<TH1D*>*> *h_time_chips;
vector<vector<TH1D*>*> *h_nclu_chips;
vector<vector<TH1D*>*> *h_clSize_chips;

vector<TH1D*> *h_bp_boards;
vector<TH1D*> *h_charge_boards;
vector<TH1D*> *h_time_boards;
vector<TH1D*> *h_nclu_boards;
vector<TH1D*> *h_clSize_boards;


TH1D* test1 = new TH1D("test1","test1",65,0,65);
TH1D* test2 = new TH1D("test2","test2",20,-10,100);

int mysize=0;
TString hname;

vmm3_clu *t;//this is the object from which we get the tree
void load_tree_objects(TString file)
{
    //Get file
    TFile *root_file = new TFile(file.Data(),"READ");
    //Load all cluster trees from file
    t  = new vmm3_clu((TTree*)root_file->Get("vmm3_clu"));
}
void setupHistos()
{
    h_bp_chips      = new vector<vector<TH1D*>*>;
    h_charge_chips  = new vector<vector<TH1D*>*>;
    h_time_chips    = new vector<vector<TH1D*>*>;
    h_nclu_chips    = new vector<vector<TH1D*>*>;
    h_clSize_chips  = new vector<vector<TH1D*>*>;

    h_bp_boards     = new vector<TH1D*>;
    h_charge_boards = new vector<TH1D*>;
    h_time_boards   = new vector<TH1D*>;
    h_nclu_boards   = new vector<TH1D*>;
    h_clSize_boards = new vector<TH1D*>;

    for(int iboard=0;iboard<nBoards;iboard++)
    {
        h_bp_chips        ->push_back(new vector<TH1D*>());
        h_charge_chips       ->push_back(new vector<TH1D*>());
        h_time_chips       ->push_back(new vector<TH1D*>());
        h_nclu_chips      ->push_back(new vector<TH1D*>());
        h_clSize_chips ->push_back(new vector<TH1D*>());

        hname = Form("Beam Profile board%i",iboard);
        h_bp_boards        ->push_back(new TH1D(hname,
                                                hname,
                                                nChipsPerBoard*64,
                                                0,
                                                nChipsPerBoard*64
                                                ));
        hname = Form("Charge board%i",iboard);
        h_charge_boards       ->push_back(new TH1D(hname,
                                                   hname,
                                                   1024,0,1024
                                                   ));
        hname = Form("Time board%i",iboard);
        h_time_boards       ->push_back(new TH1D(hname,
                                                 hname,
                                                 256,0,256
                                                 ));
        hname = Form("NoOfClusters board%i",iboard);
        h_nclu_boards      ->push_back(new TH1D(hname,
                                                hname,
                                                4096,0,4096
                                                ));
        hname = Form("ClusterSize board%i",iboard);
        h_clSize_boards ->push_back(new TH1D(hname,
                                             hname,
                                             64,0,64
                                             ));


        for(int ichip=0;ichip<nChipsPerBoard;ichip++)
        {
            hname = Form("Beam Profile vmm%i board%i",ichip,iboard);
            h_bp_chips->at(iboard)->push_back(new TH1D(hname,
                                                       hname,
                                                       64,0,64
                                                       ));
            hname = Form("Charge vmm%i board%i",ichip,iboard);
            h_charge_chips->at(iboard)->push_back(new TH1D(hname,
                                                           hname,
                                                           1024,0,1024
                                                           ));
            hname = Form("Time vmm%i board%i",ichip,iboard);
            h_time_chips->at(iboard)->push_back(new TH1D(hname,
                                                         hname,
                                                         256,0,256
                                                         ));
            hname = Form("NoOfClusters vmm%i board%i",ichip,iboard);
            h_nclu_chips->at(iboard)->push_back(new TH1D(hname,
                                                         hname,
                                                         4096,0,4096
                                                         ));
            hname = Form("ClusterSize vmm%i board%i",ichip,iboard);
            h_clSize_chips->at(iboard)->push_back(new TH1D(hname,
                                                           hname,
                                                           64,0,64
                                                           ));
        }
    }


}
void drawSomeHistos()
{
    if(drawRawChips)
    {
        TCanvas *c_bp_chips = new TCanvas("c_bp_chips","c_bp_chips",700,500);
        TCanvas *c_charge_chips = new TCanvas("c_charge_chips","c_charge_chips",700,500);
        TCanvas *c_time_chips = new TCanvas("c_time_chips","c_time_chips",700,500);
        TCanvas *c_nclu_chips = new TCanvas("c_nclu_chips","c_nclu_chips",700,500);
        TCanvas *c_clSize_chips = new TCanvas("c_clSize_chips","c_clSize_chips",700,500);


        c_bp_chips->Divide(nBoards*nChipsPerBoard);
        c_charge_chips->Divide(nBoards*nChipsPerBoard);
        c_time_chips->Divide(nBoards*nChipsPerBoard);
        c_nclu_chips->Divide(nBoards*nChipsPerBoard);
        c_clSize_chips->Divide(nBoards*nChipsPerBoard);

        int cd_index=1;
        for(int iBoard=0;iBoard<nBoards;iBoard++)
        {
            for(int iChip=0;iChip<nChipsPerBoard;iChip++)
            {
                c_bp_chips->cd(cd_index);
                h_bp_chips->at(iBoard)->at(iChip)->Draw();
                c_charge_chips->cd(cd_index);
                h_charge_chips->at(iBoard)->at(iChip)->Draw();
                c_time_chips->cd(cd_index);
                h_time_chips->at(iBoard)->at(iChip)->Draw();
                c_nclu_chips->cd(cd_index);
                h_nclu_chips->at(iBoard)->at(iChip)->Draw();
                c_clSize_chips->cd(cd_index);
                h_clSize_chips->at(iBoard)->at(iChip)->Draw();
                cd_index++;

            }//for iChip
        }//for iBoard
    }//if(drawRawChips)

    if(drawRawBoards)
    {
        TCanvas *c_bp_boards = new TCanvas("c_bp_boards","c_bp_boards",700,500);
        TCanvas *c_charge_boards = new TCanvas("c_charge_boards","c_charge_boards",700,500);
        TCanvas *c_time_boards = new TCanvas("c_time_boards","c_time_boards",700,500);
        TCanvas *c_nclu_boards = new TCanvas("c_nclu_boards","c_nclu_boards",700,500);
        TCanvas *c_clSize_boards = new TCanvas("c_clSize_boards","c_clSize_boards",700,500);


        c_bp_boards->Divide(nBoards);
        c_charge_boards->Divide(nBoards);
        c_time_boards->Divide(nBoards);
        c_nclu_boards->Divide(nBoards);
        c_clSize_boards->Divide(nBoards);

        int cd_index=1;
        for(int iBoard=0;iBoard<nBoards;iBoard++)
        {
            c_bp_boards->cd(cd_index);
            h_bp_boards->at(iBoard)->Draw();
            c_charge_boards->cd(cd_index);
            h_charge_boards->at(iBoard)->Draw();
            c_time_boards->cd(cd_index);
            h_time_boards->at(iBoard)->Draw();
            c_nclu_boards->cd(cd_index);
            h_nclu_boards->at(iBoard)->Draw();
            c_clSize_boards->cd(cd_index);
            h_clSize_boards->at(iBoard)->Draw();
            cd_index++;

        }//for iBoard
    }//if(drawRawBoards)


}//drawSomeHistos
void adjustHistos()
{

    //    histZoomFit(&h_eventSize_chips);

    for(int i=0;i<h_clSize_chips->size();i++)
    {

        //ZOOM fit for event Size
        histZoomFit(h_clSize_boards->
                    at(i));

        for(int j=0;j<h_clSize_chips->at(i)->size();j++)
        {


            //ZOOM fit for event Size
            histZoomFit(h_clSize_chips->
                        at(i)->
                        at(j));

            //make beam profile xaxis labels to go 8 by 8
            h_bp_chips->
                    at(i)->
                    at(j)->
                    GetXaxis()->
                    SetNdivisions(808,kFALSE);
            h_bp_boards->
                    at(i)->
                    GetXaxis()->
                    SetNdivisions(808,kFALSE);
        }
    }

}

void subEventLoop()
{
    //this is called inside each Entry (after calling tree->GetEntry(i)
    //loop over the event clusters
    for(int iCluster=0;iCluster< t->cl_ChipId->size();iCluster++)
    {

        h_bp_boards->
                at( t->cl_BoardId->at(iCluster) ) ->
                Fill(t->cl_clpos_strips->at(iCluster)+stripsPerChip*stripPitch*t->cl_ChipId->at(iCluster));
        /*
        h_charge_boards->
                at( t->boardId->at(iCluster) ) ->
                Fill(t->pdo->at(iCluster).at(iHit));
        h_time_boards->
                at( t->boardId->at(iCluster) ) ->
                Fill(t->tdo->at(iCluster).at(iHit));
        h_nclu_boards->
                at( t->boardId->at(iCluster) ) ->
                Fill(t->bcid->at(iCluster).at(iHit));

*/

        h_bp_chips->
                at( t->cl_BoardId->at(iCluster) ) ->
                at(t->cl_ChipId->at(iCluster))->
                Fill(t->cl_clpos_strips->at(iCluster));
/*
        h_pdo_chips->
                at( t->boardId->at(iCluster) ) ->
                at(t->chip->at(iCluster))->
                Fill(t->pdo->at(iCluster).at(iHit));
        h_tdo_chips->
                at( t->boardId->at(iCluster) ) ->
                at(t->chip->at(iCluster))->
                Fill(t->tdo->at(iCluster).at(iHit));
        h_bcid_chips->
                at( t->boardId->at(iCluster) ) ->
                at(t->chip->at(iCluster))->
                Fill(t->bcid->at(iCluster).at(iHit));
                */
    }//cluster loop



}

void run(TString file)
{
    //load file
    load_tree_objects(file);
    //after loading, init the histo vectors
    setupHistos();

    size_t nEntries = t->fChain->GetEntries();
    cout << "Entries = " << nEntries << endl;

    for(int iEntry=0;iEntry<nEntries;iEntry++)
    {
        t->GetEntry(iEntry);
        printPercentage(iEntry,nEntries);

        //foreach entry (event) we loop over the subEvents (in this case clusters)
        //chip.size = board.size = entrySize foreach event, which means

        subEventLoop();


        //eventSize
        h_clSize_boards->
                at( t->cl_BoardId->at(iEntry) ) ->
                Fill(t->cl_clsize->at(iEntry));
        h_clSize_chips->
                at( t->cl_BoardId->at(iEntry) ) ->
                at(t->cl_ChipId->at(iEntry))->
                Fill(t->cl_clsize->at(iEntry));


    }

    //    cout << "mysize = "<<mysize<<endl;

    adjustHistos();
    drawSomeHistos();
    cout <<endl;
}

