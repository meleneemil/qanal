#include "vmm.C"
#include "vmm2.C"
#include "emiltools.cpp"
//######## SETUP ################
int nBoards=2;
int nChipsPerBoard=1;

bool drawRawChips = false;
bool drawRawBoards = false;
bool drawRawTubes = true;


//bool doClusterization=true;

//TAC
float tacCounts = 255.0;//cnt
float tacDuration = 650.0;//ns
float tdoMultiplier = tacDuration/tacCounts;
//so when we have tdo counts, we make it time by:
// time = 650-tdoMultiplier*tdoCnts


//########################
vector<vector<TH1D*>*> *h_bp_chips;
vector<vector<TH1D*>*> *h_pdo_chips;
vector<vector<TH1D*>*> *h_tdo_chips;
vector<vector<TH1D*>*> *h_bcid_chips;
vector<vector<TH1D*>*> *h_eventSize_chips;

vector<TH1D*> *h_bp_boards;
vector<TH1D*> *h_pdo_boards;
vector<TH1D*> *h_tdo_boards;
vector<TH1D*> *h_bcid_boards;
vector<TH1D*> *h_eventSize_boards;

vector<vector<TH1D*>*> *h_pdo_tubes;
vector<vector<TH1D*>*> *h_time_tubes;

TH1D *h_time_allTubes0;
TH1D *h_time_allTubes1;
TH1D *h_time_allTubesBoth;

TCanvas *c_pdo_tubes0,*c_time_tubes0,*c_pdo_tubes1,*c_time_tubes1;
TCanvas *c_time_tubesCombined;


TH1D* test1 = new TH1D("test1","test1",65,0,65);
TH1D* test2 = new TH1D("test2","test2",20,-10,100);

int mysize=0;
TString hname;
TString sRunNumber;

vmm *t;//this is the object from which we get the tree

void load_tree_objects(TString file)
{
    sRunNumber = get_run_number_from_file(file);
    //Get file
    TFile *root_file = new TFile(file.Data(),"READ");
    //Load all cluster trees from file
    t  = new vmm((TTree*)root_file->Get("vmm"));
}


void setupHistos()
{

    gStyle->SetTitleFontSize(0.08);
    gStyle->SetStatH(0.5);
    gStyle->SetStatW(0.25);
    gStyle->SetStatX(0.88);
    gStyle->SetStatY(0.88);


    h_bp_chips         = new vector<vector<TH1D*>*>;
    h_pdo_chips        = new vector<vector<TH1D*>*>;
    h_tdo_chips        = new vector<vector<TH1D*>*>;
    h_bcid_chips       = new vector<vector<TH1D*>*>;
    h_eventSize_chips  = new vector<vector<TH1D*>*>;

    h_bp_boards        = new vector<TH1D*>;
    h_pdo_boards       = new vector<TH1D*>;
    h_tdo_boards       = new vector<TH1D*>;
    h_bcid_boards      = new vector<TH1D*>;
    h_eventSize_boards = new vector<TH1D*>;

    h_pdo_tubes        = new vector<vector<TH1D*>*>;
    h_time_tubes       = new vector<vector<TH1D*>*>;

    hname = Form("run"+sRunNumber+"Time all tubes board0");
    h_time_allTubes0   = new TH1D(hname,hname,650,0,650);
    hname = Form("run"+sRunNumber+"Time all tubes board1");
    h_time_allTubes1   = new TH1D(hname,hname,650,0,650);
    hname = Form("run"+sRunNumber+"Time all tubes 2boards");
    h_time_allTubesBoth   = new TH1D(hname,hname,650,0,650);

    //loop the boards
    for(int iboard=0;iboard<nBoards;iboard++)
    {

        gStyle->SetLabelSize(0.04,"X");
        gStyle->SetLabelSize(0.03,"Y");

        //board histos
        hname = Form("Beam Profile board%i",iboard);
        h_bp_boards        ->push_back(new TH1D(hname,
                                                hname,
                                                nChipsPerBoard*64,
                                                0,
                                                nChipsPerBoard*64
                                                ));
        hname = Form("PDO board%i",iboard);
        h_pdo_boards       ->push_back(new TH1D(hname,
                                                hname,
                                                1024,0,1024
                                                ));
        hname = Form("TDO board%i",iboard);
        h_tdo_boards       ->push_back(new TH1D(hname,
                                                hname,
                                                256,0,256
                                                ));
        hname = Form("BCID board%i",iboard);
        h_bcid_boards      ->push_back(new TH1D(hname,
                                                hname,
                                                4096,0,4096
                                                ));
        hname = Form("EventSize board%i",iboard);
        h_eventSize_boards ->push_back(new TH1D(hname,
                                                hname,
                                                64,0,64
                                                ));

        //chip histos

        h_bp_chips        ->push_back(new vector<TH1D*>());
        h_pdo_chips       ->push_back(new vector<TH1D*>());
        h_tdo_chips       ->push_back(new vector<TH1D*>());
        h_bcid_chips      ->push_back(new vector<TH1D*>());
        h_eventSize_chips ->push_back(new vector<TH1D*>());

        for(int ichip=0;ichip<nChipsPerBoard;ichip++)
        {
            hname = Form("Beam Profile vmm%i board%i",ichip,iboard);
            h_bp_chips->at(iboard)->push_back(new TH1D(hname,
                                                       hname,
                                                       64,0,64
                                                       ));
            hname = Form("PDO vmm%i board%i",ichip,iboard);
            h_pdo_chips->at(iboard)->push_back(new TH1D(hname,
                                                        hname,
                                                        1024,0,1024
                                                        ));
            hname = Form("TDO vmm%i board%i",ichip,iboard);
            h_tdo_chips->at(iboard)->push_back(new TH1D(hname,
                                                        hname,
                                                        256,0,256
                                                        ));
            hname = Form("BCID vmm%i board%i",ichip,iboard);
            h_bcid_chips->at(iboard)->push_back(new TH1D(hname,
                                                         hname,
                                                         4096,0,4096
                                                         ));
            hname = Form("EventSize vmm%i board%i",ichip,iboard);
            h_eventSize_chips->at(iboard)->push_back(new TH1D(hname,
                                                              hname,
                                                              64,0,64
                                                              ));
        }

        gStyle->SetLabelSize(0.09,"XY");

        //tube histos
        h_pdo_tubes       ->push_back(new vector<TH1D*>());
        h_time_tubes      ->push_back(new vector<TH1D*>());

        for(int iTube=0;iTube<25;iTube++)
        {
            hname = Form("PDO board%i tube%i",iboard,iTube);
            h_pdo_tubes->at(iboard)->push_back(new TH1D(hname,
                                                        hname,
                                                        1024,
                                                        0,
                                                        1024
                                                        ));



            hname = Form("Time board%i tube%i",iboard,iTube);
            h_time_tubes->at(iboard)->push_back(new TH1D(hname,
                                                         hname,
                                                         650,
                                                         0,
                                                         650
                                                         ));

        }

    }


}
void drawSomeHistos()
{


    if(drawRawChips)
    {
        TCanvas *c_bp_chips = new TCanvas("c_bp_chips","c_bp_chips",700,500);
        TCanvas *c_pdo_chips = new TCanvas("c_pdo_chips","c_pdo_chips",700,500);
        TCanvas *c_tdo_chips = new TCanvas("c_tdo_chips","c_tdo_chips",700,500);
        TCanvas *c_bcid_chips = new TCanvas("c_bcid_chips","c_bcid_chips",700,500);
        TCanvas *c_eventSize_chips = new TCanvas("c_eventSize_chips","c_eventSize_chips",700,500);


        c_bp_chips->Divide(nBoards*nChipsPerBoard);
        c_pdo_chips->Divide(nBoards*nChipsPerBoard);
        c_tdo_chips->Divide(nBoards*nChipsPerBoard);
        c_bcid_chips->Divide(nBoards*nChipsPerBoard);
        c_eventSize_chips->Divide(nBoards*nChipsPerBoard);

        int cd_index=1;
        for(int iBoard=0;iBoard<nBoards;iBoard++)
        {
            for(int iChip=0;iChip<nChipsPerBoard;iChip++)
            {
                c_bp_chips->cd(cd_index);
                h_bp_chips->at(iBoard)->at(iChip)->Draw();
                c_pdo_chips->cd(cd_index);
                h_pdo_chips->at(iBoard)->at(iChip)->Draw();
                c_tdo_chips->cd(cd_index);
                h_tdo_chips->at(iBoard)->at(iChip)->Draw();
                c_bcid_chips->cd(cd_index);
                h_bcid_chips->at(iBoard)->at(iChip)->Draw();
                c_eventSize_chips->cd(cd_index);
                h_eventSize_chips->at(iBoard)->at(iChip)->Draw();
                cd_index++;

            }//for iChip
        }//for iBoard
    }//if(drawRawChips)

    if(drawRawBoards)
    {
        TCanvas *c_bp_boards = new TCanvas("c_bp_boards","c_bp_boards",700,500);
        TCanvas *c_pdo_boards = new TCanvas("c_pdo_boards","c_pdo_boards",700,500);
        TCanvas *c_tdo_boards = new TCanvas("c_tdo_boards","c_tdo_boards",700,500);
        TCanvas *c_bcid_boards = new TCanvas("c_bcid_boards","c_bcid_boards",700,500);
        TCanvas *c_eventSize_boards = new TCanvas("c_eventSize_boards","c_eventSize_boards",700,500);


        c_bp_boards->Divide(nBoards);
        c_pdo_boards->Divide(nBoards);
        c_tdo_boards->Divide(nBoards);
        c_bcid_boards->Divide(nBoards);
        c_eventSize_boards->Divide(nBoards);

        int cd_index=1;
        for(int iBoard=0;iBoard<nBoards;iBoard++)
        {
            c_bp_boards->cd(cd_index);
            h_bp_boards->at(iBoard)->Draw();
            c_pdo_boards->cd(cd_index);
            h_pdo_boards->at(iBoard)->Draw();
            c_tdo_boards->cd(cd_index);
            h_tdo_boards->at(iBoard)->Draw();
            c_bcid_boards->cd(cd_index);
            h_bcid_boards->at(iBoard)->Draw();
            c_eventSize_boards->cd(cd_index);
            h_eventSize_boards->at(iBoard)->Draw();
            cd_index++;

        }//for iBoard
    }//if(drawRawBoards)

    if(drawRawTubes)
    {
        hname = Form("c_pdo_tubes board0 run"+sRunNumber);
        c_pdo_tubes0  = new TCanvas(hname,hname,700,500);
        hname = Form("c_time_tubes board0 run"+sRunNumber);
        c_time_tubes0 = new TCanvas(hname,hname,700,500);
        hname = Form("c_pdo_tubes board1 run"+sRunNumber);
        c_pdo_tubes1  = new TCanvas(hname,hname,700,500);
        hname = Form("c_time_tubes board1 run"+sRunNumber);
        c_time_tubes1 = new TCanvas(hname,hname,700,500);


        c_pdo_tubes0->Divide(4,6);
        c_time_tubes0->Divide(4,6);
        c_pdo_tubes1->Divide(4,6);
        c_time_tubes1->Divide(4,6);

        int cd_index=1;
        for(int iTube=0;iTube<24;iTube++)
        {
            h_pdo_tubes->at(0)->at(iTube)->Rebin(16);
            h_pdo_tubes->at(1)->at(iTube)->Rebin(16);
            h_time_tubes->at(0)->at(iTube)->Rebin(16);
            h_time_tubes->at(1)->at(iTube)->Rebin(16);

            c_pdo_tubes0->cd(iTube+1);
            h_pdo_tubes->at(0)->at(iTube)->Draw();
            c_time_tubes0->cd(iTube+1);
            h_time_tubes->at(0)->at(iTube)->Draw();

            c_pdo_tubes1->cd(iTube+1);
            h_pdo_tubes->at(1)->at(iTube)->Draw();
            c_time_tubes1->cd(iTube+1);
            h_time_tubes->at(1)->at(iTube)->Draw();

        }//for iBoard

        //Combine Tubes


        h_time_allTubes0->Rebin(16);
        h_time_allTubes1->Rebin(16);
        h_time_allTubesBoth->Rebin(16);
        for(int iTube=0;iTube<24;iTube++)
        {
            h_time_allTubes0->Add(h_time_tubes->at(0)->at(iTube));
            h_time_allTubes1->Add(h_time_tubes->at(1)->at(iTube));
            h_time_allTubesBoth->Add(h_time_tubes->at(0)->at(iTube));
            h_time_allTubesBoth->Add(h_time_tubes->at(1)->at(iTube));
        }

        //        h_time_allTubes0->GetXaxis()->SetTitle("[ns]");
        //        h_time_allTubes0->GetYaxis()->SetTitle("[ev/25ns]");


        hname = Form("run"+sRunNumber+" c_time_tubesCombined");
        c_time_tubesCombined  = new TCanvas(hname,hname,1000,700);
        c_time_tubesCombined->Divide(1,3);
        c_time_tubesCombined->cd(1);
        h_time_allTubes0->Draw();
        c_time_tubesCombined->cd(2);
        h_time_allTubes1->Draw();
        c_time_tubesCombined->cd(3);
        h_time_allTubesBoth->Draw();



        //        transpad(h_time_allTubes0,h_time_allTubes1,"Time AllTubes");

    }//if(drawRawTubes)

}
void adjustHistos()
{

    //    histZoomFit(&h_eventSize_chips);

    //board loop
    for(int i=0;i<h_eventSize_chips->size();i++)
    {

        //ZOOM fit for event Size
        if(drawRawBoards)
        {
            histZoomFit(h_eventSize_boards->
                        at(i));
        }
        //chip loop
        if(drawRawChips)
            for(int j=0;j<h_eventSize_chips->at(i)->size();j++)
            {


                //ZOOM fit for event Size
                histZoomFit(h_eventSize_chips->
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

        //tube loop
        if(drawRawTubes)
            for(int iTube=0;iTube<25;iTube++)
            {
                histZoomFit(
                            h_pdo_tubes->
                            at(i)->
                            at(iTube)
                            );
                if(false)
                    histZoomFit(
                                h_time_tubes->
                                at(i)->
                                at(iTube)
                                );

            }
    }

}
void saveSomeHistosAsPng()
{
    saveCanvasAsAll(c_pdo_tubes0);
    saveCanvasAsAll(c_time_tubes0);
    saveCanvasAsAll(c_pdo_tubes1);
    saveCanvasAsAll(c_time_tubes1);
    saveCanvasAsAll(c_time_tubesCombined);
}

void subEventLoop()
{
    //this is called inside each Entry (after calling tree->GetEntry(i)

    for(int iSubEvent=0;iSubEvent < t->chip->size();iSubEvent++)
    {

        //check t->boardId
        //loop over the chip's hits of this event
        for(int iHit=0;iHit<t->channel->at(iSubEvent).size();iHit++)
        {
            //                cout << t->channel->at(iChip).at(iHit) << " ";

            h_bp_boards->
                    at( t->boardId->at(iSubEvent) ) ->
                    Fill(t->channel->at(iSubEvent).at(iHit)+64*t->chip->at(iSubEvent));
            h_pdo_boards->
                    at( t->boardId->at(iSubEvent) ) ->
                    Fill(t->pdo->at(iSubEvent).at(iHit));
            h_tdo_boards->
                    at( t->boardId->at(iSubEvent) ) ->
                    Fill(t->tdo->at(iSubEvent).at(iHit));
            h_bcid_boards->
                    at( t->boardId->at(iSubEvent) ) ->
                    Fill(t->bcid->at(iSubEvent).at(iHit));



            h_bp_chips->
                    at( t->boardId->at(iSubEvent) ) ->
                    at(t->chip->at(iSubEvent))->
                    Fill(t->channel->at(iSubEvent).at(iHit));

            h_pdo_chips->
                    at( t->boardId->at(iSubEvent) ) ->
                    at(t->chip->at(iSubEvent))->
                    Fill(t->pdo->at(iSubEvent).at(iHit));
            h_tdo_chips->
                    at( t->boardId->at(iSubEvent) ) ->
                    at(t->chip->at(iSubEvent))->
                    Fill(t->tdo->at(iSubEvent).at(iHit));
            h_bcid_chips->
                    at( t->boardId->at(iSubEvent) ) ->
                    at(t->chip->at(iSubEvent))->
                    Fill(t->bcid->at(iSubEvent).at(iHit));


            h_pdo_tubes->
                    at( t->boardId->at(iSubEvent) ) ->
                    at(t->channel->at(iSubEvent).at(iHit))->
                    Fill(t->pdo->at(iSubEvent).at(iHit));

            h_time_tubes->
                    at( t->boardId->at(iSubEvent) ) ->
                    at(t->channel->at(iSubEvent).at(iHit))->
                    Fill(650-tdoMultiplier*t->tdo->at(iSubEvent).at(iHit));


        }

        //eventSize
        h_eventSize_boards->
                at( t->boardId->at(iSubEvent) ) ->
                Fill(t->channel->at(iSubEvent).size());
        h_eventSize_chips->
                at( t->boardId->at(iSubEvent) ) ->
                at(t->chip->at(iSubEvent))->
                Fill(t->channel->at(iSubEvent).size());

    }
}
void run(TString file)
{
    //load file
    load_tree_objects(file);
    //after loading, init the histo vectors
    setupHistos();

    size_t nEntries = t->fChain->GetEntries();
    cout << "run#"<<sRunNumber<<" | Entries = " << nEntries << endl;

    for(int iEntry=0;iEntry<nEntries;iEntry++)
    {
        t->GetEntry(iEntry);
        printPercentage(iEntry,nEntries);

        //foreach entry (event) we loop over the subEvents
        //chip.size = board.size = entrySize foreach event, which means
        subEventLoop();

    }

    //    cout << "mysize = "<<mysize<<endl;

    adjustHistos();
    drawSomeHistos();
    saveSomeHistosAsPng();


    cout <<endl;
}













