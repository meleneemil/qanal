#include "vmm2.C"
#include "emiltools.cpp"
//######## SETUP #########
bool drawRawChips = false;

//bool doClusterization=true;

//########################

//cluster vectors

vector<int> cl_width;
vector<int> cl_charge;
vector<int> cl_tdo;
vector<int> cl_posCentroid;
vector<vector<int>> cl_strips;

//GLOBAL cuts config
int minStripsPerChipEvent=1;
//int minStripsPerCluster=1;
int minPdo_strip=30;
int minTdo_strip=50;
bool tdoTest=true;

//########################

TGraph *tg_pdoStrip;
TGraph *tg_tdoStrip;
TGraph *tg_bcidStrip;

vector<TH1D*> *h_bp_chips;
vector<TH1D*> *h_pdo_chips;
vector<TH1D*> *h_tdo_chips;
vector<TH1D*> *h_bcid_chips;
vector<TH1D*> *h_eventSize_chips;

vector<TH2D*> *h_tdoVSpdo_chips;
TH2D *h_tdoVSpdo_chips_all;
TH2D *h_tdoVSpdo_chips_all_afterCut;

TH1D *h_tdo_min_max_diffs;

TCanvas *c_bp_chips,*c_pdo_chips,*c_tdo_chips,*c_bcid_chips,*c_evSize_chips,*c_evPdo,*c_evTdo,*c_evBcid,*c_tdoVSpdo_chips,*c_tdoVSpdo_all;

TCanvas *c_tdo_min_max_diffs;

int evSize_afterHitCuts=0;

string cmdInput="";
TString hname;
TString sRunNumber;

vmm2 *t;//this is the object from which we get the tree

void load_tree_objects(TString file)
{
    sRunNumber = get_run_number_from_file(file);
    //Get file
    TFile *root_file = new TFile(file.Data(),"READ");
    //Load all cluster trees from file
    t  = new vmm2((TTree*)root_file->Get("vmm2"));
}
void setupHistos()
{

    c_evPdo = new TCanvas("event PDO vs strip","event PDO vs strip",1000,10,900,900);
    c_evPdo->Divide(1,3);

    h_bp_chips = new vector<TH1D*>;
    h_pdo_chips = new vector<TH1D*>;
    h_tdo_chips = new vector<TH1D*>;
    h_bcid_chips = new vector<TH1D*>;
    h_eventSize_chips = new vector<TH1D*>;
    h_tdoVSpdo_chips = new vector<TH2D*>;
    hname = Form("h_tdoVSpdo_chips_all");
    h_tdoVSpdo_chips_all = new TH2D(hname,hname,1024-minPdo_strip,0+minPdo_strip,1024,256-minTdo_strip,0+minTdo_strip,256);
    hname = Form("h_tdoVSpdo_chips_all after cut");
    h_tdoVSpdo_chips_all_afterCut = new TH2D(hname,hname,1024-minPdo_strip,0+minPdo_strip,1024,256-minTdo_strip,0+minTdo_strip,256);


    hname = Form("h_tdo_min_max_diffs");
    h_tdo_min_max_diffs=new TH1D(hname,hname,185,0,185);

    for(int iChip=0;iChip<8;iChip++)
    {
        hname = Form("Beam Profile chip%i",iChip);
        h_bp_chips->push_back(new TH1D(hname,hname,64,0,64));
        hname = Form("PDO chip%i",iChip);
        h_pdo_chips->push_back(new TH1D(hname,hname,1024,0,1024));
        hname = Form("TDO chip%i",iChip);
        h_tdo_chips->push_back(new TH1D(hname,hname,256,0,256));
        hname = Form("BCID chip%i",iChip);
        h_bcid_chips->push_back(new TH1D(hname,hname,4096,0,4096));
        hname = Form("EvSize chip%i",iChip);
        h_eventSize_chips->push_back(new TH1D(hname,hname,64,0,64));

        hname = Form("pdo VS tdo chip%i",iChip);
        h_tdoVSpdo_chips->push_back( new TH2D(hname,hname,1024,0,1024,256,0,256) );
    }
}
void drawHistos()
{
    //beam profiles
    {
        hname = Form("run"+sRunNumber+" Beam Profile chips");
        c_bp_chips=new TCanvas(hname,hname,200,200,700,600);

        c_bp_chips->Divide(4,2);
        for(int iChip=0;iChip<8;iChip++)
        {
            c_bp_chips->cd(iChip+1);
            h_bp_chips->at(iChip)->Draw();

        }
    }
    //pdo
    {
        hname = Form("run"+sRunNumber+" PDO chips");
        c_pdo_chips=new TCanvas(hname,hname,200,200,700,600);

        c_pdo_chips->Divide(4,2);
        for(int iChip=0;iChip<8;iChip++)
        {
            c_pdo_chips->cd(iChip+1);
            h_pdo_chips->at(iChip)->Draw();

        }
    }
    //tdo
    {
        hname = Form("run"+sRunNumber+" TDO chips");
        c_tdo_chips=new TCanvas(hname,hname,200,200,700,600);

        c_tdo_chips->Divide(4,2);
        for(int iChip=0;iChip<8;iChip++)
        {
            c_tdo_chips->cd(iChip+1);
            h_tdo_chips->at(iChip)->Draw();

        }
    }
    //bcid
    {
        hname = Form("run"+sRunNumber+" bcid chips");
        c_bcid_chips=new TCanvas(hname,hname,200,200,700,600);

        c_bcid_chips->Divide(4,2);
        for(int iChip=0;iChip<8;iChip++)
        {
            c_bcid_chips->cd(iChip+1);
            h_bcid_chips->at(iChip)->Draw();

        }
    }
    //evSize
    {
        hname = Form("run"+sRunNumber+" evSize chips");
        c_evSize_chips=new TCanvas(hname,hname,200,200,700,600);

        c_evSize_chips->Divide(4,2);
        for(int iChip=0;iChip<8;iChip++)
        {
            c_evSize_chips->cd(iChip+1);
            h_eventSize_chips->at(iChip)->Draw();

            histZoomFit(h_eventSize_chips->at(iChip));

        }
    }
    //pdo vs tdo
    {
        hname = Form("run"+sRunNumber+" tdo VS pdo chips");
        c_tdoVSpdo_chips=new TCanvas(hname,hname,200,200,700,600);

        c_tdoVSpdo_chips->Divide(4,2);
        for(int iChip=0;iChip<8;iChip++)
        {
            c_tdoVSpdo_chips->cd(iChip+1);
            h_tdoVSpdo_chips->at(iChip)->Draw("colz");

            //histZoomFit(h_tdoVSpdo_chips->at(iChip));
        }

        hname = Form("run"+sRunNumber+" tdo VS pdo all chips");
        c_tdoVSpdo_all=new TCanvas(hname,hname,200,200,700,600);
        c_tdoVSpdo_all->Divide(2,1);
        c_tdoVSpdo_all->cd(1);
        h_tdoVSpdo_chips_all->Draw("colz");
        c_tdoVSpdo_all->cd(2);
        h_tdoVSpdo_chips_all_afterCut->Draw("colz");
    }

    //tdo diffs
    hname = Form("run"+sRunNumber+" tdo max-min all chips");
    c_tdo_min_max_diffs = new TCanvas(hname,hname,300,200,700,600);
    h_tdo_min_max_diffs->Draw();

}

void clusterize()
{
    //pdo>50

}

void subEventLoop(bool eventByEvent)
{
    //loop over the chips of the event
    for(int iSubEvent=0;iSubEvent < t->chip->size();iSubEvent++)
    {

        int chipEventSize= t->channel->at(iSubEvent).size();

        // apply cuts
        if(chipEventSize<minStripsPerChipEvent)
            continue;

        //initialize vars
        evSize_afterHitCuts=0;
        int pdos[chipEventSize]   = {0};
        int tdos[chipEventSize]   = {0};
        int bcids[chipEventSize]  = {0};
        int strips[chipEventSize] = {0};

        //        cout << "Filling chip" << t->chip->at(iSubEvent) <<endl;

        //print event by event info
        if(eventByEvent)
        {
            cout << "ev="<<t->eventFAFA<<"\ts";
            cout << "subEv="<<t->chip->at(iSubEvent) << " - \t";
        }

        for(int iHit=0;iHit<chipEventSize;iHit++)
        {

            //Apply GLOBAL strip cuts
            if(
                    t->pdo->at(iSubEvent).at(iHit) < minPdo_strip ||
                    t->tdo->at(iSubEvent).at(iHit) < minTdo_strip
                    //t->pdo->at(iSubEvent).at(iHit) % 16 == 0 ||
                    //t->tdo->at(iSubEvent).at(iHit) % 8 == 0
                    )
            {
                continue;
            }
            else
            {
                //continue;
            }

            //fill final histos
            h_bp_chips->
                    at(t->chip->
                       at(iSubEvent))->
                    Fill(t->channel->at(iSubEvent).at(iHit));

            //cout << t->channel->at(iSubEvent).at(iHit) << " ";

            h_pdo_chips->
                    at(t->chip->
                       at(iSubEvent))->
                    Fill(t->pdo->at(iSubEvent).at(iHit));

            h_tdo_chips->
                    at(t->chip->
                       at(iSubEvent))->
                    Fill(t->tdo->at(iSubEvent).at(iHit));

            h_bcid_chips->
                    at(t->chip->
                       at(iSubEvent))->
                    Fill(t->bcid->at(iSubEvent).at(iHit));

            h_tdoVSpdo_chips->
                    at(t->chip->
                       at(iSubEvent))->
                    Fill(
                        t->pdo->at(iSubEvent).at(iHit),
                        t->tdo->at(iSubEvent).at(iHit)
                        );
            h_tdoVSpdo_chips_all->
                    Fill(
                        t->pdo->at(iSubEvent).at(iHit),
                        t->tdo->at(iSubEvent).at(iHit)
                        );

            //event by event: fill TGraph arrays
            if(eventByEvent || tdoTest)
            {
                pdos[iHit] = t->pdo->at(iSubEvent).at(iHit);
                tdos[iHit] = t->tdo->at(iSubEvent).at(iHit);
                bcids[iHit] = t->bcid->at(iSubEvent).at(iHit);
                strips[iHit] = t->channel->at(iSubEvent).at(iHit);
            }


            //Apply strip cuts for other histos
            if(
                    //t->pdo->at(iSubEvent).at(iHit) < minPdo_strip ||
                    //t->tdo->at(iSubEvent).at(iHit) < minTdo_strip
                    t->pdo->at(iSubEvent).at(iHit) % 16 == 0 ||
                    t->tdo->at(iSubEvent).at(iHit) % 8 == 0
                    )
            {
                continue;
            }
            else
            {
                //                continue;
            }

            //FILL OTHER HISTOS (for after cuts comparison)
            h_tdoVSpdo_chips_all_afterCut->
                    Fill(
                        t->pdo->at(iSubEvent).at(iHit),
                        t->tdo->at(iSubEvent).at(iHit)
                        );




            //increment counter for "effective" strips in event
            evSize_afterHitCuts++;

        }//hitLoop

        //skip event if not effective strips were found
        if(evSize_afterHitCuts==0)
            continue;


        h_eventSize_chips->at(t->chip->at(iSubEvent))->Fill(evSize_afterHitCuts);
        h_tdo_min_max_diffs->Fill(
                    tdos[findMax(tdos,chipEventSize)]-
                tdos[findMin(tdos,chipEventSize)]
                );

        //event by event plots!
        if(eventByEvent)
        {
            tg_pdoStrip = new TGraph(chipEventSize,strips,pdos);
            tg_pdoStrip->SetTitle(Form("PDO Chip%i",t->chip->at(iSubEvent)));
            tg_pdoStrip->SetMarkerSize(2);
            tg_pdoStrip->SetMarkerStyle(22);
            tg_tdoStrip = new TGraph(chipEventSize,strips,tdos);
            tg_tdoStrip->SetTitle(Form("TDO Chip%i",t->chip->at(iSubEvent)));
            tg_tdoStrip->SetMarkerColor(kRed);
            tg_tdoStrip->SetMarkerSize(2);
            tg_bcidStrip = new TGraph(chipEventSize,strips,bcids);
            tg_bcidStrip->SetTitle(Form("BCID Chip%i",t->chip->at(iSubEvent)));
            tg_bcidStrip->SetMarkerColor(kRed);
            tg_bcidStrip->SetMarkerSize(2);

            cout << "evSize="<<chipEventSize;//<<endl;

            gStyle->SetLabelSize(.1,"XY");
            //            tg_pdoStrip->SetMinimum(0);
            //            tg_tdoStrip->SetMinimum(0);
            //            tg_bcidStrip->SetMinimum(0);

            c_evPdo->cd(1);
            gPad->Modified();
            gPad->Update();
            tg_pdoStrip->Draw("A*");
            c_evPdo->cd(2);
            gPad->Modified();
            gPad->Update();
            tg_tdoStrip->Draw("A*");
            c_evPdo->cd(3);
            gPad->Modified();
            gPad->Update();
            tg_bcidStrip->Draw("A*");


            c_evPdo->Resize();
            c_evPdo->Modified();
            c_evPdo->Update();

            getline(cin, cmdInput);
        }//ev by ev plots

        //if(eventNumberForDisplay[0] == '\0') continue

    }//subEvent loop
}

void run(TString file, bool eventByEvent=false)
{
    //load file
    load_tree_objects(file);

    setupHistos();

    size_t nEntries = t->fChain->GetEntries();
    cout << "run#"<<sRunNumber<<" | Entries = " << nEntries << endl;

    for(int iEntry=0;iEntry<nEntries;iEntry++)
    {
        t->GetEntry(iEntry);
        printPercentage(iEntry,nEntries);

        //foreach entry (event) we loop over the subEvents
        subEventLoop(eventByEvent);

        //        cout <<endl;

    }//entry loop
    cout <<endl;

    drawHistos();

}

