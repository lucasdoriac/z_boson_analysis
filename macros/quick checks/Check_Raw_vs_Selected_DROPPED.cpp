/*
Plots raw vs selected variables in the same canvas.

Raw: from the processed TTree, before any selection is applied.
Selected: from 'mySelectedData.root', after applying the selection criteria.

*/


//---Enumerates
enum class SampleType {
    Data,
    MC
};

enum class CollisionSystem {
    ppRef,
    PbPb
};

//---Structs
struct Dataset {
    std::string name;
    SampleType type;
    CollisionSystem system;
    int year;
    std::string treeName;
    std::string filePattern;
    std::string basePath;

    bool hasCentrality;//Or maybe is AA
    bool applyTrigger;
    ULong64_t triggerBit;
};

Dataset datasets[] = {
    
    {
        "ppRef2024_Data",
        SampleType::Data,
        CollisionSystem::ppRef,
        2024,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_ppRef2024.root",
        BasePath + "Data/ppRef2024/",
        false,
        false,
        0ULL
    },
    
    {
        "PbPb2023_Data",
        SampleType::Data,
        CollisionSystem::PbPb,
        2023,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_PbPb2023.root",
        BasePath + "Data/PbPb2023/",
        true,
        true,
        1ULL << 6 //'HLT_HIL2SingleMu7_v'
    },

    {
        "PbPb2024_Data",
        SampleType::Data,
        CollisionSystem::PbPb,
        2024,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_PbPb2024Data.root",
        BasePath + "Data/PbPb2024/",
        true,
        true,
        1ULL << 7 //'HLT_HIL2SingleMu12_v'
    },

    {
        "PbPb2025_Data",
        SampleType::Data,
        CollisionSystem::PbPb,
        2025,
        "hionia/DimuonTree",
        "HighPtMuon_PbPb2025Data.root",
        BasePath + "Data/PbPb2025/",
        true,
        true,
        1ULL << 7 //'HLT_HIL2SingleMu12_v'
    },

    {
        "PbPb2026_Data",
        SampleType::Data,
        CollisionSystem::PbPb,
        2026,
        "hionia/DimuonTree",
        "HighPtMuon_PbPb2026Data.root",
        BasePath + "Data/PbPb2026/",
        true,
        true,
        1ULL << 7 //'HLT_HIL2SingleMu12_v'
    },

    {
        "PbPb2024_MC",
        SampleType::MC,
        CollisionSystem::PbPb,
        2024,
        "hionia/myTree",
        "Oniatree_PowhegZtoMuMu_PbPb2024_*.root",
        BasePath + "MC/PbPb2024/DYto2Mu_MLL-50_TuneCP5_5p36TeV_powheg-pythia8/PowhegEmbedded_March9/260309_143939/0000/",
        false,
        false,
        0ULL
    },

    {
        "ppRef2024_MC",
        SampleType::MC,
        CollisionSystem::ppRef,
        2024,
        "hionia/myTree",
        "Oniatree_PowhegZtoMuMu_ppRef2024_*.root",
        BasePath + "MC/ppRef2024/DYToMuMu_M-50_TuneCP5_5p36TeV_powheg-pythia8/Powheg_ppRefPileup_March20/260320_125046/0000/",
        false,
        false,
        0ULL
    }
};


//List of histograms to plot raw vs selected.
std::vector<std::string> histNames = {

    //"h3D_PtMuPl_PtMuMi_Cent",
    "h1D_centrality",
    "h1D_hiHF",
    "h1D_invMass",
    "h1D_zPt",
    "h1D_zRapidity",

    "h1D_ptMuPlus",
    "h1D_ptMuMinus",

    "h1D_etaMuPlus",
    "h1D_etaMuMinus",

    "h1D_phiMuPlus",
    "h1D_phiMuMinus",

    "h1D_zVtx",

    "h1D_deltaPhiMuMu",
    "h1D_deltaRMuMu",
    "h1D_acoplanarity",

    "h1D_muonPtRelDiff",
    "h1D_muonPtDiff"

    //"h2D_muonPtRelDiff_Cent",
    //"h2D_zPt_Cent",

//    "h2D_muonPtRelDiff_zPt",
//    "h2D_muonPtRelDiff_zRapidity",
//    "h2D_zPt_zRapidity",
//    "h2D_PtMuPl_PtMuMi"
};

struct HistStruct{
    std::string name;
    TH1* rawHist;
    TH1* selectedHist;
};

std::vector<HistStruct> histStructs;

void myFunction(TFile* inputFile, Dataset& dataset);
void MakePlot(HistStruct& histStruct, Dataset& dataset);

void Check_Raw_vs_Selected(){

    gROOT->SetBatch(kTRUE);
    TFile* inputFile = new TFile("mySelectedData.root", "READ");

    myFunction(inputFile, datasets[0]);
    myFunction(inputFile, datasets[1]);
    myFunction(inputFile, datasets[2]);
    myFunction(inputFile, datasets[3]);
    myFunction(inputFile, datasets[4]);
    myFunction(inputFile, );

    inputFile->Close();
}


void myFunction(TFile* inputFile, Dataset& dataset){

    //Get raw from TTree.
    //This part is exactly the same as in ApplyGoodSelection.cpp, just without the selections.
    std::string fullPath = dataset.basePath + dataset.filePattern;

    TChain *chain = new TChain(dataset.treeName.c_str());
    chain->Add(fullPath.c_str());

    std::cout << "\n> Number of files added to TChain = " << chain->GetListOfFiles()->GetEntries() << "\n" << std::endl;

    std::cout << "> Opening files " << fullPath << "\n" << std::endl;

    std::cout << "> Running function " << __func__ << " on " << dataset.name << "\n" << std::endl;
    
    //Total number of events on Tree.
    Long64_t nEvents = chain->GetEntries();

    std::cout << "> Total number of events on tree = " << nEvents << "\n" << std::endl;

    const int MAX_DIMUON = 1000;
    const int MAX_MUON   = 1000;

    //For now, writing only a fraction of the branches i.e. writing only what i'll use.

    //Event-level variables
    Int_t Centrality;
    Float_t  zVtx;
    Float_t SumET_HF;

    chain->SetBranchAddress("zVtx", &zVtx);
    
    if(dataset.hasCentrality) {
        chain->SetBranchAddress("Centrality", &Centrality);
        chain->SetBranchAddress("SumET_HF", &SumET_HF);
    }

    //Dimuon-level variables
    Short_t Reco_Dimuon_size;

    Short_t Reco_Dimuon_sign[MAX_DIMUON];
    Short_t Reco_Dimuon_muonPlusIndex[MAX_DIMUON];
    Short_t Reco_Dimuon_muonMinusIndex[MAX_DIMUON];

    ULong64_t Reco_Dimuon_trig[MAX_DIMUON];
    Float_t Reco_Dimuon_vtxProb[MAX_DIMUON];

    std::vector<float>* Reco_Dimuon_pt = nullptr;
    std::vector<float>* Reco_Dimuon_eta = nullptr;
    std::vector<float>* Reco_Dimuon_rapidity = nullptr;
    std::vector<float>* Reco_Dimuon_phi = nullptr;
    std::vector<float>* Reco_Dimuon_invMass = nullptr;

    std::vector<float>* Reco_Dimuon_muonPtDiff = nullptr;
    std::vector<float>* Reco_Dimuon_muonPtRelDiff = nullptr;

    chain->SetBranchAddress("Reco_Dimuon_size", &Reco_Dimuon_size);

    chain->SetBranchAddress("Reco_Dimuon_sign", Reco_Dimuon_sign);
    chain->SetBranchAddress("Reco_Dimuon_muonPlusIndex", Reco_Dimuon_muonPlusIndex);
    chain->SetBranchAddress("Reco_Dimuon_muonMinusIndex", Reco_Dimuon_muonMinusIndex);

    chain->SetBranchAddress("Reco_Dimuon_trig", Reco_Dimuon_trig);
    chain->SetBranchAddress("Reco_Dimuon_vtxProb", Reco_Dimuon_vtxProb);

    chain->SetBranchAddress("Reco_Dimuon_pt", &Reco_Dimuon_pt);
    chain->SetBranchAddress("Reco_Dimuon_eta", &Reco_Dimuon_eta);
    chain->SetBranchAddress("Reco_Dimuon_rapidity", &Reco_Dimuon_rapidity);
    chain->SetBranchAddress("Reco_Dimuon_phi", &Reco_Dimuon_phi);
    chain->SetBranchAddress("Reco_Dimuon_invMass", &Reco_Dimuon_invMass);

    chain->SetBranchAddress("Reco_Dimuon_muonPtDiff", &Reco_Dimuon_muonPtDiff);
    chain->SetBranchAddress("Reco_Dimuon_muonPtRelDiff", &Reco_Dimuon_muonPtRelDiff);

    //Muon-level variables
    Short_t Reco_Muon_size;

    std::vector<float>* Reco_Muon_pt = nullptr;
    //std::vector<float>* Reco_Muon_ptErrTrk = nullptr; //I think we need to study this branch further.
    std::vector<float>* Reco_Muon_eta = nullptr;
    std::vector<float>* Reco_Muon_phi = nullptr;
    std::vector<float>* Reco_Muon_mass = nullptr;

    ULong64_t Reco_Muon_trig[MAX_MUON];
    Bool_t Reco_Muon_isTightCutBased[MAX_MUON];
    
    chain->SetBranchAddress("Reco_Muon_size", &Reco_Muon_size);

    chain->SetBranchAddress("Reco_Muon_pt", &Reco_Muon_pt);
    //chain->SetBranchAddress("Reco_Muon_ptErrTrk", &Reco_Muon_ptErrTrk);
    chain->SetBranchAddress("Reco_Muon_eta", &Reco_Muon_eta);
    chain->SetBranchAddress("Reco_Muon_phi", &Reco_Muon_phi);
    chain->SetBranchAddress("Reco_Muon_mass", &Reco_Muon_mass);

    chain->SetBranchAddress("Reco_Muon_trig", Reco_Muon_trig);
    chain->SetBranchAddress("Reco_Muon_isTightCutBased", Reco_Muon_isTightCutBased);
    

    //Raw histograms. Store them in ROOT memory.
    gROOT->cd();
    std::string rawDirName = "rawHistograms_" + whichDataset;
    TDirectory* rawDir = gDirectory->mkdir(rawDirName.c_str());
        if (!rawDir) {
            std::cerr << "Nao foi possivel criar " << rawDirName << '\n';
            return;
        }
    rawDir->cd();

    //
    TH3D* h3D_PtMuPl_PtMuMi_Cent = nullptr;
    TH1D* h1D_centrality = nullptr;
    TH1D* h1D_hiHF = nullptr;

    if(dataset.hasCentrality){
        h3D_PtMuPl_PtMuMi_Cent = new TH3D("h3D_PtMuPl_PtMuMi_Cent",
        "p_{T}^{#mu^{+}} vs p_{T}^{#mu^{-}} vs Centrality; p_{T}^{#mu^{+}} [GeV/c]; p_{T}^{#mu^{-}} [GeV/c]; Centrality [%]",
        200, 0., 200., 200, 0., 200., 200, 0., 100.);

        h1D_centrality = new TH1D("h1D_centrality",
        "Centrality of selected Z candidates;Centrality [%];N of dimuons",
        200, 0., 100.);

        h1D_hiHF = new TH1D("h1D_hiHF",
        "hiHF for selected Z candidates;hiHF;N of dimuons",
        200, 0., 8000.);
    }
    
    TH1D* h1D_invMass = new TH1D("h1D_invMass",
        "Invariant Mass of selected dimuons; M_{#mu^{+}#mu^{-}} [GeV/c^{2}]; N of dimuons",
        120, 60., 120.);

    TH1D* h1D_zPt = new TH1D("h1D_zPt",
        "Transverse momentum of selected Z candidates;p_{T}^{Z} [GeV/c];N of dimuons",
        200, 0., 200.);

    TH1D* h1D_zRapidity = new TH1D("h1D_zRapidity",
        "Rapidity of selected Z candidates;y_{Z};N of dimuons",
        48, -2.4, 2.4);

    TH1D* h1D_ptMuPlus = new TH1D("h1D_ptMuPlus",
        "Transverse momentum of selected #mu^{+};p_{T}^{#mu^{+}} [GeV/c];N of muons",
        200, 0., 200.);

    TH1D* h1D_ptMuMinus = new TH1D("h1D_ptMuMinus",
        "Transverse momentum of selected #mu^{-};p_{T}^{#mu^{-}} [GeV/c];N of muons",
        200, 0., 200.);

    TH1D* h1D_etaMuPlus = new TH1D("h1D_etaMuPlus",
        "Pseudorapidity of selected #mu^{+};#eta^{#mu^{+}};N of muons",
        48, -2.4, 2.4);

    TH1D* h1D_etaMuMinus = new TH1D("h1D_etaMuMinus",
        "Pseudorapidity of selected #mu^{-};#eta^{#mu^{-}};N of muons",
        48, -2.4, 2.4);

    TH1D* h1D_phiMuPlus = new TH1D("h1D_phiMuPlus",
        "Azimuthal angle of selected #mu^{+};#phi^{#mu^{+}} [rad];N of muons",
        64, -TMath::Pi(), TMath::Pi());

    TH1D* h1D_phiMuMinus = new TH1D("h1D_phiMuMinus",
        "Azimuthal angle of selected #mu^{-};#phi^{#mu^{-}} [rad];N of muons",
        64, -TMath::Pi(), TMath::Pi());

    TH1D* h1D_zVtx = new TH1D("h1D_zVtx",
        "Primary vertex z position for selected Z candidates;z_{vtx} [cm];N of dimuons",
        60, -15., 15.);

    //Azimuthal separation between the two daughter muons.
    TH1D* h1D_deltaPhiMuMu = new TH1D("h1D_deltaPhiMuMu",
        "Azimuthal separation of selected dimuons;#Delta#phi(#mu^{+},#mu^{-}) [rad];N of dimuons",
        64, 0., TMath::Pi());

    //Angular separation between the two daughter muons. Distance between two points in the eta-phi space.
    TH1D* h1D_deltaRMuMu = new TH1D("h1D_deltaRMuMu",
        "Angular separation of selected dimuons;#DeltaR(#mu^{+},#mu^{-});N of dimuons",
        60, 0., 6.);

    //Acoplanarity between the two daughter muons.
    TH1D* h1D_acoplanarity = new TH1D("h1D_acoplanarity",
        "Acoplanarity of selected dimuons;1 - #Delta#phi(#mu^{+},#mu^{-})/#pi;N of dimuons",
        100, 0., 1.);

    //Relative pT difference between the two daughter muons.
    TH1D* h1D_muonPtRelDiff = new TH1D("h1D_muonPtRelDiff",
        "Relative p_{T} difference of selected dimuons;(p_{T}^{#mu^{+}}-p_{T}^{#mu^{-}})/(p_{T}^{#mu^{+}}+p_{T}^{#mu^{-}});N of dimuons",
        100, -1., 1.);

    //Absolute pT difference between the two daughter muons.
    TH1D* h1D_muonPtDiff = new TH1D("h1D_muonPtDiff",
        "p_{T} difference of selected dimuons;p_{T}^{#mu^{+}} - p_{T}^{#mu^{-}} [GeV/c];N of dimuons",
        200, -100., 100.);


    //Some correlation histograms i thought could be interesting to look at:
    TH2D* h2D_muonPtRelDiff_Cent = nullptr;
    TH2D* h2D_zPt_Cent = nullptr;

    if(dataset.hasCentrality){
        //Centrality vs muonPtRelDiff 
        h2D_muonPtRelDiff_Cent = new TH2D("h2D_muonPtRelDiff_Cent",
        "Relative p_{T} difference vs Centrality;Centrality [%];(p_{T}^{#mu^{+}}-p_{T}^{#mu^{-}})/(p_{T}^{#mu^{+}}+p_{T}^{#mu^{-}})",
            200, 0., 100., 100, -1., 1.);

        //Centrality vs Z pT 
        h2D_zPt_Cent = new TH2D("h2D_zPt_Cent",
        "Z p_{T} vs Centrality;Centrality [%];p_{T}^{Z} [GeV/c]",
            200, 0., 100., 200, 0., 200.);
    }

    //Z pT vs muonPtRelDiff
    TH2D* h2D_muonPtRelDiff_zPt = new TH2D("h2D_muonPtRelDiff_zPt",
    "Relative muon p_{T} difference vs Z p_{T};p_{T}^{Z} [GeV/c];(p_{T}^{#mu^{+}}-p_{T}^{#mu^{-}})/(p_{T}^{#mu^{+}}+p_{T}^{#mu^{-}})",
        200, 0., 200., 100, -1., 1.);

    //Z rapidity vs muonPtRelDiff 
    TH2D* h2D_muonPtRelDiff_zRapidity = new TH2D("h2D_muonPtRelDiff_zRapidity",
    "Relative muon p_{T} difference vs Z rapidity;y^{Z};(p_{T}^{#mu^{+}}-p_{T}^{#mu^{-}})/(p_{T}^{#mu^{+}}+p_{T}^{#mu^{-}})",
        48, -2.4, 2.4, 100, -1., 1.);

    //Z pT vs Z rapidity
    TH2D* h2D_zPt_zRapidity = new TH2D("h2D_zPt_zRapidity",
    "Z p_{T} vs Z rapidity;p_{T}^{Z} [GeV/c];y^{Z}",
        200, 0., 200., 48, -2.4, 2.4);

    //pT(mu+) vs pT(mu-)
    TH2D* h2D_PtMuPl_PtMuMi = new TH2D("h2D_PtMuPl_PtMuMi",
    "p_{T} of muon+ vs p_{T} of muon-;p_{T}^{#mu^{+}} [GeV/c];p_{T}^{#mu^{-}} [GeV/c]",
        200, 0., 200., 200, 0., 200.);


    //Start event-level loop
    for(Long64_t i = 0; i < nEvents; ++i){//Loop through all EVENTS in the CHAIN.

        chain->GetEntry(i); //Get event i.


        for(Short_t j = 0; j < Reco_Dimuon_size; ++j){ //Loop through all reco dimuon candidates of event i.
            

            Short_t muonPlusIndex = Reco_Dimuon_muonPlusIndex[j]; //Index of antimuon in the reco muon arrays.
            Short_t muonMinusIndex = Reco_Dimuon_muonMinusIndex[j]; //Index of corresponding muon in the reco muon arrays.

            double ptplus = Reco_Muon_pt->at(muonPlusIndex); //pT of antimuon.
            double ptminus = Reco_Muon_pt->at(muonMinusIndex); //pT of corresponding muon.
            double etaplus = Reco_Muon_eta->at(muonPlusIndex); //Pseudorapidity of antimuon.
            double etaminus = Reco_Muon_eta->at(muonMinusIndex); //Pseudorapidity of corresponding muon.


            //Simple histograms
            h1D_invMass->Fill(Reco_Dimuon_invMass->at(j));
            h1D_zPt->Fill(Reco_Dimuon_pt->at(j));
            h1D_zRapidity->Fill(Reco_Dimuon_rapidity->at(j));

            h1D_ptMuPlus->Fill(ptplus);
            h1D_ptMuMinus->Fill(ptminus);
            h1D_etaMuPlus->Fill(etaplus);
            h1D_etaMuMinus->Fill(etaminus);

            double phiplus = Reco_Muon_phi->at(muonPlusIndex);
            double phiminus = Reco_Muon_phi->at(muonMinusIndex);
            h1D_phiMuPlus->Fill(phiplus);
            h1D_phiMuMinus->Fill(phiminus);

            h1D_zVtx->Fill(zVtx);

            //Azimuthal separation between the two daughter muons.
            double deltaPhiMuMu = std::abs( TVector2::Phi_mpi_pi(phiplus - phiminus) );
            h1D_deltaPhiMuMu->Fill(deltaPhiMuMu);

            //Angular separation between the two daughter muons. Distance between two points in the eta-phi space.
            double deltaEtaMuMu = etaplus - etaminus;
            double deltaRMuMu = std::sqrt(deltaEtaMuMu*deltaEtaMuMu + deltaPhiMuMu*deltaPhiMuMu);
            h1D_deltaRMuMu->Fill(deltaRMuMu);

            //Acoplanarity between the two daughter muons.
            double acoplanarity = 1.0 - ( deltaPhiMuMu / TMath::Pi() );
            h1D_acoplanarity->Fill(acoplanarity);

            //Relative pT difference between the two daughter muons.
            h1D_muonPtRelDiff->Fill(Reco_Dimuon_muonPtRelDiff->at(j));

            //Absolute pT difference between the two daughter muons.
            h1D_muonPtDiff->Fill(Reco_Dimuon_muonPtDiff->at(j));

            if (dataset.hasCentrality) {
                h3D_PtMuPl_PtMuMi_Cent->Fill(ptplus, ptminus, Centrality/2.);
                h1D_centrality->Fill(Centrality/2.);
                h1D_hiHF->Fill(SumET_HF);

                //Latest correlation histograms
                h2D_muonPtRelDiff_Cent->Fill(Centrality/2., Reco_Dimuon_muonPtRelDiff->at(j));
                h2D_zPt_Cent->Fill(Centrality/2., Reco_Dimuon_pt->at(j));
            }

            //Latest correlation histograms
            h2D_muonPtRelDiff_zPt->Fill(Reco_Dimuon_pt->at(j), Reco_Dimuon_muonPtRelDiff->at(j));
            h2D_muonPtRelDiff_zRapidity->Fill(Reco_Dimuon_rapidity->at(j), Reco_Dimuon_muonPtRelDiff->at(j));
            h2D_zPt_zRapidity->Fill(Reco_Dimuon_pt->at(j), Reco_Dimuon_rapidity->at(j));
            h2D_PtMuPl_PtMuMi->Fill(ptplus, ptminus);

        }//End of dimuon candidate loop.


        //Track progress of event loop.
        Long64_t progressStep = std::max<Long64_t>(1, nEvents / 100);
        if (i % progressStep == 0){

            int percent = static_cast<int>(100.0 * i / nEvents+0.5);

            std::cout << "\r"
                      << percent
                      << "% complete..."
                      << std::flush;
        }

    }//Exiting event-by-event loop.

    std::cout << "\n\n> Selections applied on " << dataset.name << "\n" << std::endl;
    std::cout << "> Writing histograms to output file..." << "\n" << std::endl;
    std::cout << "> Number of Z bosons selected = " << h1D_invMass->GetEntries() << std::endl;

    //Put raw histograms in the vector of HistStructs
    for(const auto& histName : histNames) {

        HistStruct histStruct;
        histStruct.name = histName;
        histStruct.rawHist = dynamic_cast<TH1*>(rawDir->Get(histName.c_str()));

        TDirectory *dir = inputFile->GetDirectory(whichDataset.c_str());
        histStruct.selectedHist = dynamic_cast<TH1*>(dir->Get(histName.c_str()));
        histStructs.push_back(histStruct);

        //Now all the raw and selected histograms are in the vector of HistStructs. We can now plot them on the same canvas.
        MakePlot(histStruct, whichDataset);
    }

    gROOT->cd();
    delete rawDir;
    delete chain;
}


void MakePlot(HistStruct& histStruct, std::string whichDataset){

    basicHistFormatting(histStruct.rawHist);
    basicHistFormatting(histStruct.selectedHist);

    TCanvas *c = new TCanvas("c", "c", 800, 600);
    basicCanvasFormatting(c);
    
    histStruct.rawHist->SetLineColor(kBlack);
    histStruct.rawHist->SetLineWidth(2);
    
    histStruct.selectedHist->SetLineColor(kOrange);
    histStruct.selectedHist->SetLineWidth(2);

    histStruct.rawHist->Draw("HIST");
    histStruct.selectedHist->Draw("HIST SAME");
    
    TLegend *leg = new TLegend(0.7, 0.78, 0.95, 0.88);
    basicLegendFormatting(leg);
    leg->AddEntry(histStruct.rawHist, "Raw", "l");
    leg->AddEntry(histStruct.selectedHist, "Selected", "l");
    leg->Draw();

    drawLatexText("#bf{CMS}", 0.12, 0.93, 0.042);
    drawLatexText("#it{Work in Progress}", 0.2, 0.93, 0.033);
    drawLatexText(whichDataset.c_str(), 0.5, 0.93, 0.033);

    c->Update();
    std::string outputName = whichDataset + "_" + histStruct.name + plot_extension;
    c->SaveAs(outputName.c_str());

    delete leg;
    delete c;
}