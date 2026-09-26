/*
Plots raw vs selected variables in the same canvas.

Raw: from the processed TTree, before any selection is applied.
Selected: from 'mySelectedData.root', after applying the selection criteria.

*/


//List of histograms to plot raw vs selected.
std::vector<std::string> histNames = {

    "h3D_PtMuPl_PtMuMi_Cent",
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
    "h1D_muonPtDiff",

    "h2D_muonPtRelDiff_Cent",
    "h2D_zPt_Cent",

    "h2D_muonPtRelDiff_zPt",
    "h2D_muonPtRelDiff_zRapidity",
    "h2D_zPt_zRapidity",
    "h2D_PtMuPl_PtMuMi"
};

struct HistStruct{
    std::string name;
    TH1* rawHist;
    TH1* selectedHist;
};

std::vector<HistStruct> histStructs;

void foo(TFile* inputFile);


void main(){

    gROOT->SetBatch(kTRUE);
    TFile* inputFile = new TFile("mySelectedData.root", "READ");

    foo(inputFile);


    inputFile->Close();
}




void foo(TFile* inputFile){

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
    

    //Histograms to be saved in the ROOT mySelectedData file.
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

            ptplus = Reco_Muon_pt->at(muonPlusIndex); //pT of antimuon.
            ptminus = Reco_Muon_pt->at(muonMinusIndex); //pT of corresponding muon.
            etaplus = Reco_Muon_eta->at(muonPlusIndex); //Pseudorapidity of antimuon.
            etaminus = Reco_Muon_eta->at(muonMinusIndex); //Pseudorapidity of corresponding muon.


            //Simple histograms
            h1D_invMass->Fill(Reco_Dimuon_invMass->at(j));
            h1D_zPt->Fill(Reco_Dimuon_pt->at(j));
            h1D_zRapidity->Fill(Reco_Dimuon_rapidity->at(j));

            h1D_ptMuPlus->Fill(ptplus);
            h1D_ptMuMinus->Fill(ptminus);
            h1D_etaMuPlus->Fill(etaplus);
            h1D_etaMuMinus->Fill(etaminus);

            phiplus = Reco_Muon_phi->at(muonPlusIndex);
            phiminus = Reco_Muon_phi->at(muonMinusIndex);
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
        histStruct.rawHist = dynamic_cast<TH1*>(gDirectory->Get(histName.c_str()));
        histStruct.selectedHist = dynamic_cast<TH1*>(inputFile->Get(histName.c_str()));
        histStructs.push_back(histStruct);
    }


    //Selected
    TDirectory *dir = inputFile->GetDirectory(whichDataset.c_str());
    TH1 *hist = 



    for(const auto& histName : selectedDistributions) {

        //Get original histograms
        TH1D* h_PbPb_original = dynamic_cast<TH1D*>(PbPb_dir->Get(histName.c_str()));
        TH1D* h_ppRef_original = dynamic_cast<TH1D*>(ppRef_dir->Get(histName.c_str()));

            if (!h_PbPb_original || !h_ppRef_original) {//Just checking if everything was found.
                std::cerr << "Error: Could not find the histogram "
                        << histName << " in the input file."
                        << std::endl;
                continue;
            }
        
        //Get histogram clones to manipulate.
        TH1D* h_PbPb = dynamic_cast<TH1D*>(h_PbPb_original->Clone(("h_PbPb_" + histName).c_str()));
        TH1D* h_ppRef = dynamic_cast<TH1D*>(h_ppRef_original->Clone(("h_ppRef_" + histName).c_str()));
        h_PbPb->SetDirectory(nullptr);
        h_ppRef->SetDirectory(nullptr);


        //Begin normalization and stuff
        h_PbPb->Scale(1.0 / h_PbPb->Integral());
        h_ppRef->Scale(1.0 / h_ppRef->Integral());

        std::string canvasName = "c_" + histName;
        TCanvas *c = new TCanvas(canvasName.c_str(), "Normalized Distributions", 800, 600);
        basicCanvasFormatting(c);
        c->SetLogy();

            if(histName == "h1D_muonPtRelDiff"){//Turn off log scale for these two distributions.
                c->SetLogy(0);
            }

        basicHistFormatting(h_PbPb);
        basicHistFormatting(h_ppRef);
        
        h_PbPb->SetMarkerStyle(21);
        h_PbPb->SetMarkerSize(0.8);
        h_PbPb->SetMarkerColor(kRed);
        h_PbPb->SetLineColor(kRed);
        
        h_ppRef->SetMarkerStyle(25);
        h_ppRef->SetMarkerSize(0.8);
        h_ppRef->SetMarkerColor(kBlack);
        h_ppRef->SetLineColor(kBlack);

        h_PbPb->GetYaxis()->SetTitle("Normalized Entries");

            if(histName == "h1D_ptMuPlus" || histName == "h1D_ptMuMinus"){
                h_PbPb->GetXaxis()->SetRangeUser(18., 100.);
            }
        
        h_PbPb->Draw("P");
        h_ppRef->Draw("P SAME");
        
        TLegend *leg = new TLegend(0.7, 0.78, 0.95, 0.88);
        basicLegendFormatting(leg);
        leg->AddEntry(h_PbPb, JointPbPb.c_str(), "p");
        leg->AddEntry(h_ppRef, "ppRef2024", "p");
        leg->Draw();

        drawLatexText("#bf{CMS}", 0.12, 0.93, 0.042);
        drawLatexText("#it{Work in Progress}", 0.2, 0.93, 0.033);
        drawLatexText(dataSamplesUsed.c_str(), 0.5, 0.93, 0.033);

        c->Update();
        std::string outputName = "Normalized_Distributions_Joined_" + histName + plot_extension;
        c->SaveAs(outputName.c_str());

        delete leg;
        delete h_PbPb;
        delete h_ppRef;
        delete c;
    }
    

}