/*
Create a few histograms of interest for the PbPb2024 and ppRef2024 datasets after applying a definite set of good selections.
*/

//---Libraries
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TString.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
#include <TVector2.h>
#include <algorithm>


//Location of datasets
//std::string BasePath = "/home/lucas/Documents/CMS/z_boson_analysis/"; //IFT
//std::string BasePath = "/home/lucasdoriac/z_boson_analysis/data/"; //Home

//Good Selection values
const double maxZvtx = 15.0;

const double minZ_Mass = 60.;
const double maxZ_Mass = 120.;
const double RapidityCutValue = 2.4;

const double EtaCutValue = 2.4;
const double ptCutValue = 20.;

ULong64_t triggerBit = 1ULL << 7; //Trigger bit for L2SingleMu12 trigger selection.


// ##############################################################################
// ##############################################################################

//---Enumerates
enum class SampleType {
    Data,
    MC
};

enum class CollisionSystem {
    PbPb2023,
    PbPb2024,
    ppRef2024
};

//---Structs
struct Dataset {
    std::string name;
    SampleType type;
    CollisionSystem system;
    std::string treeName;
    std::string filePattern;
    std::string basePath;
};

Dataset datasets[] = {
    {
        "PbPb2023_Data",
        SampleType::Data,
        CollisionSystem::PbPb2023,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_PbPb2023.root",
        BasePath + "Data/PbPb2023/"
    },

    {
        "PbPb2024_Data",
        SampleType::Data,
        CollisionSystem::PbPb2024,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_PbPb2024Data.root",
        BasePath + "Data/PbPb2024/"
    },

    {
        "ppRef2024_Data",
        SampleType::Data,
        CollisionSystem::ppRef2024,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_ppRef2024.root",
        BasePath + "Data/ppRef2024/"
    },

    {
        "PbPb2024_MC",
        SampleType::MC,
        CollisionSystem::PbPb2024,
        "hionia/myTree",
        "Oniatree_PowhegZtoMuMu_PbPb2024_*.root",
        BasePath + "MC/PbPb2024/DYto2Mu_MLL-50_TuneCP5_5p36TeV_powheg-pythia8/PowhegEmbedded_March9/260309_143939/0000/"
    },

    {
        "ppRef2024_MC",
        SampleType::MC,
        CollisionSystem::ppRef2024,
        "hionia/myTree",
        "Oniatree_PowhegZtoMuMu_ppRef2024_*.root",
        BasePath + "MC/ppRef2024/DYToMuMu_M-50_TuneCP5_5p36TeV_powheg-pythia8/Powheg_ppRefPileup_March20/260320_125046/0000/"
    }
};


//---Function declarations
void makeGoodSelection(const Dataset& dataset, TFile* outputFile);

//---Main
void ApplyGoodSelection(){

    //Currently opens PbPb2024 and ppRef2024 datasets and applies set of good selections.
    //The selected data is stored on another ROOT file containing histograms of interest.
    TFile *outputFile = new TFile("mySelectedData.root", "RECREATE");

    makeGoodSelection(datasets[1], outputFile); //PbPb2024
    makeGoodSelection(datasets[2], outputFile); //ppRef2024

    outputFile->Close();
}

void makeGoodSelection(const Dataset& dataset, TFile* outputFile){

    // Load root file.
    std::string fullPath = dataset.basePath + dataset.filePattern;

    TChain *chain = new TChain(dataset.treeName.c_str());
    chain->Add(fullPath.c_str());

    std::cout << "\n> Number of files = " << chain->GetListOfFiles()->GetEntries() << "\n" << std::endl;

    std::cout << "> Opening files " << fullPath << "\n" << std::endl;

    std::cout << "> Running function " << __func__ << " on " << dataset.name << "\n" << std::endl;
    
    //Total number of events on Tree.
    Long64_t nEvents = chain->GetEntries();

    std::cout << "> Total number of events on tree = " << nEvents << "\n" << std::endl;

    //Create a directory for this dataset in the output ROOT file.
    TDirectory* datasetDir = outputFile->mkdir(dataset.name.c_str());
    datasetDir->cd();

    const int MAX_DIMUON = 1000;
    const int MAX_MUON   = 1000;

    //For now, writing only a fraction of the branches i.e. writing only what i'll use.

    //Event-level variables
    Int_t Centrality;
    Float_t  zVtx;

    chain->SetBranchAddress("zVtx", &zVtx);
    
    if(dataset.system == CollisionSystem::PbPb2024) { //Centrality is only defined for PbPb2024 dataset.
        chain->SetBranchAddress("Centrality", &Centrality);
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

    if(dataset.system == CollisionSystem::PbPb2024){
        h3D_PtMuPl_PtMuMi_Cent = new TH3D("h3D_PtMuPl_PtMuMi_Cent",
        "p_{T}^{#mu^{+}} vs p_{T}^{#mu^{-}} vs Centrality; p_{T}^{#mu^{+}} [GeV/c]; p_{T}^{#mu^{-}} [GeV/c]; Centrality [%]",
        100, 0., 100., 100, 0., 100., 200, 0., 100.);

        h1D_centrality = new TH1D("h1D_centrality",
        "Centrality of selected Z candidates;Centrality [%];N of dimuons",
        200, 0., 100.);
    }
    
    TH1D* h1D_invMass = new TH1D("h1D_invMass",
        "Invariant Mass of selected dimuons; M_{#mu^{+}#mu^{-}} [GeV/c^{2}]; N of dimuons",
        120, 60., 120.);

    TH1D* h1D_zPt = new TH1D("h1D_zPt",
        "Transverse momentum of selected Z candidates;p_{T}^{Z} [GeV/c];N of dimuons",
        100, 0., 200.);

    TH1D* h1D_zRapidity = new TH1D("h1D_zRapidity",
        "Rapidity of selected Z candidates;y_{Z};N of dimuons",
        48, -2.4, 2.4);

    TH1D* h1D_ptMuPlus = new TH1D("h1D_ptMuPlus",
        "Transverse momentum of selected #mu^{+};p_{T}^{#mu^{+}} [GeV/c];N of muons",
        100, 0., 100.);

    TH1D* h1D_ptMuMinus = new TH1D("h1D_ptMuMinus",
        "Transverse momentum of selected #mu^{-};p_{T}^{#mu^{-}} [GeV/c];N of muons",
        100, 0., 100.);

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

    if(dataset.system == CollisionSystem::PbPb2024){//Centrality is only defined for PbPb2024 dataset.
        //muonPtRelDiff vs Centrality
        h2D_muonPtRelDiff_Cent = new TH2D("h2D_muonPtRelDiff_Cent",
        "Relative p_{T} difference vs Centrality;Centrality [%];(p_{T}^{#mu^{+}}-p_{T}^{#mu^{-}})/(p_{T}^{#mu^{+}}+p_{T}^{#mu^{-}})",
            200, 0., 100.,100, -1., 1.);

        //Z pT vs Centrality
        h2D_zPt_Cent = new TH2D("h2D_zPt_Cent",
        "Z p_{T} vs Centrality;Centrality [%];p_{T}^{Z} [GeV/c]",
            200, 0., 100., 100, 0., 200.);
    }

    //muonPtRelDiff vs Z pT
    TH2D* h2D_muonPtRelDiff_zPt = new TH2D("h2D_muonPtRelDiff_zPt",
    "Relative muon p_{T} difference vs Z p_{T};p_{T}^{Z} [GeV/c];(p_{T}^{#mu^{+}}-p_{T}^{#mu^{-}})/(p_{T}^{#mu^{+}}+p_{T}^{#mu^{-}})",
        100, 0., 200., 100, -1., 1.);

    //muonPtRelDiff vs Z rapidity
    TH2D* h2D_muonPtRelDiff_zRapidity = new TH2D("h2D_muonPtRelDiff_zRapidity",
    "Relative muon p_{T} difference vs Z rapidity;y^{Z};(p_{T}^{#mu^{+}}-p_{T}^{#mu^{-}})/(p_{T}^{#mu^{+}}+p_{T}^{#mu^{-}})",
        48, -2.4, 2.4, 100, -1., 1.);

    //Z pT vs Z rapidity
    TH2D* h2D_zPt_zRapidity = new TH2D("h2D_zPt_zRapidity",
    "Z p_{T} vs Z rapidity;p_{T}^{Z} [GeV/c];y^{Z}",
        100, 0., 200., 48, -2.4, 2.4);

    //pT(mu+) vs pT(mu-)
    TH2D* h2D_PtMuPl_PtMuMi = new TH2D("h2D_PtMuPl_PtMuMi",
    "p_{T} of muon+ vs p_{T} of muon-;p_{T}^{#mu^{+}} [GeV/c];p_{T}^{#mu^{-}} [GeV/c]",
        100, 0., 100., 100, 0., 100.);

    
    //Muon-level selection variables
    double ptplus, ptminus;
    double etaplus, etaminus;
    double phiplus, phiminus;
    bool MuPlIsTight;
    bool MuMiIsTight;

    for(Long64_t i = 0; i < nEvents; ++i){//Loop through all EVENTS in the CHAIN.

        chain->GetEntry(i); //Get event i.

        //Good event selection. No centrality selection/accepting all centralities from 0 to 100% for PbPb2024 dataset.
        bool goodVertex = (std::abs(zVtx) < maxZvtx);

        if (!goodVertex) continue;

        for(Short_t j = 0; j < Reco_Dimuon_size; ++j){ //Loop through all reco dimuon candidates of event i.
            
            //Good Z selection
            bool goodMass = (Reco_Dimuon_invMass->at(j) > minZ_Mass && Reco_Dimuon_invMass->at(j) < maxZ_Mass);
            bool goodRapidity = (std::abs(Reco_Dimuon_rapidity->at(j)) < RapidityCutValue);
            bool goodCharge = (Reco_Dimuon_sign[j] == 0); //Opposite sign muons.
            bool goodVtxProb = (Reco_Dimuon_vtxProb[j] > 0.001); //Vertex probability cut of .1% for dimuon candidates.
            bool isTriggerMatched = true;

                if (dataset.system == CollisionSystem::PbPb2024){
                    isTriggerMatched = (Reco_Dimuon_trig[j] & triggerBit);//**At least** one of the daughter muons must be matched to the trigger.
                }

            if (!goodMass) continue;
            if (!goodRapidity) continue;
            if (!goodCharge) continue;
            if (!goodVtxProb) continue;
            if (!isTriggerMatched) continue;

            //Good muon selection
            Short_t muonPlusIndex = Reco_Dimuon_muonPlusIndex[j]; //Index of antimuon in the reco muon arrays.
            Short_t muonMinusIndex = Reco_Dimuon_muonMinusIndex[j]; //Index of corresponding muon in the reco muon arrays.

            ptplus = Reco_Muon_pt->at(muonPlusIndex); //pT of antimuon.
            ptminus = Reco_Muon_pt->at(muonMinusIndex); //pT of corresponding muon.
            etaplus = Reco_Muon_eta->at(muonPlusIndex); //Pseudorapidity of antimuon.
            etaminus = Reco_Muon_eta->at(muonMinusIndex); //Pseudorapidity of corresponding muon.
            MuPlIsTight = Reco_Muon_isTightCutBased[muonPlusIndex];
            MuMiIsTight = Reco_Muon_isTightCutBased[muonMinusIndex];
            
            bool goodMuPl = (ptplus > ptCutValue)
                            && (std::abs(etaplus) < EtaCutValue)
                            && (MuPlIsTight);

            bool goodMuMi = (ptminus > ptCutValue)
                            && (std::abs(etaminus) < EtaCutValue)
                            && (MuMiIsTight);

            if (!goodMuPl || !goodMuMi) continue;
            //End of good selection for dimuon candidate j of event i.

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

            if (dataset.system == CollisionSystem::PbPb2024) {
                h3D_PtMuPl_PtMuMi_Cent->Fill(ptplus, ptminus, Centrality/2.);
                h1D_centrality->Fill(Centrality/2.);

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

    //Write everything on the current directory.
    if(dataset.system == CollisionSystem::PbPb2024){//Only write for PbPb2024 dataset.
        h3D_PtMuPl_PtMuMi_Cent->Write(); 
        h1D_centrality->Write();

        //Latest correlation histograms
        h2D_muonPtRelDiff_Cent->Write();
        h2D_zPt_Cent->Write();
    }
    h1D_invMass->Write();
    h1D_zPt->Write();
    h1D_zRapidity->Write();
    h1D_ptMuPlus->Write();
    h1D_ptMuMinus->Write();
    h1D_etaMuPlus->Write();
    h1D_etaMuMinus->Write();
    h1D_phiMuPlus->Write();
    h1D_phiMuMinus->Write();
    h1D_zVtx->Write();
    h1D_deltaPhiMuMu->Write();
    h1D_deltaRMuMu->Write();
    h1D_acoplanarity->Write();
    h1D_muonPtRelDiff->Write();
    h1D_muonPtDiff->Write();
    
    //Latest correlation histograms
    h2D_muonPtRelDiff_zPt->Write();
    h2D_muonPtRelDiff_zRapidity->Write();
    h2D_zPt_zRapidity->Write();
    h2D_PtMuPl_PtMuMi->Write();
    
    //Leave directory.
    outputFile->cd();
}