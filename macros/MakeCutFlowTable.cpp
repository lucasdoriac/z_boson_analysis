/*
cut flow table.
*/

//---Libraries
#include <TROOT.h>
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
std::string BasePath = "/home/lucas/Documents/CMS/z_boson_analysis/"; //IFT
//std::string BasePath = "/home/lucasdoriac/z_boson_analysis/data/"; //Home


//Good Selection values
const double maxZvtx = 15.0;

const double minZ_Mass = 60.;
const double maxZ_Mass = 120.;
const double RapidityCutValue = 2.4;

const double EtaCutValue = 2.4;
const double ptCutValue = 20.;


// ##############################################################################
// ##############################################################################


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


//Flow table as a struct to hold the counts for each cut stage.
struct FlowTable {
    std::string datasetName;
    Long64_t nEvents = 0;
    Long64_t nGoodZvtx = 0;
    Long64_t nRawDimuons = 0;
    Long64_t nGoodMass = 0;
    Long64_t nGoodRapidity = 0;
    Long64_t nGoodCharge = 0;
    Long64_t nGoodVtxProb = 0;
    Long64_t nTriggerMatched = 0;
    Long64_t nGoodDimuons = 0;
    Long64_t nMultiDimuonEvents = 0;
};

std::vector<FlowTable> flowTablesVector; //Vector to hold flow tables for each dataset.

void MakeTableFromTree(const Dataset& dataset);
void PrintFlowTables();

void MakeCutFlowTable(){

    gROOT->SetBatch(kTRUE);
    flowTablesVector.clear(); //Clear the vector before starting.

    MakeTableFromTree(datasets[0]); //ppRef2024
    MakeTableFromTree(datasets[1]); //PbPb2023
    MakeTableFromTree(datasets[2]); //PbPb2024
    //MakeTableFromTree(datasets[3]); //PbPb2025
    //MakeTableFromTree(datasets[4]); //PbPb2026

    //Print the flow tables.
    PrintFlowTables();
}

void PrintFlowTables(){

    std::cout << "\n\n> Cut Flow Tables Summary: (Absolute values)\n" << std::endl;

    for(const auto& flowTable : flowTablesVector){
        std::cout << "Dataset: " << flowTable.datasetName << std::endl;
        std::cout << "Total events: " << flowTable.nEvents << std::endl;
        std::cout << "Events with good zVtx: " << flowTable.nGoodZvtx << std::endl;
        std::cout << "Total dimuon candidates: " << flowTable.nRawDimuons << std::endl;
        std::cout << "Dimuons with good mass: " << flowTable.nGoodMass << std::endl;
        std::cout << "Dimuons with good rapidity: " << flowTable.nGoodRapidity << std::endl;
        std::cout << "Dimuons with good charge: " << flowTable.nGoodCharge << std::endl;
        std::cout << "Dimuons with good vertex probability: " << flowTable.nGoodVtxProb << std::endl;
        std::cout << "Dimuons with trigger matched: " << flowTable.nTriggerMatched << std::endl;
        std::cout << "Number of good dimuons: " << flowTable.nGoodDimuons << std::endl;
        std::cout << "Events with multiple dimuon candidates: " << flowTable.nMultiDimuonEvents << "\n" << std::endl;

        std::cout << "\n\n> % kept from previous cut: \n" << std::endl;
        
        std::cout << "\n\nDataset: " << flowTable.datasetName << std::endl;
        std::cout << "Total events: " << flowTable.nEvents << std::endl;
        std::cout << "Events with good zVtx: " << static_cast<double>(flowTable.nGoodZvtx)*100./flowTable.nEvents << std::endl;
        std::cout << "Total dimuon candidates: " << static_cast<double>(flowTable.nRawDimuons)*100./flowTable.nGoodZvtx << std::endl;
        std::cout << "Dimuons with good mass: " << static_cast<double>(flowTable.nGoodMass)*100./flowTable.nRawDimuons << std::endl;
        std::cout << "Dimuons with good rapidity: " << static_cast<double>(flowTable.nGoodRapidity)*100./flowTable.nGoodMass << std::endl;
        std::cout << "Dimuons with good charge: " << static_cast<double>(flowTable.nGoodCharge)*100./flowTable.nGoodRapidity << std::endl;
        std::cout << "Dimuons with good vertex probability: " << static_cast<double>(flowTable.nGoodVtxProb)*100./flowTable.nGoodCharge << std::endl;
        std::cout << "Dimuons with trigger matched: " << static_cast<double>(flowTable.nTriggerMatched)*100./flowTable.nGoodVtxProb << std::endl;
        std::cout << "Number of good dimuons: " << static_cast<double>(flowTable.nGoodDimuons)*100./flowTable.nTriggerMatched << std::endl;


        std::cout << "\n> % kept from previous cut: (Cummulative) \n" << std::endl;
        
        std::cout << "Dataset: " << flowTable.datasetName << std::endl;
        std::cout << "Total events: " << flowTable.nEvents << std::endl;
        std::cout << "Events with good zVtx: " << static_cast<double>(flowTable.nGoodZvtx)*100./flowTable.nEvents << std::endl;
        std::cout << "Total dimuon candidates: " << static_cast<double>(flowTable.nRawDimuons)*100./flowTable.nGoodZvtx << std::endl;
        std::cout << "Dimuons with good mass: " << static_cast<double>(flowTable.nGoodMass)*100./flowTable.nRawDimuons << std::endl;
        std::cout << "Dimuons with good rapidity: " << static_cast<double>(flowTable.nGoodRapidity)*100./flowTable.nRawDimuons << std::endl;
        std::cout << "Dimuons with good charge: " << static_cast<double>(flowTable.nGoodCharge)*100./flowTable.nRawDimuons << std::endl;
        std::cout << "Dimuons with good vertex probability: " << static_cast<double>(flowTable.nGoodVtxProb)*100./flowTable.nRawDimuons << std::endl;
        std::cout << "Dimuons with trigger matched: " << static_cast<double>(flowTable.nTriggerMatched)*100./flowTable.nRawDimuons << std::endl;
        std::cout << "Number of good dimuons: " << static_cast<double>(flowTable.nGoodDimuons)*100./flowTable.nRawDimuons << std::endl;
        
        std::cout << "\n\n------------------------------------------------------------\n" << std::endl;
    }

}


void MakeTableFromTree(const Dataset& dataset){

    //Load root file.
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
    
    //Flow table to hold the counts for each cut stage.
    FlowTable flowTable;
    flowTable.datasetName = dataset.name;
    flowTable.nEvents = nEvents;

    //Muon-level selection variables
    double ptplus, ptminus;
    double etaplus, etaminus;
    double phiplus, phiminus;
    bool MuPlIsTight;
    bool MuMiIsTight;

    Long64_t nMultiDimuonEvents = 0; //Tracker of n of events with more than one dimuon candidate passing the selection.

    for(Long64_t i = 0; i < nEvents; ++i){//Loop through all EVENTS in the CHAIN.

        chain->GetEntry(i); //Get event i.
        int nSelectedInEvent = 0; //Tracker of n of dimuon candidates passing the selection in a single event.

        //Good event selection. No centrality selection at this point.
        bool goodVertex = (std::abs(zVtx) < maxZvtx);

        if (!goodVertex) continue;
        flowTable.nGoodZvtx++;

        for(Short_t j = 0; j < Reco_Dimuon_size; ++j){ //Loop through all reco dimuon candidates of event i.
            
            //Good Z selection
            bool goodMass = (Reco_Dimuon_invMass->at(j) > minZ_Mass && Reco_Dimuon_invMass->at(j) < maxZ_Mass);
            bool goodRapidity = (std::abs(Reco_Dimuon_rapidity->at(j)) < RapidityCutValue);
            bool goodCharge = (Reco_Dimuon_sign[j] == 0); //Opposite sign muons.
            bool goodVtxProb = (Reco_Dimuon_vtxProb[j] > 0.001); //Vertex probability cut of .1% for dimuon candidates.
            bool isTriggerMatched = true;

                if (dataset.applyTrigger) {
                    //**At least one** of the daughter muons must be matched to the trigger.
                    isTriggerMatched = (Reco_Dimuon_trig[j] & dataset.triggerBit);
                }

            flowTable.nRawDimuons++; //Count all dimuon candidates before any selection.
            if (!goodMass) continue;
            flowTable.nGoodMass++;
            if (!goodRapidity) continue;
            flowTable.nGoodRapidity++;
            if (!goodCharge) continue;
            flowTable.nGoodCharge++;
            if (!goodVtxProb) continue;
            flowTable.nGoodVtxProb++;
            if (!isTriggerMatched) continue;
            flowTable.nTriggerMatched++;

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
            flowTable.nGoodDimuons++;
            //End of good selection for dimuon candidate j of event i.

            //Sept. 24, 2026: Temporary.
            //Keep track of the number of dimuon candidates passing the selection in this event.
            //We will have to deal with the case of more than one dimuon candidate passing the selection in a single event in the future.
            nSelectedInEvent++;

        }//End of dimuon candidate loop.

        if(nSelectedInEvent == 0) continue; //Skip to next event if no dimuon candidate passed the selection in this event.

        if(nSelectedInEvent > 1) flowTable.nMultiDimuonEvents++;

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


    std::cout << "\n\n> Made table flow for " << dataset.name << "\n" << std::endl;
    std::cout << "> Number of selected dimuon candidates = " << flowTable.nGoodDimuons << std::endl;
    std::cout << "> Number of events with multiple dimuon candidates = " << flowTable.nMultiDimuonEvents << "\n\n" << std::endl;

    //Add the flow table for this dataset to the vector.
    flowTablesVector.push_back(flowTable);

    delete chain;
}