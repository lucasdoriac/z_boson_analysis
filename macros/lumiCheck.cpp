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

using namespace std;

//---Macro settings

std::string plot_extension = ".png"; // ".png" for regular development and ".pdf" for final quality plots
std::string BasePath = "/home/lucas/Documents/CMS_analyzes/Z_boson_analysis/"; //IFT
//std::string BasePath = "/home/lucasdoriac/z_boson_analysis/data/"; //Home

const int MIN_CENTRALITY = 0.;
const int MAX_CENTRALITY = 180;
const double MAX_ZVTX = 15.0;

const double MINZ_MASS = 80.;
const double MAXZ_MASS = 100.;
const double RAPIDITYCUTVALUE = 2.4;

const float ETACUTVALUE = 2.4;
const double PTCUTVALUE = 10.;

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

void checkppRef(const Dataset& dataset);

void lumiCheck(){

    gROOT->SetBatch(kTRUE);

    checkppRef(datasets[1]);
}

void checkppRef(const Dataset& dataset){

    // Load root file.
    std::string fullPath = dataset.basePath + dataset.filePattern;

    TChain *chain = new TChain(dataset.treeName.c_str());
    chain->Add(fullPath.c_str());

    std::cout << "> Number of files = " << chain->GetListOfFiles()->GetEntries() << "\n" << std::endl;

    std::cout << "> Opening files " << fullPath << "\n" << std::endl;

    std::cout << "> Running function " << __func__ << " on " << dataset.name << "\n" << std::endl;
    
    //Total number of events on Tree.
    Long64_t nEvents = chain->GetEntries();

    UInt_t eventNb;
    UInt_t runNb;
    UInt_t LS;

    chain->SetBranchAddress("eventNb", &eventNb);
    chain->SetBranchAddress("runNb", &runNb);
    chain->SetBranchAddress("LS", &LS);

    Short_t Reco_Dimuon_size;
    Short_t Reco_Muon_size;
    chain->SetBranchAddress("Reco_Dimuon_size", &Reco_Dimuon_size);
    chain->SetBranchAddress("Reco_Muon_size", &Reco_Muon_size);

    TH2D* h_ZCandidates = new TH2D(
    "h_ZCandidates",
    "Z candidates per run and LS;LS;Run",
    200, 0, 2000,
    50, 0, 50
);

    std::set<UInt_t> runs;

for (Long64_t i = 0; i < nEvents; ++i) {

    chain->GetEntry(i);

    runs.insert(runNb);
}

std::map<UInt_t, int> runIndex;

int index = 1;

for (UInt_t run : runs) {
    runIndex[run] = index;
    index++;
}

    std::map<std::pair<UInt_t, UInt_t>, int> nZCandidates;
    
    for (Long64_t i = 0; i < nEvents; ++i) {

        chain->GetEntry(i);

        auto key = std::make_pair(runNb, LS);

        for (int j = 0; j < Reco_Dimuon_size; ++j) {

            nZCandidates[key]++;
        }

    }

    for (const auto& entry : nZCandidates) {

    UInt_t run = entry.first.first;
    UInt_t ls  = entry.first.second;

    int nZ = entry.second;

    int yBin = runIndex[run];

    h_ZCandidates->SetBinContent(
        h_ZCandidates->FindBin(ls, yBin),
        nZ
    );
}

    TCanvas *c = new TCanvas("c", "", 600, 600);

    h_ZCandidates->Draw("COLZ");
    h_ZCandidates->SetStats(0);
    c->Update();
    c->SaveAs("lumiCheckPbPb.png");

    delete chain;
}

