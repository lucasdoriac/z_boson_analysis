/*
I want to check two things in this macro:

1- Usage of classes and structs to eliminate 200 lines of hardcoding variable initialization and branch addressing.
2- Make a plot of pT asymmetry that is centrality dependendent, with and without selection cuts.

    Acho interessante colocar os 6 graficos no mesmo canvas. PbPb2023, PbPb2024, ppRef2024.

3- Add pad in the end with A(pT,PbPb2024) - A(pT,ppRef).
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

//---Macro settings

TString plot_extension = ".pdf";

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
        "/home/lucasdoriac/z_boson_analysis/Data/PbPb2023/"
    },

    {
        "PbPb2024_Data",
        SampleType::Data,
        CollisionSystem::PbPb2024,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_PbPb2024Data.root",
        "/home/lucasdoriac/z_boson_analysis/Data/PbPb2024/"
    },

    {
        "ppRef2024_Data",
        SampleType::Data,
        CollisionSystem::ppRef2024,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_ppRef2024.root",
        "/home/lucasdoriac/z_boson_analysis/Data/ppRef2024/"
    },

    {
        "PbPb2024_MC",
        SampleType::MC,
        CollisionSystem::PbPb2024,
        "hionia/myTree",
        "Oniatree_PowhegZtoMuMu_PbPb2024_*.root",
        "/home/lucasdoriac/z_boson_analysis/MC/PbPb2024/DYto2Mu_MLL-50_TuneCP5_5p36TeV_powheg-pythia8/PowhegEmbedded_March9/260309_143939/0000/"
    },

    {
        "ppRef2024_MC",
        SampleType::MC,
        CollisionSystem::ppRef2024,
        "hionia/myTree",
        "Oniatree_PowhegZtoMuMu_ppRef2024_*.root",
        "/home/lucasdoriac/z_boson_analysis/MC/ppRef2024/DYToMuMu_M-50_TuneCP5_5p36TeV_powheg-pythia8/Powheg_ppRefPileup_March20/260320_125046/0000/"
    }
};

//---Data storage

// --- PbPb2023 asymmetry
float Z_yield_PbPb2023 = 0.;
float Z_yieldErr_PbPb2023 = 0.; 
std::vector<double> mean_ApT_PbPb2023;
std::vector<double> mean_ApTErr_PbPb2023;

// --- Vectors to store PbPb2024 asymmetry
float Z_yield_PbPb2024 = 0.;
float Z_yieldErr_PbPb2024 = 0.; 
std::vector<double> mean_ApT_PbPb2024;
std::vector<double> mean_ApTErr_PbPb2024;

// --- Vectors to store ppRef asymmetry
float Z_yield_ppRef2024 = 0.;
float Z_yieldErr_ppRef2024 = 0.; 
std::vector<double> mean_ApT_ppRef;
std::vector<double> mean_ApTErr_ppRef;

//---Observables

//---Histogram formatting
void basicCanvasFormatting(TCanvas *c);
void basicHistFormatting(TH1D *hist);
void basicLegendFormatting(TLegend *leg);
void drawLatexText(TString latexText = "#bf{CMS}", double x = 0.15, double y = 0.93, double TextSize = 0.05);

//---Miscellaneous
void printTreeContents(const Dataset& dataset);

//---Still writing !!DONT USE!!
void pt_asymmetry(const Dataset& dataset, float minCentrality, float maxCentrality, bool goodSelection); //Stores the mean ApT and ApTErr in respective data vectors.
void make_plot();
//dataset in the respective centrality interval.
//The datasets will be PbPb2023, PbPb2024, ppRef2024.
//The intervals will be 0-10, 10-30, 30-80.
//With and without filters.

/*
    Dataset:
    0 → PbPb2023
    1 → PbPb2024
    2 → ppRef2024

    Centrality bin:
    0 → 0–20
    1 → 20–60
    2 → 60–160

    Selection:
    0 → no selection
    1 → good selection

index 0 -> 0-20, no selection
index 1 -> 20-60, no selection
index 2 -> 60-160, no selection
index 3 -> 0-20, good selection
index 4 -> 20-60, good selection
index 5 -> 60-160, good selection
*/

//---main()
void asymmetry_Evaluation_ppRef_vs_PbPb_v1(){

	gROOT->SetBatch(kTRUE);

    for(int i = 0; i < 3; ++i){ //Each dataset TGraph will have 6 points, 3 with goodSelection and 3 with noSelection.

        pt_asymmetry(datasets[i], 0., 20., false);
        pt_asymmetry(datasets[i], 20., 60., false);
        pt_asymmetry(datasets[i], 60., 160., false);

        pt_asymmetry(datasets[i], 0., 20., true);
        pt_asymmetry(datasets[i], 20., 60., true);
        pt_asymmetry(datasets[i], 60., 160., true);
    }

    make_plot();
}

void make_plot()
{
    const int nPoints = 3;

    // x = centrality bin centers
    double x[nPoints]    = {10., 40., 110.};
    double xErr[nPoints] = {10., 20., 50.};

    //--------------------------
    // No selection
    //--------------------------

    TGraphErrors *gPbPb2023_noSel =
        new TGraphErrors(nPoints);

    TGraphErrors *gPbPb2024_noSel =
        new TGraphErrors(nPoints);

    TGraphErrors *gppRef_noSel =
        new TGraphErrors(nPoints);

    for (int i = 0; i < nPoints; ++i) {

        gPbPb2023_noSel->SetPoint(i,
                                  x[i],
                                  mean_ApT_PbPb2023[i]);

        gPbPb2023_noSel->SetPointError(i,
                                       xErr[i],
                                       mean_ApTErr_PbPb2023[i]);

        gPbPb2024_noSel->SetPoint(i,
                                  x[i],
                                  mean_ApT_PbPb2024[i]);

        gPbPb2024_noSel->SetPointError(i,
                                       xErr[i],
                                       mean_ApTErr_PbPb2024[i]);

        gppRef_noSel->SetPoint(i,
                               x[i],
                               mean_ApT_ppRef[i]);

        gppRef_noSel->SetPointError(i,
                                    xErr[i],
                                    mean_ApTErr_ppRef[i]);
    }

    //--------------------------
    // Good selection
    //--------------------------

    TGraphErrors *gPbPb2023_sel =
        new TGraphErrors(nPoints);

    TGraphErrors *gPbPb2024_sel =
        new TGraphErrors(nPoints);

    TGraphErrors *gppRef_sel =
        new TGraphErrors(nPoints);

    for (int i = 0; i < nPoints; ++i) {

        int index = i + 3;

        gPbPb2023_sel->SetPoint(i,
                                x[i],
                                mean_ApT_PbPb2023[index]);

        gPbPb2023_sel->SetPointError(i,
                                     xErr[i],
                                     mean_ApTErr_PbPb2023[index]);

        gPbPb2024_sel->SetPoint(i,
                                x[i],
                                mean_ApT_PbPb2024[index]);

        gPbPb2024_sel->SetPointError(i,
                                     xErr[i],
                                     mean_ApTErr_PbPb2024[index]);

        gppRef_sel->SetPoint(i,
                             x[i],
                             mean_ApT_ppRef[index]);

        gppRef_sel->SetPointError(i,
                                  xErr[i],
                                  mean_ApTErr_ppRef[index]);
    }

    gPbPb2023_noSel->SetMarkerColor(kBlue);
    gPbPb2023_sel  ->SetMarkerColor(kBlue);

    gPbPb2024_noSel->SetMarkerColor(kRed);
    gPbPb2024_sel  ->SetMarkerColor(kRed);

    gppRef_noSel->SetMarkerColor(kBlack);
    gppRef_sel  ->SetMarkerColor(kBlack);

    gPbPb2023_noSel->SetMarkerStyle(20);
    gPbPb2024_noSel->SetMarkerStyle(20);
    gppRef_noSel   ->SetMarkerStyle(20);

    gPbPb2023_sel->SetMarkerStyle(21);
    gPbPb2024_sel->SetMarkerStyle(21);
    gppRef_sel   ->SetMarkerStyle(21);

    TCanvas *c = new TCanvas("c",
                         "p_{T} asymmetry",
                         900,
                         700);

    gPbPb2023_noSel->SetTitle(";Centrality (%);<#Delta p_{T}/(p_{T}^{-}+p_{T}^{+})>");

    gPbPb2023_noSel->Draw("AP");

    gPbPb2024_noSel->Draw("P SAME");
    gppRef_noSel->Draw("P SAME");

    gPbPb2023_sel->Draw("P SAME");
    gPbPb2024_sel->Draw("P SAME");
    gppRef_sel->Draw("P SAME");

    TLegend *leg = new TLegend(0.58,0.60,0.88,0.88);

    leg->AddEntry(gPbPb2023_noSel,"PbPb2023","p");
    leg->AddEntry(gPbPb2024_noSel,"PbPb2024","p");
    leg->AddEntry(gppRef_noSel,"ppRef2024","p");

    leg->AddEntry((TObject*)0,"","");

    leg->AddEntry(gPbPb2023_noSel,"#circ No selection","p");
    leg->AddEntry(gPbPb2023_sel,"#square Good selection","p");

    leg->Draw();

    c->SaveAs("ptAsymmetry_vs_Centrality.pdf");
}

void pt_asymmetry(const Dataset& dataset, float minCentrality, float maxCentrality, bool goodSelection){

    // Load root file.
    std::string fullPath = dataset.basePath + dataset.filePattern;

    TChain *chain = new TChain(dataset.treeName.c_str());
    chain->Add(fullPath.c_str());

    std::cout << "> Number of files = " << chain->GetListOfFiles()->GetEntries() << "\n" << std::endl;

    std::cout << "> Opening files " << fullPath << "\n" << std::endl;

    std::cout << "> Running function " << __func__ << " on " << dataset.name << "\n" << std::endl;
    std::cout << "> On interval " << minCentrality << " to " << maxCentrality << "\n" << std::endl;
    
    //Total number of events on Tree.
    Long64_t nEvents = chain->GetEntries();

    //Event-level variables
    UInt_t eventNb;
    UInt_t runNb;
    UInt_t LS;

    Float_t zVtx;
    Short_t nPV;
    Int_t Centrality;
    Short_t NpixelTracks;

    Int_t trigPrescale[8];
    ULong64_t HLTriggers;

    Float_t SumET_HF;
    Float_t SumET_HFplus;
    Float_t SumET_HFminus;
    Float_t SumET_HFplusEta4;
    Float_t SumET_HFminusEta4;
    Float_t SumET_ET;
    Float_t SumET_EE;
    Float_t SumET_EB;

    Int_t nEP;

    Float_t rpAng[100];
    Float_t rpSin[100];
    Float_t rpCos[100];

    Float_t rpAng_origin[100];
    Float_t rpSin_origin[100];
    Float_t rpCos_origin[100];

    chain->SetBranchAddress("eventNb", &eventNb);
    chain->SetBranchAddress("runNb", &runNb);
    chain->SetBranchAddress("LS", &LS);

    chain->SetBranchAddress("zVtx", &zVtx);
    chain->SetBranchAddress("nPV", &nPV);
    chain->SetBranchAddress("Centrality", &Centrality);
    chain->SetBranchAddress("NpixelTracks", &NpixelTracks);

    chain->SetBranchAddress("trigPrescale", trigPrescale);
    chain->SetBranchAddress("HLTriggers", &HLTriggers);

    chain->SetBranchAddress("SumET_HF", &SumET_HF);
    chain->SetBranchAddress("SumET_HFplus", &SumET_HFplus);
    chain->SetBranchAddress("SumET_HFminus", &SumET_HFminus);
    chain->SetBranchAddress("SumET_HFplusEta4", &SumET_HFplusEta4);
    chain->SetBranchAddress("SumET_HFminusEta4", &SumET_HFminusEta4);
    chain->SetBranchAddress("SumET_ET", &SumET_ET);
    chain->SetBranchAddress("SumET_EE", &SumET_EE);
    chain->SetBranchAddress("SumET_EB", &SumET_EB);

    chain->SetBranchAddress("nEP", &nEP);

    chain->SetBranchAddress("rpAng", rpAng);
    chain->SetBranchAddress("rpSin", rpSin);
    chain->SetBranchAddress("rpCos", rpCos);

    chain->SetBranchAddress("rpAng_origin", rpAng_origin);
    chain->SetBranchAddress("rpSin_origin", rpSin_origin);
    chain->SetBranchAddress("rpCos_origin", rpCos_origin);
    
    //Reconstructed dimuons
    Short_t Reco_Dimuon_size;

    Short_t Reco_Dimuon_sign[1000];
    Short_t Reco_Dimuon_muonPlusIndex[1000];
    Short_t Reco_Dimuon_muonMinusIndex[1000];

    ULong64_t Reco_Dimuon_trig[1000];

    Float_t Reco_Dimuon_ctau[1000];
    Float_t Reco_Dimuon_ctauErr[1000];
    Float_t Reco_Dimuon_cosAlpha[1000];

    Float_t Reco_Dimuon_ctau3D[1000];
    Float_t Reco_Dimuon_ctauErr3D[1000];
    Float_t Reco_Dimuon_cosAlpha3D[1000];

    Float_t Reco_Dimuon_vtxProb[1000];
    Float_t Reco_Dimuon_dca[1000];
    Float_t Reco_Dimuon_MassErr[1000];

    std::vector<float>* Reco_Dimuon_pt = nullptr;
    std::vector<float>* Reco_Dimuon_eta = nullptr;
    std::vector<float>* Reco_Dimuon_rapidity = nullptr;
    std::vector<float>* Reco_Dimuon_phi = nullptr;
    std::vector<float>* Reco_Dimuon_invMass = nullptr;

    std::vector<float>* Reco_Dimuon_muonPtDiff = nullptr;
    std::vector<float>* Reco_Dimuon_muonPtRelDiff = nullptr;

    std::vector<float>* Reco_Dimuon_vtx_xpos = nullptr;
    std::vector<float>* Reco_Dimuon_vtx_ypos = nullptr;
    std::vector<float>* Reco_Dimuon_vtx_zpos = nullptr;

    chain->SetBranchAddress("Reco_Dimuon_size", &Reco_Dimuon_size);

    chain->SetBranchAddress("Reco_Dimuon_sign", Reco_Dimuon_sign);

    chain->SetBranchAddress("Reco_Dimuon_pt", &Reco_Dimuon_pt);
    chain->SetBranchAddress("Reco_Dimuon_eta", &Reco_Dimuon_eta);
    chain->SetBranchAddress("Reco_Dimuon_rapidity", &Reco_Dimuon_rapidity);
    chain->SetBranchAddress("Reco_Dimuon_phi", &Reco_Dimuon_phi);
    chain->SetBranchAddress("Reco_Dimuon_invMass", &Reco_Dimuon_invMass);

    chain->SetBranchAddress("Reco_Dimuon_muonPtDiff", &Reco_Dimuon_muonPtDiff);
    chain->SetBranchAddress("Reco_Dimuon_muonPtRelDiff", &Reco_Dimuon_muonPtRelDiff);

    chain->SetBranchAddress("Reco_Dimuon_muonPlusIndex", Reco_Dimuon_muonPlusIndex);
    chain->SetBranchAddress("Reco_Dimuon_muonMinusIndex", Reco_Dimuon_muonMinusIndex);

    chain->SetBranchAddress("Reco_Dimuon_trig", Reco_Dimuon_trig);

    chain->SetBranchAddress("Reco_Dimuon_ctau", Reco_Dimuon_ctau);
    chain->SetBranchAddress("Reco_Dimuon_ctauErr", Reco_Dimuon_ctauErr);
    chain->SetBranchAddress("Reco_Dimuon_cosAlpha", Reco_Dimuon_cosAlpha);

    chain->SetBranchAddress("Reco_Dimuon_ctau3D", Reco_Dimuon_ctau3D);
    chain->SetBranchAddress("Reco_Dimuon_ctauErr3D", Reco_Dimuon_ctauErr3D);
    chain->SetBranchAddress("Reco_Dimuon_cosAlpha3D", Reco_Dimuon_cosAlpha3D);

    chain->SetBranchAddress("Reco_Dimuon_vtxProb", Reco_Dimuon_vtxProb);
    chain->SetBranchAddress("Reco_Dimuon_dca", Reco_Dimuon_dca);
    //chain->SetBranchAddress("Reco_Dimuon_MassErr", Reco_Dimuon_MassErr);

    chain->SetBranchAddress("Reco_Dimuon_vtx_xpos", &Reco_Dimuon_vtx_xpos);
    chain->SetBranchAddress("Reco_Dimuon_vtx_ypos", &Reco_Dimuon_vtx_ypos);
    chain->SetBranchAddress("Reco_Dimuon_vtx_zpos", &Reco_Dimuon_vtx_zpos);

    //Reconstructed muons
    Short_t Reco_Muon_size;

    std::vector<float>* Reco_Muon_pt = nullptr;
    std::vector<float>* Reco_Muon_ptErrTrk = nullptr;
    std::vector<float>* Reco_Muon_eta = nullptr;
    std::vector<float>* Reco_Muon_phi = nullptr;
    std::vector<float>* Reco_Muon_mass = nullptr;

    std::vector<float>* Reco_Muon_etaL1 = nullptr;
    std::vector<float>* Reco_Muon_phiL1 = nullptr;

    ULong64_t Reco_Muon_trig[1000];

    Bool_t Reco_Muon_isPF[1000];
    Bool_t Reco_Muon_isTracker[1000];
    Bool_t Reco_Muon_isGlobal[1000];

    Bool_t Reco_Muon_isSoftCutBased[1000];
    Bool_t Reco_Muon_isHybridSoft[1000];
    Bool_t Reco_Muon_isLooseCutBased[1000];
    Bool_t Reco_Muon_isMediumCutBased[1000];
    Bool_t Reco_Muon_isTightCutBased[1000];

    Float_t Reco_Muon_softMVAValue[1000];
    Float_t Reco_Muon_muonMVAValue[1000];

    std::vector<bool>* Reco_Muon_passesPFIsoLoose = nullptr;
    std::vector<bool>* Reco_Muon_passesPFIsoMedium = nullptr;
    std::vector<bool>* Reco_Muon_passesPFIsoTight = nullptr;
    std::vector<bool>* Reco_Muon_passesPFIsoVeryTight = nullptr;

    std::vector<float>* Reco_Muon_isoTrackSumPt = nullptr;

    std::vector<bool>* Reco_Muon_passesMultiIsoMedium = nullptr;

    std::vector<float>* Reco_Muon_HIMVAIso = nullptr;

    std::vector<bool>* Reco_Muon_HIMVAIsoWP80 = nullptr;
    std::vector<bool>* Reco_Muon_HIMVAIsoWP85 = nullptr;
    std::vector<bool>* Reco_Muon_HIMVAIsoWP90 = nullptr;
    std::vector<bool>* Reco_Muon_HIMVAIsoWP95 = nullptr;

    chain->SetBranchAddress("Reco_Muon_size", &Reco_Muon_size);

    chain->SetBranchAddress("Reco_Muon_pt", &Reco_Muon_pt);
    chain->SetBranchAddress("Reco_Muon_ptErrTrk", &Reco_Muon_ptErrTrk);
    chain->SetBranchAddress("Reco_Muon_eta", &Reco_Muon_eta);
    chain->SetBranchAddress("Reco_Muon_phi", &Reco_Muon_phi);
    chain->SetBranchAddress("Reco_Muon_mass", &Reco_Muon_mass);

    chain->SetBranchAddress("Reco_Muon_etaL1", &Reco_Muon_etaL1);
    chain->SetBranchAddress("Reco_Muon_phiL1", &Reco_Muon_phiL1);

    chain->SetBranchAddress("Reco_Muon_trig", Reco_Muon_trig);

    chain->SetBranchAddress("Reco_Muon_isPF", Reco_Muon_isPF);
    chain->SetBranchAddress("Reco_Muon_isTracker", Reco_Muon_isTracker);
    chain->SetBranchAddress("Reco_Muon_isGlobal", Reco_Muon_isGlobal);

    chain->SetBranchAddress("Reco_Muon_isSoftCutBased", Reco_Muon_isSoftCutBased);
    chain->SetBranchAddress("Reco_Muon_isHybridSoft", Reco_Muon_isHybridSoft);
    chain->SetBranchAddress("Reco_Muon_isLooseCutBased", Reco_Muon_isLooseCutBased);
    chain->SetBranchAddress("Reco_Muon_isMediumCutBased", Reco_Muon_isMediumCutBased);
    chain->SetBranchAddress("Reco_Muon_isTightCutBased", Reco_Muon_isTightCutBased);

    chain->SetBranchAddress("Reco_Muon_softMVAValue", Reco_Muon_softMVAValue);
    chain->SetBranchAddress("Reco_Muon_muonMVAValue", Reco_Muon_muonMVAValue);

    chain->SetBranchAddress("Reco_Muon_passesPFIsoLoose", &Reco_Muon_passesPFIsoLoose);
    chain->SetBranchAddress("Reco_Muon_passesPFIsoMedium", &Reco_Muon_passesPFIsoMedium);
    chain->SetBranchAddress("Reco_Muon_passesPFIsoTight", &Reco_Muon_passesPFIsoTight);
    chain->SetBranchAddress("Reco_Muon_passesPFIsoVeryTight", &Reco_Muon_passesPFIsoVeryTight);

    chain->SetBranchAddress("Reco_Muon_isoTrackSumPt", &Reco_Muon_isoTrackSumPt);

    chain->SetBranchAddress("Reco_Muon_passesMultiIsoMedium", &Reco_Muon_passesMultiIsoMedium);

    chain->SetBranchAddress("Reco_Muon_HIMVAIso", &Reco_Muon_HIMVAIso);

    chain->SetBranchAddress("Reco_Muon_HIMVAIsoWP80", &Reco_Muon_HIMVAIsoWP80);
    chain->SetBranchAddress("Reco_Muon_HIMVAIsoWP85", &Reco_Muon_HIMVAIsoWP85);
    chain->SetBranchAddress("Reco_Muon_HIMVAIsoWP90", &Reco_Muon_HIMVAIsoWP90);
    chain->SetBranchAddress("Reco_Muon_HIMVAIsoWP95", &Reco_Muon_HIMVAIsoWP95);

    //Histograms
    TH1D* hist = new TH1D("hist", "", 100, -1., 1.);

    //Event-level cut values.
    const double maxZvtx = 15.0;

    //Z-level cut values.
    const double minZ_Mass = 80.0; 
    const double maxZ_Mass = 100.0; 
    const double RapidityCutValue = 2.4;

    //Muon-level cut values.
    const float EtaCutValue = 2.4;
    const double ptCutValue = 20.0;
    bool MuPlIsTight;
    bool MuMiIsTight;

    //Observables
    TLorentzVector Z;
    double ptplus, ptminus;
    double etaplus, etaminus;
    //

    for(Long64_t i = 0; i < nEvents; ++i){ //Loop through all EVENTS in the CHAIN.

        chain->GetEntry(i); //Get event i.

        //Centrality cut out of the event selection since we want a centrality dependent histogram.
        bool goodCent = (Centrality > minCentrality && Centrality < maxCentrality);

        //Good event selection. Accept all events unless goodSelection = true.
        bool goodVertex = true;
        if (goodSelection) goodVertex = (std::abs(zVtx) < maxZvtx);

        if (!goodCent) continue;
        if (!goodVertex) continue;

        for(Short_t j = 0; j < Reco_Dimuon_size; ++j){ //Loop through all reco dimuon candidates of event i.
            
            //Good Z selection. Accept all dimuons unless goodSelection = true.
            bool goodMass = true;
            bool goodRapidity = true;

            if(goodSelection){
                goodMass = (Reco_Dimuon_invMass->at(j) > minZ_Mass && Reco_Dimuon_invMass->at(j) < maxZ_Mass);
                goodRapidity = (Reco_Dimuon_rapidity->at(j) < RapidityCutValue);
            }

            if (!goodMass) continue;
            if (!goodRapidity) continue;

            //Good muon selection. Accept all muons unless goodSelection = true.
            bool goodMuPl = true;
            bool goodMuMi = true;

            if(goodSelection){
                ptplus = Reco_Muon_pt->at(Reco_Dimuon_muonPlusIndex[j]); //pT of antimuon.
                ptminus = Reco_Muon_pt->at(Reco_Dimuon_muonMinusIndex[j]); //pT of corresponding muon.
                etaplus = Reco_Muon_eta->at(Reco_Dimuon_muonPlusIndex[j]); //Pseudorapidity of antimuon.
                etaminus = Reco_Muon_eta->at(Reco_Dimuon_muonMinusIndex[j]); //Pseudorapidity of corresponding muon.
                MuPlIsTight = Reco_Muon_isTightCutBased[Reco_Dimuon_muonPlusIndex[j]];
                MuMiIsTight = Reco_Muon_isTightCutBased[Reco_Dimuon_muonMinusIndex[j]];

                goodMuPl = (ptplus > ptCutValue) 
                                && (std::abs(etaplus) < EtaCutValue) 
                                && (MuPlIsTight);

                goodMuMi = (ptminus > ptCutValue) 
                                && (std::abs(etaminus) < EtaCutValue) 
                                && (MuMiIsTight);
            }
            
            if (!goodMuPl || !goodMuMi) continue;

            hist->Fill(Reco_Dimuon_muonPtRelDiff->at(j));
        }

        //Just to check the progress.
        if (i % (nEvents / 100) == 0){

            int percent = static_cast<int>(100.0 * i / nEvents+0.5);

            std::cout << "\r"
                      << percent
                      << "% complete..."
                      << std::flush;
        }//
    
    }//Exiting event-by-event loop.

    //End of histogram filling.
    std::cout << "\r100% complete!" << std::endl;

    //Starting histogram statistics calculation.

    //Calculate Z candidates yield and error.
    //float Z_yield, Z_yield_err;
    //Z_yield = hist->IntegralAndError(1, hist->GetNbinsX(), Z_yield_err);

    //Calculate mean of distributions (should be around ~ 0 GeV)
    float mean = hist->GetMean();
    float meanErr = hist->GetMeanError();

    if(dataset.name == "PbPb2023_Data"){
        mean_ApT_PbPb2023.push_back(mean);
        mean_ApTErr_PbPb2023.push_back(meanErr);
    }
    else if(dataset.name == "PbPb2024_Data"){
        mean_ApT_PbPb2024.push_back(mean);
        mean_ApTErr_PbPb2024.push_back(meanErr);
    }
    else if(dataset.name == "ppRef2024_Data"){
        mean_ApT_ppRef.push_back(mean);
        mean_ApTErr_ppRef.push_back(meanErr);
    }

    std::cout << "Values stored... moving on \n" << std::endl;

    //Calculate standard deviation of pT distributions (should be around ?)
    /*double std_pl = hist_ptpl->GetStdDev();
    double std_mi = hist_ptmi->GetStdDev();
    double std_pl_err = hist_ptpl->GetStdDevError();
    double std_mi_err = hist_ptmi->GetStdDevError();*/

    delete hist;
    delete chain;
}