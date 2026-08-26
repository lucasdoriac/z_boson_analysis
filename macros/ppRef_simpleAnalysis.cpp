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

TString plot_extension = ".png"; // ".png" for regular development and ".pdf" for final quality plots
std::string BasePath = "/home/lucas/Documents/CMS_analyzes/Z_boson_analysis/"; //IFT
//std::string BasePath = "/home/lucasdoriac/z_boson_analysis/data/"; //Home

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


//---Observables
TH1D* GetChargeAsymmetryHist();


//---Histogram formatting
void basicCanvasFormatting(TCanvas *c);
void basicHistFormatting(TH1D *hist);
void basicLegendFormatting(TLegend *leg);
void drawLatexText(TString latexText = "#bf{CMS}", double x = 0.15, double y = 0.93, double TextSize = 0.05);

//---Miscellaneous
void printTreeContents(const Dataset &dataset);

//---Still writing !!DONT USE!!
void invMassSpectrum(const Dataset &dataset);

//---main()
void ppRef_simpleAnalysis(){

	gROOT->SetBatch(kTRUE);

    printTreeContents(datasets[2]); // 2 -> ppRef2024.
    invMassSpectrum(datasets[2]);
}

//---Histogram formatting

void basicCanvasFormatting(TCanvas* c){
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.035);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    c->SetTickx(1);
    c->SetTicky(1);
    c->SetFillColor(0);
    c->SetFrameFillColor(0);
    c->SetFrameLineWidth(2);
}

void basicHistFormatting(TH1D* hist){
	hist->GetXaxis()->CenterTitle(false);
    hist->GetYaxis()->CenterTitle(false);
    hist->GetXaxis()->SetTitleOffset(.9);
    hist->GetYaxis()->SetTitleOffset(1.);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetTitleSize(0.055);
    hist->GetYaxis()->SetTitleSize(0.055);
    hist->GetXaxis()->SetLabelSize(0.042);
    hist->GetYaxis()->SetLabelSize(0.042);
    hist->SetTitle("");
	hist->SetStats(0);
	hist->SetMarkerStyle(20);
	//hist->SetMarkerSize(0.9);
}

void basicLegendFormatting(TLegend *leg){
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.042);
    leg->SetTextFont(42);
    leg->SetMargin(0.2);
    leg->SetEntrySeparation(0.04);
}

void drawLatexText(TString latexText, double x, double y, double TextSize){
// We can add any pT or eta selection. If no text is passed to the function the CMS Header will be drawn.
TLatex latex;
latex.SetNDC();              // For normalized coordinates
latex.SetTextSize(TextSize);
latex.SetTextFont(42);       // Helvetica
latex.SetTextAlign(11);      // Left-top aligned.
latex.DrawLatex(x, y, latexText);
}

//---Miscellaneous
void printTreeContents(const Dataset &dataset){

    std::string fullPath = dataset.basePath + dataset.filePattern;

    TChain *chain = new TChain(dataset.treeName.c_str());
    chain->Add(fullPath.c_str());

    chain->Print();
    std::cout << "Number of files = " << chain->GetListOfFiles()->GetEntries() << std::endl;

    Long64_t nEvents = chain->GetEntries();
    std::cout << "Total number of events = " << nEvents << std::endl;

    std::cout << "Dataset: " << dataset.name << std::endl;
    std::cout << "Tree: " << dataset.treeName << std::endl;
    std::cout << "Path: " << fullPath << std::endl;
}

void invMassSpectrum(const Dataset& dataset){

    // Load root file.
    std::string fullPath = dataset.basePath + dataset.filePattern;

    TChain *chain = new TChain(dataset.treeName.c_str());
    chain->Add(fullPath.c_str());

    std::cout << "> Number of files = " << chain->GetListOfFiles()->GetEntries() << "\n" << std::endl;

    std::cout << "> Opening files " << fullPath << "\n" << std::endl;

    std::cout << "> Running function " << __func__ << " on " << dataset.name << "\n" << std::endl;
    
    //Total number of events on Tree.
    Long64_t nEvents = chain->GetEntries();

    const int MAX_DIMUON = 1000;
    const int MAX_MUON   = 1000;

    // ============================================================
    // Event branches
    // ============================================================

    UInt_t   eventNb;
    UInt_t   runNb;
    UInt_t   LS;
    Float_t  zVtx;
    Short_t  nPV;

    Int_t       trigPrescale[9];
    ULong64_t   HLTriggers;


    // ============================================================
    // Dimuon branches
    // ============================================================

    Short_t Reco_Dimuon_size;

    Short_t     Reco_Dimuon_sign[MAX_DIMUON];

    std::vector<float>* Reco_Dimuon_pt       = nullptr;
    std::vector<float>* Reco_Dimuon_eta      = nullptr;
    std::vector<float>* Reco_Dimuon_rapidity = nullptr;
    std::vector<float>* Reco_Dimuon_phi      = nullptr;
    std::vector<float>* Reco_Dimuon_invMass  = nullptr;

    std::vector<float>* Reco_Dimuon_muonPtDiff    = nullptr;
    std::vector<float>* Reco_Dimuon_muonPtRelDiff = nullptr;

    Short_t Reco_Dimuon_muonPlusIndex[MAX_DIMUON];
    Short_t Reco_Dimuon_muonMinusIndex[MAX_DIMUON];

    ULong64_t Reco_Dimuon_trig[MAX_DIMUON];

    Float_t Reco_Dimuon_ctau[MAX_DIMUON];
    Float_t Reco_Dimuon_ctauErr[MAX_DIMUON];
    Float_t Reco_Dimuon_cosAlpha[MAX_DIMUON];

    Float_t Reco_Dimuon_ctau3D[MAX_DIMUON];
    Float_t Reco_Dimuon_ctauErr3D[MAX_DIMUON];
    Float_t Reco_Dimuon_cosAlpha3D[MAX_DIMUON];

    Float_t Reco_Dimuon_vtxProb[MAX_DIMUON];
    Float_t Reco_Dimuon_dca[MAX_DIMUON];
    Float_t Reco_Dimuon_invMassErr[MAX_DIMUON];

    std::vector<float>* Reco_Dimuon_vtx_xpos = nullptr;
    std::vector<float>* Reco_Dimuon_vtx_ypos = nullptr;
    std::vector<float>* Reco_Dimuon_vtx_zpos = nullptr;


    // ============================================================
    // Muon branches
    // ============================================================

    Short_t Reco_Muon_size;

    std::vector<float>* Reco_Muon_pt       = nullptr;
    std::vector<float>* Reco_Muon_ptErrTrk = nullptr;
    std::vector<float>* Reco_Muon_eta      = nullptr;
    std::vector<float>* Reco_Muon_phi      = nullptr;
    std::vector<float>* Reco_Muon_mass     = nullptr;

    std::vector<float>* Reco_Muon_etaL1 = nullptr;
    std::vector<float>* Reco_Muon_phiL1 = nullptr;

    ULong64_t Reco_Muon_trig[MAX_MUON];

    Bool_t Reco_Muon_isPF[MAX_MUON];
    Bool_t Reco_Muon_isTracker[MAX_MUON];
    Bool_t Reco_Muon_isGlobal[MAX_MUON];

    Bool_t Reco_Muon_isSoftCutBased[MAX_MUON];
    Bool_t Reco_Muon_isHybridSoft[MAX_MUON];
    Bool_t Reco_Muon_isLooseCutBased[MAX_MUON];
    Bool_t Reco_Muon_isMediumCutBased[MAX_MUON];
    Bool_t Reco_Muon_isTightCutBased[MAX_MUON];

    Float_t Reco_Muon_softMVAValue[MAX_MUON];
    Float_t Reco_Muon_muonMVAValue[MAX_MUON];

    std::vector<bool>* Reco_Muon_passesPFIsoLoose     = nullptr;
    std::vector<bool>* Reco_Muon_passesPFIsoMedium    = nullptr;
    std::vector<bool>* Reco_Muon_passesPFIsoTight     = nullptr;
    std::vector<bool>* Reco_Muon_passesPFIsoVeryTight = nullptr;

    std::vector<float>* Reco_Muon_isoTrackSumPt = nullptr;

    std::vector<bool>* Reco_Muon_passesMultiIsoMedium = nullptr;

    // ============================================================
    // Event branches
    // ============================================================

    chain->SetBranchAddress("eventNb", &eventNb);
    chain->SetBranchAddress("runNb", &runNb);
    chain->SetBranchAddress("LS", &LS);
    chain->SetBranchAddress("zVtx", &zVtx);
    chain->SetBranchAddress("nPV", &nPV);

    chain->SetBranchAddress("trigPrescale", trigPrescale);
    chain->SetBranchAddress("HLTriggers", &HLTriggers);


    // ============================================================
    // Dimuon branches
    // ============================================================

    chain->SetBranchAddress("Reco_Dimuon_size", &Reco_Dimuon_size);

    chain->SetBranchAddress("Reco_Dimuon_sign",
                            Reco_Dimuon_sign);

    chain->SetBranchAddress("Reco_Dimuon_pt",
                            &Reco_Dimuon_pt);

    chain->SetBranchAddress("Reco_Dimuon_eta",
                            &Reco_Dimuon_eta);

    chain->SetBranchAddress("Reco_Dimuon_rapidity",
                            &Reco_Dimuon_rapidity);

    chain->SetBranchAddress("Reco_Dimuon_phi",
                            &Reco_Dimuon_phi);

    chain->SetBranchAddress("Reco_Dimuon_invMass",
                            &Reco_Dimuon_invMass);

    chain->SetBranchAddress("Reco_Dimuon_muonPtDiff",
                            &Reco_Dimuon_muonPtDiff);

    chain->SetBranchAddress("Reco_Dimuon_muonPtRelDiff",
                            &Reco_Dimuon_muonPtRelDiff);

    chain->SetBranchAddress("Reco_Dimuon_muonPlusIndex",
                            Reco_Dimuon_muonPlusIndex);

    chain->SetBranchAddress("Reco_Dimuon_muonMinusIndex",
                            Reco_Dimuon_muonMinusIndex);

    chain->SetBranchAddress("Reco_Dimuon_trig",
                            Reco_Dimuon_trig);

    chain->SetBranchAddress("Reco_Dimuon_ctau",
                            Reco_Dimuon_ctau);

    chain->SetBranchAddress("Reco_Dimuon_ctauErr",
                            Reco_Dimuon_ctauErr);

    chain->SetBranchAddress("Reco_Dimuon_cosAlpha",
                            Reco_Dimuon_cosAlpha);

    chain->SetBranchAddress("Reco_Dimuon_ctau3D",
                            Reco_Dimuon_ctau3D);

    chain->SetBranchAddress("Reco_Dimuon_ctauErr3D",
                            Reco_Dimuon_ctauErr3D);

    chain->SetBranchAddress("Reco_Dimuon_cosAlpha3D",
                            Reco_Dimuon_cosAlpha3D);

    chain->SetBranchAddress("Reco_Dimuon_vtxProb",
                            Reco_Dimuon_vtxProb);

    chain->SetBranchAddress("Reco_Dimuon_dca",
                            Reco_Dimuon_dca);

    // Branch name is Reco_Dimuon_MassErr in the tree
    chain->SetBranchAddress("Reco_Dimuon_MassErr",
                            Reco_Dimuon_invMassErr);

    chain->SetBranchAddress("Reco_Dimuon_vtx_xpos",
                            &Reco_Dimuon_vtx_xpos);

    chain->SetBranchAddress("Reco_Dimuon_vtx_ypos",
                            &Reco_Dimuon_vtx_ypos);

    chain->SetBranchAddress("Reco_Dimuon_vtx_zpos",
                            &Reco_Dimuon_vtx_zpos);


    // ============================================================
    // Muon branches
    // ============================================================

    chain->SetBranchAddress("Reco_Muon_size",
                            &Reco_Muon_size);

    chain->SetBranchAddress("Reco_Muon_pt",
                            &Reco_Muon_pt);

    chain->SetBranchAddress("Reco_Muon_ptErrTrk",
                            &Reco_Muon_ptErrTrk);

    chain->SetBranchAddress("Reco_Muon_eta",
                            &Reco_Muon_eta);

    chain->SetBranchAddress("Reco_Muon_phi",
                            &Reco_Muon_phi);

    chain->SetBranchAddress("Reco_Muon_mass",
                            &Reco_Muon_mass);

    chain->SetBranchAddress("Reco_Muon_etaL1",
                            &Reco_Muon_etaL1);

    chain->SetBranchAddress("Reco_Muon_phiL1",
                            &Reco_Muon_phiL1);

    chain->SetBranchAddress("Reco_Muon_trig",
                            Reco_Muon_trig);

    chain->SetBranchAddress("Reco_Muon_isPF",
                            Reco_Muon_isPF);

    chain->SetBranchAddress("Reco_Muon_isTracker",
                            Reco_Muon_isTracker);

    chain->SetBranchAddress("Reco_Muon_isGlobal",
                            Reco_Muon_isGlobal);

    chain->SetBranchAddress("Reco_Muon_isSoftCutBased",
                            Reco_Muon_isSoftCutBased);

    chain->SetBranchAddress("Reco_Muon_isHybridSoft",
                            Reco_Muon_isHybridSoft);

    chain->SetBranchAddress("Reco_Muon_isLooseCutBased",
                            Reco_Muon_isLooseCutBased);

    chain->SetBranchAddress("Reco_Muon_isMediumCutBased",
                            Reco_Muon_isMediumCutBased);

    chain->SetBranchAddress("Reco_Muon_isTightCutBased",
                            Reco_Muon_isTightCutBased);

    chain->SetBranchAddress("Reco_Muon_softMVAValue",
                            Reco_Muon_softMVAValue);

    chain->SetBranchAddress("Reco_Muon_muonMVAValue",
                            Reco_Muon_muonMVAValue);

    chain->SetBranchAddress("Reco_Muon_passesPFIsoLoose",
                            &Reco_Muon_passesPFIsoLoose);

    chain->SetBranchAddress("Reco_Muon_passesPFIsoMedium",
                            &Reco_Muon_passesPFIsoMedium);

    chain->SetBranchAddress("Reco_Muon_passesPFIsoTight",
                            &Reco_Muon_passesPFIsoTight);

    chain->SetBranchAddress("Reco_Muon_passesPFIsoVeryTight",
                            &Reco_Muon_passesPFIsoVeryTight);

    chain->SetBranchAddress("Reco_Muon_isoTrackSumPt",
                            &Reco_Muon_isoTrackSumPt);

    chain->SetBranchAddress("Reco_Muon_passesMultiIsoMedium", &Reco_Muon_passesMultiIsoMedium);

    //Histograms
    //TH1D *h_triggerBits = new TH1D("h_triggerBits", "", 260, 0, 260);
    TH1D *h_invMass = new TH1D("h_invMass", "", 40, 70, 110);


    for(Long64_t i = 0; i < nEvents; ++i){

        chain->GetEntry(i);

        for(Short_t j = 0; j < Reco_Dimuon_size; ++j){




        }
    }


}