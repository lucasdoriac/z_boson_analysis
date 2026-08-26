/*
1- First candidate plot let's do the usual N(\mu+)(pT) and N(\mu-)(pT).
This is just the distribution of dN/dpT for the individual muon daughters and it's the distribution that
the theory paper uses to predict the shift in the peak of THESE distributions in the 30-50 GeV region.
On the bottom panel we will plot R(pT) = N(\mu+)/N(\mu-). We can do this for both PbPb and ppRef.

2- We would like to also look at the double ratio of the R(pT) distributions, i.e. R(PbPb)/R(ppRef).
For this we can plot the usual N(\mu+)(pT) and N(\mu-)(pT), as we did in the first candidate plot, but now for PbPb and ppRef in the same canvas. 
This would give four distributions in the top pannel, and in the bottom panel we can plot the double ratio R(PbPb)/R(ppRef).

Of course, we could also do that in the first candidate plot, since we are calculating R(pT) = N(\mu+)/N(\mu-) for each dataset,
i.e. we are already calculating R(PbPb) and R(ppRef). However, it is interesting at a first moment to plot both datasets separately,
and then later plot them in the same canvas.

3- Another observable option is the event-by-event relative difference in pT between the two muons. This is defined as:
    ΔpT = (pT(mu+) - pT(mu-)) / (pT(mu+) + pT(mu-)).
If we call that variable ΔpT, then we can calculate ΔpT(PbPb) and ΔpT(ppRef), event-by-event.
In an ideal scenario, one might expect ΔpT(PbPb) to be modified by the presence of the B field in the QGP,
with maximum effects on the 40-70% centrality bin (and maybe in the low pT region?).

The sweet observable would be the mean of ΔpT(PbPb), <ΔpT(PbPb)> per centrality bin,
as well as their difference with respect to the ppRef dataset, i.e. <ΔpT(PbPb)> - <ΔpT(ppRef)>.

4- The last simple observable we can try to check is the charge asymmetry of the muons in each dataset. This is defined as:
    A(X) = (N(mu+) - N(mu-)) / (N(mu+) + N(mu-)).
Notice that this is a YIELD asymmetry, NOT a pT asymmetry. It is the difference in the number of muons of each charge. (Normalized by the total number of muons?)
We can calculate A(PbPb) and A(ppRef) and plot them in the same canvas.
In the bottom panel we can plot the difference in the charge asymmetry between PbPb and ppRef, i.e. ΔA = A(PbPb) - A(ppRef).
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

//---Common functions

//---First candidate functions
void FillChargeYields(const Dataset& dataset, TH1D*& histMuPlus, TH1D*& histMuMinus);
void PlotChargeYields(const Dataset& dataset, TH1D* histMuPlus, TH1D* histMuMinus);

//---Second candidate functions

//---Third candidate functions
TH1D* MakePtRelativeDiff(const Dataset& dataset);
void PlotPtRelativeDiff(TH1D* histPbPb, TH1D* histppRef);

//---Fourth candidate functions
TH1D* MakeChargeAsymmetryHist(const Dataset& dataset);
void PlotChargeAsymmetry(TH1D* histPbPb, TH1D* histppRef);

//---Histogram formatting
void basicCanvasFormatting(TCanvas *c, TPad *pad1, TPad *pad2);
void basicHistFormatting(TH1D *hist, bool isRatio = false);
void basicLegendFormatting(TLegend *leg);
void drawLatexText(TString latexText = "#bf{CMS}", double x = 0.15, double y = 0.93, double TextSize = 0.05);

//---Main function
void Penguin(){

    gROOT->SetBatch(kTRUE);

    //First candidate plot stuff
    TH1D* histYieldMuPlus = new TH1D("histYieldMuPlus", "histYieldMuPlus", 100, 0., 100.);
    TH1D* histYieldMuMinus = new TH1D("histYieldMuMinus", "histYieldMuMinus", 100, 0., 100.);

    //PbPb2024
    FillChargeYields(datasets[1], histYieldMuPlus, histYieldMuMinus);
    PlotChargeYields(datasets[1], histYieldMuPlus, histYieldMuMinus);

    //ppRef2024
    //FillChargeYields(datasets[2], histYieldMuPlus, histYieldMuMinus);
    //PlotChargeYields(datasets[2], histYieldMuPlus, histYieldMuMinus);

    /*TH1D* histPbPb = nullptr;
    TH1D* histppRef = nullptr;

    histPbPb = MakeChargeAsymmetryHist(datasets[1]); 
    histppRef = MakeChargeAsymmetryHist(datasets[2]);

    std::cout << "PbPb entries = " << histPbPb->GetEntries() << std::endl;
    std::cout << "ppRef entries = " << histppRef->GetEntries() << std::endl;

    std::cout << "PbPb integral = " << histPbPb->Integral() << std::endl;
    std::cout << "ppRef integral = " << histppRef->Integral() << std::endl;

    PlotChargeAsymmetry(histPbPb, histppRef);

    histPbPb = MakePtRelativeDiff(datasets[1]);
    histppRef = MakePtRelativeDiff(datasets[2]);
    PlotPtRelativeDiff(histPbPb, histppRef);*/
}

void PlotChargeYields(const Dataset& dataset, TH1D* histMuPlus, TH1D* histMuMinus){

    //Plots dN/dpT distribution for a specific dataset.
    //Top pannel: dN/dpT for mu+ and mu-. (Yield distribution of daughter muons in pT bins).
    //Bottom pannel: R(pT) = N(mu+)/N(mu-). (Ratio of the yield distributions of daughter muons in pT bins).
    TCanvas *c = new TCanvas("c", "c", 800, 800);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.30, 1, 1.0);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.00, 1, 0.30);

    basicCanvasFormatting(c, pad1, pad2);
    pad1->Draw();
    pad2->Draw();

    //Top pannel
    pad1->cd();
    basicHistFormatting(histMuPlus, false);
    basicHistFormatting(histMuMinus, false);

    TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->AddEntry(histMuPlus, "#mu^{+}", "l");
    leg->AddEntry(histMuMinus, "#mu^{-}", "l");
    leg->Draw();

    drawLatexText("#bf{CMS}", 0.15, 0.93, 0.05);
    drawLatexText("#it{Work in Progress}", 0.24, 0.93, 0.042);
    //For now we are only working with PbPb2024 and ppRef2024 datasets. So we can just check the dataset.system to know which one we are plotting.
    if (dataset.system == CollisionSystem::PbPb2024) drawLatexText("PbPb 2024 (5.36 TeV)", 0.65, 0.93, 0.042);
    else drawLatexText("ppRef 2024 (5.36 TeV)", 0.65, 0.93, 0.042);

    histMuPlus->SetLineColor(kRed);
    histMuMinus->SetLineColor(kBlue);

    histMuPlus->Draw("HIST");
    histMuMinus->Draw("HIST SAME");

    //Bottom pannel
    pad2->cd();
    
    TH1D* histRatio = new TH1D("histRatio", "histRatio", histMuPlus->GetNbinsX(), histMuPlus->GetXaxis()->GetXmin(), histMuPlus->GetXaxis()->GetXmax());
    histRatio->Divide(histMuPlus, histMuMinus, 1.0, 1.0, "B");
    
    basicHistFormatting(histRatio, true);
    histRatio->SetMarkerStyle(20);
    histRatio->SetMarkerSize(0.8);
    histRatio->SetMarkerColor(kBlack);
    histRatio->SetLineColor(kBlack);
    histRatio->GetYaxis()->SetTitle("R(p_{T})");
    histRatio->GetXaxis()->SetTitle("p_{T} [GeV]");

    histRatio->Draw("P");

    //Save the canvas
    TString output = dataset.name + "_ChargeYields.png";
    c->Update();
    c->SaveAs(output);

    delete c;
}

void FillChargeYields(const Dataset& dataset, TH1D*& histMuPlus, TH1D*& histMuMinus){

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

    //For now, writing ONLY branches that are relevant to the observable we want to measure.

    //Event-level variables
    Int_t Centrality;
    Float_t  zVtx;

    chain->SetBranchAddress("zVtx", &zVtx);
    if(dataset.system == CollisionSystem::PbPb2024) {
        chain->SetBranchAddress("Centrality", &Centrality);
    }

    //Dimuon-level variables
    Short_t Reco_Dimuon_size;

    Short_t Reco_Dimuon_sign[1000];
    Short_t Reco_Dimuon_muonPlusIndex[1000];
    Short_t Reco_Dimuon_muonMinusIndex[1000];

    ULong64_t Reco_Dimuon_trig[1000];

    std::vector<float>* Reco_Dimuon_pt = nullptr;
    std::vector<float>* Reco_Dimuon_eta = nullptr;
    std::vector<float>* Reco_Dimuon_rapidity = nullptr;
    std::vector<float>* Reco_Dimuon_phi = nullptr;
    std::vector<float>* Reco_Dimuon_invMass = nullptr;

    std::vector<float>* Reco_Dimuon_muonPtDiff = nullptr;
    std::vector<float>* Reco_Dimuon_muonPtRelDiff = nullptr;

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

    //Muon-level variables
    Short_t Reco_Muon_size;

    std::vector<float>* Reco_Muon_pt = nullptr;
    //std::vector<float>* Reco_Muon_ptErrTrk = nullptr;
    std::vector<float>* Reco_Muon_eta = nullptr;
    std::vector<float>* Reco_Muon_phi = nullptr;
    std::vector<float>* Reco_Muon_mass = nullptr;

    ULong64_t Reco_Muon_trig[1000];
    Bool_t Reco_Muon_isTightCutBased[1000];
    
    chain->SetBranchAddress("Reco_Muon_size", &Reco_Muon_size);

    chain->SetBranchAddress("Reco_Muon_pt", &Reco_Muon_pt);
    //chain->SetBranchAddress("Reco_Muon_ptErrTrk", &Reco_Muon_ptErrTrk);
    chain->SetBranchAddress("Reco_Muon_eta", &Reco_Muon_eta);
    chain->SetBranchAddress("Reco_Muon_phi", &Reco_Muon_phi);
    chain->SetBranchAddress("Reco_Muon_mass", &Reco_Muon_mass);

    chain->SetBranchAddress("Reco_Muon_trig", Reco_Muon_trig);
    chain->SetBranchAddress("Reco_Muon_isTightCutBased", Reco_Muon_isTightCutBased);
    //

    //Event-level cut values.
    const int minCentrality = 0.0;
    const int maxCentrality = 180; //Cent max = maxCentrality/2.
    const double maxZvtx = 15.0;

    //Z-level cut values.
    const double minZ_Mass = 80.0; 
    const double maxZ_Mass = 100.0; 
    const double RapidityCutValue = 2.4;

    //Muon-level cut values.
    const float EtaCutValue = 2.4;
    const double ptCutValue = 10.0;
    bool MuPlIsTight;
    bool MuMiIsTight;

    //Observables
    TLorentzVector Z;
    double ptplus, ptminus;
    double etaplus, etaminus;

    //New trigger selection 'L2SingleMu12'.
    ULong64_t triggerBit = 1ULL << 7;

    //Aug.26- For now the trigger matching selection is only applied to PbPb2024 dataset, since there's an issue with ppRef2024 trigger.

    for(Long64_t i = 0; i < nEvents; ++i){//Loop through all EVENTS in the CHAIN.

        chain->GetEntry(i); //Get event i.

        //Good event selection
        bool goodVertex = (std::abs(zVtx) < maxZvtx);
        bool goodCent = true;

        if (dataset.system == CollisionSystem::PbPb2024) {
            goodCent = (Centrality > minCentrality && Centrality < maxCentrality);
        }

        if (!goodVertex) continue;
        if (!goodCent) continue;

        for(Short_t j = 0; j < Reco_Dimuon_size; ++j){ //Loop through all reco dimuon candidates of event i.
            
            //Good Z selection
            bool goodMass = (Reco_Dimuon_invMass->at(j) > minZ_Mass && Reco_Dimuon_invMass->at(j) < maxZ_Mass);
            bool goodRapidity = (std::abs(Reco_Dimuon_rapidity->at(j)) < RapidityCutValue);
            bool goodCharge = (Reco_Dimuon_sign[j] == 0);
            
            bool isTriggerMatched = true;

                if (dataset.system == CollisionSystem::PbPb2024){
                    isTriggerMatched = (Reco_Dimuon_trig[j]) & triggerBit;
                }

            if (!goodMass) continue;
            if (!goodRapidity) continue;
            if (!goodCharge) continue;
            if (!isTriggerMatched) continue;

            Short_t muonPlusIndex = Reco_Dimuon_muonPlusIndex[j];
            Short_t muonMinusIndex = Reco_Dimuon_muonMinusIndex[j];

            //Good muon selection
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
        
            //Fill histograms with yields of daughter muons in pT bins.
            histMuPlus->Fill(ptplus);
            histMuMinus->Fill(ptminus);
        }

    }//Exiting event-by-event loop.

    delete chain;
}

void PlotPtRelativeDiff(TH1D* histPbPb, TH1D* histppRef){

    TCanvas *c = new TCanvas("c", "Pt Relative Diff", 800, 600);
    //basicCanvasFormatting(c);

    basicHistFormatting(histPbPb);
    basicHistFormatting(histppRef);

    histPbPb->SetMarkerStyle(20);
    histPbPb->SetMarkerSize(0.8);
    histPbPb->SetMarkerColor(kBlue);
    histPbPb->SetLineColor(kBlue);

    histppRef->SetMarkerStyle(20);
    histppRef->SetMarkerSize(0.8);
    histppRef->SetMarkerColor(kRed);
    histppRef->SetLineColor(kRed);

    histPbPb->Draw("P");
    histppRef->Draw("P SAME");

    TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(histPbPb, "PbPb2024", "l");
    leg->AddEntry(histppRef, "ppRef2024", "l");
    leg->Draw();

    drawLatexText("#bf{CMS}", 0.15, 0.93, 0.05);
    drawLatexText("#it{Internal}", 0.24, 0.93, 0.042);
    drawLatexText("PbPb 2024 (5.36 TeV) vs ppRef 2024 (5.36 TeV)", 0.65, 0.93, 0.042);

    TString output = "PtRelDifference_PbPb2024_vs_ppRef2024.png";
    c->Update();
    c->SaveAs(output);

    delete c;
}

TH1D* MakePtRelativeDiff(const Dataset& dataset){

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

    //For now, writing ONLY branches that are relevant to the observable we want to measure.

    //Event-level variables
    Int_t Centrality;
    Float_t  zVtx;

    chain->SetBranchAddress("zVtx", &zVtx);
    if(dataset.system == CollisionSystem::PbPb2024) {
        chain->SetBranchAddress("Centrality", &Centrality);
    }

    //Dimuon-level variables
    Short_t Reco_Dimuon_size;

    Short_t Reco_Dimuon_sign[1000];
    Short_t Reco_Dimuon_muonPlusIndex[1000];
    Short_t Reco_Dimuon_muonMinusIndex[1000];

    ULong64_t Reco_Dimuon_trig[1000];

    std::vector<float>* Reco_Dimuon_pt = nullptr;
    std::vector<float>* Reco_Dimuon_eta = nullptr;
    std::vector<float>* Reco_Dimuon_rapidity = nullptr;
    std::vector<float>* Reco_Dimuon_phi = nullptr;
    std::vector<float>* Reco_Dimuon_invMass = nullptr;

    std::vector<float>* Reco_Dimuon_muonPtDiff = nullptr;
    std::vector<float>* Reco_Dimuon_muonPtRelDiff = nullptr;

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

    //Muon-level variables
    Short_t Reco_Muon_size;

    std::vector<float>* Reco_Muon_pt = nullptr;
    //std::vector<float>* Reco_Muon_ptErrTrk = nullptr;
    std::vector<float>* Reco_Muon_eta = nullptr;
    std::vector<float>* Reco_Muon_phi = nullptr;
    std::vector<float>* Reco_Muon_mass = nullptr;

    ULong64_t Reco_Muon_trig[1000];
    Bool_t Reco_Muon_isTightCutBased[1000];
    
    chain->SetBranchAddress("Reco_Muon_size", &Reco_Muon_size);

    chain->SetBranchAddress("Reco_Muon_pt", &Reco_Muon_pt);
    //chain->SetBranchAddress("Reco_Muon_ptErrTrk", &Reco_Muon_ptErrTrk);
    chain->SetBranchAddress("Reco_Muon_eta", &Reco_Muon_eta);
    chain->SetBranchAddress("Reco_Muon_phi", &Reco_Muon_phi);
    chain->SetBranchAddress("Reco_Muon_mass", &Reco_Muon_mass);

    chain->SetBranchAddress("Reco_Muon_trig", Reco_Muon_trig);
    chain->SetBranchAddress("Reco_Muon_isTightCutBased", Reco_Muon_isTightCutBased);
    //

    TH1D* hist_Dimuon_muonPtRelDiff = new TH1D("hist_Dimuon_muonPtRelDiff", "", 20, -1., 1.);

    //Event-level cut values.
    const int minCentrality = 0.0;
    const int maxCentrality = 180; //Cent max = maxCentrality/2.
    const double maxZvtx = 15.0;

    //Z-level cut values.
    const double minZ_Mass = 80.0; 
    const double maxZ_Mass = 100.0; 
    const double RapidityCutValue = 2.4;

    //Muon-level cut values.
    const float EtaCutValue = 2.4;
    const double ptCutValue = 10.0;
    bool MuPlIsTight;
    bool MuMiIsTight;

    //Observables
    TLorentzVector Z;
    double ptplus, ptminus;
    double etaplus, etaminus;

    //New trigger selection 'L2SingleMu12'
    ULong64_t triggerBit = 1ULL << 7;

    for(Long64_t i = 0; i < nEvents; ++i){//Loop through all EVENTS in the CHAIN.

        chain->GetEntry(i); //Get event i.

        //Good event selection
        bool goodVertex = (std::abs(zVtx) < maxZvtx);
        bool goodCent = true;

        if (dataset.system == CollisionSystem::PbPb2024) {
            goodCent = (Centrality > minCentrality && Centrality < maxCentrality);
        }

        if (!goodVertex) continue;
        if (!goodCent) continue;

        for(Short_t j = 0; j < Reco_Dimuon_size; ++j){ //Loop through all reco dimuon candidates of event i.
            
            //Good Z selection
            bool goodMass = (Reco_Dimuon_invMass->at(j) > minZ_Mass && Reco_Dimuon_invMass->at(j) < maxZ_Mass);
            bool goodRapidity = (std::abs(Reco_Dimuon_rapidity->at(j)) < RapidityCutValue);
            bool goodCharge = (Reco_Dimuon_sign[j] == 0);
            
            bool isTriggerMatched = true;

                if (dataset.system == CollisionSystem::PbPb2024){
                    isTriggerMatched = (Reco_Dimuon_trig[j]) & triggerBit;
                }

            if (!goodMass) continue;
            if (!goodRapidity) continue;
            if (!goodCharge) continue;
            if (!isTriggerMatched) continue;

            Short_t muonPlusIndex = Reco_Dimuon_muonPlusIndex[j];
            Short_t muonMinusIndex = Reco_Dimuon_muonMinusIndex[j];

            //Good muon selection
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
        
            //Fill histograms for muon pT distributions
            hist_Dimuon_muonPtRelDiff->Fill(Reco_Dimuon_muonPtRelDiff->at(j));
        }
    }//Exiting event-by-event loop.

    return hist_Dimuon_muonPtRelDiff;
}

TH1D* MakeChargeAsymmetryHist(const Dataset& dataset){

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

    //For now, writing ONLY branches that are relevant to the observable we want to measure.

    //Event-level variables
    Int_t Centrality;
    Float_t  zVtx;

    chain->SetBranchAddress("zVtx", &zVtx);
    if(dataset.system == CollisionSystem::PbPb2024) {
        chain->SetBranchAddress("Centrality", &Centrality);
    }

    //Dimuon-level variables
    Short_t Reco_Dimuon_size;

    Short_t Reco_Dimuon_sign[1000];
    Short_t Reco_Dimuon_muonPlusIndex[1000];
    Short_t Reco_Dimuon_muonMinusIndex[1000];

    ULong64_t Reco_Dimuon_trig[1000];

    std::vector<float>* Reco_Dimuon_pt = nullptr;
    std::vector<float>* Reco_Dimuon_eta = nullptr;
    std::vector<float>* Reco_Dimuon_rapidity = nullptr;
    std::vector<float>* Reco_Dimuon_phi = nullptr;
    std::vector<float>* Reco_Dimuon_invMass = nullptr;

    std::vector<float>* Reco_Dimuon_muonPtDiff = nullptr;
    std::vector<float>* Reco_Dimuon_muonPtRelDiff = nullptr;

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

    //Muon-level variables
    Short_t Reco_Muon_size;

    std::vector<float>* Reco_Muon_pt = nullptr;
    //std::vector<float>* Reco_Muon_ptErrTrk = nullptr;
    std::vector<float>* Reco_Muon_eta = nullptr;
    std::vector<float>* Reco_Muon_phi = nullptr;
    std::vector<float>* Reco_Muon_mass = nullptr;

    ULong64_t Reco_Muon_trig[1000];
    Bool_t Reco_Muon_isTightCutBased[1000];
    
    chain->SetBranchAddress("Reco_Muon_size", &Reco_Muon_size);

    chain->SetBranchAddress("Reco_Muon_pt", &Reco_Muon_pt);
    //chain->SetBranchAddress("Reco_Muon_ptErrTrk", &Reco_Muon_ptErrTrk);
    chain->SetBranchAddress("Reco_Muon_eta", &Reco_Muon_eta);
    chain->SetBranchAddress("Reco_Muon_phi", &Reco_Muon_phi);
    chain->SetBranchAddress("Reco_Muon_mass", &Reco_Muon_mass);

    chain->SetBranchAddress("Reco_Muon_trig", Reco_Muon_trig);
    chain->SetBranchAddress("Reco_Muon_isTightCutBased", Reco_Muon_isTightCutBased);
    //

    TH1D* histMuPl = new TH1D("histMuPl", "", 100, 0., 100.);
    TH1D* histMuMi = new TH1D("histMuMi", "", 100, 0., 100.);

    //Event-level cut values.
    const int minCentrality = 0.0;
    const int maxCentrality = 180; //Cent max = maxCentrality/2.
    const double maxZvtx = 15.0;

    //Z-level cut values.
    const double minZ_Mass = 80.0; 
    const double maxZ_Mass = 100.0; 
    const double RapidityCutValue = 2.4;

    //Muon-level cut values.
    const float EtaCutValue = 2.4;
    const double ptCutValue = 10.0;
    bool MuPlIsTight;
    bool MuMiIsTight;

    //Observables
    TLorentzVector Z;
    double ptplus, ptminus;
    double etaplus, etaminus;

    //New trigger selection 'L2SingleMu12'
    ULong64_t triggerBit = 1ULL << 7;

    for(Long64_t i = 0; i < nEvents; ++i){//Loop through all EVENTS in the CHAIN.

        chain->GetEntry(i); //Get event i.

        //Good event selection
        bool goodVertex = (std::abs(zVtx) < maxZvtx);
        bool goodCent = true;

        if (dataset.system == CollisionSystem::PbPb2024) {
            goodCent = (Centrality > minCentrality && Centrality < maxCentrality);
        }

        if (!goodVertex) continue;
        if (!goodCent) continue;

        for(Short_t j = 0; j < Reco_Dimuon_size; ++j){ //Loop through all reco dimuon candidates of event i.
            
            //Good Z selection
            bool goodMass = (Reco_Dimuon_invMass->at(j) > minZ_Mass && Reco_Dimuon_invMass->at(j) < maxZ_Mass);
            bool goodRapidity = (std::abs(Reco_Dimuon_rapidity->at(j)) < RapidityCutValue);
            bool goodCharge = (Reco_Dimuon_sign[j] == 0);
            
            bool isTriggerMatched = true;

                if (dataset.system == CollisionSystem::PbPb2024){
                    isTriggerMatched = (Reco_Dimuon_trig[j]) & triggerBit;
                }

            if (!goodMass) continue;
            if (!goodRapidity) continue;
            if (!goodCharge) continue;
            if (!isTriggerMatched) continue;

            Short_t muonPlusIndex = Reco_Dimuon_muonPlusIndex[j];
            Short_t muonMinusIndex = Reco_Dimuon_muonMinusIndex[j];

            //Good muon selection
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
        
            //Fill histograms for muon pT distributions
            histMuPl->Fill(ptplus);
            histMuMi->Fill(ptminus);
        }
    }//Exiting event-by-event loop.

    // Create a new histogram for charge asymmetry
    TString histName = "hist_ChargeAsymmetry_" + dataset.name;
    TH1D* hist_ChargeAsymmetry = new TH1D(histName, "Charge Asymmetry",
                                            histMuPl->GetNbinsX(), 
                                            histMuPl->GetXaxis()->GetXmin(), 
                                            histMuPl->GetXaxis()->GetXmax());

    hist_ChargeAsymmetry->SetTitle("Charge Asymmetry");
    hist_ChargeAsymmetry->GetYaxis()->SetTitle("Asymmetry (N_{+} - N_{-}) / (N_{+} + N_{-})");
    hist_ChargeAsymmetry->GetXaxis()->SetTitle(histMuPl->GetXaxis()->GetTitle());

    // Calculate the asymmetry for each bin
    for(int i = 1; i <= histMuPl->GetNbinsX(); ++i){
        
        double N_plus = histMuPl->GetBinContent(i);
        double N_minus = histMuMi->GetBinContent(i);
        double asymmetry = 0.0;
        
        if ((N_plus + N_minus) != 0) asymmetry = (N_plus - N_minus) / (N_plus + N_minus);
        
        hist_ChargeAsymmetry->SetBinContent(i, asymmetry);
    }
    
    delete histMuPl;
    delete histMuMi;
    delete chain;

    return hist_ChargeAsymmetry;
}

void PlotChargeAsymmetry(TH1D* histPbPb, TH1D* histppRef){

    //Plots charge asymmetry histograms of PbPb2024 and ppRef2024 datasets.
    //The plot will be in the same canvas, with the PbPb2024 histogram in blue and the ppRef2024 histogram in red.
    //The canvas will have a second pad with the ratio of the two histograms.

    //Still missing pads.

    TCanvas *c = new TCanvas("c", "Charge Asymmetry", 800, 600);
    //basicCanvasFormatting(c);

    basicHistFormatting(histPbPb);
    basicHistFormatting(histppRef);

    histPbPb->SetMarkerStyle(20);
    histPbPb->SetMarkerSize(0.8);
    histPbPb->SetMarkerColor(kBlue);
    histPbPb->SetLineColor(kBlue);

    histppRef->SetMarkerStyle(20);
    histppRef->SetMarkerSize(0.8);
    histppRef->SetMarkerColor(kRed);
    histppRef->SetLineColor(kRed);

    histPbPb->Draw("P");
    histppRef->Draw("P SAME");

    TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(histPbPb, "PbPb2024", "l");
    leg->AddEntry(histppRef, "ppRef2024", "l");
    leg->Draw();

    drawLatexText("#bf{CMS}", 0.15, 0.93, 0.05);
    drawLatexText("#it{Internal}", 0.24, 0.93, 0.042);
    drawLatexText("PbPb 2024 (5.36 TeV) vs ppRef 2024 (5.36 TeV)", 0.65, 0.93, 0.042);

    TString output = "ChargeAsymmetry_PbPb2024_vs_ppRef2024.png";
    c->Update();
    c->SaveAs(output);

    delete c;
}

void basicCanvasFormatting(TCanvas* c, TPad* pad1, TPad* pad2){
    // Canvas
    c->SetFillColor(0);
    c->SetFrameFillColor(0);
    c->SetTickx(1);
    c->SetTicky(1);

    // Top pad
    pad1->SetLeftMargin(0.12);
    pad1->SetRightMargin(0.035);
    pad1->SetTopMargin(0.08);
    pad1->SetBottomMargin(0.02);

    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad1->SetFillColor(0);
    pad1->SetFrameFillColor(0);
    pad1->SetFrameLineWidth(2);

    // Bottom pad
    pad2->SetLeftMargin(0.12);
    pad2->SetRightMargin(0.035);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.30);

    pad2->SetTickx(1);
    pad2->SetTicky(1);
    pad2->SetFillColor(0);
    pad2->SetFrameFillColor(0);
    pad2->SetFrameLineWidth(2);
}

void basicHistFormatting(TH1D* hist, bool isRatio = false){

    hist->SetTitle("");
    hist->SetStats(0);

    hist->GetXaxis()->CenterTitle(false);
    hist->GetYaxis()->CenterTitle(true);

    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);

    if (!isRatio) {

        // Main plot
        hist->GetXaxis()->SetTitleSize(0.055);
        hist->GetYaxis()->SetTitleSize(0.055);

        hist->GetXaxis()->SetLabelSize(0.0);  // Hide x labels
        hist->GetXaxis()->SetTitleSize(0.0);  // Hide x title

        hist->GetYaxis()->SetLabelSize(0.042);

        hist->GetYaxis()->SetTitleOffset(1.0);

    } else {

        // Ratio plot
        hist->GetXaxis()->SetTitleSize(0.12);
        hist->GetYaxis()->SetTitleSize(0.10);

        hist->GetXaxis()->SetLabelSize(0.10);
        hist->GetYaxis()->SetLabelSize(0.08);

        hist->GetXaxis()->SetTitleOffset(1.1);
        hist->GetYaxis()->SetTitleOffset(0.25);
    }

    hist->SetMarkerStyle(20);
}

void basicLegendFormatting(TLegend* leg){

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->SetTextFont(42);
    leg->SetTextSize(0.042);

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