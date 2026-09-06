/*
Calculates the observable
    ΔpT = (pT(mu+) - pT(mu-)) / (pT(mu+) + pT(mu-))
for both PbPb2024 and ppRef2024 samples, and plots the difference 
    <ΔpT(PbPb)> - <ΔpT(ppRef)>
as a function of centrality.
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

std::string plot_extension = ".pdf"; // ".png" for regular development and ".pdf" for final quality plots
std::string BasePath = "/home/lucas/Documents/CMS_analyzes/Z_boson_analysis/"; //IFT
//std::string BasePath = "/home/lucasdoriac/z_boson_analysis/data/"; //Home

const double MAX_ZVTX = 15.0;

const double MINZ_MASS = 80.;
const double MAXZ_MASS = 100.;
const double RAPIDITYCUTVALUE = 2.4;

const float ETACUTVALUE = 2.4;
const double PTCUTVALUE = 20.;

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

//Centrality bins
std::vector<std::pair<int,int>> centralityBins = {
    {0, 10},
    {10, 20},
    {20, 30},
    {30, 100},
};

/*std::vector<std::pair<int, int>> centralityBins = {
    {0, 10},
    {10, 30},
    {30, 50},
    {50, 100},
};*/

//Reference values.
double ppRef_DeltaPtMean;
double ppRef_DeltaPtMeanError;

//Vectors to store PbPb2024 values for each centrality bin.
std::vector<double> PbPb_DeltaPtMeanVector;
std::vector<double> PbPb_DeltaPtMeanErrorVector;

//Vectors to store final observable for each centrality bin.
std::vector<double> DeltaPtMeanVector;
std::vector<double> DeltaPtMeanErrorVector;

std::pair<double, double> CalculatePtRelDiff(const Dataset& dataset, std::pair<int, int> centralityBins);
void MakePtRelDiffPlot();
void drawLatexText(TString latexText = "#bf{CMS}", double x = 0.15, double y = 0.93, double TextSize = 0.05);


void PtRelDiff(){

    //gROOT->SetBatch(kTRUE); //Run in batch mode, without opening canvases.

    //Get value from ppRef first
    std::tie(ppRef_DeltaPtMean, ppRef_DeltaPtMeanError) = CalculatePtRelDiff(datasets[2], {0,0}); //ppRef2024

    for(int cBin = 0; cBin < centralityBins.size(); ++cBin){

        std::cout << "Centrality bin: " << centralityBins[cBin].first << "-" << centralityBins[cBin].second << std::endl;

        auto values = CalculatePtRelDiff(datasets[1], centralityBins[cBin]); //PbPb2024

        double PbPb_DeltaPtMean = values.first;
        double PbPb_DeltaPtMeanError = values.second;

        PbPb_DeltaPtMeanVector.push_back(PbPb_DeltaPtMean);
        PbPb_DeltaPtMeanErrorVector.push_back(PbPb_DeltaPtMeanError);

        double DeltaPtMeanDiff = PbPb_DeltaPtMean - ppRef_DeltaPtMean;
        double DeltaPtMeanDiffError = std::sqrt(std::pow(PbPb_DeltaPtMeanError, 2) + std::pow(ppRef_DeltaPtMeanError, 2));

        cout << "DeltaPtMeanDiff = " << DeltaPtMeanDiff << " +/- " << DeltaPtMeanDiffError << endl;

        DeltaPtMeanVector.push_back(DeltaPtMeanDiff);
        DeltaPtMeanErrorVector.push_back(DeltaPtMeanDiffError);
    }

    MakePtRelDiffPlot();


    // ============================================================
    // Print results
    // ============================================================

    std::cout << "\n========================================" << std::endl;
    std::cout << "            PtRelDiff results" << std::endl;
    std::cout << "========================================" << std::endl;

    std::cout << "ppRef2024:" << std::endl;
    std::cout << "  <Delta pT> = "
              << ppRef_DeltaPtMean
              << " +/- "
              << ppRef_DeltaPtMeanError
              << std::endl;

    std::cout << "\nPbPb2024:" << std::endl;

    for(int cBin = 0; cBin < centralityBins.size(); ++cBin){

        std::cout
            << "  Centrality "
            << centralityBins[cBin].first
            << "-"
            << centralityBins[cBin].second
            << "%:"
            << std::endl;

        std::cout
            << "    PbPb <Delta pT> = "
            << PbPb_DeltaPtMeanVector[cBin]
            << " +/- "
            << PbPb_DeltaPtMeanErrorVector[cBin]
            << std::endl;

        std::cout
            << "    PbPb - ppRef   = "
            << DeltaPtMeanVector[cBin]
            << " +/- "
            << DeltaPtMeanErrorVector[cBin]
            << std::endl;
    }

    std::cout << "========================================\n"
              << std::endl;

}

void MakePtRelDiffPlot(){

    // Number of centrality bins
    const int nBins = centralityBins.size();

    // Vectors containing the x-coordinates and their errors
    std::vector<double> x(nBins);
    std::vector<double> xError(nBins);

    // Fill x coordinates
    for(int i = 0; i < nBins; ++i){

        // Center of the centrality bin
        x[i] = (centralityBins[i].first + centralityBins[i].second) / 2.0;

        // Half-width of the centrality bin
        xError[i] = (centralityBins[i].second - centralityBins[i].first) / 2.0;
    }

    TGraphErrors* graph = new TGraphErrors(nBins, x.data(), DeltaPtMeanVector.data(), xError.data(), DeltaPtMeanErrorVector.data());

    TGraphErrors* graphPbPb = new TGraphErrors(
        nBins,
        x.data(),
        PbPb_DeltaPtMeanVector.data(),
        xError.data(),
        PbPb_DeltaPtMeanErrorVector.data()
    );

    TCanvas* canvas = new TCanvas("canvas_PtRelDiff", "pT relative difference vs centrality", 800, 700);
    canvas->SetLeftMargin(0.14);
    canvas->SetRightMargin(0.035);
    canvas->SetBottomMargin(0.12);
    canvas->SetTopMargin(0.08);
    canvas->SetTickx(1);
    canvas->SetTicky(1);
    canvas->SetFillColor(0);
    canvas->SetFrameFillColor(0);

    graph->GetXaxis()->CenterTitle(false);
    graph->GetYaxis()->CenterTitle(true);
    graph->GetXaxis()->SetTitleOffset(1.);
    graph->GetYaxis()->SetTitleOffset(1.3);
    graph->GetXaxis()->SetTitleFont(42);
    graph->GetYaxis()->SetTitleFont(42);
    graph->GetXaxis()->SetLabelFont(42);
    graph->GetYaxis()->SetLabelFont(42);
    graph->GetXaxis()->SetTitleSize(0.04);
    graph->GetYaxis()->SetTitleSize(0.04);
    graph->GetXaxis()->SetLabelSize(0.035);
    graph->GetYaxis()->SetLabelSize(0.035);

    graph->SetTitle("");
    graph->SetMarkerStyle(25);
    graph->SetMarkerSize(1.1);
    graph->SetMarkerColor(kRed+1);
    graph->SetLineColor(kRed+1);
    graph->SetLineWidth(2);

    graph->GetXaxis()->SetTitle("Centrality (%)");
    graph->GetYaxis()->SetTitle("#LT #Delta p_{T} #GT");

    graph->GetXaxis()->SetRangeUser(0, 100);
    graph->Draw("AP");

    //PbPb2024 mean
    graphPbPb->SetMarkerStyle(20);
    graphPbPb->SetMarkerSize(1.1);
    graphPbPb->SetMarkerColor(kBlack);
    graphPbPb->SetLineColor(kBlack);
    graphPbPb->SetLineWidth(2);

    graphPbPb->Draw("P SAME");

    //ppRef horizontal line
    double xMin = 0.0;
    double xMax = 100.0;

    TLine* ppRefLine = new TLine(
        xMin,
        ppRef_DeltaPtMean,
        xMax,
        ppRef_DeltaPtMean
    );

    ppRefLine->SetLineColor(kBlack);
    ppRefLine->SetLineStyle(2);
    ppRefLine->SetLineWidth(2);

    ppRefLine->Draw("SAME");

    TLegend* legend = new TLegend(0.4, 0.75, 0.6, 0.85);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.025);
    legend->AddEntry(graph, "#LT #Delta p_{T} #GT", "p");
    legend->AddEntry(graphPbPb, "#LT #Delta p_{T} #GT (PbPb2024)", "p");
    legend->AddEntry(ppRefLine, "#LT #Delta p_{T} #GT (ppRef2024)", "l");
    legend->Draw();

    drawLatexText("#bf{CMS}", 0.15, 0.93, 0.038);
    drawLatexText("#it{Work in Progress}", 0.23, 0.93, 0.03);
    drawLatexText("PbPb 2024, ppRef 2024 (5.36 TeV)", 0.55, 0.93, 0.03);

    canvas->Update();
    TString output = "DeltaPtRelDiff_PbPb2024_vs_ppRef2024" + plot_extension;
    canvas->SaveAs(output);

}

std::pair<double, double> CalculatePtRelDiff(const Dataset& dataset, std::pair<int, int> centralityBins){

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

    TString CentString = std::to_string(centralityBins.first) + "_" + std::to_string(centralityBins.second);
    TString histName = "hist_muonPtRelDiff_" + TString(dataset.name.c_str()) + "_" + CentString;
    TH1D* hist_muonPtRelDiff = new TH1D(histName, "", 200, -1., 1.);

    //Event-level cut values.
    const int minCentrality = 2*centralityBins.first; //Cent min = minCentrality/2.
    const int maxCentrality = 2*centralityBins.second; //Cent max = maxCentrality/2.
    const double maxZvtx = MAX_ZVTX;

    //Z-level cut values.
    const double minZ_Mass = MINZ_MASS; 
    const double maxZ_Mass = MAXZ_MASS; 
    const double RapidityCutValue = RAPIDITYCUTVALUE;

    //Muon-level cut values.
    const float EtaCutValue = ETACUTVALUE;
    const double ptCutValue = PTCUTVALUE;
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
            goodCent = (Centrality >= minCentrality && Centrality < maxCentrality);
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
            hist_muonPtRelDiff->Fill(Reco_Dimuon_muonPtRelDiff->at(j));
        }

    }//Exiting event-by-event loop.


    //Get mean and mean error of the histogram
    double mean = hist_muonPtRelDiff->GetMean();
    double meanError = hist_muonPtRelDiff->GetMeanError();

    return std::make_pair(mean, meanError);
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