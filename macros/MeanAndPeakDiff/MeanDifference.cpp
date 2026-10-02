/*
Mean and peak difference as function of centrality bin.
Needs to be directly from the TREE because any projected histogram loses their statistics from the filling time.

**Important**: If we want the skewness of this distribution as well we need to modify this code
to calculate the skewness directly from the Tree.
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
#include "../headers/basicFormatting.h"


//---Macro settings
std::string plot_extension = ".pdf"; // ".png" for regular development and ".pdf" for final quality plots
std::string BasePath = "/home/lucas/Documents/CMS/z_boson_analysis/"; //IFT
//std::string BasePath = "/home/lucasdoriac/z_boson_analysis/data/"; //Home

std::string whichDataset = "PbPb_Run3_Data"; // "PbPb_Run3_Data", "PbPb2023_Data", "PbPb2024_Data".
std::string JointPbPb = "PbPb2023-2026"; //"PbPb2023+2024", "PbPb2023", "PbPb2024".
std::string dataSamplesUsed = "PbPb 2023-2026, ppRef 2024 (5.36 TeV)";


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


//Set of centrality bins for PbPb2024 data. We can decide to change the centrality bins later if we want to.
std::vector<std::pair<double, double>> CentralitySet = {
    {0., 10.},
    {10., 20.},
    {20., 30.},
    {30., 100.},
    {0., 100.}
};


//Second proposed set of centrality bins for PbPb2024 data.
/*std::vector<std::pair<double, double>> CentralitySet = {
    {0., 10.},
    {10., 30.},
    {30., 50.},
    {50., 100.},
    {0., 100.}
};*/


//Vectors to store mean and peak
std::vector<std::pair<double, double>> PbPbMean;
std::vector<std::pair<double, double>> ppRefMean;

std::vector<std::pair<double, double>> PbPbPeak;
std::vector<std::pair<double, double>> ppRefPeak;

//Vector to save data for TGraphErrors at the end.
std::vector<std::pair<double, double>> MeanDiffAndError_PbPb_vs_ppRef;

//Vectors to save peak difference and error.
std::vector<std::pair<double, double>> PeakDiffAndError_PbPb_vs_ppRef;


void FillPtHistograms(const Dataset& dataset, double lowCent, double highCent, TH1D* h_MuPl, TH1D* h_MuMi);
void GetMeanAndPeakDifference(const Dataset& dataset, TH1D* h_MuPl, TH1D* h_MuMi);
void PlotMeanDifference();
void PlotPeakDifference();


//---Main function()
void MeanDifference(){

    gROOT->SetBatch(kTRUE);
    
    //Clear vectors to avoid contamination from previous runs.
    PbPbMean.clear();
    ppRefMean.clear();
    PbPbPeak.clear();
    ppRefPeak.clear();
    MeanDiffAndError_PbPb_vs_ppRef.clear();
    PeakDiffAndError_PbPb_vs_ppRef.clear();


    //Histograms to calculate statistics from.
    TH1D* h_MuPl = new TH1D("h_MuPl", "Muon Plus pT; pT [GeV]; Entries", 200, 0., 200.);
    TH1D* h_MuMi = new TH1D("h_MuMi", "Muon Minus pT; pT [GeV]; Entries", 200, 0., 200.);
    

    //Get reference values.
    FillPtHistograms(datasets[0], 0., 100., h_MuPl, h_MuMi); //ppRef2024 dataset
    GetMeanAndPeakDifference(datasets[0], h_MuPl, h_MuMi);


    //For PbPb datasets loop over centrality bins. 
    for(const auto& cBin : CentralitySet){

        //Clear histograms before each centrality call..
        h_MuPl->Reset();
        h_MuMi->Reset();

        std::cout << "\nCalculating Mean Diff for centrality bin: " << cBin.first << "-" << cBin.second << "%\n";

        FillPtHistograms(datasets[1], cBin.first, cBin.second, h_MuPl, h_MuMi); //PbPb2023 dataset
        FillPtHistograms(datasets[2], cBin.first, cBin.second, h_MuPl, h_MuMi); //PbPb2024 dataset
        FillPtHistograms(datasets[3], cBin.first, cBin.second, h_MuPl, h_MuMi); //PbPb2025 dataset
        FillPtHistograms(datasets[4], cBin.first, cBin.second, h_MuPl, h_MuMi); //PbPb2026 dataset

        //Histograms contain PbPb2023+2024 at this point.
        //datasets[1] is passed only to identify this as a PbPb result.
        GetMeanAndPeakDifference(datasets[1], h_MuPl, h_MuMi);
    }

    PlotMeanDifference();
    PlotPeakDifference();

    delete h_MuPl;
    delete h_MuMi;
}

void PlotPeakDifference(){

    //Calculate relative difference to reference sample:
    for(size_t i = 0; i < PbPbPeak.size(); ++i){
        double peakDiff = PbPbPeak[i].first - ppRefPeak[0].first;
        double peakError = std::sqrt(std::pow(PbPbPeak[i].second, 2) + std::pow(ppRefPeak[0].second, 2));
        PeakDiffAndError_PbPb_vs_ppRef.push_back(std::make_pair(peakDiff, peakError));
    }

    int nPoints = PeakDiffAndError_PbPb_vs_ppRef.size();
    if(nPoints != CentralitySet.size()){
        std::cerr << "Error: Number of points in DeltaPtAndError does not match number of centrality bins." << std::endl;
        return;
    }

    std::vector<double> xValues(nPoints);
    std::vector<double> yValues(nPoints);
    std::vector<double> xErrors(nPoints);
    std::vector<double> yErrors(nPoints);

    for(int i = 0; i < nPoints; ++i){
        xValues[i] = i+1;
        xErrors[i] = 0.;
        yValues[i] = PeakDiffAndError_PbPb_vs_ppRef[i].first;
        yErrors[i] = PeakDiffAndError_PbPb_vs_ppRef[i].second;
    }


    //TGraphErrors
    TGraphErrors* graph = new TGraphErrors(nPoints, xValues.data(), yValues.data(), xErrors.data(), yErrors.data());
    basicGraphFormatting(graph);

    //Plot
    TCanvas* c = new TCanvas("c", "Peak diff vs Centrality", 800, 600);
    basicCanvasFormatting(c);

    //Frame TH1 helper to set the x-axis labels for centrality bins.
    TH1D* frame = new TH1D("frame_peak", "", nPoints, 0.5, nPoints + 0.5);
    basicHistFormatting(frame);
    for(int i = 0; i < nPoints; ++i){
        std::string label = Form("%.0f-%.0f%%", CentralitySet[i].first, CentralitySet[i].second);
        frame->GetXaxis()->SetBinLabel(i + 1,label.c_str());
    }

    //Give range information to new frame histogram:
    double yMin = yValues[0] - yErrors[0];
    double yMax = yValues[0] + yErrors[0];
    for(int i = 1; i < nPoints; ++i){
        yMin = std::min(yMin,yValues[i] - yErrors[i]);
        yMax = std::max(yMax,yValues[i] + yErrors[i]);
    }
    double yRange = yMax - yMin;
    yMin -= 0.20 * yRange;
    yMax += 0.20 * yRange;
    
    //Make sure zero is visible:
    yMin = std::min(yMin, 0.0);
    yMax = std::max(yMax, 0.0);
    frame->SetMinimum(yMin);
    frame->SetMaximum(yMax);
    //

    frame->GetXaxis()->SetTickLength(0.0);
    frame->GetXaxis()->SetTitle("Centrality bin");
    frame->GetYaxis()->SetTitle("#Delta p_{T,peak}^{PbPb} - #Delta #bar{p_{T,peak}}^{ppRef}");
    frame->GetYaxis()->CenterTitle(true);
    frame->GetYaxis()->SetTitleOffset(1.4);

    //Draw only the axis frame.
    frame->Draw("AXIS");

    //Format graph.
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.9);
    graph->SetMarkerColorAlpha(kRed+1, 1.);
    graph->SetLineColorAlpha(kRed-7, 0.8);
    graph->SetLineWidth(2);

    //Draw
    graph->Draw("P SAME");

    //Grey line at y=0
    TLine* line = new TLine(0.5, 0.0, 0.5+nPoints, 0.0);
    line->SetLineColor(kGray);
    line->SetLineStyle(7);
    line->SetLineWidth(2);
    line->Draw();

    drawLatexText("#bf{CMS}", 0.14, 0.93, 0.04);
    drawLatexText("#it{Internal}", 0.21, 0.93, 0.03);
    drawLatexText(dataSamplesUsed.c_str(), 0.55, 0.93, 0.03);
    
    //Plot specifications
    drawLatexText("p_{T}^{#mu} > 20 GeV, |#eta^{#mu}| < 2.4", 0.8, 0.85, 0.03);
    drawLatexText("60 < M_{#mu#mu} < 120 GeV", 0.8, 0.8, 0.03);

    c->Update();
    std::string outputName = "PeakDiff_vs_Centrality_FromTree" + plot_extension;
    c->SaveAs(outputName.c_str());

    delete frame;
    delete line;
    delete graph;
    delete c;
}

void PlotMeanDifference() {

    //Calculate relative difference to reference sample:
    for(size_t i = 0; i < PbPbMean.size(); ++i){
        double meanDiff = PbPbMean[i].first - ppRefMean[0].first;
        double meanError = std::sqrt(std::pow(PbPbMean[i].second, 2) + std::pow(ppRefMean[0].second, 2));
        MeanDiffAndError_PbPb_vs_ppRef.push_back(std::make_pair(meanDiff, meanError));
    }

    int nPoints = MeanDiffAndError_PbPb_vs_ppRef.size();
    if(nPoints != CentralitySet.size()){
        std::cerr << "Error: Number of points in DeltaPtAndError does not match number of centrality bins." << std::endl;
        return;
    }

    std::vector<double> xValues(nPoints);
    std::vector<double> yValues(nPoints);
    std::vector<double> xErrors(nPoints);
    std::vector<double> yErrors(nPoints);

    for(int i = 0; i < nPoints; ++i){
        xValues[i] = i+1;
        xErrors[i] = 0.;
        yValues[i] = MeanDiffAndError_PbPb_vs_ppRef[i].first;
        yErrors[i] = MeanDiffAndError_PbPb_vs_ppRef[i].second;
    }


    //TGraphErrors
    TGraphErrors* graph = new TGraphErrors(nPoints, xValues.data(), yValues.data(), xErrors.data(), yErrors.data());
    basicGraphFormatting(graph);

    //Plot
    TCanvas* c = new TCanvas("c", "Mean diff vs Centrality", 800, 600);
    basicCanvasFormatting(c);

    //Frame TH1 helper to set the x-axis labels for centrality bins.
    TH1D* frame = new TH1D("frame_mean", "", nPoints, 0.5, nPoints + 0.5);
    basicHistFormatting(frame);
    for(int i = 0; i < nPoints; ++i){
        std::string label = Form("%.0f-%.0f%%", CentralitySet[i].first, CentralitySet[i].second);
        frame->GetXaxis()->SetBinLabel(i + 1,label.c_str());
    }

    //Give range information to new frame histogram:
    double yMin = yValues[0] - yErrors[0];
    double yMax = yValues[0] + yErrors[0];
    for(int i = 1; i < nPoints; ++i){
        yMin = std::min(yMin,yValues[i] - yErrors[i]);
        yMax = std::max(yMax,yValues[i] + yErrors[i]);
    }
    double yRange = yMax - yMin;
    yMin -= 0.20 * yRange;
    yMax += 0.20 * yRange;
    
    //Make sure zero is visible:
    yMin = std::min(yMin, 0.0);
    yMax = std::max(yMax, 0.0);
    frame->SetMinimum(-0.28);
    frame->SetMaximum(0.28);
    //

    frame->GetXaxis()->SetTickLength(0.0);
    frame->GetXaxis()->SetTitle("Centrality bin");
    frame->GetYaxis()->SetTitle("#Delta #bar{p_{T}}^{PbPb} - #Delta #bar{p_{T}}^{ppRef}");
    frame->GetYaxis()->CenterTitle(true);
    frame->GetYaxis()->SetTitleOffset(1.3);
    frame->GetYaxis()->SetLabelSize(0.03);

    //Draw only the axis frame.
    frame->Draw("AXIS");

    //Format graph.
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.9);
    graph->SetMarkerColorAlpha(kRed+1, 1.);
    graph->SetLineColorAlpha(kRed-7, 0.8);
    graph->SetLineWidth(2);

    //Grey line at y=0
    TLine* line = new TLine(0.5, 0.0, 0.5+nPoints, 0.0);
    line->SetLineColor(kGray);
    line->SetLineStyle(7);
    line->SetLineWidth(2);
    line->Draw();

    //Also plot ppRef line for reference.
    TLine* ppRefLine = new TLine(0.5, ppRefMean[0].first, 0.5+nPoints, ppRefMean[0].first);
    ppRefLine->SetLineColor(kCyan-3);
    ppRefLine->SetLineStyle(1);
    ppRefLine->SetLineWidth(2);
    ppRefLine->Draw();
    //Draw the uncertainties also as a TBox
    TBox* ppRefBand = new TBox(0.5, ppRefMean[0].first - ppRefMean[0].second, 0.5+nPoints, ppRefMean[0].first + ppRefMean[0].second);
    ppRefBand->SetFillStyle(1001);
    ppRefBand->SetFillColorAlpha(kCyan-3, 0.25);
    ppRefBand->SetLineWidth(0);
    ppRefBand->Draw();

    //Draw graph points for last.
    graph->Draw("P SAME");

    TLegend* legend = new TLegend(0.71, 0.8, 0.88, 0.9);
    basicLegendFormatting(legend);
    legend->AddEntry(graph, "PbPb 2023-2026", "p");
    legend->AddEntry(ppRefLine, "ppRef 2024", "l");
    legend->Draw();

    drawLatexText("#bf{CMS}", 0.13, 0.93, 0.042);
    drawLatexText("#it{Internal}", 0.2, 0.93, 0.034);
    drawLatexText(dataSamplesUsed.c_str(), 0.55, 0.93, 0.03);
    
    //Plot specifications
    drawLatexText("p_{T}^{#mu} > 20 GeV, |#eta^{#mu}| < 2.4", 0.2, 0.85, 0.03);
    drawLatexText("60 < M_{#mu#mu} < 120 GeV", 0.2, 0.8, 0.03);

    c->Update();
    std::string outputName = "MeanDiff_vs_Centrality_FromTree" + plot_extension;
    c->SaveAs(outputName.c_str());

    delete ppRefLine;
    delete ppRefBand;
    delete frame;
    delete line;
    delete graph;
    delete c;
}


void FillPtHistograms(const Dataset& dataset, double lowCent, double highCent, TH1D* h_MuPl, TH1D* h_MuMi) {

    // Load root file.
    std::string fullPath = dataset.basePath + dataset.filePattern;

    TChain *chain = new TChain(dataset.treeName.c_str());
    chain->Add(fullPath.c_str());

    std::cout << "> Number of files added to TChain = " << chain->GetListOfFiles()->GetEntries() << "\n" << std::endl;
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

    if(dataset.hasCentrality) {//PbPb2023, PbPb2024.
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


    //dN/dpT histograms for MuPl and MuMi.
    //TH1D* h_MuPl = new TH1D("h_MuPl", "Muon Plus pT; pT [GeV]; Entries", 100, 0., 100.);
    //TH1D* h_MuMi = new TH1D("h_MuMi", "Muon Minus pT; pT [GeV]; Entries", 100, 0., 100.);


    //Muon-level selection variables
    double ptplus, ptminus;
    double etaplus, etaminus;
    bool MuPlIsTight;
    bool MuMiIsTight;
    
    //Centrality interval passed as argument to the function.
    float minCentrality = 2.*lowCent;
    float maxCentrality = 2.*highCent;
    //

    for(Long64_t i = 0; i < nEvents; ++i){//Loop through all EVENTS in the CHAIN.

        chain->GetEntry(i); //Get event i.

        //Good event selection
        bool goodVertex = (std::abs(zVtx) < maxZvtx);
        bool goodCent = true;

        if (dataset.hasCentrality) {
            goodCent = (Centrality >= minCentrality && Centrality < maxCentrality);
        }

        if (!goodVertex) continue;
        if (!goodCent) continue;

        for(Short_t j = 0; j < Reco_Dimuon_size; ++j){ //Loop through all reco dimuon candidates of event i.
            
            //Good Z selection
            bool goodMass = (Reco_Dimuon_invMass->at(j) > minZ_Mass && Reco_Dimuon_invMass->at(j) < maxZ_Mass);
            bool goodRapidity = (std::abs(Reco_Dimuon_rapidity->at(j)) < RapidityCutValue);
            bool goodCharge = (Reco_Dimuon_sign[j] == 0); //Opposite sign muons.
            bool goodVtxProb = (Reco_Dimuon_vtxProb[j] > 0.001); //Vertex probability cut of .1% for dimuon candidates.
            bool isTriggerMatched = true;

                if (dataset.applyTrigger) {
                    //**At least** one of the daughter muons must be matched to the trigger.
                    isTriggerMatched = (Reco_Dimuon_trig[j] & dataset.triggerBit);//Corresponding triggerBit to that dataset.
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
            
            bool goodMuPl = (ptplus > ptCutValue && ptplus < 200.)
                            && (std::abs(etaplus) < EtaCutValue)
                            && (MuPlIsTight);

            bool goodMuMi = (ptminus > ptCutValue && ptminus < 200.)
                            && (std::abs(etaminus) < EtaCutValue)
                            && (MuMiIsTight);

            if (!goodMuPl || !goodMuMi) continue;
            //End of good selection for dimuon candidate j of event i.
            
            //Fill each histogram with respective muon pT.
            h_MuPl->Fill(ptplus);
            h_MuMi->Fill(ptminus);

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

    //Check the n of dimuons for this centrality bin
    std::cout << "\n> Number of dimuons in centrality bin " << lowCent << "-" << highCent << "%: " << h_MuPl->GetEntries() << std::endl;

    delete chain;
}


void GetMeanAndPeakDifference(const Dataset& dataset, TH1D* h_MuPl, TH1D* h_MuMi){

    //pT(mu+) statistics
    //Mean and error
    double MuPl_mean = h_MuPl->GetMean();
    double MuPl_meanError = h_MuPl->GetMeanError();

    //Variance and RMS
    double MuPl_variance = h_MuPl->GetStdDev();
    double MuPl_RMS = h_MuPl->GetRMS();

    //Peak and "error" (bin width)
    double MuPl_peak = h_MuPl->GetBinCenter(h_MuPl->GetMaximumBin());
    double MuPl_peakError = h_MuPl->GetBinWidth(h_MuPl->GetMaximumBin());

    //pT(mu-) statistics
    //Mean and error
    double MuMi_mean = h_MuMi->GetMean();
    double MuMi_meanError = h_MuMi->GetMeanError();

    //Variance and RMS
    double MuMi_variance = h_MuMi->GetStdDev();
    double MuMi_RMS = h_MuMi->GetRMS();

    //Peak and "error" (bin width)
    double MuMi_peak = h_MuMi->GetBinCenter(h_MuMi->GetMaximumBin());
    double MuMi_peakError = h_MuMi->GetBinWidth(h_MuMi->GetMaximumBin());
    
    //Calculate mean difference and error
    double meanDifference = MuPl_mean - MuMi_mean;
    double meanDifferenceError = std::sqrt(std::abs(MuPl_meanError * MuPl_meanError - MuMi_meanError * MuMi_meanError));

    //**Nothing** with variance and RMS for now.

    //Peak difference and error
    double peakDifference = MuPl_peak - MuMi_peak;
    double peakDifferenceError = std::sqrt(std::abs(MuPl_peakError * MuPl_peakError - MuMi_peakError * MuMi_peakError));


    //Save results in respective vectors.
    if(dataset.hasCentrality) { //PbPb datasets
        PbPbMean.push_back(std::make_pair(meanDifference, meanDifferenceError));
        PbPbPeak.push_back(std::make_pair(peakDifference, peakDifferenceError));
    }

    else if(dataset.system == CollisionSystem::ppRef){ //ppRef.
        ppRefMean.push_back(std::make_pair(meanDifference, meanDifferenceError));
        ppRefPeak.push_back(std::make_pair(peakDifference, peakDifferenceError));
    }

}