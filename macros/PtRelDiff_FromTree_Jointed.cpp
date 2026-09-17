/*
Makes PtRelDiff (relative to ppRef) vs centrality bin without using histograms.
Calculates directly from the Tree.
Calculates mean and skewness of the PtRelDiff distribution for each centrality bin.
Analyzes both PbPb2023 and PbPb2024 datasets.
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
#include <tuple>
#include <iomanip>
#include "../headers/basicFormatting.h"


//---Macro settings
std::string plot_extension = ".pdf"; // ".png" for regular development and ".pdf" for final quality plots
std::string BasePath = "/home/lucas/Documents/CMS/z_boson_analysis/"; //IFT
//std::string BasePath = "/home/lucasdoriac/z_boson_analysis/data/"; //Home
std::string dataSamplesUsed = "PbPb 2023+2024, ppRef 2024 (5.36 TeV)";


//Good selection threshold values
const double MAX_ZVTX = 15.0;

const double MINZ_MASS = 60.;
const double MAXZ_MASS = 120.;
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

    bool hasCentrality;//Or maybe is AA
    bool applyTrigger;
    ULong64_t triggerBit;
};

Dataset datasets[] = {
    {
        "PbPb2023_Data",
        SampleType::Data,
        CollisionSystem::PbPb2023,
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
        CollisionSystem::PbPb2024,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_PbPb2024Data.root",
        BasePath + "Data/PbPb2024/",
        true,
        true,
        1ULL << 7 //'HLT_HIL2SingleMu12_v'
    },

    {
        "ppRef2024_Data",
        SampleType::Data,
        CollisionSystem::ppRef2024,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_ppRef2024.root",
        BasePath + "Data/ppRef2024/",
        false,
        false,
        0ULL
    },

    {
        "PbPb2024_MC",
        SampleType::MC,
        CollisionSystem::PbPb2024,
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
        CollisionSystem::ppRef2024,
        "hionia/myTree",
        "Oniatree_PowhegZtoMuMu_ppRef2024_*.root",
        BasePath + "MC/ppRef2024/DYToMuMu_M-50_TuneCP5_5p36TeV_powheg-pythia8/Powheg_ppRefPileup_March20/260320_125046/0000/",
        false,
        false,
        0ULL
    }
};

//Centrality bins for PbPb2024 data
std::vector<std::pair<double, double>> CentralitySet = {
    {0., 10.},
    {10., 20.},
    {20., 30.},
    {30., 100.},
    {0., 100.}
};

/*std::vector<std::pair<double, double>> CentralitySet = {
    {0., 10.},
    {10., 30.},
    {30., 50.},
    {50., 100.},
    {0., 100.}
};*/

//A struct is better to organize results in this case.
struct PtRelDiffResult {
    std::string centralityBinStr; // For PbPb2024 dataset
    Long64_t Zcount;
    double mean;
    double meanError;
    double variance;
    double skewness;
    double skewnessError;
};


//Vectors to store the results of the PtRelDiff calculation for ppRef and PbPb datasets.
std::vector<PtRelDiffResult> ppRefResults; //dummy string, n, mean, meanError, variance, skewness, skewnessError.
std::vector<PtRelDiffResult> PbPbResults; //centralityBinStr, n, mean, meanError, variance, skewness, skewnessError.
std::vector<PtRelDiffResult> FinalResults; //centralityBinStr, n, mean, meanError, variance, skewness, skewnessError.

//---Function declarations
void CalculatePtRelativeDiff(const Dataset& dataset, double lowCent, double highCent, std::vector<double>& ptRelDiffValues);
void FillResultStructFromVector(bool hasCentrality, double lowCent, double highCent, const std::vector<double>& ptRelDiffValues);
std::tuple<Long64_t, double, double, double, double, double> GetStatistics(const std::vector<double>& values);
void PlotPtRelativeDiff();
void PlotSkewnessRelativeDiff();


//---Main function
void PtRelDiff_FromTree_Jointed(){

    gROOT->SetBatch(kTRUE);
    ppRefResults.clear();
    PbPbResults.clear();
    FinalResults.clear();

    //Calculate PtRelDiff for ppRef2024 dataset.
    std::vector<double> ptRelDiffValues_ppRef;
    CalculatePtRelativeDiff(datasets[2], 0., 100., ptRelDiffValues_ppRef); //ppRef2024 dataset.
    FillResultStructFromVector(false, 0., 100., ptRelDiffValues_ppRef); //ppRef2024 dataset.

    //Now for PbPb2024 dataset. Loop over centrality bins. 
    for(const auto& centralityBin : CentralitySet){

        std::cout << "\nCalculating PtRelDiff for centrality bin: " << centralityBin.first << "-" << centralityBin.second << "%\n";

        //Vector to calculate statistics from.
        std::vector<double> ptRelDiffValues;
        CalculatePtRelativeDiff(datasets[0], centralityBin.first, centralityBin.second, ptRelDiffValues); //PbPb2023 dataset
        CalculatePtRelativeDiff(datasets[1], centralityBin.first, centralityBin.second, ptRelDiffValues); //PbPb2024 dataset

        //Fill a struct that contains statistics for each centrality.
        FillResultStructFromVector(true, centralityBin.first, centralityBin.second, ptRelDiffValues);
    }

    //Plot the final results.
    //One small caution: The final results are filled within the PlotPtRelativeDiff() function.
    PlotPtRelativeDiff();
    PlotSkewnessRelativeDiff();
}

void PlotSkewnessRelativeDiff(){
    //Values for TGraphErrors are just skewness diff and skewness error diff.
    int nPoints = FinalResults.size();
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
        yValues[i] = FinalResults[i].skewness;
        yErrors[i] = FinalResults[i].skewnessError;
    }


    //TGraphErrors
    TGraphErrors* graph = new TGraphErrors(nPoints, xValues.data(), yValues.data(), xErrors.data(), yErrors.data());
    basicGraphFormatting(graph);

    //Plot
    TCanvas* c = new TCanvas("c", "Skewness diff vs Centrality", 800, 600);
    basicCanvasFormatting(c);
    c->SetLeftMargin(0.13);

    //Frame TH1 helper to set the x-axis labels for centrality bins.
    TH1D* frame = new TH1D("frame_skewness","",nPoints,0.5,nPoints + 0.5);
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
    frame->GetYaxis()->SetTitle("#Delta #gamma_{i}");
    frame->GetYaxis()->CenterTitle(true);
    frame->GetYaxis()->SetTitleOffset(1.4);

    //Draw only the axis frame.
    frame->Draw("AXIS");

    //Format graph.
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.9);
    graph->SetMarkerColorAlpha(kGreen+1, 1.);
    graph->SetLineColorAlpha(kGreen-7, 0.8);
    graph->SetLineWidth(2);

    //Axes configurations
    //graph->GetXaxis()->SetLimits(0, 100);
    //graph->GetXaxis()->SetTitle("Centrality (%)");
    //graph->GetYaxis()->SetTitle("#LT#Delta p_{T}^{rel}#GT_{PbPb} - #LT#Delta p_{T}^{rel}#GT_{ppRef}");
    //graph->GetYaxis()->CenterTitle(true);
    //graph->GetYaxis()->SetTitleOffset(1.4);

    graph->Draw("P SAME");

    //Grey line at y=0
    TLine* line = new TLine(0.5, 0.0, 0.5+nPoints, 0.0);
    line->SetLineColor(kGray);
    line->SetLineStyle(7);
    line->SetLineWidth(2);
    line->Draw();

    drawLatexText("#bf{CMS}", 0.14, 0.93, 0.04);
    drawLatexText("#it{Work in Progress}", 0.21, 0.93, 0.03);
    drawLatexText(dataSamplesUsed.c_str(), 0.55, 0.93, 0.03);
    
    //Plot specifications
    drawLatexText("p_{T}^{#mu} > 20 GeV, |#eta^{#mu}| < 2.4", 0.2, 0.8, 0.03);
    drawLatexText("60 < M_{#mu#mu} < 120 GeV", 0.2, 0.75, 0.03);

    c->Update();
    std::string outputName = "SkewnessDiff_vs_Centrality_FROMTREE" + plot_extension;
    c->SaveAs(outputName.c_str());

    delete frame;
    delete line;
    delete graph;
    delete c;
}

void PlotPtRelativeDiff(){

    //Before the actual plot we have to take the difference of PbPb values from ppRef values.
    const auto& ppRefResult = ppRefResults[0];
    PtRelDiffResult finalResult;

    //Fill finalResult vector.
    for (size_t i = 0; i < PbPbResults.size(); ++i) {
        const auto& PbPbResult = PbPbResults[i];
        finalResult.centralityBinStr = PbPbResult.centralityBinStr;
        finalResult.Zcount = PbPbResult.Zcount;
        finalResult.mean = PbPbResult.mean - ppRefResult.mean;
        finalResult.meanError = std::sqrt(std::pow(PbPbResult.meanError, 2) + std::pow(ppRefResult.meanError, 2));
        finalResult.variance = finalResult.meanError * finalResult.meanError;
        finalResult.skewness = PbPbResult.skewness - ppRefResult.skewness;
        finalResult.skewnessError = std::sqrt(std::pow(PbPbResult.skewnessError, 2) + std::pow(ppRefResult.skewnessError, 2));
        FinalResults.push_back(finalResult);
    }

    //Values for TGraphErrors are just mean diff and mean error diff so
    int nPoints = FinalResults.size();
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
        yValues[i] = FinalResults[i].mean;
        yErrors[i] = FinalResults[i].meanError;
    }


    //TGraphErrors
    TGraphErrors* graph = new TGraphErrors(nPoints, xValues.data(), yValues.data(), xErrors.data(), yErrors.data());
    basicGraphFormatting(graph);

    //Plot
    TCanvas* c = new TCanvas("c", "Pt Relative Difference FROMTREE vs Centrality", 800, 600);
    basicCanvasFormatting(c);
    c->SetLeftMargin(0.13);

    //Frame TH1 helper to set the x-axis labels for centrality bins.
    TH1D* frame = new TH1D("frame_mean","",nPoints,0.5,nPoints + 0.5);
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
    frame->GetYaxis()->SetTitle("#Delta #LT p^{rel}_{T} #GT^{i}");
    frame->GetYaxis()->CenterTitle(true);
    frame->GetYaxis()->SetTitleOffset(1.4);

    //Draw only the axis frame.
    frame->Draw("AXIS");

    //Format graph.
    graph->SetMarkerStyle(21);
    graph->SetMarkerSize(0.9);
    graph->SetMarkerColorAlpha(kRed+1, 1.);
    graph->SetLineColorAlpha(kRed-7, 0.8);
    graph->SetLineWidth(2);

    //Draw graph on top of frame.
    graph->Draw("P SAME");

    //Grey line at y=0
    TLine* line = new TLine(0.5, 0.0, 0.5+nPoints, 0.0);
    line->SetLineColor(kGray);
    line->SetLineStyle(7);
    line->SetLineWidth(2);
    line->Draw();

    drawLatexText("#bf{CMS}", 0.14, 0.93, 0.04);
    drawLatexText("#it{Work in Progress}", 0.21, 0.93, 0.03);
    drawLatexText(dataSamplesUsed.c_str(), 0.55, 0.93, 0.03);
    
    //Plot specifications
    drawLatexText("p_{T}^{#mu} > 20 GeV, |#eta^{#mu}| < 2.4", 0.2, 0.8, 0.03);
    drawLatexText("60 < M_{#mu#mu} < 120 GeV", 0.2, 0.75, 0.03);

    c->Update();
    std::string outputName = "DeltaPtRelDiff_vs_Centrality_FROMTREE" + plot_extension;
    c->SaveAs(outputName.c_str());

    delete frame;
    delete line;
    delete graph;
    delete c;
}

void CalculatePtRelativeDiff(const Dataset& dataset, double lowCent, double highCent, std::vector<double>& ptRelDiffValues) {

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


    //Event-level cut values.
    double maxZvtx = MAX_ZVTX;

    //Z-level cut values.
    double minZ_Mass = MINZ_MASS; 
    double maxZ_Mass = MAXZ_MASS; 
    double RapidityCutValue = RAPIDITYCUTVALUE;

    //Muon-level cut values.
    float EtaCutValue = ETACUTVALUE;
    double ptCutValue = PTCUTVALUE;
    bool MuPlIsTight;
    bool MuMiIsTight;

    //Observables
    TLorentzVector Z;
    double ptplus, ptminus;
    double etaplus, etaminus;

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
            bool goodCharge = (Reco_Dimuon_sign[j] == 0);
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
            
            bool goodMuPl = (ptplus > ptCutValue)
                            && (std::abs(etaplus) < EtaCutValue)
                            && (MuPlIsTight);

            bool goodMuMi = (ptminus > ptCutValue)
                            && (std::abs(etaminus) < EtaCutValue)
                            && (MuMiIsTight);

            if (!goodMuPl || !goodMuMi) continue;
            //End of good selection for dimuon candidate j of event i.
            
            //Store the PtRelDiff value.
            ptRelDiffValues.push_back(Reco_Dimuon_muonPtRelDiff->at(j));
        
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


    delete chain;
}


void FillResultStructFromVector(bool hasCentrality, double lowCent, double highCent, const std::vector<double>& ptRelDiffValues){
    
    //Get statistics from the vector of PtRelDiff values.
    auto [n, mean, meanError, variance, skewness, skewnessError] = GetStatistics(ptRelDiffValues);

    PtRelDiffResult result;
    if(hasCentrality) {
        //Centrality to string
        std::string centString = std::to_string(static_cast<int>(lowCent)) + "-" + std::to_string(static_cast<int>(highCent));

        result.centralityBinStr = centString;
        result.Zcount = n;
        result.mean = mean;
        result.meanError = meanError;
        result.variance = variance;
        result.skewness = skewness;
        result.skewnessError = skewnessError;

        PbPbResults.push_back(result);
    }

    else {
        //Store ppRef results on its own vector.
        std::string dummyString = "ppRef";

        result.centralityBinStr = dummyString;
        result.Zcount = n;
        result.mean = mean;
        result.meanError = meanError;
        result.variance = variance;
        result.skewness = skewness;
        result.skewnessError = skewnessError;

        ppRefResults.push_back(result);
    }
}

std::tuple<Long64_t, double, double, double, double, double> GetStatistics(const std::vector<double>& values) {

    const Long64_t n = values.size();

    double mean = 0.0;
    for(double x : values)
    {
        mean += x;
    }
    mean = mean/static_cast<double>(n);

    //Test for 0-100% cent bin:
    //good: -0.0026304
    //raw: -0.0241666

    double M2 = 0.0;
    double M3 = 0.0;
    for(double x : values){

        double delta = x - mean;
        M2 += delta * delta;//Central moment of order 2 (variance).
        M3 += delta * delta * delta;//Central moment of order 3 (skewness).
    }

    double variance = M2/static_cast<double>(n-1);
    //The variance error depends on the fourth central moment, which we are not calculating here. For now, we will NOT calculate the variance error.

    //The calculation of skewness is a bit complicated.
    //It is defined as
    // g1 = sqrt(n) * M3 / M2^(3/2)
    //For now im calculating the skewness error in the same way as ROOT.
    //However this is a simplification of the calculation and if we decide to keep calculating the skewness
    //in the future we can think of a better way to estimate its error.
    double skewness = std::sqrt(static_cast<double>(n))*M3/std::pow(M2, 1.5);
    double skewnessError = std::sqrt(6.0 / static_cast<double>(n));

    double meanError = std::sqrt(variance/static_cast<double>(n));//Std dev. = sqrt(variance), Std dev. of the mean = std dev./sqrt(n).

    return std::make_tuple(n, mean, meanError, variance, skewness, skewnessError);
}