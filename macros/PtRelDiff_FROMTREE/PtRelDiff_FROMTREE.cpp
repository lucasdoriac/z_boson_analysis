
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
    {50., 100.}
};*/

//A struct is better to organize results in this case.
struct PtRelDiffResult {
    std::string centralityBinStr; // For PbPb2024 dataset
    Long64_t Zcount;
    double mean;
    double meanError;
    double variance;
    double skewness;
};


//Vectors to store the results of the PtRelDiff calculation for ppRef and PbPb datasets.
std::vector<PtRelDiffResult> ppRefResults; //dummy string, n, mean, meanError, variance, skewness.
std::vector<PtRelDiffResult> PbPbResults; //centralityBinStr, n, mean, meanError, variance, skewness.
std::vector<PtRelDiffResult> FinalResults; //centralityBinStr, n, mean, meanError, variance, skewness.

//---Function declarations
void CalculatePtRelativeDiff(const Dataset& dataset, double lowCent, double highCent);
void PlotPtRelativeDiff();
void PrintStatistics();
std::tuple<Long64_t, double, double, double, double> GetStatistics(const std::vector<double>& values);

//---Main function
void PtRelDiff_FROMTREE(){

    gROOT->SetBatch(kTRUE);
    ppRefResults.clear();
    PbPbResults.clear();
    FinalResults.clear();

    //Calculate PtRelDiff for ppRef2024 dataset.
    CalculatePtRelativeDiff(datasets[2], 0., 100.); //ppRef2024 dataset.

    //Now for PbPb2024 dataset. Loop over centrality bins. 
    for(const auto& centralityBin : CentralitySet){
        std::cout << "\nCalculating PtRelDiff_PbPb for centrality bin: " << centralityBin.first << "-" << centralityBin.second << "%\n";
        CalculatePtRelativeDiff(datasets[1], centralityBin.first, centralityBin.second); //PbPb2024 dataset
    }

    //Plot the final results.
    PlotPtRelativeDiff();

    //Print statistics for comparison with binned approach.
    PrintStatistics();
}

void PrintStatistics(){

    std::cout << std::scientific << std::setprecision(10);

    // ============================================================
    // ppRef
    // ============================================================

    std::cout << "\n\n";
    std::cout << "============================================================\n";
    std::cout << "ppRef 2024 - FROM TREE\n";
    std::cout << "============================================================\n";

    if(ppRefResults.empty()){
        std::cout << "No ppRef results available.\n";
    }
    else{

        const auto& result = ppRefResults[0];

        std::cout
            << "Z count     = " << result.Zcount    << "\n"
            << "Mean        = " << result.mean      << "\n"
            << "Mean error  = " << result.meanError << "\n"
            << "Variance    = " << result.variance  << "\n"
            << "Skewness    = " << result.skewness  << "\n";
    }


    // ============================================================
    // PbPb
    // ============================================================

    std::cout << "\n";
    std::cout << "============================================================\n";
    std::cout << "PbPb 2024 - FROM TREE\n";
    std::cout << "============================================================\n";

    std::cout
        << std::left
        << std::setw(14) << "Centrality"
        << std::setw(14) << "Z count"
        << std::setw(20) << "Mean"
        << std::setw(20) << "Mean error"
        << std::setw(20) << "Variance"
        << std::setw(20) << "Skewness"
        << "\n";

    std::cout << std::string(108, '-') << "\n";

    for(const auto& result : PbPbResults){

        std::cout
            << std::left
            << std::setw(14) << result.centralityBinStr
            << std::setw(14) << result.Zcount
            << std::setw(20) << result.mean
            << std::setw(20) << result.meanError
            << std::setw(20) << result.variance
            << std::setw(20) << result.skewness
            << "\n";
    }


    // ============================================================
    // Final observable
    // ============================================================

    std::cout << "\n";
    std::cout << "============================================================\n";
    std::cout << "FINAL OBSERVABLE - FROM TREE\n";
    std::cout << "<PtRelDiff>_PbPb - <PtRelDiff>_ppRef\n";
    std::cout << "============================================================\n";

    std::cout
        << std::left
        << std::setw(14) << "Centrality"
        << std::setw(14) << "Z count"
        << std::setw(20) << "Delta mean"
        << std::setw(20) << "Mean error"
        << std::setw(20) << "Var(delta)"
        << std::setw(20) << "Delta skew."
        << "\n";

    std::cout << std::string(108, '-') << "\n";

    for(const auto& result : FinalResults){

        std::cout
            << std::left
            << std::setw(14) << result.centralityBinStr
            << std::setw(14) << result.Zcount
            << std::setw(20) << result.mean
            << std::setw(20) << result.meanError
            << std::setw(20) << result.variance
            << std::setw(20) << result.skewness
            << "\n";
    }

    std::cout << "\n";

}

void PlotPtRelativeDiff(){

    //Before the actual plot we have to take the difference of PbPb values from ppRef values.
    const auto& ppRefResult = ppRefResults[0];
    PtRelDiffResult finalResult;

    for (size_t i = 0; i < PbPbResults.size(); ++i) {
        const auto& PbPbResult = PbPbResults[i];
        finalResult.centralityBinStr = PbPbResult.centralityBinStr;
        finalResult.Zcount = PbPbResult.Zcount;
        finalResult.mean = PbPbResult.mean - ppRefResult.mean;
        finalResult.meanError = std::sqrt(std::pow(PbPbResult.meanError, 2) + std::pow(ppRefResult.meanError, 2));
        finalResult.variance = finalResult.meanError * finalResult.meanError;
        finalResult.skewness = PbPbResult.skewness - ppRefResult.skewness;
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
        xErrors[i] = 0.
        yValues[i] = FinalResults[i].mean;
        yErrors[i] = FinalResults[i].meanError;
    }

    /*std::vector<double> xValues, yValues, xErrors, yErrors;
    int idx = 0;
    for (const auto& result : FinalResults) {
        double lowCent = CentralitySet[idx].first;
        double highCent = CentralitySet[idx].second;
        double centralityBinCenter = (lowCent + highCent) / 2.0;

        xValues.push_back(centralityBinCenter);
        yValues.push_back(result.mean);
        xErrors.push_back((highCent - lowCent) / 2.0); // Half-width of the centrality bin
        yErrors.push_back(result.meanError);

        ++idx;
    }*/

    //TGraphErrors
    //TGraphErrors* graph = new TGraphErrors(xValues.size(), xValues.data(), yValues.data(), xErrors.data(), yErrors.data());
    TGraphErrors* graph = new TGraphErrors(nPoints, xValues.data(), yValues.data(), xErrors.data(), yErrors.data());
    basicGraphFormatting(graph);

    //Plot
    TCanvas* c = new TCanvas("c", "Pt Relative Difference FROMTREE vs Centrality", 800, 600);
    basicCanvasFormatting(c);
    c->SetLeftMargin(0.13);

    //Frame TH1 helper to set the x-axis labels for centrality bins.
    TH1D* frame = new TH1D("frame","",nPoints,0.5,nPoints + 0.5);
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
    frame->GetYaxis()->SetTitle("#Delta p^{rel}_{T,PbPb} - #Delta p^{rel}_{T,ppRef}");
    frame->GetYaxis()->CenterTitle(true);
    frame->GetYaxis()->SetTitleOffset(1.5);

    //Draw only the axis frame.
    frame->Draw("AXIS");

    //Format graph.
    graph->SetMarkerStyle(21);
    graph->SetMarkerSize(0.9);
    graph->SetMarkerColorAlpha(kRed+1, 1.);
    graph->SetLineColorAlpha(kRed-7, 0.8);
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

    drawLatexText("#bf{CMS}", 0.12, 0.93, 0.042);
    drawLatexText("#it{Work in Progress}", 0.2, 0.93, 0.033);
    drawLatexText("PbPb 2024, ppRef 2024 (5.36 TeV)", 0.6, 0.93, 0.033);
    
    //Plot specifications
    drawLatexText("p_{T}^{#mu} > 20 GeV, |#eta^{#mu}| < 2.4", 0.2, 0.8, 0.03);
    drawLatexText("60 < M_{#mu #mu} < 120 GeV", 0.2, 0.75, 0.03);

    c->Update();
    std::string outputName = "DeltaPtRelDiff_vs_Centrality_FROMTREE" + plot_extension;
    c->SaveAs(outputName.c_str());

    delete graph;
    delete c;
}

void CalculatePtRelativeDiff(const Dataset& dataset, double lowCent, double highCent){

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

    if(dataset.system == CollisionSystem::PbPb2024) {//Centrality is only defined for PbPb2024 dataset.
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

    //Trigger selection: 'L2SingleMu12'.
    ULong64_t triggerBit = 1ULL << 7;

    //Centrality interval for PbPb2024 dataset. For ppRef2024 dataset, this is ignored.
    const float minCentrality = 2.*lowCent;
    const float maxCentrality = 2.*highCent;

    //Vector to calculate statistics from.
    std::vector<double> ptRelDiffValues;


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

    //Calculate statistics from the collected PtRelDiff values.
    auto [n, mean, meanError, variance, skewness] = GetStatistics(ptRelDiffValues);

    PtRelDiffResult result;
    if(dataset.system == CollisionSystem::PbPb2024){
        //Centrality to string
        std::string centString = std::to_string(static_cast<int>(lowCent)) + "-" + std::to_string(static_cast<int>(highCent));

        result.centralityBinStr = centString;
        result.Zcount = n;
        result.mean = mean;
        result.meanError = meanError;
        result.variance = variance;
        result.skewness = skewness;

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
        ppRefResults.push_back(result);
    }

    delete chain;
}

std::tuple<Long64_t, double, double, double, double> GetStatistics(const std::vector<double>& values) {

    const Long64_t n = values.size();

    //Mean (M1)
    double mean = 0.0;

    for(double x : values)
    {
        mean += x;
    }

    mean = mean/static_cast<double>(n);


    //Central moments (M2, M3)
    double M2 = 0.0;
    double M3 = 0.0;

    for(double x : values){

        double delta = x - mean;

        M2 += delta * delta;
        M3 += delta * delta * delta;
    }

    //Variance (M2/n)
    double variance = M2/static_cast<double>(n); // std dev = \sqrt(variance).

    // Moment skewness (M3/(M2^(3/2)))
    double skewness = std::sqrt(static_cast<double>(n))*M3/std::pow(M2, 1.5);

    //Mean error (sqrt(variance/n))
    double meanError = std::sqrt(variance/static_cast<double>(n));

    return std::make_tuple(n, mean, meanError, variance, skewness);
}