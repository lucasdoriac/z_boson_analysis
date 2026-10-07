/*
Search for asymmetry in ptplus - ptminus, event by event.
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
#include <stdexcept>
#include "../headers/basicFormatting.h"



//Location of datasets
std::string BasePath = "/home/lucas/Documents/CMS/z_boson_analysis/"; //IFT
//std::string BasePath = "/home/lucasdoriac/z_boson_analysis/data/"; //Home

std::string plot_extension = ".pdf"; // ".png" for regular development and ".pdf" for final quality plots
std::string dataSamplesUsed = "PbPb 2023-2026, ppRef 2024 (5.36 TeV)";

double rho = 1.;

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


//Centrality intervals for PbPb datasets.
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


//Struct for PtDiffCounts
struct PtDiffCount {
    int nPtDiffPos = 0;
    int nPtDiffNeg = 0;
    int zero = 0;
};

PtDiffCount counts_ppRef;
double PtDiffAsymmetry_ppRef = 0.;
double PtDiffAsymmetryError_ppRef = 0.;


//---Function declarations
void FillPtDiffHistogram(const Dataset& dataset, double lowCent, double highCent, PtDiffCount& counts);
std::pair<double, double> CalculateSignAsymmetry(const PtDiffCount& counts);
void PlotPtDiffAsymmetry(const std::vector<std::pair<double, double>>& asymmetryResults);


//---Main()
void PtDiffAsymmetry(){

    gROOT->SetBatch(kTRUE);

    //First on ppRef
    FillPtDiffHistogram(datasets[0], 0., 100., counts_ppRef);
    auto ppResult = CalculateSignAsymmetry(counts_ppRef);
    PtDiffAsymmetry_ppRef = ppResult.first;
    PtDiffAsymmetryError_ppRef = ppResult.second;
    //
    std::cout << "\n PtDiff Asymmetry in ppRef2024: " << PtDiffAsymmetry_ppRef << " +/- " << PtDiffAsymmetryError_ppRef << std::endl;


    //Run on PbPb datasets for each centrality interval.
    std::vector<std::pair<double, double>> asymmetryResults;
    for(const auto& cBin : CentralitySet) {

        std::cout << "\n Processing centrality interval: " << cBin.first << "-" << cBin.second << "%" << std::endl;

        PtDiffCount counts;
        FillPtDiffHistogram(datasets[1], cBin.first, cBin.second, counts);
        FillPtDiffHistogram(datasets[2], cBin.first, cBin.second, counts);
        FillPtDiffHistogram(datasets[3], cBin.first, cBin.second, counts);
        FillPtDiffHistogram(datasets[4], cBin.first, cBin.second, counts);
     
        asymmetryResults.push_back(CalculateSignAsymmetry(counts));
    }

    PlotPtDiffAsymmetry(asymmetryResults);
}


void PlotPtDiffAsymmetry(const std::vector<std::pair<double, double>>& asymmetryResults) {

    int nPoints = asymmetryResults.size();
    if(nPoints != CentralitySet.size()){
        std::cerr << "Error: Number of points in asymmetryResults does not match number of centrality bins." << std::endl;
        return;
    }

    std::vector<double> xValues(nPoints);
    std::vector<double> yValues(nPoints);
    std::vector<double> xErrors(nPoints);
    std::vector<double> yErrors(nPoints);

    for(int i = 0; i < nPoints; ++i){
        xValues[i] = i+1;
        xErrors[i] = 0.;
        yValues[i] = asymmetryResults[i].first - PtDiffAsymmetry_ppRef; //Subtract ppRef asymmetry
        yErrors[i] = std::hypot(asymmetryResults[i].second, PtDiffAsymmetryError_ppRef);
    }


    //TGraphErrors
    TGraphErrors* graph = new TGraphErrors(nPoints, xValues.data(), yValues.data(), xErrors.data(), yErrors.data());
    basicGraphFormatting(graph);

    TCanvas* c = new TCanvas("c", "PtDiff Asymmetry", 800, 600);
    basicCanvasFormatting(c);
    c->SetLeftMargin(0.13);

    //Frame TH1 helper to set the x-axis labels for centrality bins.
    TH1D* frame = new TH1D("frame_ptdiff","",nPoints,0.5,nPoints + 0.5);
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
    frame->SetMinimum(-0.012);
    frame->SetMaximum(0.012);
    //

    frame->GetXaxis()->SetTickLength(0.0);
    frame->GetXaxis()->SetTitle("Centrality bin");
    frame->GetYaxis()->SetTitle("PtDiff Asymmetry");
    frame->GetYaxis()->CenterTitle(true);
    frame->GetYaxis()->SetTitleOffset(1.4);

    //Draw only the axis frame.
    frame->Draw("AXIS");

    //Format graph.
    graph->SetMarkerStyle(21);
    graph->SetMarkerSize(1.);
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
    drawLatexText("#it{Internal}", 0.21, 0.93, 0.033);
    drawLatexText(dataSamplesUsed.c_str(), 0.55, 0.93, 0.03);
    
    //Plot specifications
    drawLatexText("p_{T}^{#mu} > 20 GeV, |#eta^{#mu}| < 2.4", 0.7, 0.8, 0.03);
    drawLatexText("60 < M_{#mu#mu} < 120 GeV", 0.7, 0.75, 0.03);

    c->Update();
    std::string outputName = "PtDiff_vs_Centrality" + plot_extension;
    c->SaveAs(outputName.c_str());

    delete frame;
    delete line;
    delete graph;
    delete c;
}


std::pair<double, double> CalculateSignAsymmetry(const PtDiffCount& counts){

    double N_plus = static_cast<double>(counts.nPtDiffPos);
    double N_minus = static_cast<double>(counts.nPtDiffNeg);
    double total = N_plus + N_minus;

    double N_plusError = std::sqrt(N_plus); //Poisson error for N_plus
    double N_minusError = std::sqrt(N_minus); //Poisson error for N_minus

    if (total <= 0.) {//Protect against division by zero.
        std::cerr << "Error: Total count of events is zero or negative. Cannot calculate asymmetry." << std::endl;
        return {0., 0.};
    }

    double asymmetry = (N_plus - N_minus)/total;
    double asymmetryError = 2. * std::abs(N_minus * N_plusError - N_plus * N_minusError) / (total * total);
    //double asymmetryError = std::sqrt((1. - asymmetry * asymmetry) / total);

    std::cout << "A = " << asymmetry << " +/- " << asymmetryError << std::endl;

    return {asymmetry, asymmetryError};
}

void FillPtDiffHistogram(const Dataset& dataset, double lowCent, double highCent, PtDiffCount& counts){

    // Load root file.
    std::string fullPath = dataset.basePath + dataset.filePattern;

    TChain *chain = new TChain(dataset.treeName.c_str());
    chain->Add(fullPath.c_str());

    std::cout << "> Number of files on TChain = " << chain->GetListOfFiles()->GetEntries() << "\n" << std::endl;
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


    //Muon variables
    double ptplus, ptminus;
    double etaplus, etaminus;
    bool MuPlIsTight;
    bool MuMiIsTight;

    //Centrality interval passed as argument to the function.
    float minCentrality = 2.*lowCent;
    float maxCentrality = 2.*highCent;


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
            
            bool goodMuPl = (ptplus > ptCutValue && ptplus < 200.)
                            && (std::abs(etaplus) < EtaCutValue)
                            && (MuPlIsTight);

            bool goodMuMi = (ptminus > ptCutValue && ptminus < 200.)
                            && (std::abs(etaminus) < EtaCutValue)
                            && (MuMiIsTight);

            if (!goodMuPl || !goodMuMi) continue;
            //End of good selection for dimuon candidate j of event i.
            
            //Calculate PtDiff and increase the count.
            double PtDiff = Reco_Dimuon_muonPtDiff->at(j); //ptplus - ptminus

            if (PtDiff > 0.) ++counts.nPtDiffPos;
            else if (PtDiff < 0.) ++counts.nPtDiffNeg;
            else ++counts.zero;

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