/*
Peak difference directly from the TTree.
I was using pT histograms with bin width of 1 GeV until now, but the peak difference is not very precise with that.
Since we are using Run 3 data now and we increased the statistics,
we can increase the granularity of the pT histograms to 0.4 GeV and get a more precise peak difference (I hope).
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


//Vectors to store peak shift and measurement error for each centrality bin.
std::vector<std::pair<double, double>> ppRefValues;
std::vector<std::pair<double, double>> PbPbValues;
std::vector<std::pair<double, double>> PeakDiffAndError;


//---Function declarations
void FillPtHistograms(const Dataset& dataset, double lowCent, double highCent, TH1D* h_MuPl, TH1D* h_MuMi);
void GetPeakDifference(const Dataset& dataset, TH1D* h_MuPl, TH1D* h_MuMi);
void PlotPeakShift();

//---Main()
void NewPeakDiff(){

    gROOT->SetBatch(kTRUE);

    //Histograms to calculate peak shift from.
    TH1D* h_MuPl = new TH1D("h_MuPl", "Muon Plus pT; pT [GeV]; Entries", 100, 20., 70.); //0.5 GeV bin width.
    TH1D* h_MuMi = new TH1D("h_MuMi", "Muon Minus pT; pT [GeV]; Entries", 100, 20., 70.);
    FillPtHistograms(datasets[0], 0., 100., h_MuPl, h_MuMi); //ppRef2024 dataset
    GetPeakDifference(datasets[0], h_MuPl, h_MuMi);
    

    for(const auto& cBin : CentralitySet){

        //Clear histograms before each centrality call..
        h_MuPl->Reset();
        h_MuMi->Reset();

        std::cout << "\nCalculating Peak Diff for centrality bin: " << cBin.first << "-" << cBin.second << "%\n";

        FillPtHistograms(datasets[1], cBin.first, cBin.second, h_MuPl, h_MuMi); //PbPb2023 dataset
        FillPtHistograms(datasets[2], cBin.first, cBin.second, h_MuPl, h_MuMi); //PbPb2024 dataset
        FillPtHistograms(datasets[3], cBin.first, cBin.second, h_MuPl, h_MuMi); //PbPb2025 dataset
        FillPtHistograms(datasets[4], cBin.first, cBin.second, h_MuPl, h_MuMi); //PbPb2026 dataset

        //Sanity check: print number of selected dimuons.
        std::cout << "\n> Number of selected dimuons in this centrality bin: " << h_MuPl->Integral() << std::endl;

        //Another sanity check: number of muons in the bin 40 of the pT histogram.
        int bin40 = h_MuPl->FindBin(40.0);
        std::cout << "> Number of selected muons(+) in the pT bin 40-40.4 GeV: " << h_MuPl->GetBinContent(bin40) << std::endl;
        std::cout << "> Number of selected muons(-) in the pT bin 40-40.4 GeV: " << h_MuMi->GetBinContent(bin40) << std::endl;

        //Histograms contain all Run 3 PbPb data at this point for the selected centrality.
        GetPeakDifference(datasets[1], h_MuPl, h_MuMi);
    }

    PlotPeakShift();

    delete h_MuPl;
    delete h_MuMi;
}

void PlotPeakShift(){

    int nPoints = PbPbValues.size();
    std::vector<double> xValues(nPoints);
    std::vector<double> yValues(nPoints);
    std::vector<double> xErrors(nPoints);
    std::vector<double> yErrors(nPoints);

    for(int i = 0; i < nPoints; ++i){
        xValues[i] = i+1;
        xErrors[i] = 0.;
        yValues[i] = PbPbValues[i].first - ppRefValues[0].first;

        std::cout << "\n PbPb peak diff: " << PbPbValues[i].first << " GeV, ppRef peak diff: " << ppRefValues[0].first << " GeV" << std::endl;
        std::cout << " Peak diff (PbPb - ppRef): " << yValues[i] << " GeV" << std::endl;
        //Add errors in quadrature.
        yErrors[i] = std::sqrt(std::pow(PbPbValues[i].second, 2) + std::pow(ppRefValues[0].second, 2));
    }

    //TGraphErrors
    TGraphErrors* graph = new TGraphErrors(nPoints, xValues.data(), yValues.data(), xErrors.data(), yErrors.data());
    basicGraphFormatting(graph);

    //Plot
    TCanvas* c = new TCanvas("c", "Peak diff vs Centrality", 800, 600);
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
    frame->SetMinimum(yMin);
    frame->SetMaximum(yMax);
    //

    frame->GetXaxis()->SetTickLength(0.0);
    frame->GetXaxis()->SetTitle("Centrality bin");
    frame->GetYaxis()->SetTitle("#Delta peak p_{T} (PbPb - ppRef) [GeV]");
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

    //Draw graph points for last.
    graph->Draw("P SAME");

    drawLatexText("#bf{CMS}", 0.13, 0.93, 0.042);
    drawLatexText("#it{Internal}", 0.2, 0.93, 0.034);
    drawLatexText(dataSamplesUsed.c_str(), 0.55, 0.93, 0.03);
    
    //Plot specifications
    drawLatexText("p_{T}^{#mu} > 20 GeV, |#eta^{#mu}| < 2.4", 0.6, 0.85, 0.03);
    drawLatexText("60 < M_{#mu#mu} < 120 GeV", 0.6, 0.8, 0.03);

    c->Update();
    std::string outputName = "PeakDiff_vs_Centrality" + plot_extension;
    c->SaveAs(outputName.c_str());

    delete frame;
    delete line;
    delete graph;
    delete c;
}


void FillPtHistograms(const Dataset& dataset, double lowCent, double highCent, TH1D* h_MuPl, TH1D* h_MuMi){

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

    if(dataset.hasCentrality) {//PbPb.
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
    //TH1D* h_MuPl = new TH1D("h_MuPl", "Muon Plus pT; pT [GeV]; Entries", 100, 20., 60.); //0.4 GeV bin width.
    //TH1D* h_MuMi = new TH1D("h_MuMi", "Muon Minus pT; pT [GeV]; Entries", 100, 20., 60.);


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

    delete chain;
}


void GetPeakDifference(const Dataset& dataset, TH1D* h_MuPl, TH1D* h_MuMi){

    //Find peak position of each histogram and calculate the difference.
    double peakMuPl = h_MuPl->GetBinCenter(h_MuPl->GetMaximumBin());
    double peakMuMi = h_MuMi->GetBinCenter(h_MuMi->GetMaximumBin());
    
    double peakDiff = peakMuPl - peakMuMi;

    //Measurement error set to be the bin width of the histograms (0.4 GeV).
    double peakError = h_MuPl->GetBinWidth(1);

    std::cout << "\n> Peak MuPl: " << peakMuPl << " GeV" << std::endl;
    std::cout << "> Peak MuMi: " << peakMuMi << " GeV" << std::endl;
    std::cout << "> Peak Difference (MuPl - MuMi): " << peakDiff << " GeV" << std::endl;

    //Save results in respective vectors for ppRef and PbPb datasets.
    if(dataset.system == CollisionSystem::ppRef){
        ppRefValues.push_back(std::make_pair(peakDiff, peakError));
    } 
    
    else if (dataset.system == CollisionSystem::PbPb){
        PbPbValues.push_back(std::make_pair(peakDiff, peakError));
    }
}