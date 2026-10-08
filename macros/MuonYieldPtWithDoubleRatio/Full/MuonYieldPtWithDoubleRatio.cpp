/*
Mini macro to plot pT(\mu+) and pT(\mu-) as asked by Cesar on the Z boson analysis gDoc.
Double ratio.
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
#include <tuple>
#include <utility>
#include <TROOT.h>
#include <TPad.h>
#include <TLine.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
#include <TVector2.h>
#include <algorithm>
#include <TF1.h>
#include "../headers/basicFormatting.h"


//---Macro settings
std::string BasePath = "/home/lucas/Documents/CMS/z_boson_analysis/"; //IFT
//std::string BasePath = "/home/lucasdoriac/z_boson_analysis/data/"; //Home

std::string plot_extension = ".pdf"; // ".png" for regular development and ".pdf" for final quality plots
std::string whichDataset = "PbPb_Run3_Data"; // "PbPb2023_2024_Data", "PbPb2023_Data", "PbPb2024_Data".
std::string JointPbPb = "PbPb2023-2026"; //"PbPb2023+2024", "PbPb2023", "PbPb2024".
std::string dataSamplesUsed = "PbPb 2023-2026, ppRef 2024 (5.36 TeV)"; //"PbPb 2023+2024, ppRef 2024 (5.36 TeV)", "PbPb 2023, ppRef 2024 (5.36 TeV)", "PbPb 2024, ppRef 2024 (5.36 TeV)".

double delta = 1e-6;
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


std::vector<std::pair<double, double>> CentralityBinsSet = {
    {0., 10.},
    {10., 20.},
    {20., 30.},
    {30., 100.},
    {0., 100.}
};

/*std::vector<std::pair<double, double>> CentralityBinsSet = {
    {0., 10.},
    {10., 30.},
    {30., 50.},
    {50., 100.},
    {0., 100.}
};*/


//ppRef histograms
TH1D* h1D_PtMuPl_ppRef = nullptr;
TH1D* h1D_PtMuMi_ppRef = nullptr;


//---Function declarations
void FillPtHistograms(const Dataset& dataset, float lowCent, float highCent, TH1D* h1D_PtMuPl, TH1D* h1D_PtMuMi);
void PlotPtDistributionsWithRatio(double lowCent, double highCent, TH1D* h_PtMuPl, TH1D* h_PtMuMi, 
                                    TH1D* h_ppPlus, TH1D* h_ppMinus, CollisionSystem system);


//---Main()
void MuonYieldPtWithDoubleRatio(){

    gROOT->SetBatch(kTRUE);

    //Plot the ppRef just for reference.
    h1D_PtMuPl_ppRef = new TH1D("h1D_PtMuPl_ppRef", "Muon Plus pT; pT [GeV]; Entries", 200, 0., 200.);
    h1D_PtMuMi_ppRef = new TH1D("h1D_PtMuMi_ppRef", "Muon Minus pT; pT [GeV]; Entries", 200, 0., 200.);
    FillPtHistograms(datasets[0], 0., 100., h1D_PtMuPl_ppRef, h1D_PtMuMi_ppRef); //ppRef

    //For PbPb datasets loop over centrality bins. 
    for(const auto& cBin : CentralityBinsSet){

        //One histogram for each distribution: dN/dpT(mu+) and dN/dpT(mu-)
        TH1D* h1D_PtMuPl = new TH1D("h1D_PtMuPl", "Muon Plus pT; pT [GeV]; Entries", 200, 0., 200.);
        TH1D* h1D_PtMuMi = new TH1D("h1D_PtMuMi", "Muon Minus pT; pT [GeV]; Entries", 200, 0., 200.);

        FillPtHistograms(datasets[1], cBin.first, cBin.second, h1D_PtMuPl, h1D_PtMuMi); //PbPb2023 dataset
        FillPtHistograms(datasets[2], cBin.first, cBin.second, h1D_PtMuPl, h1D_PtMuMi); //PbPb2024 dataset
        FillPtHistograms(datasets[3], cBin.first, cBin.second, h1D_PtMuPl, h1D_PtMuMi); //PbPb2025 dataset
        FillPtHistograms(datasets[4], cBin.first, cBin.second, h1D_PtMuPl, h1D_PtMuMi); //PbPb2026 dataset

        //Check the n of dimuons for this centrality bin
        std::cout << "\n> Number of dimuons in centrality bin " 
        << cBin.first << "-" << cBin.second << "%: " << h1D_PtMuPl->GetEntries() << std::endl;

        TH1D* ppPlusCopy  = static_cast<TH1D*>(h1D_PtMuPl_ppRef->Clone("ppPlusCopy"));
        TH1D* ppMinusCopy = static_cast<TH1D*>(h1D_PtMuMi_ppRef->Clone("ppMinusCopy"));
        ppPlusCopy->SetDirectory(nullptr);
        ppMinusCopy->SetDirectory(nullptr);

        //At this point the histograms have all candidates in each centrality bin.
        PlotPtDistributionsWithRatio(cBin.first, cBin.second, h1D_PtMuPl, h1D_PtMuMi, ppPlusCopy, ppMinusCopy, CollisionSystem::PbPb);

        delete h1D_PtMuPl;
        delete h1D_PtMuMi;
        delete ppPlusCopy;
        delete ppMinusCopy;
    }

    delete h1D_PtMuPl_ppRef;
    delete h1D_PtMuMi_ppRef;
}


//---Function definitions
void PlotPtDistributionsWithRatio(double lowCent, double highCent, TH1D* h_PtMuPl, TH1D* h_PtMuMi, TH1D* h_ppPlus, TH1D* h_ppMinus, CollisionSystem system){

    //Calculate Z count and its error for the given centrality range.
    double Zcount = 0.;
    double ZcountError = 0.;
    Zcount = h_PtMuPl->IntegralAndError(1, h_PtMuPl->GetNbinsX(), ZcountError);

    //Create double ratio histogram before normalization
    TH1D* histRatio = new TH1D("histRatio", "histRatio", h_PtMuPl->GetNbinsX(), h_PtMuPl->GetXaxis()->GetXmin(), h_PtMuPl->GetXaxis()->GetXmax());
    //histRatio->Divide(h_PtMuPl, h_PtMuMi, 1.0, 1.0, "B");//Old histRatio definition. Binomial error propagation.

    int ptfirstbin = histRatio->FindBin(ptCutValue + delta);
    //Calculate the ratio and its error for each bin, taking into account the COMPLETE correlation between the two histograms.
    for(int i = ptfirstbin; i <= 100; ++i){

        double Nplus_PbPb  = h_PtMuPl->GetBinContent(i);
        double Nminus_PbPb = h_PtMuMi->GetBinContent(i);
        double sigmaPlus_PbPb  = h_PtMuPl->GetBinError(i);
        double sigmaMinus_PbPb = h_PtMuMi->GetBinError(i);

        double Nplus_ppRef  = h_ppPlus->GetBinContent(i);
        double Nminus_ppRef = h_ppMinus->GetBinContent(i);
        double sigmaPlus_ppRef  = h_ppPlus->GetBinError(i);
        double sigmaMinus_ppRef = h_ppMinus->GetBinError(i);

        if(Nminus_PbPb <= 0.0 || Nplus_PbPb <= 0.0 || Nminus_ppRef <= 0.0 || Nplus_ppRef <= 0.0){
            std::cout << "Attention: Bin content in pT bin " << i << " is zero or negative." << std::endl;
            histRatio->SetBinContent(i, 0.);
            histRatio->SetBinError(i, 0.);
            continue;
        }

        double ratio_PbPb = Nplus_PbPb / Nminus_PbPb;
        double ratio_ppRef = Nplus_ppRef / Nminus_ppRef;
        if(ratio_ppRef <= 0.0){
            std::cout << "Attention: Ratio in pT bin " << i << " is zero or negative." << std::endl;
            histRatio->SetBinContent(i, 0.);
            histRatio->SetBinError(i, 0.);
            continue;
        }
        double DoubleRatio = ratio_PbPb / ratio_ppRef;

        //Statistical uncertainty assuming COMPLETE correlation between the two histograms. From Lara's paper.
        double ratio_PbPb_Error = ratio_PbPb * std::sqrt(sigmaPlus_PbPb*sigmaPlus_PbPb/(Nplus_PbPb*Nplus_PbPb)
            + sigmaMinus_PbPb*sigmaMinus_PbPb/(Nminus_PbPb*Nminus_PbPb)
            - 2.*rho*sigmaPlus_PbPb*sigmaMinus_PbPb/(Nplus_PbPb*Nminus_PbPb));

        double ratio_ppRef_Error = ratio_ppRef * std::sqrt(sigmaPlus_ppRef*sigmaPlus_ppRef/(Nplus_ppRef*Nplus_ppRef)
            + sigmaMinus_ppRef*sigmaMinus_ppRef/(Nminus_ppRef*Nminus_ppRef)
            - 2.*rho*sigmaPlus_ppRef*sigmaMinus_ppRef/(Nplus_ppRef*Nminus_ppRef));

        //Add them in quadrature to get the error of the double ratio.
        double DoubleRatioError = DoubleRatio * std::hypot(ratio_PbPb_Error / ratio_PbPb,ratio_ppRef_Error / ratio_ppRef);


        histRatio->SetBinContent(i, DoubleRatio);
        histRatio->SetBinError(i, DoubleRatioError);
    }

    //Create canvas for plotting
    TCanvas *c = new TCanvas("c", "Muon pT Distributions with Ratio", 800, 800);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
    basicPaddedCanvasFormatting(c, pad1, pad2);
    pad1->Draw();
    pad2->Draw();

    //Top pad: pT distributions of mu+ and mu- from both PbPb and ppRef datasets.
    pad1->cd();
    basicPaddedHistFormatting(h_PtMuPl, false);//The false means its not the ratio.
    basicPaddedHistFormatting(h_PtMuMi, false);
    basicPaddedHistFormatting(h_ppPlus, false);
    basicPaddedHistFormatting(h_ppMinus, false);

    //Normalize histograms to unity for comparison
    h_PtMuPl->Scale(1.0 / h_PtMuPl->Integral());
    h_PtMuMi->Scale(1.0 / h_PtMuMi->Integral());
    h_ppPlus->Scale(1.0 / h_ppPlus->Integral());
    h_ppMinus->Scale(1.0 / h_ppMinus->Integral());


    h_ppPlus->SetLineColorAlpha(kRed-9, 0.8);
    h_ppPlus->SetMarkerStyle(20);
    h_ppPlus->SetMarkerSize(0.65);
    h_ppPlus->SetMarkerColorAlpha(kRed, 0.8);

    h_ppMinus->SetLineColorAlpha(kBlue-9, 0.8);
    h_ppMinus->SetMarkerStyle(20);
    h_ppMinus->SetMarkerSize(0.65);
    h_ppMinus->SetMarkerColorAlpha(kBlue, 0.8);

    h_ppPlus->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    h_ppPlus->GetYaxis()->SetTitle("N of muons [GeV/c]^{-1}");
    h_ppPlus->GetXaxis()->SetRangeUser(18., 100.);
    h_ppPlus->GetYaxis()->SetTitleOffset(1.);

    
    //Filling histogram. No border.
    auto* fillPl = static_cast<TH1*>(h_PtMuPl->Clone("h_fill"));
    fillPl->SetDirectory(nullptr);
    fillPl->SetFillStyle(1001);
    fillPl->SetFillColorAlpha(kRed-7, 0.55);
    fillPl->SetLineColorAlpha(kRed-7, 0.0);

    auto* fillMi = static_cast<TH1*>(h_PtMuMi->Clone("h_fill"));
    fillMi->SetDirectory(nullptr);
    fillMi->SetFillStyle(1001);
    fillMi->SetFillColorAlpha(kBlue-9, 0.5);
    fillMi->SetLineColorAlpha(kBlue-9, 0.0);


    //Draw
    h_ppPlus->Draw("E1");
    h_ppMinus->Draw("E1 SAME");
    fillPl->Draw("HIST ][ SAME");
    fillMi->Draw("HIST ][ SAME");


    //Add Z count and error on the plot
    drawLatexText(Form("Z count: %.0f #pm %.0f", Zcount, ZcountError), 0.7, 0.35, 0.03);

    //Selections and cuts
    drawLatexText(Form("p_{T} > %.0f GeV, |#eta| < %.1f", ptCutValue, EtaCutValue), 0.7, 0.3, 0.03);
    drawLatexText(Form("|y| < %.1f", RapidityCutValue), 0.7, 0.25, 0.03);
    drawLatexText(Form("%.0f < M_{#mu#mu} < %.0f GeV", minZ_Mass, maxZ_Mass), 0.7, 0.2, 0.03);
    if(system == CollisionSystem::PbPb){//Only draw centrality string for PbPb collision system.
        drawLatexText(Form("Centrality: %.0f - %.0f %%", lowCent, highCent), 0.7, 0.15, 0.03);
    }

    TLegend* leg = new TLegend(0.75, 0.65, 0.97, 0.85);
    basicLegendFormatting(leg);
    leg->AddEntry(fillPl, "p_{T}(#mu^{+})", "f");
    leg->AddEntry(fillMi, "p_{T}(#mu^{-})", "f");
    leg->AddEntry(h_ppPlus, "p_{T}(#mu^{+}) pp", "lep");
    leg->AddEntry(h_ppMinus, "p_{T}(#mu^{-}) pp", "lep");
    leg->Draw();
    pad1->Update();

    //
    //Bottom pad. Double ratio of pT distributions of mu+ and mu-.
    pad2->cd();
    basicPaddedHistFormatting(histRatio, true); //true=is ratio histogram.
    
    histRatio->SetMarkerStyle(24);
    histRatio->SetMarkerSize(0.8);
    histRatio->SetMarkerColor(kBlack);
    histRatio->SetLineColor(kBlack);
    histRatio->SetLineWidth(1);

    histRatio->GetYaxis()->SetTitle("Double Ratio");
    histRatio->GetXaxis()->SetTitle("p_{T} [GeV]");
    histRatio->GetXaxis()->SetRangeUser(18., 100.);
    histRatio->GetYaxis()->SetRangeUser(0., 2.);

    histRatio->Draw("E1");

    //A horizontal line to represent the null hypothesis of no difference between mu+ and mu- pT distributions.
    TLine *line = new TLine(18., 1.0, 100., 1.0);
    line->SetLineColor(kMagenta+2);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");

    TLegend *leg2 = new TLegend(0.2, 0.78, 0.32, 0.96);
    basicLegendFormatting(leg2);
    leg2->SetTextSize(0.058);
    leg2->AddEntry(line, "R = 1 (null hypothesis)", "l");
    leg2->Draw();
    pad2->Update();

    //Perform a Chi2 test to see if the double ratio is compatible with 1.
    //auto [ndf, chi2, pValue] = MakeChi2Test(histRatio);

    //Back to the canvas.
    c->cd();
    drawLatexText("#bf{CMS}", 0.11, 0.95, 0.04);
    drawLatexText("#it{Internal}", 0.2, 0.95, 0.03);
    drawLatexText(dataSamplesUsed.c_str(), 0.5, 0.95, 0.026);
    //drawLatexText(Form("#chi^{2}/ndf = %.2f/%d", chi2, ndf), 0.67, 0.72, 0.022);
    //drawLatexText(Form("p-value = %.2f", pValue), 0.67, 0.69, 0.022);

    //Save
    std::string centString = Form("_Cent%.0f-%.0f", lowCent, highCent);
    std::string output = "MuonPtMuPlMuMiHistWithDoubleRatio" + centString + plot_extension;
    c->Update();
    c->SaveAs(output.c_str());


    delete fillPl;
    delete fillMi;
    delete leg;
    delete leg2;
    delete histRatio;
    delete line;
    delete c;
}


void FillPtHistograms(const Dataset& dataset, float lowCent, float highCent, TH1D* h1D_PtMuPl, TH1D* h1D_PtMuMi){

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
            h1D_PtMuPl->Fill(ptplus);
            h1D_PtMuMi->Fill(ptminus);

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