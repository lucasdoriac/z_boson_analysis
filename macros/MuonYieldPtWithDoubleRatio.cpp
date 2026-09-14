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


//Location of datasets
//std::string BasePath = "/home/lucas/Documents/CMS/z_boson_analysis/"; //IFT
std::string BasePath = "/home/lucasdoriac/z_boson_analysis/data/"; //Home


//---Macro settings
std::string plot_extension = ".pdf"; // ".png" for regular development and ".pdf" for final quality plots
double delta = 1e-6; //Small value to avoid binning issues when projecting histograms.

double MINZ_MASS = 60.;
double MAXZ_MASS = 120.;
double RAPIDITYCUTVALUE = 2.4;
double ETACUTVALUE = 2.4;
double PTCUTVALUE = 20.;

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

//Set of centrality bins for PbPb2024 data. We can decide to change the centrality bins later if we want to.
std::vector<std::pair<double, double>> CentralityBinsSet = {
    {0., 10.},
    {10., 20.},
    {20., 30.},
    {30., 100.}
};

//Second proposed set of centrality bins for PbPb2024 data.
/*
std::vector<std::pair<double, double>> CentralityBinsSet = {
    {0., 10.},
    {10., 30.},
    {30., 50.},
    {50., 100.}
};
//*/

void MuonPtMuPlMuMiHistWithDoubleRatio(TFile* inputFile, const Dataset& dataset, const Dataset& ref_dataset, double lowCent = 0., double highCent = 100.);

void MuonYieldPtWithDoubleRatio(){

    gROOT->SetBatch(kTRUE);
    TFile* inputFile = new TFile("mySelectedData.root", "READ");

    for(const auto& cBin : CentralityBinsSet){
        double lowCent = cBin.first;
        double highCent = cBin.second;
        std::cout << "> Processing centrality bin: " << lowCent << " - " << highCent << std::endl;

        MuonPtMuPlMuMiHistWithDoubleRatio(inputFile, datasets[1], datasets[2], lowCent, highCent); //PbPb2024 and ppRef2024.
    }
    
    MuonPtMuPlMuMiHistWithDoubleRatio(inputFile, datasets[1], datasets[2], 0., 100.); //PbPb2024 and ppRef2024, full centrality range.
}

void MuonPtMuPlMuMiHistWithDoubleRatio(TFile* inputFile, const Dataset& dataset, const Dataset& ref_dataset, double lowCent, double highCent){

    //Identify the dataset and get the corresponding directory.
    std::string dirName = dataset.name;
    std::string ref_dir = ref_dataset.name;

    //Get directory
    TDirectory *dir = inputFile->GetDirectory(dirName.c_str());
    TDirectory *refDir = inputFile->GetDirectory(ref_dir.c_str());

    //Get histograms
    TH1D* h_ogpl = nullptr;
    TH1D* h_ogmi = nullptr;
    TH1D* h_refpl = nullptr;
    TH1D* h_refmi = nullptr;

    //We need to get the 3D histogram and project it onto the pT(mu+) vs pT(mu-) plane for the given centrality range.
    TH3D* h3D_PtMuPl_PtMuMi_Cent = dynamic_cast<TH3D*>(dir->Get("h3D_PtMuPl_PtMuMi_Cent"));
    TH3D* h3D_PtMuPl_PtMuMi_Cent_clone = dynamic_cast<TH3D*>(h3D_PtMuPl_PtMuMi_Cent->Clone("h3D_PtMuPl_PtMuMi_Cent_clone"));
    
    /*
    z-axis -> centrality
    y-axis -> pT of mu-
    x-axis -> pT of mu+
    */

    //Select centrality range and project onto pT(mu+) vs pT(mu-) plane.
    int binLow = h3D_PtMuPl_PtMuMi_Cent_clone->GetZaxis()->FindBin(lowCent + delta);
    int binHigh = h3D_PtMuPl_PtMuMi_Cent_clone->GetZaxis()->FindBin(highCent - delta);
    std::cout << "> Centrality range: " << lowCent << " - " << highCent << std::endl;
    h3D_PtMuPl_PtMuMi_Cent_clone->GetZaxis()->SetRange(binLow, binHigh);
    TH2D* h2D_PtMuPl_PtMuMi = dynamic_cast<TH2D*>(h3D_PtMuPl_PtMuMi_Cent_clone->Project3D("yx"));
    h_ogpl = h2D_PtMuPl_PtMuMi->ProjectionX("h_ogpl",1,h2D_PtMuPl_PtMuMi->GetNbinsY(),"e");
    h_ogmi = h2D_PtMuPl_PtMuMi->ProjectionY("h_ogmi",1,h2D_PtMuPl_PtMuMi->GetNbinsX(),"e");

    //Reference histograms.
    h_refpl = dynamic_cast<TH1D*>(refDir->Get("h1D_ptMuPlus"));
    h_refmi = dynamic_cast<TH1D*>(refDir->Get("h1D_ptMuMinus"));

    //Clone histograms
    TH1D* h_PtMuPl = dynamic_cast<TH1D*>(h_ogpl->Clone("h_PtMuPl"));
    TH1D* h_PtMuMi = dynamic_cast<TH1D*>(h_ogmi->Clone("h_PtMuMi"));
    TH1D* h_RefPl = dynamic_cast<TH1D*>(h_refpl->Clone("h_RefPl"));
    TH1D* h_RefMi = dynamic_cast<TH1D*>(h_refmi->Clone("h_RefMi"));

    //
    TCanvas *c = new TCanvas("c", "c", 800, 800);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.4, 1, 1.0);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.4);
    basicPaddedCanvasFormatting(c, pad1, pad2);
    pad1->Draw();
    pad2->Draw();

    //Top pad. pT distributions of mu+ and mu- for both datasets.
    pad1->cd();
    basicPaddedHistFormatting(h_PtMuPl, false);
    basicPaddedHistFormatting(h_PtMuMi, false);
    basicPaddedHistFormatting(h_RefPl, false);
    basicPaddedHistFormatting(h_RefMi, false);

    //Normalize histograms to unit area for comparison.
    h_PtMuPl->Scale(1.0 / h_PtMuPl->Integral());
    h_PtMuMi->Scale(1.0 / h_PtMuMi->Integral());
    h_RefPl->Scale(1.0 / h_RefPl->Integral());
    h_RefMi->Scale(1.0 / h_RefMi->Integral());

    //Reference histograms with points
    h_RefPl->SetFillStyle(0);
    h_RefPl->SetMarkerStyle(22);
    h_RefPl->SetMarkerSize(0.75);
    h_RefPl->SetMarkerColorAlpha(kRed, 1.);

    h_RefMi->SetFillStyle(0);
    h_RefMi->SetMarkerStyle(24);
    h_RefMi->SetMarkerSize(0.75);
    h_RefMi->SetMarkerColorAlpha(kBlue, 1.);

    h_RefPl->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    h_RefPl->GetYaxis()->SetTitle("Normalized Yield");
    h_RefPl->GetXaxis()->SetRangeUser(18., 100.);
    h_RefPl->GetYaxis()->SetTitleOffset(1.);

    h_RefPl->Draw("HIST P");
    h_RefMi->Draw("HIST P SAME");

    //Filling histogram with PbPb values. No border.
    auto* fillPl = static_cast<TH1*>(h_PtMuPl->Clone("h_fill"));
    fillPl->SetDirectory(nullptr);
    fillPl->SetFillStyle(1001);
    fillPl->SetFillColorAlpha(kRed-10, 0.6);
    fillPl->SetLineColorAlpha(kRed-10, 0.0);
    fillPl->Draw("HIST ][ SAME");

    auto* fillMi = static_cast<TH1*>(h_PtMuMi->Clone("h_fill"));
    fillMi->SetDirectory(nullptr);
    fillMi->SetFillStyle(1001);
    fillMi->SetFillColorAlpha(kBlue-10, 0.6);
    fillMi->SetLineColorAlpha(kBlue-10, 0.0);
    fillMi->Draw("HIST ][ SAME");

    //Selections and cuts
    drawLatexText(Form("p_{T} > %.0f GeV, |#eta| < %.1f", PTCUTVALUE, ETACUTVALUE), 0.7, 0.3, 0.03);
    drawLatexText(Form("|y| < %.1f", RAPIDITYCUTVALUE), 0.7, 0.25, 0.03);
    drawLatexText(Form("%.0f < M_{#mu#mu} < %.0f GeV", MINZ_MASS, MAXZ_MASS), 0.7, 0.2, 0.03);
    drawLatexText(Form("Centrality: %.0f - %.0f %%", lowCent, highCent), 0.7, 0.15, 0.03);

    TLegend* leg = new TLegend(0.72, 0.7, 0.95, 0.85);
    basicLegendFormatting(leg);
    leg->AddEntry(fillPl, "p_{T}(#mu^{+})", "f");
    leg->AddEntry(fillMi, "p_{T}(#mu^{-})", "f");
    leg->AddEntry(h_RefPl, "p_{T}(#mu^{+}) ppRef", "p");
    leg->AddEntry(h_RefMi, "p_{T}(#mu^{-}) ppRef", "p");
    leg->Draw();
    pad1->Update();

    //Bottom pad. Ratio of pT distributions of mu+ and mu-.
    pad2->cd();

    TH1D* histRatio = new TH1D("histRatio", "histRatio", h_PtMuPl->GetNbinsX(), h_PtMuPl->GetXaxis()->GetXmin(), h_PtMuPl->GetXaxis()->GetXmax());

    int ptfirstbin = histRatio->FindBin(PTCUTVALUE + delta);
    //Calculate the ratio and its error for each bin, taking into account the COMPLETE correlation between the two histograms.
    for(int i = ptfirstbin; i <= h_PtMuPl->GetNbinsX(); ++i){

        //Calculate reference ratio
        double NplusRef = h_RefPl->GetBinContent(i);
        double NminusRef = h_RefMi->GetBinContent(i);

        double sigmaPlusRef  = h_RefPl->GetBinError(i);
        double sigmaMinusRef = h_RefMi->GetBinError(i);
        if(NminusRef <= 0.0 || NplusRef <= 0.0){
            //std::cout << "Warning: zero content in reference bin " << i << std::endl;
            histRatio->SetBinContent(i, 0.0);
            histRatio->SetBinError(i, 0.0);
            continue;
        }

        double refRatio = NplusRef / NminusRef;
        double refRatioError_completeCorr = refRatio * std::abs(sigmaPlusRef/NplusRef - sigmaMinusRef/NminusRef);

        //PbPb ratio now
        double NplusPbPb  = h_PtMuPl->GetBinContent(i);
        double NminusPbPb = h_PtMuMi->GetBinContent(i);

        double sigmaPlusPbPb  = h_PtMuPl->GetBinError(i);
        double sigmaMinusPbPb = h_PtMuMi->GetBinError(i);
        if(NminusPbPb <= 0.0 || NplusPbPb <= 0.0){
            //std::cout << "Warning: zero content in PbPb bin " << i << std::endl;
            histRatio->SetBinContent(i, 0.0);
            histRatio->SetBinError(i, 0.0);
            continue;
        }

        double PbPbratio = NplusPbPb / NminusPbPb;
        double PbPbratioError_completeCorr = PbPbratio * std::abs(sigmaPlusPbPb/NplusPbPb - sigmaMinusPbPb/NminusPbPb);

        double doubleRatio = PbPbratio / refRatio;
        double doubleRatioIndError = doubleRatio * std::sqrt(
            std::pow(PbPbratioError_completeCorr/PbPbratio, 2) + std::pow(refRatioError_completeCorr/refRatio, 2));

        histRatio->SetBinContent(i, doubleRatio);
        histRatio->SetBinError(i, doubleRatioIndError);
    }

    basicPaddedHistFormatting(histRatio, true);
    
    histRatio->SetMarkerStyle(24);
    histRatio->SetMarkerSize(0.8);
    histRatio->SetMarkerColor(kBlack);
    histRatio->SetLineColor(kBlack);

    histRatio->GetYaxis()->SetTitle("R_{PbPb}/R_{ppRef}");
    histRatio->GetXaxis()->SetTitle("p_{T} [GeV]");
    histRatio->GetXaxis()->SetRangeUser(18., 100.);
    histRatio->GetYaxis()->SetRangeUser(0., 2.);

    //Because i made the pad 2 larger i have to re-configure the histRatio.
    histRatio->GetXaxis()->SetTitleOffset(1.1);
    histRatio->GetYaxis()->SetTitleOffset(0.55);
    histRatio->GetXaxis()->SetTitleSize(0.072);
    histRatio->GetYaxis()->SetTitleSize(0.068);
    histRatio->GetXaxis()->SetLabelSize(0.055);
    histRatio->GetYaxis()->SetLabelSize(0.05);

    histRatio->Draw("P");

    //A horizontal line to represent the null hypothesis of no difference between mu+ and mu- pT distributions.
    TLine *line = new TLine(18., 1.0, histRatio->GetXaxis()->GetXmax(), 1.0);
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

    //Back to the canvas.
    c->cd();
    drawLatexText("#bf{CMS}", 0.11, 0.96, 0.035);
    drawLatexText("#it{Work in Progress}", 0.2, 0.96, 0.025);
    drawLatexText("PbPb 2024, ppRef 2024", 0.7, 0.96, 0.025);
    
    //Save
    std::string centString = Form("_Cent%.0f-%.0f", lowCent, highCent);
    std::string output = "MuonPtMuPlMuMiHistWithDoubleRatio" + centString + plot_extension;
    c->Update();
    c->SaveAs(output.c_str());

    delete h_PtMuPl;
    delete h_PtMuMi;
    delete histRatio;
    delete line;
    delete c;
}