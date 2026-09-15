/*
Mini macro to plot pT(\mu+) and pT(\mu-) as asked by Cesar on the Z boson analysis gDoc.
Bottom pad is difference.
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

void MuonPtMuPlMuMiHistWithSingleDiff(TFile* inputFile, const Dataset& dataset, double lowCent = 0., double highCent = 100.);

void MuonYieldPtWithDifference(){

    gROOT->SetBatch(kTRUE);
    TFile* inputFile = new TFile("mySelectedData.root", "READ");

    for(const auto& cBin : CentralityBinsSet){
        double lowCent = cBin.first;
        double highCent = cBin.second;

        std::cout << "> Processing centrality bin: " << lowCent << " - " << highCent << std::endl;

        //centrality bins are defined and applied only for PbPb2024.
        MuonPtMuPlMuMiHistWithSingleDiff(inputFile, datasets[1], lowCent, highCent); //PbPb2024
    }
    
    MuonPtMuPlMuMiHistWithSingleDiff(inputFile, datasets[1], 0., 100.); //PbPb2024
    MuonPtMuPlMuMiHistWithSingleDiff(inputFile, datasets[2], 0., 100.); //ppRef2024
}

void MuonPtMuPlMuMiHistWithSingleDiff(TFile* inputFile, const Dataset& dataset, double lowCent, double highCent){

    //Identify the dataset and get the corresponding directory.
    std::string dirName = dataset.name;

    //Get directory
    TDirectory *dir = inputFile->GetDirectory(dirName.c_str());

    //Get histograms
    TH1D* h_ogpl = nullptr;
    TH1D* h_ogmi = nullptr;

    if(dataset.system == CollisionSystem::PbPb2024){
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
    }

    else if(dataset.system == CollisionSystem::ppRef2024){
        //For ppRef2024, we have individual histograms for mu+ and mu-.
        h_ogpl = dynamic_cast<TH1D*>(dir->Get("h1D_ptMuPlus"));
        h_ogmi = dynamic_cast<TH1D*>(dir->Get("h1D_ptMuMinus"));
    }

    //Clone histograms
    TH1D* h_PtMuPl = dynamic_cast<TH1D*>(h_ogpl->Clone("h_PtMuPl"));
    TH1D* h_PtMuMi = dynamic_cast<TH1D*>(h_ogmi->Clone("h_PtMuMi"));

    //PbPb2024 canvas
    TCanvas *c = new TCanvas("c", "c", 800, 800);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.30, 1, 1.0);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.00, 1, 0.30);
    basicPaddedCanvasFormatting(c, pad1, pad2);
    pad1->Draw();
    pad2->Draw();

    //Top pad. pT distributions of mu+ and mu-.
    pad1->cd();
    basicPaddedHistFormatting(h_PtMuPl, false);
    basicPaddedHistFormatting(h_PtMuMi, false);

    h_PtMuPl->SetFillStyle(0);
    h_PtMuPl->SetLineWidth(1);
    h_PtMuPl->SetLineColorAlpha(kRed+1, 1.);
    h_PtMuPl->SetMarkerStyle(20);
    h_PtMuPl->SetMarkerSize(0.75);
    h_PtMuPl->SetMarkerColorAlpha(kRed+1, 1.);

    h_PtMuMi->SetFillStyle(0);
    h_PtMuMi->SetLineWidth(1);
    h_PtMuMi->SetLineColorAlpha(kBlue+1, 1.);
    h_PtMuMi->SetMarkerStyle(20);
    h_PtMuMi->SetMarkerSize(0.75);
    h_PtMuMi->SetMarkerColorAlpha(kBlue+1, 1.);

    h_PtMuPl->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    h_PtMuPl->GetYaxis()->SetTitle("N of muons [GeV/c]^{-1}");
    h_PtMuPl->GetXaxis()->SetRangeUser(18., 100.);
    h_PtMuPl->GetYaxis()->SetTitleOffset(1.);

    h_PtMuPl->Draw("E1");
    h_PtMuMi->Draw("E1 SAME");

    //Filling histogram. No border.
    auto* fillPl = static_cast<TH1*>(h_PtMuPl->Clone("h_fill"));
    fillPl->SetDirectory(nullptr);
    fillPl->SetFillStyle(1001);
    fillPl->SetFillColorAlpha(kRed-10, 0.5);
    fillPl->SetLineColorAlpha(kRed-10, 0.0);
    fillPl->Draw("HIST ][ SAME");

    auto* fillMi = static_cast<TH1*>(h_PtMuMi->Clone("h_fill"));
    fillMi->SetDirectory(nullptr);
    fillMi->SetFillStyle(1001);
    fillMi->SetFillColorAlpha(kBlue-10, 0.5);
    fillMi->SetLineColorAlpha(kBlue-10, 0.0);
    fillMi->Draw("HIST ][ SAME");

    //Calculate Z count and its error for the given centrality range.
    double Zcount = 0.0;
    double ZcountError = 0.0;
    for(int i = 1; i <= h_PtMuPl->GetNbinsX(); ++i){
        Zcount += h_PtMuPl->GetBinContent(i);
        ZcountError += std::pow(h_PtMuPl->GetBinError(i), 2);
    }
    //ZcountError = std::sqrt(ZcountError);
    Zcount = h_PtMuPl->IntegralAndError(1, h_PtMuPl->GetNbinsX(), ZcountError);
    //Add Z count and error on the plot
    drawLatexText(Form("Z count: %.0f #pm %.0f", Zcount, ZcountError), 0.7, 0.35, 0.03);

    //Selections and cuts
    drawLatexText(Form("p_{T} > %.0f GeV, |#eta| < %.1f", PTCUTVALUE, ETACUTVALUE), 0.7, 0.3, 0.03);
    drawLatexText(Form("|y| < %.1f", RAPIDITYCUTVALUE), 0.7, 0.25, 0.03);
    drawLatexText(Form("%.0f < M_{#mu#mu} < %.0f GeV", MINZ_MASS, MAXZ_MASS), 0.7, 0.2, 0.03);
    drawLatexText(Form("Centrality: %.0f - %.0f %%", lowCent, highCent), 0.7, 0.15, 0.03);

    TLegend* leg = new TLegend(0.75, 0.75, 0.94, 0.86);
    basicLegendFormatting(leg);
    leg->AddEntry(h_PtMuPl, "p_{T}(#mu^{+})", "l");
    leg->AddEntry(h_PtMuMi, "p_{T}(#mu^{-})", "l");
    leg->Draw();
    pad1->Update();

    //Bottom pad. Difference of pT distributions of mu+ and mu-.
    pad2->cd();

    TH1D* histDiff = new TH1D("histDiff","histDiff",h_PtMuPl->GetNbinsX(),h_PtMuPl->GetXaxis()->GetXmin(), 
                                h_PtMuPl->GetXaxis()->GetXmax());
    histDiff->SetDirectory(nullptr);

    int ptfirstbin = histDiff->FindBin(PTCUTVALUE + delta);
    //Calculate the ratio and its error for each bin, taking into account the COMPLETE correlation between the two histograms.
    for(int i = ptfirstbin; i <= h_PtMuPl->GetNbinsX(); ++i){

        double contentPl = h_PtMuPl->GetBinContent(i);
        double contentMi = h_PtMuMi->GetBinContent(i);
        double errorPl = h_PtMuPl->GetBinError(i);
        double errorMi = h_PtMuMi->GetBinError(i);

        //Difference
        double diff = contentPl - contentMi;
        histDiff->SetBinContent(i, diff);

        //Error propagation for difference, considering complete correlation.
        double diffError = std::abs(errorPl - errorMi);
        histDiff->SetBinError(i, diffError);
    }

    basicPaddedHistFormatting(histDiff, true);
    
    histDiff->SetMarkerStyle(24);
    histDiff->SetMarkerSize(0.8);
    histDiff->SetMarkerColor(kBlack);
    histDiff->SetLineColor(kBlack);

    histDiff->GetYaxis()->SetTitle("Difference");
    histDiff->GetXaxis()->SetTitle("p_{T} [GeV]");
    histDiff->GetXaxis()->SetRangeUser(18., 100.);

    histDiff->Draw("P");

    //A horizontal line to represent the null hypothesis of no difference between mu+ and mu- pT distributions.
    TLine *line = new TLine(18.,0.0,histDiff->GetXaxis()->GetXmax(),0.0);
    line->SetLineColor(kMagenta+2);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");

    if(dataset.system == CollisionSystem::PbPb2024){
        TLegend *leg2 = new TLegend(0.72, 0.78, 0.84, 0.96);
        basicLegendFormatting(leg2);
        leg2->SetTextSize(0.058);
        leg2->AddEntry(line, "A = 0 (null hypothesis)", "l");
        leg2->Draw();
        pad2->Update();
    }
    else if(dataset.system == CollisionSystem::ppRef2024){
        TLegend *leg2 = new TLegend(0.72, 0.78, 0.84, 0.96);
        basicLegendFormatting(leg2);
        leg2->SetTextSize(0.058);
        leg2->AddEntry(line, "A = 0 (null hypothesis)", "l");
        leg2->Draw();
        pad2->Update();
    }

    //Back to the canvas.
    c->cd();
    drawLatexText("#bf{CMS}", 0.11, 0.95, 0.04);
    drawLatexText("#it{Work in Progress}", 0.2, 0.95, 0.026);
    if(dataset.system == CollisionSystem::PbPb2024) drawLatexText("PbPb 2024 (5.36 TeV)", 0.7, 0.95, 0.026);
    else if(dataset.system == CollisionSystem::ppRef2024) drawLatexText("ppRef 2024 (5.36 TeV)", 0.7, 0.95, 0.026);

    //Save
    std::string centString = Form("_Cent%.0f-%.0f", lowCent, highCent);
    std::string output = dataset.name + "_MuonPtMuPlMuMiHistWithSingleDiff" + centString + plot_extension;
    c->Update();
    c->SaveAs(output.c_str());

    delete h_PtMuPl;
    delete h_PtMuMi;
    delete histDiff;
    delete line;
    delete c;
}