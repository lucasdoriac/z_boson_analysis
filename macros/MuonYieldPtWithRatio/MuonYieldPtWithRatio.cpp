/*
Mini macro to plot pT(\mu+) and pT(\mu-) in top pad.
Ratio of N(\mu+)/N(\mu-) per pT bin in bottom pad.
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
std::string whichDataset = "PbPb2023_2024_Data"; // "PbPb2023_2024_Data", "PbPb2023_Data", "PbPb2024_Data".
std::string JointPbPb = "PbPb2023+2024"; //"PbPb2023+2024", "PbPb2023", "PbPb2024".
std::string dataSamplesUsed = "PbPb 2023+2024, ppRef 2024 (5.36 TeV)"; //"PbPb 2023+2024, ppRef 2024 (5.36 TeV)", "PbPb 2023, ppRef 2024 (5.36 TeV)", "PbPb 2024, ppRef 2024 (5.36 TeV)".
double delta = 1e-6; //Small value to avoid binning issues when projecting histograms.


double MINZ_MASS = 60.;
double MAXZ_MASS = 120.;
double RAPIDITYCUTVALUE = 2.4;
double ETACUTVALUE = 2.4;
double PTCUTVALUE = 20.;


// ##############################################################################
// ##############################################################################


//Set of centrality bins for PbPb2024 data.
/*std::vector<std::pair<double, double>> CentralityBinsSet = {
    {0., 10.},
    {10., 20.},
    {20., 30.},
    {30., 100.},
    {0., 100.}
};*/

//Second proposed set of centrality bins for PbPb2024 data.
std::vector<std::pair<double, double>> CentralityBinsSet = {
    {0., 10.},
    {10., 30.},
    {30., 50.},
    {50., 100.},
    {0., 100.}
};
//

void MuonPtMuPlMuMiHistWithSingleRatio(TFile* inputFile, std::string datasetName, double lowCent = 0., double highCent = 100.);
std::tuple<int, double, double> MakeChi2Test(const TH1D* hist);

void MuonYieldPtWithRatio(){

    gROOT->SetBatch(kTRUE);
    TFile* inputFile = new TFile("mySelectedData.root", "READ");

    for(const auto& cBin : CentralityBinsSet){
        double lowCent = cBin.first;
        double highCent = cBin.second;
        std::cout << "> Processing centrality bin: " << lowCent << " - " << highCent << std::endl;

        MuonPtMuPlMuMiHistWithSingleRatio(inputFile, whichDataset, lowCent, highCent); //PbPb2024
    }
    
    MuonPtMuPlMuMiHistWithSingleRatio(inputFile, "ppRef2024_Data", 0., 100.); //ppRef2024
    inputFile->Close();
}

void MuonPtMuPlMuMiHistWithSingleRatio(TFile* inputFile, std::string datasetName, double lowCent, double highCent){

    //Get directory
    TDirectory *dir = inputFile->GetDirectory(datasetName.c_str());

    //Get histograms
    TH1D* h_ogpl = nullptr;
    TH1D* h_ogmi = nullptr;

    if(datasetName != "ppRef2024_Data"){
        //For PbPb we want to separate the data into centrality bins.
        //Thus get the 3D histogram and project it onto a TH2 pT(mu+) vs pT(mu-) for the given centrality range.
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
        std::cout << "> Project in centrality range: " << lowCent << " - " << highCent << std::endl;
        h3D_PtMuPl_PtMuMi_Cent_clone->GetZaxis()->SetRange(binLow, binHigh);
        TH2D* h2D_PtMuPl_PtMuMi = dynamic_cast<TH2D*>(h3D_PtMuPl_PtMuMi_Cent_clone->Project3D("yx"));

        h_ogpl = h2D_PtMuPl_PtMuMi->ProjectionX("h_ogpl",1,h2D_PtMuPl_PtMuMi->GetNbinsY(),"e");
        h_ogmi = h2D_PtMuPl_PtMuMi->ProjectionY("h_ogmi",1,h2D_PtMuPl_PtMuMi->GetNbinsX(),"e");
    }

    else if(datasetName == "ppRef2024_Data"){
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
    h_PtMuPl->SetLineColorAlpha(kRed+1, 0.8);
    h_PtMuPl->SetMarkerStyle(20);
    h_PtMuPl->SetMarkerSize(0.6);
    h_PtMuPl->SetMarkerColorAlpha(kRed+1, 1.);

    h_PtMuMi->SetFillStyle(0);
    h_PtMuMi->SetLineWidth(1);
    h_PtMuMi->SetLineColorAlpha(kBlue+1, 0.8);
    h_PtMuMi->SetMarkerStyle(20);
    h_PtMuMi->SetMarkerSize(0.6);
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
    fillPl->SetFillColorAlpha(kRed-10, 0.6);
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

    //Bottom pad. Ratio of pT distributions of mu+ and mu-.
    pad2->cd();

    TH1D* histRatio = new TH1D("histRatio", "histRatio", h_PtMuPl->GetNbinsX(), h_PtMuPl->GetXaxis()->GetXmin(), h_PtMuPl->GetXaxis()->GetXmax());
    //histRatio->Divide(h_PtMuPl, h_PtMuMi, 1.0, 1.0, "B");//Old histRatio definition. Binomial error propagation.

    int ptfirstbin = histRatio->FindBin(PTCUTVALUE + delta);
    //Calculate the ratio and its error for each bin, taking into account the COMPLETE correlation between the two histograms.
    for(int i = ptfirstbin; i <= h_PtMuPl->GetNbinsX(); ++i){

        double Nplus  = h_PtMuPl->GetBinContent(i);
        double Nminus = h_PtMuMi->GetBinContent(i);

        double sigmaPlus  = h_PtMuPl->GetBinError(i);
        double sigmaMinus = h_PtMuMi->GetBinError(i);

        if(Nminus <= 0.0 || Nplus <= 0.0){
            std::cout << "Warning: zero content in bin " << i << std::endl;
            histRatio->SetBinContent(i, 0.0);
            histRatio->SetBinError(i, 0.0);
            continue;
        }

        double ratio = Nplus / Nminus;
        //Statistical uncertainty assuming COMPLETE correlation between the two histograms.
        double ratioError_completeCorr = ratio * std::abs(sigmaPlus/Nplus - sigmaMinus/Nminus);

        histRatio->SetBinContent(i, ratio);
        histRatio->SetBinError(i, ratioError_completeCorr);
    }

    basicPaddedHistFormatting(histRatio, true);
    
    histRatio->SetMarkerStyle(24);
    histRatio->SetMarkerSize(0.8);
    histRatio->SetMarkerColor(kBlack);
    histRatio->SetLineColor(kBlack);

    histRatio->GetYaxis()->SetTitle("Ratio");
    histRatio->GetXaxis()->SetTitle("p_{T} [GeV]");
    histRatio->GetXaxis()->SetRangeUser(18., 100.);
    histRatio->GetYaxis()->SetRangeUser(0., 2.);

    histRatio->Draw("P");

    //A horizontal line to represent the null hypothesis of no difference between mu+ and mu- pT distributions.
    TLine *line = new TLine(18., 1.0, histRatio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineColor(kMagenta+2);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");

    if(datasetName != "ppRef2024_Data"){
        TLegend *leg2 = new TLegend(0.2, 0.78, 0.32, 0.96);
        basicLegendFormatting(leg2);
        leg2->SetTextSize(0.058);
        leg2->AddEntry(line, "R = 1 (null hypothesis)", "l");
        leg2->Draw();
        pad2->Update();
    }
    else if(datasetName == "ppRef2024_Data"){
        TLegend *leg2 = new TLegend(0.15, 0.42, 0.27, 0.6);
        basicLegendFormatting(leg2);
        leg2->SetTextSize(0.058);
        leg2->AddEntry(line, "R = 1 (null hypothesis)", "l");
        leg2->Draw();
        pad2->Update();
    }

    //Perform a Chi2 test to see if the double ratio is compatible with 1.
    auto [ndf, chi2, pValue] = MakeChi2Test(histRatio);

    //Back to the canvas.
    c->cd();
    drawLatexText("#bf{CMS}", 0.11, 0.95, 0.04);
    drawLatexText("#it{Work in Progress}", 0.2, 0.95, 0.026);
    if(datasetName != "ppRef2024_Data") drawLatexText("PbPb 2023+2024 (5.36 TeV)", 0.6, 0.95, 0.026);
    else if(datasetName == "ppRef2024_Data") drawLatexText("ppRef 2024 (5.36 TeV)", 0.7, 0.95, 0.026);
    drawLatexText(Form("#chi^{2}/ndf = %.2f/%d", chi2, ndf), 0.67, 0.72, 0.022);
    drawLatexText(Form("p-value = %.2f", pValue), 0.67, 0.69, 0.022);

    //Save
    std::string centString = Form("_Cent%.0f-%.0f", lowCent, highCent);
    std::string output = datasetName + "_MuonPtMuPlMuMiHistWithSingleRatio" + centString + plot_extension;
    c->Update();
    c->SaveAs(output.c_str());

    delete h_PtMuPl;
    delete h_PtMuMi;
    delete histRatio;
    delete line;
    delete c;
}

std::tuple<int, double, double> MakeChi2Test(const TH1D* hist){

    double ptMin = 40.;
    double ptMax = 65.;

    TF1* nullHypothesis = new TF1("nullHypothesis","1.",ptMin,ptMax);

    double chi2 = hist->Chisquare(nullHypothesis, "R");

    int ndf = 0;
    //Count the number of bins with non-zero content and within the pT range of interest to determine the degrees of freedom.
    for(int i = 1; i <= hist->GetNbinsX(); ++i){
        double x = hist->GetBinCenter(i);
        if(x < ptMin || x > ptMax) continue;
        if(hist->GetBinContent(i) == 0.) continue;
        if(hist->GetBinError(i) <= 0.) continue;
        ndf++;
    }

    //Calculate the p-value from the Chi2 and ndf.
    double pValue = TMath::Prob(chi2, ndf);
    
    std::cout
    << "Null hypothesis: R = 1" << std::endl
    << "chi2 = " << chi2 << std::endl
    << "ndf = " << ndf << std::endl
    << "chi2/ndf = " << chi2 / ndf << std::endl
    << "p-value = " << pValue << std::endl;

    delete nullHypothesis;
    return {ndf, chi2, pValue};
}