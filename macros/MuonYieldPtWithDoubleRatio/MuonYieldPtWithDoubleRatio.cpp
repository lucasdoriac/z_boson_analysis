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
std::string plot_extension = ".pdf"; // ".png" for regular development and ".pdf" for final quality plots
std::string whichDataset = "PbPb2024_Data"; // "PbPb2023_2024_Data", "PbPb2023_Data", "PbPb2024_Data". 
std::string JointPbPb = "PbPb2024"; //"PbPb2023+2024", "PbPb2023", "PbPb2024".
std::string dataSamplesUsed = "PbPb 2024, ppRef 2024 (5.36 TeV)"; //"PbPb 2023+2024, ppRef 2024 (5.36 TeV)", 
                                                                       //"PbPb 2023, ppRef 2024 (5.36 TeV)", 
                                                                       //"PbPb 2024, ppRef 2024 (5.36 TeV)".
double delta = 1e-6;


double MINZ_MASS = 60.;
double MAXZ_MASS = 120.;
double RAPIDITYCUTVALUE = 2.4;
double ETACUTVALUE = 2.4;
double PTCUTVALUE = 20.;


// ##############################################################################
// ##############################################################################


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


void MuonPtMuPlMuMiHistWithDoubleRatio(TFile* inputFile, std::string whichDataset, std::string refString, double lowCent = 0., double highCent = 100.);
std::tuple<int, double, double> MakeChi2Test(const TH1D* histRatio);


void MuonYieldPtWithDoubleRatio(){

    gROOT->SetBatch(kTRUE);
    TFile* inputFile = new TFile("mySelectedData.root", "READ");

    for(const auto& cBin : CentralityBinsSet){
        double lowCent = cBin.first;
        double highCent = cBin.second;
        std::cout << "> Processing centrality bin: " << lowCent << " - " << highCent << std::endl;

        MuonPtMuPlMuMiHistWithDoubleRatio(inputFile, whichDataset, "ppRef2024_Data", lowCent, highCent); //PbPb2024 and ppRef2024.
    }
    
}

void MuonPtMuPlMuMiHistWithDoubleRatio(TFile* inputFile, std::string whichDataset, std::string refString, double lowCent, double highCent){

    std::string centString = Form("_Cent%.0f-%.0f", lowCent, highCent);
    //Identify the dataset and get the corresponding directory.
    std::string dirName = whichDataset;
    std::string ref_dir = refString;

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
    h_ogpl = h2D_PtMuPl_PtMuMi->ProjectionX(("h_ogpl" + centString).c_str(),1,h2D_PtMuPl_PtMuMi->GetNbinsY(),"e");
    h_ogmi = h2D_PtMuPl_PtMuMi->ProjectionY(("h_ogmi" + centString).c_str(),1,h2D_PtMuPl_PtMuMi->GetNbinsX(),"e");

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

    //We need to manipulate clones of the histograms since we normalize them before doing the ratio.
    //By normalizing only 'plot versions' of the histograms, we avoid that.
    //Clone histograms for plotting.
    TH1D* h_PtMuPl_plot = dynamic_cast<TH1D*>(h_PtMuPl->Clone("h_PtMuPl_plot"));
    TH1D* h_PtMuMi_plot = dynamic_cast<TH1D*>(h_PtMuMi->Clone("h_PtMuMi_plot"));
    TH1D* h_RefPl_plot = dynamic_cast<TH1D*>(h_RefPl->Clone("h_RefPl_plot"));
    TH1D* h_RefMi_plot = dynamic_cast<TH1D*>(h_RefMi->Clone("h_RefMi_plot"));

    basicPaddedHistFormatting(h_PtMuPl_plot, false);
    basicPaddedHistFormatting(h_PtMuMi_plot, false);
    basicPaddedHistFormatting(h_RefPl_plot, false);
    basicPaddedHistFormatting(h_RefMi_plot, false);

    //Normalize histograms to unit area for comparison.
    h_PtMuPl_plot->Scale(1.0 / h_PtMuPl_plot->Integral());
    h_PtMuMi_plot->Scale(1.0 / h_PtMuMi_plot->Integral());
    h_RefPl_plot->Scale(1.0 / h_RefPl_plot->Integral());
    h_RefMi_plot->Scale(1.0 / h_RefMi_plot->Integral());

    //Reference histograms with points
    h_RefPl_plot->SetFillStyle(0);
    h_RefPl_plot->SetMarkerStyle(22);
    h_RefPl_plot->SetMarkerSize(0.75);
    h_RefPl_plot->SetMarkerColorAlpha(kRed, 1.);

    h_RefMi_plot->SetFillStyle(0);
    h_RefMi_plot->SetMarkerStyle(24);
    h_RefMi_plot->SetMarkerSize(0.75);
    h_RefMi_plot->SetMarkerColorAlpha(kBlue, 1.);

    h_RefPl_plot->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    h_RefPl_plot->GetYaxis()->SetTitle("Normalized Yield");
    h_RefPl_plot->GetXaxis()->SetRangeUser(18., 100.);
    h_RefPl_plot->GetYaxis()->SetTitleOffset(1.);

    //Since we plot h_RefPl_plot first, we set the maximum on it:
    /*double maxY = std::max({
        h_PtMuPl_plot->GetMaximum(),
        h_PtMuMi_plot->GetMaximum(),
        h_RefPl_plot->GetMaximum(),
        h_RefMi_plot->GetMaximum()
    });
    h_RefPl_plot->SetMaximum(1.25 * maxY);*/

    //Draw.
    //h_RefPl_plot->Draw("HIST P");
    //h_RefMi_plot->Draw("HIST P SAME");

    //Filling histogram with PbPb values. No border.
    auto* fillPl = static_cast<TH1*>(h_PtMuPl_plot->Clone(("h_fillPl" + centString).c_str()));
    fillPl->SetDirectory(nullptr);
    fillPl->SetFillStyle(1001);
    fillPl->SetFillColorAlpha(kRed-10, 0.6);
    fillPl->SetLineColorAlpha(kRed-10, 0.0);
    //fillPl->Draw("HIST ][ SAME");

    auto* fillMi = static_cast<TH1*>(h_PtMuMi_plot->Clone(("h_fillMi" + centString).c_str()));
    fillMi->SetDirectory(nullptr);
    fillMi->SetFillStyle(1001);
    fillMi->SetFillColorAlpha(kBlue-10, 0.6);
    fillMi->SetLineColorAlpha(kBlue-10, 0.0);
    //fillMi->Draw("HIST ][ SAME");


    //Testing new draw style. Axis first, then fill, then points.
    h_RefPl_plot->Draw("AXIS");

    fillPl->Draw("HIST ][ SAME");
    fillMi->Draw("HIST ][ SAME");

    h_RefPl_plot->Draw("P SAME");
    h_RefMi_plot->Draw("P SAME");


    //Selections and cuts
    drawLatexText(Form("p_{T} > %.0f GeV, |#eta| < %.1f", PTCUTVALUE, ETACUTVALUE), 0.7, 0.3, 0.03);
    drawLatexText(Form("|y| < %.1f", RAPIDITYCUTVALUE), 0.7, 0.25, 0.03);
    drawLatexText(Form("%.0f < M_{#mu#mu} < %.0f GeV", MINZ_MASS, MAXZ_MASS), 0.7, 0.2, 0.03);
    drawLatexText(Form("Centrality: %.0f - %.0f %%", lowCent, highCent), 0.7, 0.15, 0.03);

    TLegend* leg = new TLegend(0.72, 0.7, 0.95, 0.85);
    basicLegendFormatting(leg);
    leg->AddEntry(fillPl, "p_{T}(#mu^{+})", "f");
    leg->AddEntry(fillMi, "p_{T}(#mu^{-})", "f");
    leg->AddEntry(h_RefPl_plot, "p_{T}(#mu^{+}) ppRef", "p");
    leg->AddEntry(h_RefMi_plot, "p_{T}(#mu^{-}) ppRef", "p");
    leg->Draw();
    pad1->Update();

    //Bottom pad. Ratio of pT distributions of mu+ and mu-.
    //Bottom pad histograms are calculated using the original histograms, not the normalized ones.
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
            std::pow(PbPbratioError_completeCorr/PbPbratio, 2) + std::pow(refRatioError_completeCorr/refRatio, 2));//Independent error propagation.

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

    //Perform a Chi2 test to see if the double ratio is compatible with 1.0.
    auto [ndf, chi2, pValue] = MakeChi2Test(histRatio);

    //Back to the canvas.
    c->cd();
    drawLatexText("#bf{CMS}", 0.11, 0.96, 0.035);
    drawLatexText("#it{Work in Progress}", 0.2, 0.96, 0.025);
    drawLatexText(dataSamplesUsed.c_str(), 0.7, 0.96, 0.025);
    drawLatexText(Form("#chi^{2}/ndf = %.2f/%d", chi2, ndf), 0.67, 0.72, 0.022);
    drawLatexText(Form("p-value = %.2f", pValue), 0.67, 0.69, 0.022);

    //Save
    std::string output = "MuonPtMuPlMuMiHistWithDoubleRatio" + centString + plot_extension;
    c->Update();
    c->SaveAs(output.c_str());

    //Clean owned memory
    delete fillPl;
    delete fillMi;
    delete h_PtMuPl_plot;
    delete h_PtMuMi_plot;
    delete h_RefPl_plot;
    delete h_RefMi_plot;
    delete h_RefPl;
    delete h_RefMi;
    delete h_ogpl;
    delete h_ogmi;
    delete h2D_PtMuPl_PtMuMi;
    delete h3D_PtMuPl_PtMuMi_Cent_clone;
    delete h_PtMuPl;
    delete h_PtMuMi;
    delete histRatio;
    delete line;
    delete leg;
    delete leg2;
    delete c;
}

std::tuple<int, double, double> MakeChi2Test(const TH1D* histRatio){

    double ptMin = 40.;
    double ptMax = 65.;

    TF1* nullHypothesis = new TF1("nullHypothesis","1.0",ptMin,ptMax);

    double chi2 = histRatio->Chisquare(nullHypothesis, "R");

    int ndf = 0;
    //Count the number of bins with non-zero content and within the pT range of interest to determine the degrees of freedom.
    for(int i = 1; i <= histRatio->GetNbinsX(); ++i){
        if(histRatio->GetBinContent(i) == 0.0)
            continue;
        double x = histRatio->GetBinCenter(i);
        if(x < ptMin || x > ptMax)
            continue;
        if(histRatio->GetBinError(i) <= 0.0)
            continue;
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