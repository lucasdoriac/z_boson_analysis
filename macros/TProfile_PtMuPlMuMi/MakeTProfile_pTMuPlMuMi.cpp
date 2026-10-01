//Make TProfile from mySelectedData.root

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
std::string whichDataset = "PbPb2023_2024_Data"; // "PbPb2023_2024_Data", "PbPb2023_Data", "PbPb2024_Data". 
std::string JointPbPb = "PbPb2023+2024"; //"PbPb2023+2024", "PbPb2023", "PbPb2024".
std::string dataSamplesUsed = "PbPb 2023+2024, ppRef 2024 (5.36 TeV)"; //"PbPb 2023+2024, ppRef 2024 (5.36 TeV)", "PbPb 2023, ppRef 2024 (5.36 TeV)", "PbPb 2024, ppRef 2024 (5.36 TeV)".
double delta = 1e-6; //Small value to avoid binning issues when projecting histograms.
double ptMin = 20.; //Minimum pT value for the histograms.


// ##############################################################################
// ##############################################################################

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


void MakeTProfile_PbPb(TFile* inputFile, double lowCent, double highCent);


//---Main()
void MakeTProfile_pTMuPlMuMi(){

    gROOT->SetBatch(kTRUE);
    TFile* inputFile = new TFile("mySelectedData.root", "READ");

    for (const auto& bin : CentralityBinsSet) {
        MakeTProfile_PbPb(inputFile, bin.first, bin.second);
    }
    inputFile->Close();
}


void MakeTProfile_PbPb(TFile* inputFile, double lowCent, double highCent){

    std::string centString = std::to_string(static_cast<int>(lowCent)) + "-" + std::to_string(static_cast<int>(highCent));

    TDirectory* dir = inputFile->GetDirectory(whichDataset.c_str());
    std::string histName = "h3D_PtMuPl_PtMuMi_Cent";

    //Get original TH3 histogram
    TH3D* h_original = dynamic_cast<TH3D*>(dir->Get(histName.c_str()));

    //Clone the histogram to manipulate.
    TH3D* h_PbPb = dynamic_cast<TH3D*>(h_original->Clone(("h_PbPb_" + centString).c_str()));

    /*
    z-axis: Centrality
    y-axis: PtMuMi
    x-axis: PtMuPl
    */

    //Get centrality range.
    int binLow  = h_PbPb->GetXaxis()->FindBin(lowCent + delta);
    int binHigh = h_PbPb->GetXaxis()->FindBin(highCent - delta);

    //Select centrality range and project onto pT(mu+) vs pT(mu-) plane.
    h_PbPb->GetZaxis()->SetRange(binLow, binHigh);
    TH2D* h2D_PtMuPl_PtMuMi = dynamic_cast<TH2D*>(h_PbPb->Project3D("yx"));
    
    //Now h2D_PtMuPl_PtMuMi is a 2D histogram with PtMuPl on the x-axis and PtMuMi on the y-axis.
    //Now we can create a TProfile to get the mean of PtMuMi for each bin of PtMuPl. 
    TProfile* h_PbPb_TProfile = h2D_PtMuPl_PtMuMi->ProfileX(("h_PbPb_TProfile_" + centString).c_str());

    TCanvas* c = new TCanvas(("c_PbPb_" + centString).c_str(), ("c_PbPb_" + centString).c_str(), 800, 600);
    basicCanvasFormatting(c);
    basicHistFormatting(h2D_PtMuPl_PtMuMi);
    c->SetLeftMargin(0.1);
    c->SetRightMargin(0.14);


    //TH2 formatting
    h2D_PtMuPl_PtMuMi->GetXaxis()->SetTitle("p_{T}^{#mu^{+}} [GeV]");
    h2D_PtMuPl_PtMuMi->GetYaxis()->SetTitle("p_{T}^{#mu^{-}} [GeV]");
    h2D_PtMuPl_PtMuMi->GetYaxis()->SetTitleOffset(1.);
    h2D_PtMuPl_PtMuMi->GetXaxis()->CenterTitle(true);
    h2D_PtMuPl_PtMuMi->GetYaxis()->CenterTitle(true);
    h2D_PtMuPl_PtMuMi->GetZaxis()->SetTitle("Dimuons");
    h2D_PtMuPl_PtMuMi->GetZaxis()->SetLabelSize(0.04);
    h2D_PtMuPl_PtMuMi->GetZaxis()->SetLabelOffset(0.01);
    h2D_PtMuPl_PtMuMi->GetZaxis()->SetTitleOffset(1.2);
    h2D_PtMuPl_PtMuMi->GetZaxis()->SetTitleSize(0.033);

    //TProfile formatting
    h_PbPb_TProfile->SetLineColorAlpha(kGray+3, 0.8);
    h_PbPb_TProfile->SetMarkerColorAlpha(kBlack, 1.);
    h_PbPb_TProfile->SetMarkerStyle(20);
    h_PbPb_TProfile->SetMarkerSize(0.6);
    h_PbPb_TProfile->SetLineWidth(1);

    //Range
    h2D_PtMuPl_PtMuMi->GetXaxis()->SetRangeUser(ptMin, 100.);
    h2D_PtMuPl_PtMuMi->GetYaxis()->SetRangeUser(ptMin, 100.);
    h_PbPb_TProfile->GetXaxis()->SetRangeUser(ptMin, 100.);
    h_PbPb_TProfile->GetYaxis()->SetRangeUser(ptMin, 100.);
    h2D_PtMuPl_PtMuMi->SetContour(100);

    //Draw
    h2D_PtMuPl_PtMuMi->Draw("COLZ");
    h_PbPb_TProfile->Draw("E1 SAME");

    //Reference y=x diagonal line
    TLine* diagonal = new TLine(ptMin, ptMin,  h2D_PtMuPl_PtMuMi->GetXaxis()->GetXmax(), h2D_PtMuPl_PtMuMi->GetXaxis()->GetXmax());

    diagonal->SetLineStyle(7);
    diagonal->SetLineWidth(2);
    diagonal->SetLineColor(kRed);
    diagonal->Draw("SAME");

    drawLatexText("#bf{CMS}", 0.12, 0.93, 0.04);
    drawLatexText("#it{Work in Progress}", 0.19, 0.93, 0.03);
    drawLatexText(dataSamplesUsed.c_str(), 0.48, 0.93, 0.03);
    
    //Plot specifications
    drawLatexText("p_{T}^{#mu} > 20 GeV, |#eta^{#mu}| < 2.4", 0.6, 0.87, 0.025);
    drawLatexText("60 < M_{#mu#mu} < 120 GeV", 0.6, 0.82, 0.025);


    c->Update();
    std::string outputFileName = "TProfile_PtMuPl_PtMuMi_Cent_" + centString + plot_extension;
    c->SaveAs(outputFileName.c_str());

    delete h_PbPb;
    delete h2D_PtMuPl_PtMuMi;
    delete h_PbPb_TProfile;
    delete c;
    delete diagonal;
}