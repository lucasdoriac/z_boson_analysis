/*
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
std::string plot_extension = ".png"; // ".png" for regular development and ".pdf" for final quality plots
double delta = 1e-6; //Small value to avoid binning issues when projecting histograms.


// ##############################################################################
// ##############################################################################


//Set of centrality bins for PbPb2024 data. We can decide to change the centrality bins later if we want to.
std::vector<std::pair<double, double>> CentralityBinsSet = {
    {0., 10.},
    {10., 20.},
    {20., 30.},
    {30., 100.}
};

//Second proposed set of centrality bins for PbPb2024 data.
/*std::vector<std::pair<double, double>> CentralityBinsSet = {
    {0., 10.},
    {10., 30.},
    {30., 50.},
    {50., 100.}
};*/

//Vector to save data for TGraphErrors at the end.
std::vector<std::pair<double, double>> MeanDifferenceAndError_PbPb_vs_ppRef;


//---Function declarations
std::pair<double, double> Get_ppRefValues(TFile* inputFile);
void CalculateMeanDifference_PbPb_vs_ppRef(TFile* inputFile, double lowCent, double highCent);
void PlotMeanDifference_PbPb_vs_ppRef();


//---Main()
void MuonYieldAsymmetry(){

    gROOT->SetBatch(kTRUE);
    TFile* inputFile = new TFile("mySelectedData.root", "READ");

    //Get ppRef values for mean difference between pT of mu+ and mu- and respective error of measurement.
    std::pair<double, double> ppRefValues = Get_ppRefValues(inputFile);

    for(const auto& centBin : CentralityBinsSet) {
        double lowCent = centBin.first;
        double highCent = centBin.second;
        CalculateMeanDifference_PbPb_vs_ppRef(inputFile, lowCent, highCent);
    }

    PlotMeanDifference_PbPb_vs_ppRef();
    inputFile->Close();
}


void PlotMeanDifference_PbPb_vs_ppRef(){

    int nPoints = MeanDifferenceAndError_PbPb_vs_ppRef.size();
    if(nPoints != CentralityBinsSet.size()){
            std::cerr << "Error: Number of points in DeltaPtAndError does not match number of centrality bins." << std::endl;
            return;
    }

    std::vector<double> xValues(nPoints);
    std::vector<double> yValues(nPoints);
    std::vector<double> xErrors(nPoints);
    std::vector<double> yErrors(nPoints);

    //Get values that will be plotted with TGraphErrors.
    for(int i = 0; i < nPoints; i++){

        xValues[i] = (CentralityBinsSet[i].first + CentralityBinsSet[i].second) / 2.0;//Centrality bin center for now.
        xErrors[i] = (CentralityBinsSet[i].second - CentralityBinsSet[i].first) / 2.0;//Half-width of the centrality bin.
        yValues[i] = MeanDifferenceAndError_PbPb_vs_ppRef[i].first;
        yErrors[i] = MeanDifferenceAndError_PbPb_vs_ppRef[i].second;
    }

    //Set the TGraphErrors
    TGraphErrors* graph = new TGraphErrors(nPoints, xValues.data(), yValues.data(), xErrors.data(), yErrors.data());
    basicGraphFormatting(graph);

    //Create canvas.
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    basicCanvasFormatting(c);

    //Format graph.
    graph->SetMarkerStyle(21);
    graph->SetMarkerSize(1.0);
    graph->SetMarkerColor(kRed+1);
    graph->SetLineColor(kRed+1);

    //Axes configurations
    graph->GetXaxis()->SetLimits(0, 100);
    graph->GetXaxis()->SetTitle("Centrality (%)");
    graph->GetYaxis()->SetTitle("#LT p_{T}^{#mu^{+}} - p_{T}^{#mu^{-}} #GT [GeV/c]");
    graph->GetYaxis()->CenterTitle(true);
    graph->GetYaxis()->SetTitleOffset(1.3);
    
    //Draw
    graph->Draw("AP SAME");

    //Grey line at y=0
    TLine* line = new TLine(graph->GetXaxis()->GetXmin(), 0.0, graph->GetXaxis()->GetXmax(), 0.0);
    line->SetLineColor(kGray);
    line->SetLineStyle(7);
    line->SetLineWidth(2);
    line->Draw();

    drawLatexText("#bf{CMS}", 0.12, 0.93, 0.042);
    drawLatexText("#it{Work in Progress}", 0.2, 0.93, 0.033);
    drawLatexText("PbPb 2024, ppRef 2024 (5.36 TeV)", 0.6, 0.93, 0.033);
    
    //Plot specifications
    drawLatexText("p_{T}^{#mu} > 20 GeV, |#eta^{#mu}| < 2.4", 0.2, 0.75, 0.03);
    drawLatexText("60 < M_{#mu #mu} < 120 GeV", 0.2, 0.7, 0.03);

    c->Update();
    std::string outputName = "MeanDifference_vs_Centrality" + plot_extension;
    c->SaveAs(outputName.c_str());

    delete graph;
    delete c;
}

void CalculateMeanDifference_PbPb_vs_ppRef(TFile* inputFile, double lowCent, double highCent){

    //Get directory
    TDirectory *PbPb_dir = inputFile->GetDirectory("PbPb2024_Data");

    //Original histogram
    TH3D* h3D_PtMuPl_PtMumi_Cent = dynamic_cast<TH3D*>(PbPb_dir->Get("h3D_PtMuPl_PtMumi_Cent"));
    if(!h3D_PtMuPl_PtMumi_Cent) {
        std::cerr << "Error: Could not find the histogram h3D_PtMuPl_PtMumi_Cent in the input file." << std::endl;
        return;
    }

    /*
    z-axis -> centrality
    y-axis -> pT of mu-
    x-axis -> pT of mu+
    */

    //Make TH2D projection for the specified centrality range.
    int binLow = h3D_PtMuPl_PtMumi_Cent->GetZaxis()->FindBin(lowCent + delta);
    int binHigh = h3D_PtMuPl_PtMumi_Cent->GetZaxis()->FindBin(highCent - delta);

    //Select centrality range and project onto pT(mu+) vs pT(mu-) plane.
    h3D_PtMuPl_PtMumi_Cent->GetZaxis()->SetRange(binLow, binHigh);
    TH2D* h2D_PtMuPl_PtMumi = dynamic_cast<TH2D*>(h3D_PtMuPl_PtMumi_Cent->Project3D("xy"));

    //Calculate mean pT of mu+ and mu- in this centrality range.
    double MuPl_mean = h2D_PtMuPl_PtMumi->GetMean(1); //1- x-axis
    double MuPl_mean_error = h2D_PtMuPl_PtMumi->GetMeanError(1);

    double MuMi_mean = h2D_PtMuPl_PtMumi->GetMean(2); //2- y-axis
    double MuMi_mean_error = h2D_PtMuPl_PtMumi->GetMeanError(2);

    //Calculate the difference and its error.
    double MeanDifference = MuPl_mean - MuMi_mean;
    double MeanDifferenceError = std::sqrt(std::pow(MuPl_mean_error, 2) + std::pow(MuMi_mean_error, 2)); //NEEDS REVIEW. MEASUREMENT MAY BE CORRELATED.

    //Store the mean difference and its error for this centrality bin.
    MeanDifferenceAndError_PbPb_vs_ppRef.push_back(std::make_pair(MeanDifference, MeanDifferenceError));

    delete h2D_PtMuPl_PtMumi;
}


std::pair<double, double> Get_ppRefValues(TFile* inputFile){

    //Get directory
    TDirectory *ppRef_dir = inputFile->GetDirectory("ppRef2024_Data");

    //Original histograms
    TH1D* og_plus = dynamic_cast<TH1D*>(ppRef_dir->Get("h1D_ptMuPlus"));
    TH1D* og_minus = dynamic_cast<TH1D*>(ppRef_dir->Get("h1D_ptMuMinus"));

    //Get clone to manipulate
    TH1D* h_MuPl = dynamic_cast<TH1D*>(og_plus->Clone("h_MuPl"));
    TH1D* h_MuMi = dynamic_cast<TH1D*>(og_minus->Clone("h_MuMi"));
    h_MuPl->SetDirectory(nullptr);
    h_MuMi->SetDirectory(nullptr);

    //Statistics calculation
    //Get mean of MuPl
    double MuPl_mean = h_MuPl->GetMean();
    double MuPl_mean_error = h_MuPl->GetMeanError();

    //Get mean of MuMi
    double MuMi_mean = h_MuMi->GetMean();
    double MuMi_mean_error = h_MuMi->GetMeanError();

    //What we want is their difference.
    double MeanDiff = MuPl_mean - MuMi_mean;
    double MeanDiffError = std::sqrt(std::pow(MuPl_mean_error, 2) + std::pow(MuMi_mean_error, 2)); //NEEDS REVIEW. MEASUREMENT MAY BE CORRELATED.

    //We may want to return their peak values also. Review that later.
    delete h_MuPl;
    delete h_MuMi;

    return std::make_pair(MeanDiff, MeanDiffError);
}