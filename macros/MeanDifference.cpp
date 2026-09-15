/*
Mean and peak difference as function of centrality bin.
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

//Vectors to save peak difference and error.
std::vector<std::pair<double, double>> PeakDiffAndError_ppRef;
std::vector<std::pair<double, double>> PeakDifferenceAndError_PbPb_vs_ppRef;

//---Function declarations
std::pair<double, double> Get_ppRefValues(TFile* inputFile);
void CalculateMeanDifference_PbPb_vs_ppRef(TFile* inputFile, double lowCent, double highCent, std::pair<double, double> ppRefValues);
void PlotMeanDifference_PbPb_vs_ppRef();
void CalculatePeakDifference_PbPb_vs_ppRef(TFile* inputFile, double lowCent, double highCent, std::pair<double, double> ppRefValues);
void PlotPeakDifference_PbPb_vs_ppRef();

//---Main()
void MeanDifference(){

    gROOT->SetBatch(kTRUE);

    MeanDifferenceAndError_PbPb_vs_ppRef.clear();
    PeakDiffAndError_ppRef.clear();
    PeakDifferenceAndError_PbPb_vs_ppRef.clear();

    TFile* inputFile = new TFile("mySelectedData.root", "READ");

    //Get ppRef values for mean difference between pT of mu+ and mu- and respective error of measurement.
    //Peak difference and error will also be saved in the global vector PeakDiffAndError_ppRef.
    std::pair<double, double> ppRefValues = Get_ppRefValues(inputFile);

    for(const auto& centBin : CentralityBinsSet) {
        double lowCent = centBin.first;
        double highCent = centBin.second;
        CalculateMeanDifference_PbPb_vs_ppRef(inputFile, lowCent, highCent, ppRefValues);
        CalculatePeakDifference_PbPb_vs_ppRef(inputFile, lowCent, highCent, PeakDiffAndError_ppRef[0]); //Using the peak difference and error from ppRef as reference.
    }

    PlotMeanDifference_PbPb_vs_ppRef();
    PlotPeakDifference_PbPb_vs_ppRef();
    inputFile->Close();
}

void PlotPeakDifference_PbPb_vs_ppRef(){

    int nPoints = PeakDifferenceAndError_PbPb_vs_ppRef.size();
    if(nPoints != CentralityBinsSet.size()){
            std::cerr << "Error: Number of points in PeakDifferenceAndError does not match number of centrality bins." << std::endl;
            return;
    }

    std::vector<double> xValues(nPoints);
    std::vector<double> yValues(nPoints);
    std::vector<double> xErrors(nPoints);
    std::vector<double> yErrors(nPoints);

    for(int i = 0; i < nPoints; ++i){
        xValues[i] = i + 1;
        xErrors[i] = 0.;
        yValues[i] = PeakDifferenceAndError_PbPb_vs_ppRef[i].first;
        yErrors[i] = PeakDifferenceAndError_PbPb_vs_ppRef[i].second;
    }

    //This block is outdated because i opted to plot with cBins as x-coordinates.
    //Get values that will be plotted with TGraphErrors.
    /*for(int i = 0; i < nPoints; i++){

        xValues[i] = (CentralityBinsSet[i].first + CentralityBinsSet[i].second) / 2.0;//Centrality bin center for now.
        xErrors[i] = (CentralityBinsSet[i].second - CentralityBinsSet[i].first) / 2.0;//Half-width of the centrality bin.
        yValues[i] = PeakDifferenceAndError_PbPb_vs_ppRef[i].first;
        yErrors[i] = PeakDifferenceAndError_PbPb_vs_ppRef[i].second;
    }*/

    //Set the TGraphErrors
    TGraphErrors* graph = new TGraphErrors(nPoints, xValues.data(), yValues.data(), xErrors.data(), yErrors.data());
    basicGraphFormatting(graph);

    //Create canvas.
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    basicCanvasFormatting(c);

    //Frame TH1 helper to set the x-axis labels for centrality bins.
    TH1D* frame = new TH1D("frame","",nPoints,0.5,nPoints + 0.5);
    basicHistFormatting(frame);
    for(int i = 0; i < nPoints; ++i){
        std::string label = Form("%.0f-%.0f%%", CentralityBinsSet[i].first, CentralityBinsSet[i].second);
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
    yMax = std::max(yMax+1.5, 0.0);
    frame->SetMinimum(yMin);
    frame->SetMaximum(yMax);
    //

    frame->GetXaxis()->SetTickLength(0.0);
    frame->GetXaxis()->SetTitle("Centrality bin");
    frame->GetYaxis()->SetTitle("#Delta #bar{p}^{PbPb}_{T} - #Delta #bar{p}^{ppRef}_{T} [GeV/c]");

    //Draw only the axis frame.
    frame->Draw("AXIS");

    //Format graph.
    graph->SetMarkerStyle(21);
    graph->SetMarkerSize(0.9);
    graph->SetMarkerColorAlpha(kRed+1, 1.);
    graph->SetLineColorAlpha(kRed-7, 0.8);
    graph->SetLineWidth(2);

    //Axes configurations. Obsolete.
    //graph->GetXaxis()->SetLimits(0, 100);
    //graph->GetXaxis()->SetTitle("Centrality (%)");
    //graph->GetYaxis()->SetTitle("#Delta p^{PbPb}_{T,peak} - #Delta p^{ppRef}_{T,peak} [GeV/c]");
    //graph->GetYaxis()->SetTitleOffset(1.3);
    
    //Draw
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
    drawLatexText("p_{T}^{#mu} > 20 GeV, |#eta^{#mu}| < 2.4", 0.22, 0.3, 0.03);
    drawLatexText("60 < M_{#mu#mu} < 120 GeV", 0.22, 0.25, 0.03);

    c->Update();
    std::string outputName = "PeakDifference_vs_Centrality" + plot_extension;
    c->SaveAs(outputName.c_str());

    delete frame;
    delete line;
    delete graph;
    delete c;
}

void CalculatePeakDifference_PbPb_vs_ppRef(TFile* inputFile, double lowCent, double highCent, std::pair<double, double> ppRefValues){

    //Get directory
    TDirectory *PbPb_dir = inputFile->GetDirectory("PbPb2024_Data");

    //Original histogram
    TH3D* h3D_PtMuPl_PtMumi_Cent = dynamic_cast<TH3D*>(PbPb_dir->Get("h3D_PtMuPl_PtMuMi_Cent"));
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
    TH2D* h2D_PtMuPl_PtMumi = dynamic_cast<TH2D*>(h3D_PtMuPl_PtMumi_Cent->Project3D("yx"));

    //Project further into TH1 histograms
    TH1D* h1D_PtMuPl = dynamic_cast<TH1D*>(h2D_PtMuPl_PtMumi->ProjectionX("h1D_PtMuPl"));
    TH1D* h1D_PtMuMi = dynamic_cast<TH1D*>(h2D_PtMuPl_PtMumi->ProjectionY("h1D_PtMuMi"));

    //Calculate peak pT of mu+ and mu- in this centrality range.
    double MuPl_peak = h1D_PtMuPl->GetBinCenter(h1D_PtMuPl->GetMaximumBin());
    double MuPl_peak_error = h1D_PtMuPl->GetBinWidth(h1D_PtMuPl->GetMaximumBin()); //For now error is taken to be the width of the bin with maximum content.
    double MuMi_peak = h1D_PtMuMi->GetBinCenter(h1D_PtMuMi->GetMaximumBin());
    double MuMi_peak_error = h1D_PtMuMi->GetBinWidth(h1D_PtMuMi->GetMaximumBin());


    //Calculate the difference and its error in PbPb sample.
    double PeakDifferencePbPb = MuPl_peak - MuMi_peak;
    double PeakDifferenceErrorPbPb = std::sqrt(MuPl_peak_error * MuPl_peak_error + MuMi_peak_error * MuMi_peak_error);


    //Subtract from reference values.
    double PeakDifference = PeakDifferencePbPb - ppRefValues.first;
    //Now with independent error propagation:
    double PeakDifferenceError = std::sqrt(std::pow(PeakDifferenceErrorPbPb, 2) + std::pow(ppRefValues.second, 2));

    //Store the peak difference and its error for this centrality bin.
    PeakDifferenceAndError_PbPb_vs_ppRef.push_back(std::make_pair(PeakDifference, PeakDifferenceError));

    delete h2D_PtMuPl_PtMumi;
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

    for(int i = 0; i < nPoints; ++i){
        xValues[i] = i + 1;
        xErrors[i] = 0.;
        yValues[i] = MeanDifferenceAndError_PbPb_vs_ppRef[i].first;
        yErrors[i] = MeanDifferenceAndError_PbPb_vs_ppRef[i].second;
    }

    //This block is outdated because i opted to plot with cBins as x-coordinates.
    //Get values that will be plotted with TGraphErrors.
    /*for(int i = 0; i < nPoints; i++){

        xValues[i] = (CentralityBinsSet[i].first + CentralityBinsSet[i].second) / 2.0;//Centrality bin center for now.
        xErrors[i] = (CentralityBinsSet[i].second - CentralityBinsSet[i].first) / 2.0;//Half-width of the centrality bin.
        yValues[i] = MeanDifferenceAndError_PbPb_vs_ppRef[i].first;
        yErrors[i] = MeanDifferenceAndError_PbPb_vs_ppRef[i].second;
    }*/

    //Set the TGraphErrors
    TGraphErrors* graph = new TGraphErrors(nPoints, xValues.data(), yValues.data(), xErrors.data(), yErrors.data());
    basicGraphFormatting(graph);

    //Create canvas.
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    basicCanvasFormatting(c);

    //Frame TH1 helper to set the x-axis labels for centrality bins.
    TH1D* frame = new TH1D("frame","",nPoints,0.5,nPoints + 0.5);
    basicHistFormatting(frame);
    for(int i = 0; i < nPoints; ++i){
        std::string label = Form("%.0f-%.0f%%", CentralityBinsSet[i].first, CentralityBinsSet[i].second);
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
    frame->GetYaxis()->SetTitle("#Delta #bar{p}^{PbPb}_{T} - #Delta #bar{p}^{ppRef}_{T} [GeV/c]");

    //Draw only the axis frame.
    frame->Draw("AXIS");

    //Format graph.
    graph->SetMarkerStyle(21);
    graph->SetMarkerSize(0.9);
    graph->SetMarkerColorAlpha(kRed+1, 1.);
    graph->SetLineColorAlpha(kRed-7, 0.8);
    graph->SetLineWidth(2);

    //Axes configurations. Obsolete.
    //graph->GetXaxis()->SetLimits(0, 100);
    //graph->GetXaxis()->SetTitle("Centrality (%)");
    //graph->GetYaxis()->SetTitle("#Delta #bar{p}^{PbPb}_{T} - #Delta #bar{p}^{ppRef}_{T} [GeV/c]");
    //graph->GetYaxis()->SetTitleOffset(1.3);
    
    //Draw
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
    drawLatexText("p_{T}^{#mu} > 20 GeV, |#eta^{#mu}| < 2.4", 0.68, 0.3, 0.03);
    drawLatexText("60 < M_{#mu#mu} < 120 GeV", 0.68, 0.25, 0.03);

    c->Update();
    std::string outputName = "MeanDifference_vs_Centrality" + plot_extension;
    c->SaveAs(outputName.c_str());

    delete graph;
    delete frame;
    delete line;
    delete c;
}

void CalculateMeanDifference_PbPb_vs_ppRef(TFile* inputFile, double lowCent, double highCent, std::pair<double, double> ppRefValues){

    //Get directory
    TDirectory *PbPb_dir = inputFile->GetDirectory("PbPb2024_Data");

    //Original histogram
    TH3D* h3D_PtMuPl_PtMumi_Cent = dynamic_cast<TH3D*>(PbPb_dir->Get("h3D_PtMuPl_PtMuMi_Cent"));
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
    TH2D* h2D_PtMuPl_PtMumi = dynamic_cast<TH2D*>(h3D_PtMuPl_PtMumi_Cent->Project3D("yx"));

    //Calculate mean pT of mu+ and mu- in this centrality range.
    double MuPl_mean = h2D_PtMuPl_PtMumi->GetMean(1); //1- x-axis
    double MuPl_mean_error = h2D_PtMuPl_PtMumi->GetMeanError(1);

    double MuMi_mean = h2D_PtMuPl_PtMumi->GetMean(2); //2- y-axis
    double MuMi_mean_error = h2D_PtMuPl_PtMumi->GetMeanError(2);

    //Calculate the difference and its error in PbPb sample.
    double MeanDifferencePbPb = MuPl_mean - MuMi_mean;
    //Error here using completely correlated formula from Lara's reference.
    double MeanDifferenceErrorPbPb = std::sqrt( std::abs(MuPl_mean_error * MuPl_mean_error - MuMi_mean_error * MuMi_mean_error));

    //Subtract from reference values.
    double MeanDifference = MeanDifferencePbPb - ppRefValues.first;
    //Now with independent error propagation:
    double MeanDifferenceError = std::sqrt(std::pow(MeanDifferenceErrorPbPb,2) + std::pow(ppRefValues.second,2)); //NEEDS REVIEW. MEASUREMENT MAY BE CORRELATED.

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
    double MeanDiffError = std::sqrt( std::abs(MuPl_mean_error * MuPl_mean_error - MuMi_mean_error * MuMi_mean_error));//From Lara's reference. Completely correlated.

    //Get also peak value
    double MuPl_peak = h_MuPl->GetBinCenter(h_MuPl->GetMaximumBin());
    double MuMi_peak = h_MuMi->GetBinCenter(h_MuMi->GetMaximumBin());
    double MuPl_peak_error = h_MuPl->GetBinWidth(h_MuPl->GetMaximumBin()); //For now error is taken to be the width of the bin with maximum content.
    double MuMi_peak_error = h_MuMi->GetBinWidth(h_MuMi->GetMaximumBin());
    double PeakDiff = MuPl_peak - MuMi_peak;
    double PeakDiffError = std::sqrt(MuPl_peak_error * MuPl_peak_error + MuMi_peak_error * MuMi_peak_error);//Completely correlated.

    //We may want to return their peak values also. Review that later.
    delete h_MuPl;
    delete h_MuMi;

    //Saving the peak difference and error for later use.
    PeakDiffAndError_ppRef.push_back(std::make_pair(PeakDiff, PeakDiffError));

    return std::make_pair(MeanDiff, MeanDiffError);
}