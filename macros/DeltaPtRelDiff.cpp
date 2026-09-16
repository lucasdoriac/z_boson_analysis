/*
This macro currently does three things. First, it makes the normalized distribution of 'PtRelDiff'
for PbPb2024 and ppRef2024 for the full centrality range (0-100%).

Second, it makes the normalized distribution of 'PtRelDiff' for PbPb2024 and ppRef2024 for each centrality bin defined in the vector 'CentralityBinsSet'.

Third, it calculates the difference between the mean of 'PtRelDiff' for PbPb2024 relative to ppRef2024 for each centrality bin
and stores the results in a vector of pairs.

The most important result is the difference between the mean of 'PtRelDiff' for PbPb2024 relative to ppRef2024.
That is,

    DeltaPt = <PtRelDiff>_{PbPb2024} - <PtRelDiff>_{ppRef2024}.

The ppRef value comes from a single distribution of DeltaPtRelDiff calculated on the ppRef2024 sample.
The PbPb value is calculated at each **centrality bin** before it is subtracted from the ppRef value.
Therefore, DeltaPt is a function of centrality bin, as it was our goal.

DeltaPt values are recorded in a vector of pairs, where the first element is the nominal value and the second element is the measurement error.
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
#include <iomanip>
#include "../headers/basicFormatting.h"


//---Macro settings
std::string plot_extension = ".pdf"; // ".png" for regular development and ".pdf" for final quality plots
double delta = 1e-6; //Small value to avoid binning issues when projecting histograms.


// ##############################################################################
// ##############################################################################


//Vectors to store final data cent dependent.
std::vector<std::pair<double, double>> DeltaPtAndError;


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


//---Function declarations
void MakeNormalizedDistPtRelDiff_PbPb_vs_ppRef_0_100(TFile* inputFile);
void MakeNormalizedDistPtRelDiff_PbPb_vs_ppRef_CentBin(TFile* inputFile, double lowCent, double highCent);
void PlotDeltaPt_vsCentralityBin();

//---Main()
void DeltaPtRelDiff(){

    DeltaPtAndError.clear();

    gROOT->SetBatch(kTRUE);
    TFile* inputFile = new TFile("mySelectedData.root", "READ");

    MakeNormalizedDistPtRelDiff_PbPb_vs_ppRef_0_100(inputFile);

    for(const auto& centBin : CentralityBinsSet){
        double lowCent = centBin.first;
        double highCent = centBin.second;
        std::cout << "Processing Centrality bin: " << lowCent << "-" << highCent << "%" << std::endl;
    
        MakeNormalizedDistPtRelDiff_PbPb_vs_ppRef_CentBin(inputFile, lowCent, highCent);
    }

    PlotDeltaPt_vsCentralityBin();

    inputFile->Close();
}

void PlotDeltaPt_vsCentralityBin(){

    int nPoints = DeltaPtAndError.size();//Number of centrality bins. Needs to be equal to 4.
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
        yValues[i] = DeltaPtAndError[i].first;
        yErrors[i] = DeltaPtAndError[i].second;
    }

    //This block is outdated because i opted to plot with cBins as x-coordinates.
    //Get values that will be plotted with TGraphErrors.
    /*for(int i = 0; i < nPoints; i++){

        xValues[i] = (CentralityBinsSet[i].first + CentralityBinsSet[i].second) / 2.0;//Centrality bin center for now.
        xErrors[i] = (CentralityBinsSet[i].second - CentralityBinsSet[i].first) / 2.0;//Half-width of the centrality bin.
        yValues[i] = DeltaPtAndError[i].first;
        yErrors[i] = DeltaPtAndError[i].second;
    }*/

    //Set the TGraphErrors
    TGraphErrors* graph = new TGraphErrors(nPoints, xValues.data(), yValues.data(), xErrors.data(), yErrors.data());
    basicGraphFormatting(graph);

    //Create canvas.
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    basicCanvasFormatting(c);
    c->SetLeftMargin(0.14);

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
    frame->GetYaxis()->SetTitle("#Delta p^{rel}_{T,PbPb} - #Delta p^{rel}_{T,ppRef}");
    frame->GetYaxis()->CenterTitle(true);
    frame->GetYaxis()->SetTitleOffset(1.5);

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
    //graph->GetYaxis()->SetTitle("#Delta p_{T} = #LT #Delta p_{T}#GT_{PbPb} - #LT #Delta p_{T}#GT_{ppRef}");
    //graph->GetYaxis()->CenterTitle(true);
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
    drawLatexText("p_{T}^{#mu} > 20 GeV, |#eta^{#mu}| < 2.4", 0.2, 0.75, 0.03);
    drawLatexText("60 < M_{#mu #mu} < 120 GeV", 0.2, 0.7, 0.03);

    c->Update();
    std::string outputName = "DeltaPt_vs_Centrality" + plot_extension;
    c->SaveAs(outputName.c_str());

    delete frame;
    delete line;
    delete graph;
    delete c;
}

void MakeNormalizedDistPtRelDiff_PbPb_vs_ppRef_CentBin(TFile* inputFile, double lowCent, double highCent){

    //Get directories
    TDirectory *PbPb_dir = inputFile->GetDirectory("PbPb2023_2024_Data");
    TDirectory *ppRef_dir = inputFile->GetDirectory("ppRef2024_Data");

    std::string histName = "h2D_muonPtRelDiff_Cent";
    std::string histNameRef = "h1D_muonPtRelDiff";

    //Get original TH2 histograms
    TH2D* h_PbPb_original = dynamic_cast<TH2D*>(PbPb_dir->Get(histName.c_str()));
    TH1D* h_ppRef_original = dynamic_cast<TH1D*>(ppRef_dir->Get(histNameRef.c_str()));

        if (!h_PbPb_original || !h_ppRef_original) {//Just checking if everything was found.
            std::cerr << "Error: Could not find the histogram "
                    << histName << " in the input file."
                    << std::endl;
            return;
        }

    //Get centrality range.
    int binLow = h_PbPb_original->GetXaxis()->FindBin(lowCent + delta);
    int binHigh = h_PbPb_original->GetXaxis()->FindBin(highCent - delta);

    std::string centString = std::to_string(static_cast<int>(lowCent)) + "-" + std::to_string(static_cast<int>(highCent));

    //Project the 2D histograms onto the Y-axis for the specified centrality range.
    TH1D* h_PbPb = h_PbPb_original->ProjectionY(("h_PbPb_" + centString).c_str(), binLow, binHigh);
    TH1D* h_ppRef = dynamic_cast<TH1D*>(h_ppRef_original->Clone("h_ppRef_Comparison"));

    h_PbPb->SetDirectory(nullptr);
    h_ppRef->SetDirectory(nullptr);

    //Statistics calculation
    double PbPb_mean = h_PbPb->GetMean();
    double PbPb_mean_error = h_PbPb->GetMeanError();

    double ppRef_mean = h_ppRef->GetMean();
    double ppRef_mean_error = h_ppRef->GetMeanError();
    
    //**Get DeltaPt and its error for this centrality bin and store it in corresponding vector of pairs.**
    double DeltaPt = PbPb_mean - ppRef_mean;
    double DeltaPtError = std::sqrt(std::pow(PbPb_mean_error, 2) + std::pow(ppRef_mean_error, 2));
    DeltaPtAndError.push_back(std::make_pair(DeltaPt, DeltaPtError));


    //Extra statistics calculation to compare with 'from tree'. I want to specially check the skewness.
    //The mean and mean error are already being calculated.
    //All thats left is to calculate the variance/std dev and skewness.
    //From PbPb
    double PbPb_stddev = h_PbPb->GetStdDev();
    double PbPb_variance = PbPb_stddev * PbPb_stddev;
    double PbPb_skewness = h_PbPb->GetSkewness();

    //From ppRef
    double ppRef_stddev = h_ppRef->GetStdDev();
    double ppRef_variance = ppRef_stddev * ppRef_stddev;
    double ppRef_skewness = h_ppRef->GetSkewness();

    //For the final observable
    double DeltaPtVariance = DeltaPtError * DeltaPtError;
    double DeltaSkewness = PbPb_skewness - ppRef_skewness;

    //Print statistics to terminal for comparison with FROMTREE.
    std::cout << std::scientific << std::setprecision(10);

    std::cout << "\n";
    std::cout << "============================================================\n";
    std::cout << "Centrality " << centString << "% - FROM HISTOGRAM\n";
    std::cout << "============================================================\n";

    std::cout << "\nPbPb:\n";
    std::cout << "Mean        = " << PbPb_mean       << "\n";
    std::cout << "Mean error  = " << PbPb_mean_error << "\n";
    std::cout << "Variance    = " << PbPb_variance   << "\n";
    std::cout << "Std dev     = " << PbPb_stddev     << "\n";
    std::cout << "Skewness    = " << PbPb_skewness   << "\n";

    std::cout << "\nppRef:\n";
    std::cout << "Mean        = " << ppRef_mean       << "\n";
    std::cout << "Mean error  = " << ppRef_mean_error << "\n";
    std::cout << "Variance    = " << ppRef_variance   << "\n";
    std::cout << "Std dev     = " << ppRef_stddev     << "\n";
    std::cout << "Skewness    = " << ppRef_skewness   << "\n";

    std::cout << "\nFinal observable:\n";
    std::cout << "Delta mean  = " << DeltaPt         << "\n";
    std::cout << "Mean error  = " << DeltaPtError    << "\n";
    std::cout << "Var(delta)  = " << DeltaPtVariance << "\n";
    std::cout << "Delta skew. = " << DeltaSkewness   << "\n";

    //Begin normalization and stuff
    h_PbPb->Scale(1.0 / h_PbPb->Integral());
    h_ppRef->Scale(1.0 / h_ppRef->Integral());

    //
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    basicCanvasFormatting(c);

    basicHistFormatting(h_PbPb);
    basicHistFormatting(h_ppRef);
    
    h_PbPb->SetMarkerStyle(21);
    h_PbPb->SetMarkerSize(0.8);
    h_PbPb->SetMarkerColor(kRed);
    h_PbPb->SetLineColor(kRed);
    
    h_ppRef->SetMarkerStyle(25);
    h_ppRef->SetMarkerSize(0.8);
    h_ppRef->SetMarkerColor(kBlack);
    h_ppRef->SetLineColor(kBlack);

    h_PbPb->GetYaxis()->SetTitle("Normalized Entries");
    h_PbPb->GetXaxis()->SetTitle("#Delta p^{rel}_{T}");
    
    h_PbPb->Draw("P");
    h_ppRef->Draw("P SAME");
    
    TLegend *leg = new TLegend(0.74, 0.76, 0.94, 0.86);
    basicLegendFormatting(leg);
    leg->AddEntry(h_PbPb, "PbPb2024", "p");
    leg->AddEntry(h_ppRef, "ppRef2024", "p");
    leg->Draw();

    drawLatexText("#bf{CMS}", 0.12, 0.93, 0.042);
    drawLatexText("#it{Work in Progress}", 0.2, 0.93, 0.033);
    drawLatexText("PbPb 2024, ppRef 2024 (5.36 TeV)", 0.6, 0.93, 0.033);

    //Plot specifications
    drawLatexText("p_{T}^{#mu} > 20 GeV, |#eta^{#mu}| < 2.4", 0.17, 0.6, 0.03);
    drawLatexText("Cent." + centString, 0.17, 0.55, 0.03);

    //Add statistics to the plot
    std::string PbPb_mean_text = "#LT #Delta p^{rel}_{T} #GT_{PbPb} = " + std::to_string(PbPb_mean) + " #pm " + std::to_string(PbPb_mean_error);
    std::string ppRef_mean_text = "#LT #Delta p^{rel}_{T} #GT_{ppRef} = " + std::to_string(ppRef_mean) + " #pm " + std::to_string(ppRef_mean_error);
    drawLatexText(ppRef_mean_text, 0.17, 0.78, 0.03);
    drawLatexText(PbPb_mean_text, 0.17, 0.83, 0.03);

    c->Update();
    std::string outputName = "Normalized_" + histName + centString + "_PbPb_vs_ppRef" + plot_extension;
    c->SaveAs(outputName.c_str());

    delete leg;
    delete h_PbPb;
    delete h_ppRef;
    delete c;
}


void MakeNormalizedDistPtRelDiff_PbPb_vs_ppRef_0_100(TFile* inputFile){

    //Get directories
    TDirectory *PbPb_dir = inputFile->GetDirectory("PbPb2023_2024_Data");
    TDirectory *ppRef_dir = inputFile->GetDirectory("ppRef2024_Data");

    std::string histName = "h1D_muonPtRelDiff";

    //Get original histograms
    TH1D* h_PbPb_original = dynamic_cast<TH1D*>(PbPb_dir->Get(histName.c_str()));
    TH1D* h_ppRef_original = dynamic_cast<TH1D*>(ppRef_dir->Get(histName.c_str()));

        if (!h_PbPb_original || !h_ppRef_original) {//Just checking if everything was found.
            std::cerr << "Error: Could not find the histogram "
                    << histName << " in the input file."
                    << std::endl;
            return;
        }

    //Get histogram clones to manipulate.
    TH1D* h_PbPb = dynamic_cast<TH1D*>(h_PbPb_original->Clone(("h_PbPb_" + histName).c_str()));
    TH1D* h_ppRef = dynamic_cast<TH1D*>(h_ppRef_original->Clone(("h_ppRef_" + histName).c_str()));
    h_PbPb->SetDirectory(nullptr);
    h_ppRef->SetDirectory(nullptr);

    
    //Statistics calculation
    //Mean and error
    double PbPb_mean = h_PbPb->GetMean();
    double PbPb_mean_error = h_PbPb->GetMeanError();
    double ppRef_mean = h_ppRef->GetMean();
    double ppRef_mean_error = h_ppRef->GetMeanError();

    //Std dev and error
    double PbPb_stddev = h_PbPb->GetStdDev();
    double PbPb_stddev_error = h_PbPb->GetStdDevError();
    double PbPb_variance = PbPb_stddev * PbPb_stddev;
    double ppRef_stddev = h_ppRef->GetStdDev();
    double ppRef_variance = ppRef_stddev * ppRef_stddev;
    double ppRef_stddev_error = h_ppRef->GetStdDevError();

    //Skewness and error
    double PbPb_skewness = h_PbPb->GetSkewness();
    double PbPb_skewness_error = h_PbPb->GetSkewness(11);
    double ppRef_skewness = h_ppRef->GetSkewness();
    double ppRef_skewness_error = h_ppRef->GetSkewness(11);

    double DeltaPt = PbPb_mean - ppRef_mean;
    double DeltaPtError = std::sqrt(std::pow(PbPb_mean_error, 2) + std::pow(ppRef_mean_error, 2));
    double DeltaPtVariance = DeltaPtError * DeltaPtError;
    double DeltaSkewness = PbPb_skewness - ppRef_skewness;
    
    //Print statistics to terminal for comparison with FROMTREE.
    std::cout << std::scientific << std::setprecision(10);

    std::cout << "\n";
    std::cout << "============================================================\n";
    std::cout << "Centrality " << "0-100%" << "% - FROM HISTOGRAM\n";
    std::cout << "============================================================\n";

    std::cout << "\nPbPb:\n";
    std::cout << "Mean        = " << PbPb_mean       << "\n";
    std::cout << "Mean error  = " << PbPb_mean_error << "\n";
    std::cout << "Variance    = " << PbPb_variance   << "\n";
    std::cout << "Std dev     = " << PbPb_stddev     << "\n";
    std::cout << "Skewness    = " << PbPb_skewness   << "\n";

    std::cout << "\nppRef:\n";
    std::cout << "Mean        = " << ppRef_mean       << "\n";
    std::cout << "Mean error  = " << ppRef_mean_error << "\n";
    std::cout << "Variance    = " << ppRef_variance   << "\n";
    std::cout << "Std dev     = " << ppRef_stddev     << "\n";
    std::cout << "Skewness    = " << ppRef_skewness   << "\n";

    std::cout << "\nFinal observable:\n";
    std::cout << "Delta mean  = " << DeltaPt         << "\n";
    std::cout << "Mean error  = " << DeltaPtError    << "\n";
    std::cout << "Var(delta)  = " << DeltaPtVariance << "\n";
    std::cout << "Delta skew. = " << DeltaSkewness   << "\n";

    //Begin normalization and stuff
    h_PbPb->Scale(1.0 / h_PbPb->Integral());
    h_ppRef->Scale(1.0 / h_ppRef->Integral());

    TCanvas* c = new TCanvas("c", "c", 800, 600);
    basicCanvasFormatting(c);

    basicHistFormatting(h_PbPb);
    basicHistFormatting(h_ppRef);
    
    h_PbPb->SetMarkerStyle(21);
    h_PbPb->SetMarkerSize(0.8);
    h_PbPb->SetMarkerColor(kRed);
    h_PbPb->SetLineColor(kRed);
    
    h_ppRef->SetMarkerStyle(25);
    h_ppRef->SetMarkerSize(0.8);
    h_ppRef->SetMarkerColor(kBlack);
    h_ppRef->SetLineColor(kBlack);

    h_PbPb->GetYaxis()->SetTitle("Normalized Entries");
    h_PbPb->GetXaxis()->SetTitle("#Delta p^{rel}_{T}");
    
    h_PbPb->Draw("P");
    h_ppRef->Draw("P SAME");
    
    TLegend *leg = new TLegend(0.74, 0.76, 0.94, 0.86);
    basicLegendFormatting(leg);
    leg->AddEntry(h_PbPb, "PbPb2024", "p");
    leg->AddEntry(h_ppRef, "ppRef2024", "p");
    leg->Draw();

    drawLatexText("#bf{CMS}", 0.12, 0.93, 0.042);
    drawLatexText("#it{Work in Progress}", 0.2, 0.93, 0.033);
    drawLatexText("PbPb 2024, ppRef 2024 (5.36 TeV)", 0.6, 0.93, 0.033);

    //Plot specifications
    drawLatexText("p_{T}^{#mu} > 20 GeV, |#eta^{#mu}| < 2.4", 0.17, 0.6, 0.03);
    drawLatexText("Cent. 0-100%", 0.17, 0.55, 0.03);

    //Add statistics to the plot
    std::string PbPb_mean_text = "#LT #Delta p^{rel}_{T} #GT_{PbPb} = " + std::to_string(PbPb_mean) + " #pm " + std::to_string(PbPb_mean_error);
    std::string ppRef_mean_text = "#LT #Delta p^{rel}_{T} #GT_{ppRef} = " + std::to_string(ppRef_mean) + " #pm " + std::to_string(ppRef_mean_error);
    drawLatexText(ppRef_mean_text, 0.17, 0.78, 0.03);
    drawLatexText(PbPb_mean_text, 0.17, 0.83, 0.03);

    c->Update();
    std::string outputName = "Normalized_Joined_" + histName + "_PbPb_vs_ppRef" + plot_extension;
    c->SaveAs(outputName.c_str());

    delete leg;
    delete h_PbPb;
    delete h_ppRef;
    delete c;
}