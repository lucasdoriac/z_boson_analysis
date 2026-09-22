/*
A = (N_{mu+} - N_{mu-}) / (N_{mu+} + N_{mu-}).
From 'mySelectedData.root'.
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
std::string dataSamplesUsed = "PbPb 2023+2024, ppRef 2024 (5.36 TeV)";
double delta = 1e-6; //Small value to avoid binning issues when projecting histograms.

double rho = 1.0; //Correlation coefficient for error propagation in asymmetry calculation.
double fitMin = 40.; //Minimum x-value for the fit range of the asymmetry histogram.
double fitMax = 65.; //Maximum x-value for the fit range of the asymmetry histogram.


// ##############################################################################
// ##############################################################################



//Set of centrality bins for PbPb2024 data. We can decide to change the centrality bins later if we want to.
std::vector<std::pair<double, double>> CentralityBinsSet = {
    {0., 10.},
    {10., 20.},
    {20., 30.},
    {30., 100.},
    {0., 100.}
};

//Second proposed set of centrality bins for PbPb2024 data.
/*std::vector<std::pair<double, double>> CentralityBinsSet = {
    {0., 10.},
    {10., 30.},
    {30., 50.},
    {50., 100.},
    {0., 100.}
};*/

//Vector of histograms
std::vector<TH1D*> AsymmetryHists;

//Vectors to store the mean asymmetry and its error for each centrality bin.
std::vector<double> meanAsymmetryValues;
std::vector<double> meanAsymmetryErrors;


//---Function declarations
void MakeAsymmetryHist_PbPb(TFile* inputFile, double lowCent, double highCent);
void PlotAsymmetryHists();
void PlotAsymmetry_vs_Centrality();

//---Main()
void MuonYieldAsymmetry(){

    //Clear vector
    AsymmetryHists.clear();

    gROOT->SetBatch(kTRUE);
    TFile* inputFile = new TFile("mySelectedData.root", "READ");


    for(const auto& centralityBin : CentralityBinsSet){
        double lowCent = centralityBin.first;
        double highCent = centralityBin.second;

        MakeAsymmetryHist_PbPb(inputFile, lowCent, highCent);
    }

    PlotAsymmetryHists();
    PlotAsymmetry_vs_Centrality();

    inputFile->Close();
}

void PlotAsymmetryHists(){

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-8);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);

    for(size_t i = 0; i < AsymmetryHists.size(); ++i){

        //Since the fit is showing to be too sensitive to the value of rho, i'll use the simple mean for now.
        /*std::string fitCentString = Form("pol0Fit_%.0f-%.0f", CentralityBinsSet[i].first, CentralityBinsSet[i].second);
        //Part one: Make fits and fill vectors with fit results.
        TF1* fitFunc = new TF1(fitCentString.c_str(), "[0]", fitMin, fitMax);
        AsymmetryHists[i]->Fit(fitFunc, "NO R EX0 Q", "", fitMin, fitMax);
        AsymmetryHists[i]->Fit(fitFunc, "NO R EX0 Q", "", fitMin, fitMax);
        TFitResultPtr fitResult = AsymmetryHists[i]->Fit(fitFunc, "NO R EX0 M S", "", fitMin, fitMax);

        //Get results to terminal
        double chi2 = fitResult->Chi2();
        int ndf = fitResult->Ndf();
        double pValue = TMath::Prob(chi2, ndf);
        std::cout << "Centrality bin: " << CentralityBinsSet[i].first << "-" << CentralityBinsSet[i].second << "%, chi2: " << chi2
                    << ", ndf: " << ndf << ", p-value: " << pValue << std::endl;

        double meanDeltaA = fitResult->Parameter(0);
        double meanError = fitResult->ParError(0);*/

        //Simple mean of y-coordinates of the histogram within the fit range.
        //Find bin of fitMin and fitMax
        int binMin = AsymmetryHists[i]->FindBin(fitMin + delta);
        int binMax = AsymmetryHists[i]->FindBin(fitMax - delta);

        double y = 0;
        double error = 0;
        double n = 0;
        for(int bin = binMin; bin <= binMax; ++bin){

            //Print the x values of the bins being used for the mean calculation.
            std::cout << "Bin " << bin << " x: " << AsymmetryHists[i]->GetBinCenter(bin) << std::endl;

            y+= AsymmetryHists[i]->GetBinContent(bin);
            error+= std::pow(AsymmetryHists[i]->GetBinError(bin),2); //Variance sum.
            n++;
        }

        std::cout << "Total of bins used = " << n << std::endl;
        double meanDeltaA = y/n;
        double meanError = std::sqrt(error) / n;

        //Store results for this centrality bin in the vectors
        meanAsymmetryValues.push_back(meanDeltaA);
        meanAsymmetryErrors.push_back(meanError);

        //Part two: Plot the asymmetry histograms with fit for each centrality bin.
        TCanvas* c = new TCanvas("c", "Muon Yield Asymmetry", 800, 600);
        basicCanvasFormatting(c);
        c->SetLeftMargin(0.1);

        basicHistFormatting(AsymmetryHists[i]);
        AsymmetryHists[i]->SetMarkerStyle(20);
        AsymmetryHists[i]->SetMarkerSize(0.8);
        AsymmetryHists[i]->SetMarkerColorAlpha(kRed+1, 1.);
        AsymmetryHists[i]->SetLineColorAlpha(kRed-7, 0.8);

        AsymmetryHists[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
        AsymmetryHists[i]->GetYaxis()->SetTitle("A(p_{T})");
        AsymmetryHists[i]->GetYaxis()->CenterTitle(true);
        AsymmetryHists[i]->GetYaxis()->SetTitleOffset(1.1);
        AsymmetryHists[i]->GetXaxis()->SetRangeUser(18., 100.);

        AsymmetryHists[i]->GetXaxis()->SetRangeUser(38., 70.);
        AsymmetryHists[i]->GetYaxis()->SetRangeUser(-0.2, 0.2);


        AsymmetryHists[i]->Draw("E1");

        //Draw fit
        //fitFunc->SetLineColor(kCyan+2);
        //fitFunc->SetLineWidth(2);
        //fitFunc->SetLineStyle(4);
        //fitFunc->Draw("SAME");

        //Grey line at y=0
        TLine* line = new TLine(38., 0.0, 70., 0.0);
        line->SetLineColor(kGray);
        line->SetLineStyle(7);
        line->SetLineWidth(1);
        line->Draw("SAME");

        drawLatexText("#bf{CMS}", 0.11, 0.93, 0.042);
        drawLatexText("#it{Work in Progress}", 0.18, 0.93, 0.033);
        drawLatexText(dataSamplesUsed.c_str(), 0.55, 0.93, 0.033);

        //Plot specifications
        drawLatexText("p_{T}^{#mu} > 20 GeV, |#eta^{#mu}| < 2.4", 0.2, 0.28, 0.03);
        drawLatexText("60 < M_{#mu#mu} < 120 GeV", 0.2, 0.24, 0.03);

        c->Update();
        std::string outputName = Form("Asymmetry_hist_cent_%g_%g", CentralityBinsSet[i].first, CentralityBinsSet[i].second) + plot_extension;
        c->SaveAs(outputName.c_str());

        delete c;
    }

}


void PlotAsymmetry_vs_Centrality(){

    //Make the TGraphErrors for mean asymmetry vs centrality.
    int nPoints = meanAsymmetryValues.size();
    if(nPoints != CentralityBinsSet.size()){
            std::cerr << "Error: Number of points in PeakDifferenceAndError does not match number of centrality bins." << std::endl;
            return;
    }

    std::vector<double> xValues(nPoints);
    std::vector<double> xErrors(nPoints);
    std::vector<double> yValues(nPoints);
    std::vector<double> yErrors(nPoints);

    for(int i = 0; i < nPoints; ++i){
        xValues[i] = i + 1;
        xErrors[i] = 0.;
        yValues[i] = meanAsymmetryValues[i];
        yErrors[i] = meanAsymmetryErrors[i];
    }
    

    TGraphErrors* graph = new TGraphErrors(nPoints, xValues.data(), yValues.data(), xErrors.data(), yErrors.data());
    basicGraphFormatting(graph);


    TCanvas* c = new TCanvas("c", "Muon Yield Asymmetry vs Centrality", 800, 600);
    basicCanvasFormatting(c);
    c->SetLeftMargin(0.13);

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
    frame->GetYaxis()->SetTitle("#LT #Delta A #GT (fit)");
    frame->GetYaxis()->CenterTitle(true);
    frame->GetYaxis()->SetTitleOffset(1.4);

    //Draw only the axis frame.
    frame->Draw("AXIS");

    //Format graph.
    graph->SetMarkerStyle(21);
    graph->SetMarkerSize(0.9);
    graph->SetMarkerColorAlpha(kRed+1, 1.);
    graph->SetLineColorAlpha(kRed-7, 0.8);
    graph->SetLineWidth(2);

    //Draw
    graph->Draw("P SAME");

    //Grey line at y=0
    TLine* line = new TLine(0.5, 0.0, 0.5+nPoints, 0.0);
    line->SetLineColor(kGray);
    line->SetLineStyle(7);
    line->SetLineWidth(2);
    line->Draw();

    drawLatexText("#bf{CMS}", 0.13, 0.93, 0.042);
    drawLatexText("#it{Work in Progress}", 0.2, 0.93, 0.033);
    drawLatexText(dataSamplesUsed.c_str(), 0.55, 0.93, 0.033);
    
    //Plot specifications
    drawLatexText("p_{T}^{#mu} > 20 GeV, |#eta^{#mu}| < 2.4", 0.22, 0.22, 0.03);
    drawLatexText("60 < M_{#mu#mu} < 120 GeV", 0.22, 0.17, 0.03);

    c->Update();
    std::string outputName = "Asymmetry_vs_Centrality" + plot_extension;
    c->SaveAs(outputName.c_str());

    delete frame;
    delete line;
    delete graph;
    delete c;
}

void MakeAsymmetryHist_PbPb(TFile* inputFile, double lowCent, double highCent){

    //**ppRef:**

    //Get directory
    TDirectory *ppRef_dir = inputFile->GetDirectory("ppRef2024_Data");

    //In this case we don't have a TH3D histogram so we can just use the individual TH1 ones.
    //Originals
    TH1D* h1D_PtMuPl = dynamic_cast<TH1D*>(ppRef_dir->Get("h1D_ptMuPlus"));
    TH1D* h1D_PtMuMi = dynamic_cast<TH1D*>(ppRef_dir->Get("h1D_ptMuMinus"));

    //Now fill asymmetry histogram
    TH1D* h1D_Asymmetry_ppRef = new TH1D("h1D_Asymmetry_ppRef", "Muon Yield Asymmetry; p_{T} [GeV/c]; A(p_{T})", 100, 0., 100.);
    for(int i = 1; i <= h1D_Asymmetry_ppRef->GetNbinsX(); ++i){
        //Get bin contents for mu+ and mu-
        double N_plus = h1D_PtMuPl->GetBinContent(i);
        double N_minus = h1D_PtMuMi->GetBinContent(i);

        //Get the poisson errors
        double N_plusError = h1D_PtMuPl->GetBinError(i);
        double N_minusError = h1D_PtMuMi->GetBinError(i);

        if(N_plus + N_minus == 0){//Just to avoid division by zero
            h1D_Asymmetry_ppRef->SetBinContent(i, 0.); 
            h1D_Asymmetry_ppRef->SetBinError(i, 0.);
            continue;
        }

        double denominator = N_plus + N_minus;
        double asymmetry = (N_plus - N_minus) / denominator;

        //Error propagation for asymmetry calculation. Lara's eq. for ratio.
        //double asymmetryError = (N_plus/N_minus)*std::sqrt((N_plusError/N_plus)*(N_plusError/N_plus) + (N_minusError/N_minus)*(N_minusError/N_minus) - 2.*(N_plusError*N_plusError)/(N_plus*N_minus));
        double CovTerm = 2.*N_plus*N_minus*N_plusError*N_minusError*rho;
        double asymmetryError = 2./(denominator*denominator)*std::sqrt(N_minus*N_minus*N_plusError*N_plusError + N_plus*N_plus*N_minusError*N_minusError - CovTerm);


        h1D_Asymmetry_ppRef->SetBinContent(i, asymmetry);
        h1D_Asymmetry_ppRef->SetBinError(i, asymmetryError);
    }
    //End of ppRef part.

    //**PbPb:**
    //Get directory
    TDirectory *PbPb_dir = inputFile->GetDirectory(whichDataset.c_str());
    TH3D* h_original = dynamic_cast<TH3D*>(PbPb_dir->Get("h3D_PtMuPl_PtMuMi_Cent"));

    /*
    z-axis -> centrality
    y-axis -> pT of mu-
    x-axis -> pT of mu+
    */

    //Clone
    TH3D* h3D_PtMuPl_PtMuMi_Cent = dynamic_cast<TH3D*>(h_original->Clone("h3D_PtMuPl_PtMuMi_Cent"));

    //Select centrality range and project onto pT(mu+) vs pT(mu-) plane.
    int binLow = h3D_PtMuPl_PtMuMi_Cent->GetZaxis()->FindBin(lowCent + delta);
    int binHigh = h3D_PtMuPl_PtMuMi_Cent->GetZaxis()->FindBin(highCent - delta);

    std::cout << "> Centrality range: " << lowCent << " - " << highCent << std::endl;
    h3D_PtMuPl_PtMuMi_Cent->GetZaxis()->SetRange(binLow, binHigh);

    TH2D* h2D_PtMuPl_PtMuMi = dynamic_cast<TH2D*>(h3D_PtMuPl_PtMuMi_Cent->Project3D("yx"));

    //Project into TH1Ds to calculate asymmetry histogram later
    h1D_PtMuPl = dynamic_cast<TH1D*>(h2D_PtMuPl_PtMuMi->ProjectionX("h1D_PtMuPl"));
    h1D_PtMuMi = dynamic_cast<TH1D*>(h2D_PtMuPl_PtMuMi->ProjectionY("h1D_PtMuMi"));

    //Now fill asymmetry histogram
    TH1D* h1D_AsymmetryPbPb = new TH1D("h1D_AsymmetryPbPb", "Muon Yield Asymmetry; p_{T} [GeV/c]; A(p_{T})", 100, 0., 100.);
    for(int i = 1; i <= h1D_AsymmetryPbPb->GetNbinsX(); ++i){
        //Get bin contents for mu+ and mu-
        double N_plus = h1D_PtMuPl->GetBinContent(i);
        double N_minus = h1D_PtMuMi->GetBinContent(i);

        //Get the poisson errors
        double N_plusError = h1D_PtMuPl->GetBinError(i);
        double N_minusError = h1D_PtMuMi->GetBinError(i);

        if(N_plus + N_minus == 0){
            h1D_AsymmetryPbPb->SetBinContent(i, 0.); //Avoid division by zero
            h1D_AsymmetryPbPb->SetBinError(i, 0.);
            continue;
        }

        double denominator = N_plus + N_minus;
        double asymmetry = (N_plus - N_minus) / denominator;

        //Error propagation for asymmetry calculation. Lara's eq. for ratio.
        //double asymmetryError = (N_plus/N_minus)*std::sqrt((N_plusError/N_plus)*(N_plusError/N_plus) + (N_minusError/N_minus)*(N_minusError/N_minus) - 2.*(N_plusError*N_plusError)/(N_plus*N_minus));
        double CovTerm = 2.*N_plus*N_minus*N_plusError*N_minusError*rho;
        double asymmetryError = 2./(denominator*denominator)*std::sqrt(N_minus*N_minus*N_plusError*N_plusError + N_plus*N_plus*N_minusError*N_minusError - CovTerm);

        h1D_AsymmetryPbPb->SetBinContent(i, asymmetry);
        h1D_AsymmetryPbPb->SetBinError(i, asymmetryError);
    }
    //End of PbPb part.

    //Make relative asymmetry histogram
    std::string histName = Form("h1D_RelAsym_cent_%g_%g", lowCent, highCent);
    TH1D* h1D_RelativeAsymmetry = new TH1D(histName.c_str(), "Muon Yield Asymmetry Difference; p_{T} [GeV/c]; #Delta A(p_{T})", 100, 0., 100.);
    

    for(int i = 1; i <= h1D_AsymmetryPbPb->GetNbinsX(); ++i){

        //Set bin contents with the difference between PbPb and ppRef asymmetry histograms.
        double asymmetry_ppRef = h1D_Asymmetry_ppRef->GetBinContent(i);
        double asymmetry_PbPb = h1D_AsymmetryPbPb->GetBinContent(i);
        double relative_asymmetry = (asymmetry_PbPb - asymmetry_ppRef);
        h1D_RelativeAsymmetry->SetBinContent(i, relative_asymmetry);

        //Independent error propagation for the relative asymmetry histogram:
        double error_ppRef = h1D_Asymmetry_ppRef->GetBinError(i);
        double error_PbPb = h1D_AsymmetryPbPb->GetBinError(i);
        double relative_asymmetry_error = std::sqrt(error_ppRef*error_ppRef + error_PbPb*error_PbPb);
        h1D_RelativeAsymmetry->SetBinError(i, relative_asymmetry_error);
    }

    //Store the relative asymmetry histogram in the vector.
    AsymmetryHists.push_back(h1D_RelativeAsymmetry);

    
    delete h1D_Asymmetry_ppRef;
    delete h1D_AsymmetryPbPb;
    delete h1D_PtMuPl;
    delete h1D_PtMuMi;
    delete h2D_PtMuPl_PtMuMi;
    delete h3D_PtMuPl_PtMuMi_Cent;
    delete h_original;
}