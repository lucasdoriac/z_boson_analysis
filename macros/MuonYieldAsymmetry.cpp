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
double fitMin = 40.; //Minimum x-value for the fit range of the asymmetry histogram.
double fitMax = 60.; //Maximum x-value for the fit range of the asymmetry histogram.

// ##############################################################################
// ##############################################################################


//Set of centrality bins for PbPb2024 data. We can decide to change the centrality bins later if we want to.
/*std::vector<std::pair<double, double>> CentralityBinsSet = {
    {0., 10.},
    {10., 20.},
    {20., 30.},
    {30., 100.}
};*/

//Second proposed set of centrality bins for PbPb2024 data.
std::vector<std::pair<double, double>> CentralityBinsSet = {
    {0., 10.},
    {10., 30.},
    {30., 50.},
    {50., 100.}
};

//Vector of histograms
std::vector<TH1D*> AsymmetryHists;


//---Function declarations
void MakeAsymmetryHist_PbPb(TFile* inputFile, double lowCent, double highCent);
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

    //Now we can plot the histograms in AsymmetryHists vector in different ways.
    //The one i'm choosing for the current version of this macro is to plot ONLY the constant b of the FIT
    //that will be made on the asymmetry histogram.
    PlotAsymmetry_vs_Centrality();

    inputFile->Close();
}

void PlotAsymmetry_vs_Centrality(){

    //We have different sections for this function.

    //1- Make the FIT for each centrality bin.
    //Vectors to store the mean asymmetry and its error for each centrality bin.
    std::vector<double> meanAsymmetryValues;
    std::vector<double> meanAsymmetryErrors;

    for(size_t i = 0; i < AsymmetryHists.size(); ++i){

        TH1D* hist = AsymmetryHists[i];

        //Recover centrality bin.
        double lowCent = CentralityBinsSet[i].first;
        double highCent = CentralityBinsSet[i].second;

        TF1* constantFit = new TF1(Form("constantFit_%d_%d", (int)lowCent, (int)highCent),"pol0",fitMin, fitMax);
        hist->Fit(constantFit, "RQ"); //Quiet fit

        //Get the fitted constant value and its error.
        //This will be the final scalar representing the asymmetry for centrality bin.
        double meanDeltaA = constantFit->GetParameter(0);
        double meanError  = constantFit->GetParError(0);

        meanAsymmetryValues.push_back(meanDeltaA);
        meanAsymmetryErrors.push_back(meanError);
    }

    //2- Make the TGraphErrors for mean asymmetry vs centrality.
    int nPoints = meanAsymmetryValues.size();
    std::vector<double> xValues(nPoints);
    std::vector<double> xErrors(nPoints);
    std::vector<double> yValues(nPoints);
    std::vector<double> yErrors(nPoints);

    for(int i = 0; i < nPoints; ++i){
        xValues[i] = (CentralityBinsSet[i].first + CentralityBinsSet[i].second) / 2.0; //Centrality bin center
        xErrors[i] = (CentralityBinsSet[i].second - CentralityBinsSet[i].first) / 2.0; //Half-width of the centrality bin
        yValues[i] = meanAsymmetryValues[i];
        yErrors[i] = 0.;

        std::cout << "> Centrality bin: " << CentralityBinsSet[i].first << " - " << CentralityBinsSet[i].second
                  << ", Mean Asymmetry: " << yValues[i] << " ± " << yErrors[i] << std::endl;
    }


    TGraphErrors* graph = new TGraphErrors(nPoints, xValues.data(), yValues.data(), xErrors.data(), yErrors.data());
    basicGraphFormatting(graph);

    TCanvas* c = new TCanvas("c", "Muon Yield Asymmetry vs Centrality", 800, 600);
    basicCanvasFormatting(c);

    graph->SetMarkerStyle(21);
    graph->SetMarkerSize(1.2);
    graph->Draw("AP");

    c->Update();
    std::string outputName = "MuonYieldAsymmetry_vs_Centrality" + plot_extension;
    c->SaveAs(outputName.c_str());

    delete graph;
    delete c;
}

void MakeAsymmetryHist_PbPb(TFile* inputFile, double lowCent, double highCent){

    //ppRef part:
    //Get directory
    TDirectory *ppRef_dir = inputFile->GetDirectory("ppRef2024_Data");

    //In this case we don't have a TH3D histogram so we can just use the individual TH1 ones.
    //Originals
    TH1D* h1 = dynamic_cast<TH1D*>(ppRef_dir->Get("h1D_ptMuPlus"));
    TH1D* h2 = dynamic_cast<TH1D*>(ppRef_dir->Get("h1D_ptMuMinus"));

    //Clone
    TH1D* h1D_PtMuPl = dynamic_cast<TH1D*>(h1->Clone("h1D_PtMuPl"));
    TH1D* h1D_PtMuMi = dynamic_cast<TH1D*>(h2->Clone("h1D_PtMuMi"));

    //Now fill asymmetry histogram
    TH1D* h1D_Asymmetry_ppRef = new TH1D("h1D_Asymmetry_ppRef", "Muon Yield Asymmetry; p_{T} [GeV/c]; A(p_{T})", 100, 0., 100.);
    for(int i = 1; i <= h1D_PtMuPl->GetNbinsX(); ++i){
        double N_plus = h1D_PtMuPl->GetBinContent(i);
        double N_minus = h1D_PtMuMi->GetBinContent(i);
        if(N_plus + N_minus == 0){
            h1D_Asymmetry_ppRef->SetBinContent(i, 0.); //Avoid division by zero
            continue;
        }
        double asymmetry = (N_plus - N_minus) / (N_plus + N_minus);
        h1D_Asymmetry_ppRef->SetBinContent(i, asymmetry);
    }
    //End of ppRef part.

    //PbPb part:
    //Get directory
    TDirectory *PbPb_dir = inputFile->GetDirectory("PbPb2024_Data");
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
    for(int i = 1; i <= h1D_PtMuPl->GetNbinsX(); ++i){
        double N_plus = h1D_PtMuPl->GetBinContent(i);
        double N_minus = h1D_PtMuMi->GetBinContent(i);
        if(N_plus + N_minus == 0){
            h1D_AsymmetryPbPb->SetBinContent(i, 0.); //Avoid division by zero
            continue;
        }
        double asymmetry = (N_plus - N_minus) / (N_plus + N_minus);
        h1D_AsymmetryPbPb->SetBinContent(i, asymmetry);
    }
    //End of PbPb part.

    //Make relative asymmetry histogram
    const TString histName = Form("h1D_RelativeAsymmetry_cent_%g_%g", lowCent, highCent);

    TH1D* h1D_RelativeAsymmetry = new TH1D(histName.Data(),
        "Muon Yield Asymmetry Difference; p_{T} [GeV/c]; #Delta A(p_{T})",
        100, 0., 100.
    );
    
    for(int i = 1; i <= h1D_AsymmetryPbPb->GetNbinsX(); ++i){
        double asymmetry_ppRef = h1D_Asymmetry_ppRef->GetBinContent(i);
        double asymmetry_PbPb = h1D_AsymmetryPbPb->GetBinContent(i);
        double relative_asymmetry = (asymmetry_PbPb - asymmetry_ppRef);
        h1D_RelativeAsymmetry->SetBinContent(i, relative_asymmetry);
    }

    //Store the relative asymmetry histogram in the vector.
    h1D_RelativeAsymmetry->SetDirectory(nullptr);
    AsymmetryHists.push_back(h1D_RelativeAsymmetry);

    delete h1D_Asymmetry_ppRef;
    delete h1D_AsymmetryPbPb;
    delete h1D_PtMuPl;
    delete h1D_PtMuMi;
    delete h2D_PtMuPl_PtMuMi;
    delete h3D_PtMuPl_PtMuMi_Cent;
}