//From mySelectedData.root

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
std::string plot_extension = ".png"; // ".png" for regular development and ".pdf" for final quality plots
std::string dataSamplesUsed = "PbPb 2023+2024, ppRef 2024 (5.36 TeV)";


// ##############################################################################
// ##############################################################################

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


void MakeTProfile_PbPb(TFile* inputFile, double lowCent, double highCent);


//---Main()
void MakeTProfile_pTMuPlMuMi(){

    TFile* inputFile = new TFile("mySelectedData.root", "READ");

    for (const auto& bin : CentralityBinsSet) {
        MakeTProfile_PbPb(inputFile, bin.first, bin.second);
    }
    inputFile->Close();
}

void MakeTProfile_PbPb(TFile* inputFile, double lowCent, double highCent){

    std::string centString = std::to_string(static_cast<int>(lowCent)) + "-" + std::to_string(static_cast<int>(highCent));

    TDirectory* PbPb_dir = inputFile->GetDirectory("PbPb2023_2024_Data");
    std::string histName = "h3D_PtMuPl_PtMuMi_Cent";

    //Get original TH3 histogram
    TH3D* h_PbPb_original = dynamic_cast<TH3D*>(PbPb_dir->Get(histName.c_str()));

    //Clone the histogram to manipulate.
    TH3D* h_PbPb = dynamic_cast<TH3D*>(h_PbPb_original->Clone(("h_PbPb_" + centString).c_str()));

    /*
    z-axis: Centrality
    y-axis: PtMuMi
    x-axis: PtMuPl
    */

    //Get centrality range and project.
    int binLow = h_PbPb_original->GetXaxis()->FindBin(lowCent + delta);
    int binHigh = h_PbPb_original->GetXaxis()->FindBin(highCent - delta);
    TH2D* h_PbPb_2D = h_PbPb->ProjectionY("h_PbPb_2D", binLow, binHigh);

    //Now we already have a 2D histogram with PtMuPl on the x-axis and PtMuMi on the y-axis.
    //Now we can create a TProfile to get the mean of PtMuMi for each bin of PtMuPl.
    TProfile* h_PbPb_TProfile = h_PbPb_2D->ProfileX("h_PbPb_TProfile");

    h_PbPb->SetDirectory(nullptr);
    h_ppRef->SetDirectory(nullptr);

}