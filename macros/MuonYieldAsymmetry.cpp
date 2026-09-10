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


//---Function declarations
void myFunction(TFile* inputFile);

//---Main()
void MuonYieldAsymmetry(){

    gROOT->SetBatch(kTRUE);

    TFile* inputFile = new TFile("mySelectedData.root", "READ");

    myFunction(inputFile);

    inputFile->Close();
}


void myFunction(TFile* inputFile){

    //Get directories
    TDirectory *PbPb_dir = inputFile->GetDirectory("PbPb2024_Data");
    TDirectory *ppRef_dir = inputFile->GetDirectory("ppRef2024_Data");

    

}