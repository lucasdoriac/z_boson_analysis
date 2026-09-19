/*
Check distribution of events in each centrality set.
*/

//---Libraries
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TH1.h>
#include <TString.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstring>
#include <vector>
#include "../headers/basicFormatting.h"

//---Macro settings
std::string plot_extension = ".pdf"; // ".png" for regular development and ".pdf" for final quality plots
std::string whichDataset = "PbPb2023_2024_Data"; // "PbPb2023_2024_Data", "PbPb2023_Data", "PbPb2024_Data". 
std::string JointPbPb = "PbPb2023+2024"; //"PbPb2023+2024", "PbPb2023", "PbPb2024".
std::string dataSamplesUsed = "PbPb 2023+2024, ppRef 2024 (5.36 TeV)"; //"PbPb 2023+2024, ppRef 2024 (5.36 TeV)", "PbPb 2023, ppRef 2024 (5.36 TeV)", "PbPb 2024, ppRef 2024 (5.36 TeV)".
double delta = 1e-6;

// ##############################################################################
// ##############################################################################

std::vector<std::pair<double, double>> cBinsSet1 = {
    {0., 10.},
    {10., 20.},
    {20., 30.},
    {30., 100.},
};

std::vector<std::pair<double, double>> cBinsSet2 = {
    {0., 10.},
    {10., 30.},
    {30., 50.},
    {50., 100.},
};


std::vector<std::vector<double, double, double>> GetEventsOnCBin(TFile* inputFile, std::vector<std::pair<double, double>> cBinsSet);
void MakeTGraph(std::vector<std::vector<double, double, double>> fromSet1,
                std::vector<std::vector<double, double, double>> fromSet2);

//---Main()
void CheckDistributionOfEvents(){

    gROOT->SetBatch(kTRUE);
    TFile *inputFile = new TFile("mySelectedData.root", "READ");

    auto fromSet1 = GetEventsOnCBin(inputFile, cBinsSet1);
    auto fromSet2 = GetEventsOnCBin(inputFile, cBinsSet2);

    MakeTGraph(fromSet1, fromSet2);

    inputFile->Close();
}

void MakeTGraph(std::vector<std::vector<double, double, double>> fromSet1,
                std::vector<std::vector<double, double, double>> fromSet2){

    TGraph *graph1 = new TGraph(fromSet1.size());
    TGraph *graph2 = new TGraph(fromSet2.size());

    for (size_t i = 0; i < fromSet1.size(); i++){
        double cMin = fromSet1[i][0];
        double cMax = fromSet1[i][1];
        double nEvents = fromSet1[i][2];

        graph1->SetPoint(i, (cMin + cMax) / 2., nEvents);
    }

    for (size_t i = 0; i < fromSet2.size(); i++){
        double cMin = fromSet2[i][0];
        double cMax = fromSet2[i][1];
        double nEvents = fromSet2[i][2];

        graph2->SetPoint(i, (cMin + cMax) / 2., nEvents);
    }


    
}

std::vector<std::vector<double, double, double>> GetEventsOnCBin(TFile* inputFile, std::vector<std::pair<double, double>> cBinsSet){

    //Get directory
    TDirectory *dir = (TDirectory*)inputFile->Get("PbPb2023_2024_Data");

    //Get n of dimuons vs centrality hist
    TH1D *h1D_centrality = (TH1D*)dir->Get("h1D_centrality");

    std::vector<std::vector<double, double, double>> values;

    //Loop over the centrality bins
    for (size_t i = 0; i < cBinsSet.size(); i++){
        double cMin = cBinsSet[i].first;
        double cMax = cBinsSet[i].second;

        //Get the bin numbers corresponding to the centrality range
        int binMin = h1D_centrality->FindBin(cMin + delta);
        int binMax = h1D_centrality->FindBin(cMax - delta);

        //Calculate the number of events in this centrality range
        double nEventsInRange = h1D_centrality->Integral(binMin, binMax);

        values.push_back({cMin, cMax, nEventsInRange});

        //Print the results
        std::cout << "Centrality range: " << cBinsSet[i].first << "% - " << cBinsSet[i].second << "%, Number of events: " << nEventsInRange << std::endl;
    }

    return values;
}