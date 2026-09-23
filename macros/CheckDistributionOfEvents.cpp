/*
Check distribution of events in each centrality set.
*/

//---Libraries
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TH1.h>
#include <TString.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <array>
#include <algorithm>
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


std::vector<std::array<double,3>> GetEventsOnCBin(TFile* inputFile, std::vector<std::pair<double, double>> cBinsSet);
void MakeTGraph(std::vector<std::array<double,3>> fromSet1,
                std::vector<std::array<double,3>> fromSet2);
TGraph* FormatGraph(TGraph* graph, Color_t color);


//---Main()
void CheckDistributionOfEvents(){

    gROOT->SetBatch(kTRUE);
    TFile *inputFile = new TFile("mySelectedData.root", "READ");

    auto fromSet1 = GetEventsOnCBin(inputFile, cBinsSet1);
    auto fromSet2 = GetEventsOnCBin(inputFile, cBinsSet2);

    MakeTGraph(fromSet1, fromSet2);

    inputFile->Close();
}


void MakeTGraph(std::vector<std::array<double,3>> fromSet1,
                std::vector<std::array<double,3>> fromSet2){

    TGraph *graph_set1;
    TGraph *graph_set2;

    int nPoints = fromSet1.size();
    if(nPoints != cBinsSet1.size()){
            std::cerr << "Error: Number of points in PeakDifferenceAndError does not match number of centrality bins." << std::endl;
            return;
    }

    std::vector<double> xValues(nPoints);
    std::vector<double> yValues(nPoints);

    for(int i = 0; i < nPoints; ++i){
        xValues[i] = i + 1;
        yValues[i] = fromSet1[i][2]; //Number of events in this centrality range
    }
    // Create the TGraph for Set1.
    TGraph *graph = new TGraph(nPoints, xValues.data(), yValues.data());
    graph_set1 = FormatGraph(graph, kBlue);

    for(int i = 0; i < nPoints; ++i){
        xValues[i] = i + 1;
        yValues[i] = fromSet2[i][2]; //Number of events in this centrality range
    }
    // Create the TGraph for Set2.
    graph = new TGraph(nPoints, xValues.data(), yValues.data());
    graph_set2 = FormatGraph(graph, kRed);

    //Create canvas.
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    basicCanvasFormatting(c);

    //Frame TH1 helper to set the x-axis labels for centrality bins.
    TH1D* frame = new TH1D("frame","",nPoints,0.5,nPoints + 0.5);
    basicHistFormatting(frame);
    for(int i = 0; i < nPoints; ++i){
        frame->GetXaxis()->SetBinLabel(i+1, Form("%d", i+1));
    }

    double yMax = 0.;
    for(const auto& value : fromSet1) yMax = std::max(yMax, value[2]);
    for(const auto& value : fromSet2) yMax = std::max(yMax, value[2]);
    frame->SetMinimum(0.);
    frame->SetMaximum(1.20 * yMax);
    //

    frame->GetXaxis()->SetTickLength(0.0);
    frame->GetXaxis()->SetTitle("cBin index");
    frame->GetYaxis()->SetTitle("Candidates");
    frame->GetYaxis()->SetTitleOffset(1.4);

    //Draw only the axis frame.
    frame->Draw("AXIS");


    //Draw graphs
    graph_set1->Draw("LP SAME");
    graph_set2->Draw("LP SAME");

    TLegend* leg = new TLegend(0.62, 0.72, 0.88, 0.84);
    basicLegendFormatting(leg);

    leg->AddEntry(graph_set1, "Centrality set 1", "lp");
    leg->AddEntry(graph_set2, "Centrality set 2", "lp");

    leg->Draw();

    drawLatexText("#bf{CMS}", 0.12, 0.93, 0.042);
    drawLatexText("#it{Internal}", 0.2, 0.93, 0.033);
    drawLatexText(dataSamplesUsed.c_str(), 0.6, 0.93, 0.033);
    
    //Plot specifications
    drawLatexText("p_{T}^{#mu} > 20 GeV, |#eta^{#mu}| < 2.4", 0.22, 0.3, 0.03);
    drawLatexText("60 < M_{#mu#mu} < 120 GeV", 0.22, 0.25, 0.03);

    c->Update();
    std::string outputName = "DistributionOfEvents_vs_Centrality" + plot_extension;
    c->SaveAs(outputName.c_str());

    delete frame;
    delete graph_set1;
    delete graph_set2;
    delete leg;
    delete c;
}


std::vector<std::array<double,3>> GetEventsOnCBin(TFile* inputFile, std::vector<std::pair<double, double>> cBinsSet){

    //Get directory
    TDirectory *dir = (TDirectory*)inputFile->Get(whichDataset.c_str());

    //Get n of dimuons vs centrality hist
    TH1D *h1D_centrality = (TH1D*)dir->Get("h1D_centrality");

    std::vector<std::array<double,3>> values;

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


TGraph* FormatGraph(TGraph* graph, Color_t color){

    //Basic formatting
    graph->SetStats(0);
    graph->SetTitle("");
    graph->GetXaxis()->CenterTitle(false);
    graph->GetYaxis()->CenterTitle(false);
    graph->GetXaxis()->SetTitleOffset(1.1);
    graph->GetYaxis()->SetTitleOffset(1.2);
    graph->GetXaxis()->SetTitleFont(42);
    graph->GetYaxis()->SetTitleFont(42);
    graph->GetXaxis()->SetLabelFont(42);
    graph->GetYaxis()->SetLabelFont(42);
    graph->GetXaxis()->SetTitleSize(0.042);
    graph->GetYaxis()->SetTitleSize(0.042);

    //Style graph.
    graph->SetMarkerStyle(21);
    graph->SetMarkerSize(0.9);
    graph->SetMarkerColorAlpha(color+1, 1.);
    graph->SetLineColorAlpha(color-7, 0.8);
    graph->SetLineWidth(2);
    return graph;
}