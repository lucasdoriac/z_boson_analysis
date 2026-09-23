/*
Check distribution of events in each dataset. One canvas.
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
double delta = 1e-6;


// ##############################################################################
// ##############################################################################


std::vector<std::pair<double, double>> cBinsSet1 = {
    {0., 10.},
    {10., 20.},
    {20., 30.},
    {30., 100.},
};


std::map<std::string, std::vector<double>> results;


std::vector<double> GetEventsOnCBin(TFile* inputFile, std::vector<std::pair<double, double>> cBinsSet, std::string whichDataset);
void MakeTGraph(const std::map<std::string, std::vector<double>>& results);
TGraph* FormatGraph(TGraph* graph, Color_t color);


//---Main()
void CheckEvent_vs_Cent(){

    gROOT->SetBatch(kTRUE);
    TFile *inputFile = new TFile("mySelectedData.root", "READ");

    results.clear();
    results["PbPb2023"] = GetEventsOnCBin(inputFile, cBinsSet1, "PbPb2023_Data");
    results["PbPb2024"] = GetEventsOnCBin(inputFile, cBinsSet1, "PbPb2024_Data");
    results["PbPb2025"] = GetEventsOnCBin(inputFile, cBinsSet1, "PbPb2025_Data");
    results["PbPb2026"] = GetEventsOnCBin(inputFile, cBinsSet1, "PbPb2026_Data");

    MakeTGraph(results);

    inputFile->Close();
}


void MakeTGraph(const std::map<std::string, std::vector<double>>& results){

    const int nPoints = cBinsSet1.size();

    std::vector<Color_t> colors = {kBlue + 1, kRed + 1, kGreen + 2, kYellow + 1};

    double yMax = 0.0;
    for (const auto& [dataset, counts] : results) {
        if (static_cast<int>(counts.size()) != nPoints) {
            std::cerr << "Número incorreto de bins em " << dataset << '\n';
            return;
        }
        for (double count : counts)
            yMax = std::max(yMax, count);
    }

    TCanvas canvas("c_events_cent", "Events vs centrality", 800, 600);
    basicCanvasFormatting(&canvas);

    TH1D* frame = new TH1D("frame_events_cent", "", nPoints, 0.5, nPoints + 0.5);
    basicHistFormatting(frame);

    for (int i = 0; i < nPoints; ++i) {
        frame->GetXaxis()->SetBinLabel(i + 1, Form("%.0f-%.0f%%", cBinsSet1[i].first, cBinsSet1[i].second));
    }

    frame->SetMinimum(0.0);
    frame->SetMaximum(yMax > 0.0 ? 1.2 * yMax : 1.0);
    frame->GetXaxis()->SetTitle("Centralidade");
    frame->GetYaxis()->SetTitle("Candidatos");
    frame->GetYaxis()->SetTitleOffset(1.4);
    frame->Draw();

    TLegend* leg = new TLegend(0.62, 0.68, 0.88, 0.88);
    basicLegendFormatting(leg);

    std::vector<std::unique_ptr<TGraph>> graphs;
    int datasetIndex = 0;

    for (const auto& [dataset, counts] : results) {
        std::vector<double> x(nPoints);
        for (int i = 0; i < nPoints; ++i)
            x[i] = i + 1;

        auto graph = std::make_unique<TGraph>(
            nPoints, x.data(), counts.data()
        );

        FormatGraph(graph.get(), colors[datasetIndex % colors.size()]);
        graph->Draw("LP SAME");
        legend.AddEntry(graph.get(), dataset.c_str(), "lp");

        graphs.push_back(std::move(graph));
        ++datasetIndex;
    }

    legend.Draw();
    canvas.SaveAs(
        ("DistributionOfEvents_vs_Centrality" + plot_extension).c_str()
    );
}


std::vector<double> GetEventsOnCBin(TFile* inputFile, std::vector<std::pair<double, double>> cBinsSet, std::string whichDataset){

    TDirectory* dir = inputFile->GetDirectory(whichDataset.c_str());

    TH1D *h1D_centrality = (TH1D*)dir->Get("h1D_centrality");

    std::vector<double> counts;

    //Loop over the centrality bins
    for (size_t i = 0; i < cBinsSet.size(); i++){
        double cMin = cBinsSet[i].first;
        double cMax = cBinsSet[i].second;

        //Get the bin numbers corresponding to the centrality range
        int binMin = h1D_centrality->FindBin(cMin + delta);
        int binMax = h1D_centrality->FindBin(cMax - delta);

        //Calculate the number of events in this centrality range
        double nEventsInRange = h1D_centrality->Integral(binMin, binMax);

        counts.push_back(nEventsInRange);
    }

    return counts;
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