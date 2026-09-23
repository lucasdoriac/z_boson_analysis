/*
Check centrality distribution
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


TH1D* getCentFromDataset(std::string datasetName, TFile* inputFile);
void plotCent(TH1D* h1D_cent2023, TH1D* h1D_cent2024, TH1D* h1D_cent2025, TH1D* h1D_cent2026);


void Check_cent(){

    gROOT->SetBatch(kTRUE);
    TFile* inputFile = TFile::Open("mySelectedData.root", "READ");

    TH1D* h1D_cent2023 = getCentFromDataset("PbPb2023_Data", inputFile);
    TH1D* h1D_cent2024 = getCentFromDataset("PbPb2024_Data", inputFile);
    TH1D* h1D_cent2025 = getCentFromDataset("PbPb2025_Data", inputFile);
    TH1D* h1D_cent2026 = getCentFromDataset("PbPb2026_Data", inputFile);

    plotCent(h1D_cent2023, h1D_cent2024, h1D_cent2025, h1D_cent2026);

    inputFile->Close();
}

void plotCent(TH1D* h1D_cent2023, TH1D* h1D_cent2024, TH1D* h1D_cent2025, TH1D* h1D_cent2026){

    TCanvas* c1 = new TCanvas("c1", "Centrality Distribution", 800, 600);
    basicCanvasFormatting(c1);
    c1->SetLogy();

    basicHistFormatting(h1D_cent2023);
    basicHistFormatting(h1D_cent2024);
    basicHistFormatting(h1D_cent2025);
    basicHistFormatting(h1D_cent2026);

    //Normalize histograms to unity
    h1D_cent2023->Scale(1.0 / h1D_cent2023->Integral());
    h1D_cent2024->Scale(1.0 / h1D_cent2024->Integral());
    h1D_cent2025->Scale(1.0 / h1D_cent2025->Integral());
    h1D_cent2026->Scale(1.0 / h1D_cent2026->Integral());

    h1D_cent2023->SetLineColor(kRed);
    h1D_cent2024->SetLineColor(kBlue);
    h1D_cent2025->SetLineColor(kGreen+1);
    h1D_cent2026->SetLineColor(kYellow+1);

    h1D_cent2023->SetLineWidth(2);
    h1D_cent2024->SetLineWidth(2);
    h1D_cent2025->SetLineWidth(2);
    h1D_cent2026->SetLineWidth(2);

    h1D_cent2026->GetYaxis()->SetTitle("N of dimuons");
    h1D_cent2026->GetXaxis()->SetTitle("Centrality [%]");
    h1D_cent2026->GetXaxis()->SetRangeUser(0., 10.);

    h1D_cent2026->Draw("HIST");
    h1D_cent2025->Draw("HIST SAME");
    h1D_cent2024->Draw("HIST SAME");
    h1D_cent2023->Draw("HIST SAME");

    TLegend* legend = new TLegend(0.8, 0.7, 0.9, 0.9);
    basicLegendFormatting(legend);
    legend->AddEntry(h1D_cent2023, "PbPb 2023", "l");
    legend->AddEntry(h1D_cent2024, "PbPb 2024", "l");
    legend->AddEntry(h1D_cent2025, "PbPb 2025", "l");
    legend->AddEntry(h1D_cent2026, "PbPb 2026", "l");
    legend->Draw();

    drawLatexText("#bf{CMS}", 0.12, 0.93, 0.042);
    drawLatexText("#it{Internal}", 0.2, 0.93, 0.033);
    drawLatexText("Run 3, ppRef 2024 (5.36 TeV)", 0.6, 0.93, 0.033);
    
    c1->Update();
    c1->SaveAs("centrality_distribution.png");

    delete c1;
    delete legend;
    delete h1D_cent2023;
    delete h1D_cent2024;
    delete h1D_cent2025;
    delete h1D_cent2026;
}


TH1D* getCentFromDataset(std::string datasetName, TFile* inputFile){

    TDirectory* dir = inputFile->GetDirectory(datasetName.c_str());
    TH1D* h1D_cent = dynamic_cast<TH1D*>(dir->Get("h1D_centrality"));

    return h1D_cent;
}