/*
Check apparent issue in centrality distribution in 2025 and 2026 PbPb data. Compare to 2023 and 2024 PbPb data.
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


// ##############################################################################
// ##############################################################################



//---Function declarations
void CheckTH1D_Centrality(TFile* inputFile);
void CheckTH1D_hiHF(TFile* inputFile);
void CountEventsInCBin(TFile* inputFile);


//---Main()
void Check_centrality(){

    gROOT->SetBatch(kTRUE);
    TFile* inputFile = new TFile("mySelectedData.root", "READ");

    CheckTH1D_Centrality(inputFile);
    CheckTH1D_hiHF(inputFile);
    CountEventsInCBin(inputFile);

    inputFile->Close();
    delete inputFile;
}


//---Function definitions
void CheckTH1D_hiHF(TFile* inputFile){

    //Get 2023 and 2024 PbPb data hiHF distributions from mySelectedData.root
    TH1D* h1D_hiHF_2023 = dynamic_cast<TH1D*>(inputFile->Get("PbPb2023_Data/h1D_hiHF"));
    TH1D* h1D_hiHF_2024 = dynamic_cast<TH1D*>(inputFile->Get("PbPb2024_Data/h1D_hiHF"));
    TH1D* h1D_hiHF_2025 = dynamic_cast<TH1D*>(inputFile->Get("PbPb2025_Data/h1D_hiHF"));
    TH1D* h1D_hiHF_2026 = dynamic_cast<TH1D*>(inputFile->Get("PbPb2026_Data/h1D_hiHF"));

    TCanvas *c = new TCanvas("c", "c", 800, 600);
    basicCanvasFormatting(c);
    c->SetLeftMargin(0.11);

    basicHistFormatting(h1D_hiHF_2023);
    basicHistFormatting(h1D_hiHF_2024);
    basicHistFormatting(h1D_hiHF_2025);
    basicHistFormatting(h1D_hiHF_2026);

    h1D_hiHF_2023->Scale(1.0 / h1D_hiHF_2023->Integral());
    h1D_hiHF_2024->Scale(1.0 / h1D_hiHF_2024->Integral());
    h1D_hiHF_2025->Scale(1.0 / h1D_hiHF_2025->Integral());
    h1D_hiHF_2026->Scale(1.0 / h1D_hiHF_2026->Integral());

    double yMax = std::max({
        h1D_hiHF_2023->GetMaximum(),
        h1D_hiHF_2024->GetMaximum(),
        h1D_hiHF_2025->GetMaximum(),
        h1D_hiHF_2026->GetMaximum()
    });
    yMax *= 1.2;

    h1D_hiHF_2023->SetLineColor(kRed);
    h1D_hiHF_2024->SetLineColor(kBlue);
    h1D_hiHF_2025->SetLineColor(kGreen);
    h1D_hiHF_2026->SetLineColor(kOrange);

    h1D_hiHF_2023->SetLineWidth(2);
    h1D_hiHF_2024->SetLineWidth(2);
    h1D_hiHF_2025->SetLineWidth(2);
    h1D_hiHF_2026->SetLineWidth(2);

    h1D_hiHF_2026->GetYaxis()->SetTitle("Entries (Norm.)");
    h1D_hiHF_2026->GetXaxis()->SetTitle("hiHF");
    h1D_hiHF_2026->GetYaxis()->SetRangeUser(0., yMax);

    h1D_hiHF_2026->Draw("HIST");
    h1D_hiHF_2024->Draw("HIST SAME");
    h1D_hiHF_2025->Draw("HIST SAME");
    h1D_hiHF_2023->Draw("HIST SAME");

    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    basicLegendFormatting(legend);
    legend->AddEntry(h1D_hiHF_2023, "PbPb 2023", "l");
    legend->AddEntry(h1D_hiHF_2024, "PbPb 2024", "l");
    legend->AddEntry(h1D_hiHF_2025, "PbPb 2025", "l");
    legend->AddEntry(h1D_hiHF_2026, "PbPb 2026", "l");
    legend->Draw();

    drawLatexText("#bf{CMS}", 0.12, 0.93, 0.042);
    drawLatexText("#it{Internal}", 0.2, 0.93, 0.033);
    drawLatexText("PbPb 2023-2026 (5.36 TeV)", 0.66, 0.93, 0.033);

    c->Update();
    std::string outputName = "hiHF_check" + plot_extension;
    c->SaveAs(outputName.c_str());

    delete legend;
    delete c;
    delete h1D_hiHF_2023;
    delete h1D_hiHF_2024;
    delete h1D_hiHF_2025;
    delete h1D_hiHF_2026;
}



void CheckTH1D_Centrality(TFile* inputFile){

    //Get 2023 and 2024 PbPb data centrality distributions from mySelectedData.root
    TH1D* h1D_centrality_2023 = dynamic_cast<TH1D*>(inputFile->Get("PbPb2023_Data/h1D_centrality"));
    TH1D* h1D_centrality_2024 = dynamic_cast<TH1D*>(inputFile->Get("PbPb2024_Data/h1D_centrality"));
    TH1D* h1D_centrality_2025 = dynamic_cast<TH1D*>(inputFile->Get("PbPb2025_Data/h1D_centrality"));
    TH1D* h1D_centrality_2026 = dynamic_cast<TH1D*>(inputFile->Get("PbPb2026_Data/h1D_centrality"));

    TCanvas *c = new TCanvas("c", "c", 800, 600);
    basicCanvasFormatting(c);
    c->SetLeftMargin(0.1);
    c->SetLogy();

    basicHistFormatting(h1D_centrality_2023);
    basicHistFormatting(h1D_centrality_2024);
    basicHistFormatting(h1D_centrality_2025);
    basicHistFormatting(h1D_centrality_2026);

    h1D_centrality_2023->SetLineColor(kRed);
    h1D_centrality_2024->SetLineColor(kBlue);
    h1D_centrality_2025->SetLineColor(kGreen);
    h1D_centrality_2026->SetLineColor(kOrange);

    h1D_centrality_2023->SetLineWidth(2);
    h1D_centrality_2024->SetLineWidth(2);
    h1D_centrality_2025->SetLineWidth(2);
    h1D_centrality_2026->SetLineWidth(2);

    h1D_centrality_2026->GetYaxis()->SetTitle("Entries");
    h1D_centrality_2026->GetXaxis()->SetTitle("Centrality (%)");
    h1D_centrality_2026->GetXaxis()->SetRangeUser(0., 10.);
    h1D_centrality_2026->GetYaxis()->SetTitleOffset(1.);

    h1D_centrality_2026->Draw("HIST");
    h1D_centrality_2024->Draw("HIST SAME");
    h1D_centrality_2025->Draw("HIST SAME");
    h1D_centrality_2023->Draw("HIST SAME");

    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    basicLegendFormatting(legend);
    legend->AddEntry(h1D_centrality_2023, "PbPb 2023", "l");
    legend->AddEntry(h1D_centrality_2024, "PbPb 2024", "l");
    legend->AddEntry(h1D_centrality_2025, "PbPb 2025", "l");
    legend->AddEntry(h1D_centrality_2026, "PbPb 2026", "l");
    legend->Draw();

    drawLatexText("#bf{CMS}", 0.11, 0.93, 0.042);
    drawLatexText("#it{Internal}", 0.18, 0.93, 0.033);
    drawLatexText("PbPb 2023-2026 (5.36 TeV)", 0.66, 0.93, 0.033);
    drawLatexText("bin width = 0.5", 0.2, 0.2, 0.033);

    c->Update();
    std::string outputName = "Centrality_check" + plot_extension;
    c->SaveAs(outputName.c_str());

    delete legend;
    delete c;
    delete h1D_centrality_2023;
    delete h1D_centrality_2024;
    delete h1D_centrality_2025;
    delete h1D_centrality_2026;
}

void CountEventsInCBin(TFile* inputFile){

    //Get 2023 and 2024 PbPb data centrality distributions from mySelectedData.root
    TH1D* h1D_centrality_2023 = dynamic_cast<TH1D*>(inputFile->Get("PbPb2023_Data/h1D_eventsVsCentrality"));
    TH1D* h1D_centrality_2024 = dynamic_cast<TH1D*>(inputFile->Get("PbPb2024_Data/h1D_eventsVsCentrality"));
    TH1D* h1D_centrality_2025 = dynamic_cast<TH1D*>(inputFile->Get("PbPb2025_Data/h1D_eventsVsCentrality"));
    TH1D* h1D_centrality_2026 = dynamic_cast<TH1D*>(inputFile->Get("PbPb2026_Data/h1D_eventsVsCentrality"));

    //I want to count the number of events in the first CBin (0 to 0.5% centrality) for each year and print it to the console.
    int firstBinContent_2023 = h1D_centrality_2023->GetBinContent(1);
    int firstBinContent_2024 = h1D_centrality_2024->GetBinContent(1);
    int firstBinContent_2025 = h1D_centrality_2025->GetBinContent(1);
    int firstBinContent_2026 = h1D_centrality_2026->GetBinContent(1);

    std::cout << "Number of events in 0-0.5% centrality bin:" << std::endl;
    std::cout << "PbPb 2023: " << firstBinContent_2023 << std::endl;
    std::cout << "PbPb 2024: " << firstBinContent_2024 << std::endl;
    std::cout << "PbPb 2025: " << firstBinContent_2025 << std::endl;
    std::cout << "PbPb 2026: " << firstBinContent_2026 << std::endl;

    int totalEvents_2023 = h1D_centrality_2023->Integral();
    int totalEvents_2024 = h1D_centrality_2024->Integral();
    int totalEvents_2025 = h1D_centrality_2025->Integral();
    int totalEvents_2026 = h1D_centrality_2026->Integral();

    std::cout << "Total number of events:" << std::endl;
    std::cout << "PbPb 2023: " << totalEvents_2023 << std::endl;
    std::cout << "PbPb 2024: " << totalEvents_2024 << std::endl;
    std::cout << "PbPb 2025: " << totalEvents_2025 << std::endl;
    std::cout << "PbPb 2026: " << totalEvents_2026 << std::endl;

    std::cout << "Fraction of events in 0-0.5% centrality bin:" << std::endl;
    std::cout << "PbPb 2023: " << static_cast<double>(firstBinContent_2023*100) / totalEvents_2023 << std::endl;
    std::cout << "PbPb 2024: " << static_cast<double>(firstBinContent_2024*100) / totalEvents_2024 << std::endl;
    std::cout << "PbPb 2025: " << static_cast<double>(firstBinContent_2025*100) / totalEvents_2025 << std::endl;
    std::cout << "PbPb 2026: " << static_cast<double>(firstBinContent_2026*100) / totalEvents_2026 << std::endl;
    

    delete h1D_centrality_2023;
    delete h1D_centrality_2024;
    delete h1D_centrality_2025;
    delete h1D_centrality_2026;
}