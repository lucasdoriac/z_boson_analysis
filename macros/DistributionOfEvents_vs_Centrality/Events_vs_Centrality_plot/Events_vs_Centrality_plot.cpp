/*
Makes plot of events vs centrality.
From 'mySelectedData.root'.
Draw lines to delimit and illustrate centrality bins.
*/

//---Libraries
#include <TROOT.h>
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
std::string BasePath = "/home/lucas/Documents/CMS/z_boson_analysis/"; //IFT
//std::string BasePath = "/home/lucasdoriac/z_boson_analysis/data/"; //Home

std::string whichDataset = "PbPb_Run3_Data"; // "PbPb_Run3_Data", "PbPb2023_Data", "PbPb2024_Data".
std::string JointPbPb = "PbPb2023-2026"; //"PbPb2023+2024", "PbPb2023", "PbPb2024".
std::string dataSamplesUsed = "PbPb 2023-2026, ppRef 2024 (5.36 TeV)";



//Set of centrality bins for PbPb2024 data. We can decide to change the centrality bins later if we want to.
std::vector<std::pair<double, double>> CentralitySet = {
    {0., 10.},
    {10., 20.},
    {20., 30.},
    {30., 100.},
    {0., 100.}
};


void makePlot(TFile* inputFile);


void Events_vs_Centrality_plot(){

    gROOT->SetBatch(kTRUE);

    TFile* inputFile = TFile::Open("mySelectedData.root", "READ");
    makePlot(inputFile);

    inputFile->Close();
    delete inputFile;
}



void makePlot(TFile* inputFile){

    TDirectory* dir = (TDirectory*)inputFile->Get("PbPb_Run3_Data");
    TH1D* h_Centrality = (TH1D*)dir->Get("h1D_eventsVsCentrality");

    TCanvas* c1 = new TCanvas("c1", "Events vs Centrality", 800, 600);
    basicCanvasFormatting(c1);
    c1->SetLeftMargin(0.11);
    c1->SetFrameLineWidth(2);
    //c1->SetLogy();

    //Format centrality histogram
    basicHistFormatting(h_Centrality);
    h_Centrality->SetLineColor(kCyan-3);
    h_Centrality->SetLineWidth(2);
    
    //Set range and titles
    h_Centrality->GetXaxis()->SetRangeUser(0., 100.);
    h_Centrality->GetYaxis()->SetRangeUser(1., h_Centrality->GetMaximum()*1.2);
    h_Centrality->GetYaxis()->SetTitle("Events");
    h_Centrality->GetXaxis()->SetTitle("Centrality (%)");
    h_Centrality->GetYaxis()->SetTitleOffset(1.3);

    //Three lines delimiting centrality bins: 0-10%, 10-20%, 20-30%, 30-100%.
    TLine* line1 = new TLine(10., 1., 10., h_Centrality->GetMaximum());
    line1->SetLineColor(kRed);
    line1->SetLineStyle(7);
    line1->SetLineWidth(2);

    TLine* line2 = new TLine(20., 1., 20., h_Centrality->GetMaximum());
    line2->SetLineColor(kRed);
    line2->SetLineStyle(7);
    line2->SetLineWidth(2);

    TLine* line3 = new TLine(30., 1., 30., h_Centrality->GetMaximum());
    line3->SetLineColor(kRed);
    line3->SetLineStyle(7);
    line3->SetLineWidth(2);

    //Draw
    h_Centrality->Draw("HIST");
    line1->Draw("SAME");
    line2->Draw("SAME");
    line3->Draw("SAME");

    //Draw labels to the lines
    //drawLatexVerticalText("0-10%", 0.11, 0.8, 0.028);
    drawLatexVerticalText("#bf{10-20%}", 0.215, 0.75, 0.025);
    drawLatexVerticalText("#bf{20-30%}", 0.3, 0.75, 0.025);
    drawLatexVerticalText("#bf{30-100%}", 0.387, 0.75, 0.025);

    //CMS labels and dataset information
    drawLatexText("#bf{CMS}", 0.11, 0.93, 0.04);
    drawLatexText("#it{Internal}", 0.18, 0.93, 0.033);
    drawLatexText(dataSamplesUsed.c_str(), 0.6, 0.93, 0.03);

    //Save
    std::string output = "Events_vs_Centrality_plot" + plot_extension;
    c1->Update();
    c1->SaveAs(output.c_str());


    delete c1;
    delete line1;
    delete line2;
    delete line3;
    delete h_Centrality;
    delete dir;
}