/*
Creates normalized distributions for selected histograms from the input ROOT file and saves them as plots.
*/

//---Libraries
#include <TFile.h>
#include <TDirectory.h>
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

//Selected distributions. Array of histogram names.
std::vector<std::string> selectedDistributions = {
    "h1D_zPt",
    "h1D_zRapidity",
    "h1D_ptMuPlus",
    "h1D_ptMuMinus",
    "h1D_muonPtRelDiff",
    "h1D_acoplanarity"
};


//---Function declarations
void CheckSelected(TFile* inputFile);

//---Main()
void CheckNormalizedDistributions(){

    gROOT->SetBatch(kTRUE);
    TFile *inputFile = new TFile("mySelectedData.root", "READ");

    CheckSelected(inputFile);
    inputFile->Close();
}

void CheckSelected(TFile* inputFile){

    //Get directories
    TDirectory *PbPb_dir = inputFile->GetDirectory("PbPb2023_2024_Data");
    TDirectory *ppRef_dir = inputFile->GetDirectory("ppRef2024_Data");

    for(const auto& histName : selectedDistributions) {

        //Get original histograms
        TH1D* h_PbPb_original = dynamic_cast<TH1D*>(PbPb_dir->Get(histName.c_str()));
        TH1D* h_ppRef_original = dynamic_cast<TH1D*>(ppRef_dir->Get(histName.c_str()));

            if (!h_PbPb_original || !h_ppRef_original) {//Just checking if everything was found.
                std::cerr << "Error: Could not find the histogram "
                        << histName << " in the input file."
                        << std::endl;
                continue;
            }
        
        //Get histogram clones to manipulate.
        TH1D* h_PbPb = dynamic_cast<TH1D*>(h_PbPb_original->Clone(("h_PbPb_" + histName).c_str()));
        TH1D* h_ppRef = dynamic_cast<TH1D*>(h_ppRef_original->Clone(("h_ppRef_" + histName).c_str()));
        h_PbPb->SetDirectory(nullptr);
        h_ppRef->SetDirectory(nullptr);


        //Begin normalization and stuff
        h_PbPb->Scale(1.0 / h_PbPb->Integral());
        h_ppRef->Scale(1.0 / h_ppRef->Integral());

        std::string canvasName = "c_" + histName;
        TCanvas *c = new TCanvas(canvasName.c_str(), "Normalized Distributions", 800, 600);
        basicCanvasFormatting(c);
        c->SetLogy();

            if(histName == "h1D_muonPtRelDiff"){//Turn off log scale for these two distributions.
                c->SetLogy(0);
            }

        basicHistFormatting(h_PbPb);
        basicHistFormatting(h_ppRef);
        
        h_PbPb->SetMarkerStyle(21);
        h_PbPb->SetMarkerSize(0.8);
        h_PbPb->SetMarkerColor(kRed);
        h_PbPb->SetLineColor(kRed);
        
        h_ppRef->SetMarkerStyle(25);
        h_ppRef->SetMarkerSize(0.8);
        h_ppRef->SetMarkerColor(kBlack);
        h_ppRef->SetLineColor(kBlack);

        h_PbPb->GetYaxis()->SetTitle("Normalized Entries");

            if(histName == "h1D_ptMuPlus" || histName == "h1D_ptMuMinus"){
                h_PbPb->GetXaxis()->SetRangeUser(18., 100.);
            }
        
        h_PbPb->Draw("P");
        h_ppRef->Draw("P SAME");
        
        TLegend *leg = new TLegend(0.7, 0.78, 0.95, 0.88);
        basicLegendFormatting(leg);
        leg->AddEntry(h_PbPb, "PbPb2023+2024", "p");
        leg->AddEntry(h_ppRef, "ppRef2024", "p");
        leg->Draw();

        drawLatexText("#bf{CMS}", 0.12, 0.93, 0.042);
        drawLatexText("#it{Work in Progress}", 0.2, 0.93, 0.033);
        drawLatexText("PbPb 2023+2024, ppRef 2024 (5.36 TeV)", 0.5, 0.93, 0.033);

        c->Update();
        std::string outputName = "Normalized_Distributions_Joined_" + histName + plot_extension;
        c->SaveAs(outputName.c_str());

        delete leg;
        delete h_PbPb;
        delete h_ppRef;
        delete c;
    }

}