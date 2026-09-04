/*
1- First candidate plot let's do the usual N(\mu+)(pT) and N(\mu-)(pT).
This is just the distribution of dN/dpT for the individual muon daughters and it's the distribution that
the theory paper uses to predict the shift in the peak of THESE distributions in the 30-50 GeV region.
On the bottom panel we will plot R(pT) = N(\mu+)/N(\mu-). We can do this for both PbPb and ppRef.

2- We would like to also look at the double ratio of the R(pT) distributions, i.e. R(PbPb)/R(ppRef).
For this we can plot the usual N(\mu+)(pT) and N(\mu-)(pT), as we did in the first candidate plot, but now for PbPb and ppRef in the same canvas. 
This would give four distributions in the top pannel, and in the bottom panel we can plot the double ratio R(PbPb)/R(ppRef).

Of course, we could also do that in the first candidate plot, since we are calculating R(pT) = N(\mu+)/N(\mu-) for each dataset,
i.e. we are already calculating R(PbPb) and R(ppRef). However, it is interesting at a first moment to plot both datasets separately,
and then later plot them in the same canvas.

3- Another observable option is the event-by-event relative difference in pT between the two muons. This is defined as:
    ΔpT = (pT(mu+) - pT(mu-)) / (pT(mu+) + pT(mu-)).
If we call that variable ΔpT, then we can calculate ΔpT(PbPb) and ΔpT(ppRef), event-by-event.
In an ideal scenario, one might expect ΔpT(PbPb) to be modified by the presence of the B field in the QGP,
with maximum effects on the 40-70% centrality bin (and maybe in the low pT region?).

The sweet observable would be the mean of ΔpT(PbPb), <ΔpT(PbPb)> per centrality bin,
as well as their difference with respect to the ppRef dataset, i.e. <ΔpT(PbPb)> - <ΔpT(ppRef)>.

4- The last simple observable we can try to check is the charge asymmetry of the muons in each dataset. This is defined as:
    A(X) = (N(mu+) - N(mu-)) / (N(mu+) + N(mu-)).
Notice that this is a YIELD asymmetry, NOT a pT asymmetry. It is the difference in the number of muons of each charge.
(Normalized by the total number of muons?)
We can calculate A(PbPb) and A(ppRef) and plot them in the same canvas.
In the bottom panel we can plot the difference in the charge asymmetry between PbPb and ppRef, i.e. ΔA = A(PbPb) - A(ppRef).

-------------------------------------
After doing the four candidate plots, now we move to another plot that gives a centrality-dependent observable.

5- 
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

using namespace std;

//---Macro settings

std::string plot_extension = ".png"; // ".png" for regular development and ".pdf" for final quality plots
std::string BasePath = "/home/lucas/Documents/CMS_analyzes/Z_boson_analysis/"; //IFT
//std::string BasePath = "/home/lucasdoriac/z_boson_analysis/data/"; //Home

const int MIN_CENTRALITY = 0.;
const int MAX_CENTRALITY = 180;
const double MAX_ZVTX = 15.0;

const double MINZ_MASS = 80.;
const double MAXZ_MASS = 100.;
const double RAPIDITYCUTVALUE = 2.4;

const float ETACUTVALUE = 2.4;
const double PTCUTVALUE = 10.;

// ##############################################################################
// ##############################################################################

//---Enumerates
enum class SampleType {
    Data,
    MC
};

enum class CollisionSystem {
    PbPb2023,
    PbPb2024,
    ppRef2024
};

//---Structs
struct Dataset {
    std::string name;
    SampleType type;
    CollisionSystem system;
    std::string treeName;
    std::string filePattern;
    std::string basePath;
};

Dataset datasets[] = {
    {
        "PbPb2023_Data",
        SampleType::Data,
        CollisionSystem::PbPb2023,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_PbPb2023.root",
        BasePath + "Data/PbPb2023/"
    },

    {
        "PbPb2024_Data",
        SampleType::Data,
        CollisionSystem::PbPb2024,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_PbPb2024Data.root",
        BasePath + "Data/PbPb2024/"
    },

    {
        "ppRef2024_Data",
        SampleType::Data,
        CollisionSystem::ppRef2024,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_ppRef2024.root",
        BasePath + "Data/ppRef2024/"
    },

    {
        "PbPb2024_MC",
        SampleType::MC,
        CollisionSystem::PbPb2024,
        "hionia/myTree",
        "Oniatree_PowhegZtoMuMu_PbPb2024_*.root",
        BasePath + "MC/PbPb2024/DYto2Mu_MLL-50_TuneCP5_5p36TeV_powheg-pythia8/PowhegEmbedded_March9/260309_143939/0000/"
    },

    {
        "ppRef2024_MC",
        SampleType::MC,
        CollisionSystem::ppRef2024,
        "hionia/myTree",
        "Oniatree_PowhegZtoMuMu_ppRef2024_*.root",
        BasePath + "MC/ppRef2024/DYToMuMu_M-50_TuneCP5_5p36TeV_powheg-pythia8/Powheg_ppRefPileup_March20/260320_125046/0000/"
    }
};

//For the second candidate plot, we will need to create a new vector to hold the histograms for each dataset.
//This will allow us to easily access the histograms for each dataset and plot them in the same canvas.
struct TmpHistStruct {
    TH1D* histMuPlus;
    TH1D* histMuMinus;
    TH1D* histRatio;
};


//---Common functions
void FillChargeYields(const Dataset& dataset, TH1D*& histMuPlus, TH1D*& histMuMinus);

//---First candidate functions
void PlotChargeYields(const Dataset& dataset, TH1D* histMuPlus, TH1D* histMuMinus);

//---Second candidate functions
void PlotDoubleRatio(const Dataset& datasetPbPb, const Dataset& datasetppRef, TH1D* histPbPbMuPlus, TH1D* histPbPbMuMinus, TH1D* histppRefMuPlus, TH1D* histppRefMuMinus);

//---Third candidate functions
TH1D* MakePtRelativeDiff(const Dataset& dataset);
void PlotPtRelativeDiff(TH1D* histPbPb, TH1D* histppRef);

//---Fourth candidate functions
TH1D* MakeChargeAsymmetryHist(const Dataset& dataset);
void PlotChargeAsymmetry(TH1D* histPbPb, TH1D* histppRef);

//---Histogram formatting
void basicCanvasFormatting(TCanvas *c, TPad *pad1, TPad *pad2);
void basicHistFormatting(TH1D *hist, bool isRatio = false);
void basicLegendFormatting(TLegend *leg);
void drawLatexText(TString latexText = "#bf{CMS}", double x = 0.15, double y = 0.93, double TextSize = 0.05);
void PrintBinInfo(TH1D* hist);

//---Main function
void Penguin(){

    gROOT->SetBatch(kTRUE);

    //First candidate plot stuff
    TH1D* histYieldMuPlus = new TH1D("histYieldMuPlus", "histYieldMuPlus", 100, 0., 100.);
    TH1D* histYieldMuMinus = new TH1D("histYieldMuMinus", "histYieldMuMinus", 100, 0., 100.);

    //PbPb2024
    FillChargeYields(datasets[1], histYieldMuPlus, histYieldMuMinus);
    PlotChargeYields(datasets[1], histYieldMuPlus, histYieldMuMinus);

    //Empty histograms before reusing them.
    histYieldMuPlus->Reset();
    histYieldMuMinus->Reset();

    //ppRef2024
    FillChargeYields(datasets[2], histYieldMuPlus, histYieldMuMinus);
    PlotChargeYields(datasets[2], histYieldMuPlus, histYieldMuMinus);
    //---

    //Second candidate plot stuff
    TH1D* histPbPbMuPlus = new TH1D("histPbPbMuPlus", "histPbPbMuPlus", 100, 0., 100.);
    TH1D* histPbPbMuMinus = new TH1D("histPbPbMuMinus", "histPbPbMuMinus", 100, 0., 100.);
    TH1D* histppRefMuPlus = new TH1D("histppRefMuPlus", "histppRefMuPlus", 100, 0., 100.);
    TH1D* histppRefMuMinus = new TH1D("histppRefMuMinus", "histppRefMuMinus", 100, 0., 100.);

    //Fill the four histograms.
    FillChargeYields(datasets[1], histPbPbMuPlus, histPbPbMuMinus);
    FillChargeYields(datasets[2], histppRefMuPlus, histppRefMuMinus);
    //Just one plot function for the four distributions plus the double ratio on bottom pad.
    PlotDoubleRatio(datasets[1], datasets[2], histPbPbMuPlus, histPbPbMuMinus, histppRefMuPlus, histppRefMuMinus);
    //---

    //Third candidate plot stuff
    TH1D* histPbPbPtRelDiff = MakePtRelativeDiff(datasets[1]);
    TH1D* histppRefPtRelDiff = MakePtRelativeDiff(datasets[2]);
    PlotPtRelativeDiff(histPbPbPtRelDiff, histppRefPtRelDiff);
    //---

    //Fourth candidate plot stuff
    TH1D* histPbPbChargeAsymmetry = MakeChargeAsymmetryHist(datasets[1]);
    TH1D* histppRefChargeAsymmetry = MakeChargeAsymmetryHist(datasets[2]);
    PlotChargeAsymmetry(histPbPbChargeAsymmetry, histppRefChargeAsymmetry);
    //---

    //Clean up
    /*delete histYieldMuPlus;
    delete histYieldMuMinus;
    delete histPbPbMuPlus;
    delete histPbPbMuMinus;
    delete histppRefMuPlus;
    delete histppRefMuMinus;
    delete histPbPbPtRelDiff;
    delete histppRefPtRelDiff;
    delete histPbPbChargeAsymmetry;
    delete histppRefChargeAsymmetry;
    delete histPbPb;
    delete histppRef;*/
}

void PlotDoubleRatio(const Dataset& datasetPbPb, const Dataset& datasetppRef, TH1D* histPbPbMuPlus, TH1D* histPbPbMuMinus, TH1D* histppRefMuPlus, TH1D* histppRefMuMinus){

    //Plots double ratio for a specific dataset.
    TCanvas *c = new TCanvas("c", "c", 800, 800);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.30, 1, 1.0);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.00, 1, 0.30);

    basicCanvasFormatting(c, pad1, pad2);
    //Since we are normalizing it the left margin needs to be a little higher
    pad1->SetLeftMargin(0.12);
    pad2->SetLeftMargin(0.12);
    pad1->Draw();
    pad2->Draw();

    //Top pannel
    pad1->cd();
    basicHistFormatting(histPbPbMuPlus, false);
    basicHistFormatting(histPbPbMuMinus, false);
    basicHistFormatting(histppRefMuPlus, false);
    basicHistFormatting(histppRefMuMinus, false);

    //We don't need to calculate statistics since we already plotted these distributions before.

    //Histogram formatting
    histPbPbMuPlus->SetMarkerStyle(24);
    histPbPbMuPlus->SetMarkerSize(0.8);
    histPbPbMuPlus->SetMarkerColor(kRed);
    histPbPbMuPlus->SetLineColor(kRed);

    histPbPbMuMinus->SetMarkerStyle(24);
    histPbPbMuMinus->SetMarkerSize(0.8);
    histPbPbMuMinus->SetMarkerColor(kBlue);
    histPbPbMuMinus->SetLineColor(kBlue);

    histppRefMuPlus->SetMarkerStyle(25);
    histppRefMuPlus->SetMarkerSize(0.8);
    histppRefMuPlus->SetMarkerColor(kBlack);
    histppRefMuPlus->SetLineColor(kBlack);

    histppRefMuMinus->SetMarkerStyle(25);
    histppRefMuMinus->SetMarkerSize(0.8);
    histppRefMuMinus->SetMarkerColor(kGray+1);
    histppRefMuMinus->SetLineColor(kGray+1);


    histPbPbMuPlus->GetYaxis()->SetTitle("1/N dN/dp_{T} [GeV]^{-1}");
    histPbPbMuPlus->GetYaxis()->SetTitleOffset(1.1);

    //These need to be normalized since we are comparing two different collision systems with very different Z yields.
    histPbPbMuPlus->Scale(1.0/histPbPbMuPlus->Integral());
    histPbPbMuMinus->Scale(1.0/histPbPbMuMinus->Integral());
    histppRefMuPlus->Scale(1.0/histppRefMuPlus->Integral());
    histppRefMuMinus->Scale(1.0/histppRefMuMinus->Integral());

    histPbPbMuPlus->Draw("P");
    histPbPbMuMinus->Draw("P SAME");
    histppRefMuPlus->Draw("P SAME");
    histppRefMuMinus->Draw("P SAME");

    //Selections and cuts
    drawLatexText(Form("p_{T} > %.0f GeV, |#eta| < %.1f", PTCUTVALUE, ETACUTVALUE), 0.7, 0.35, 0.03);
    drawLatexText(Form("|y| < %.1f", RAPIDITYCUTVALUE), 0.7, 0.25, 0.03);
    drawLatexText(Form("%.0f < M_{#mu#mu} < %.0f GeV", MINZ_MASS, MAXZ_MASS), 0.7, 0.2, 0.03);
    drawLatexText(Form("%d < Cent. < %.0f%%", MIN_CENTRALITY, MAX_CENTRALITY / 2.0), 0.7, 0.3, 0.03);

    TLegend *leg = new TLegend(0.2, 0.58, 0.35, 0.85);
    basicLegendFormatting(leg);
    leg->AddEntry(histPbPbMuPlus, "PbPb #mu^{+}", "p");
    leg->AddEntry(histPbPbMuMinus, "PbPb #mu^{-}", "p");
    leg->AddEntry(histppRefMuPlus, "pp #mu^{+}", "p");
    leg->AddEntry(histppRefMuMinus, "pp #mu^{-}", "p");
    leg->Draw();

    pad1->Update();

    //Bottom pannel
    pad2->cd();
    TH1D* histPbPbRatio = new TH1D("histPbPbRatio", "histPbPbRatio", histPbPbMuPlus->GetNbinsX(), histPbPbMuPlus->GetXaxis()->GetXmin(), histPbPbMuPlus->GetXaxis()->GetXmax());
    TH1D* histppRefRatio = new TH1D("histppRefRatio", "histppRefRatio", histppRefMuPlus->GetNbinsX(), histppRefMuPlus->GetXaxis()->GetXmin(), histppRefMuPlus->GetXaxis()->GetXmax());
    histPbPbRatio->Divide(histPbPbMuPlus, histPbPbMuMinus, 1.0, 1.0, "B");
    histppRefRatio->Divide(histppRefMuPlus, histppRefMuMinus, 1.0, 1.0, "B");

    TH1D* histDoubleRatio = new TH1D("histDoubleRatio", "histDoubleRatio", histPbPbMuPlus->GetNbinsX(), histPbPbMuPlus->GetXaxis()->GetXmin(), histPbPbMuPlus->GetXaxis()->GetXmax());
    histDoubleRatio->Divide(histPbPbRatio, histppRefRatio, 1.0, 1.0, "B");

    basicHistFormatting(histDoubleRatio, true);
    histDoubleRatio->SetMarkerStyle(24);
    histDoubleRatio->SetMarkerSize(0.8);
    histDoubleRatio->SetMarkerColor(kBlack);
    histDoubleRatio->SetLineColor(kBlack);
    histDoubleRatio->GetYaxis()->SetTitle("R_{PbPb}/R_{ppRef}");
    histDoubleRatio->GetXaxis()->SetTitle("p_{T} [GeV]");
    histDoubleRatio->Draw("P");

    //A horizontal line just to guide the eye.
    TLine *line = new TLine(histDoubleRatio->GetXaxis()->GetXmin(), 1.0, histDoubleRatio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineColor(kMagenta+2);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");

    pad2->Update();

    //Back to the canvas
    c->cd();
    drawLatexText("#bf{CMS}", 0.12, 0.95, 0.038);
    drawLatexText("#it{Work in Progress}", 0.2, 0.95, 0.03);
    drawLatexText("PbPb 2024, ppRef 2024 (5.36 TeV)", 0.52, 0.95, 0.03);

    //Save
    TString output = datasetPbPb.name + "_" + datasetppRef.name + "_DoubleRatio" + plot_extension;
    c->Update();
    c->SaveAs(output);

    delete histPbPbRatio;
    delete histppRefRatio;
    delete histDoubleRatio;
    delete line;
    delete c;
}

void PlotChargeYields(const Dataset& dataset, TH1D* histMuPlus, TH1D* histMuMinus){

    //Plots dN/dpT distribution for a specific dataset.
    //Top pannel: dN/dpT for mu+ and mu-. (Yield distribution of daughter muons in pT bins).
    //Bottom pannel: R(pT) = N(mu+)/N(mu-). (Ratio of the yield distributions of daughter muons in pT bins).
    TCanvas *c = new TCanvas("c", "c", 800, 800);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.30, 1, 1.0);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.00, 1, 0.30);

    basicCanvasFormatting(c, pad1, pad2);
    pad1->Draw();
    pad2->Draw();

    //Top pannel
    pad1->cd();
    basicHistFormatting(histMuPlus, false);
    basicHistFormatting(histMuMinus, false);

    //Calculating histogram statistics
    //Z yield
    double muPlusEntries = histMuPlus->GetEntries();
    double muMinusEntries = histMuMinus->GetEntries();
    if(muPlusEntries != muMinusEntries) std::cout << ">>> WARNING: Number of mu+ and mu- entries are different! <<<" << std::endl;
    
    double Z_yield_err = 0.0;
    double Z_yield = histMuPlus->IntegralAndError(1, histMuPlus->GetNbinsX(), Z_yield_err);

    //Peak and mean of dN/dpT distribution
    double muPlusPeak = histMuPlus->GetBinCenter(histMuPlus->GetMaximumBin());
    double muMinusPeak = histMuMinus->GetBinCenter(histMuMinus->GetMaximumBin());
    double Peak_diff = muPlusPeak - muMinusPeak;

    double muPlusMean = histMuPlus->GetMean();
    double muPlusMean_err = histMuPlus->GetMeanError();
    double muMinusMean = histMuMinus->GetMean();
    double muMinusMean_err = histMuMinus->GetMeanError();

    //Histogram formatting
    //histMuPlus->SetLineWidth(2);
    histMuPlus->SetMarkerStyle(25);
    histMuPlus->SetMarkerSize(0.8);
    histMuPlus->SetMarkerColor(kRed);
    histMuPlus->SetLineColor(kRed);

    //histMuMinus->SetLineWidth(2);
    histMuMinus->SetMarkerStyle(25);
    histMuMinus->SetMarkerSize(0.8);
    histMuMinus->SetMarkerColor(kBlue);
    histMuMinus->SetLineColor(kBlue);

    //Scaling them before plotting
    histMuPlus->Scale(1e-2);
    histMuMinus->Scale(1e-2);

    //If you want to normalize the histograms just uncomment the following. Just notice you have to unscale and comment the LatexText of 10^2.
    //histMuPlus->Scale(1.0/histMuPlus->Integral());
    //histMuMinus->Scale(1.0/histMuMinus->Integral());
    histMuPlus->GetYaxis()->SetTitle("dN/dp_{T} [GeV]^{-1}");

    histMuPlus->Draw("P");
    histMuMinus->Draw("P SAME");

    //Write histogram statistics
    TLegend *info = new TLegend(0.62, 0.62, 0.88, 0.86);
    info->SetBorderSize(0);
    info->SetFillStyle(0);
    info->SetTextFont(42);
    info->SetTextSize(0.035);
    info->SetMargin(0.0);
    info->SetEntrySeparation(0.02);
    info->AddEntry((TObject*)nullptr,
                Form("Z yield = %.1f #pm %.1f", Z_yield, Z_yield_err), "");
    info->AddEntry((TObject*)nullptr,
                Form("Peak diff. = %.2f GeV", Peak_diff), "");
    info->AddEntry((TObject*)nullptr,
                Form("<p_{T}>_{#mu^{+}} = %.2f #pm %.2f GeV",
                        muPlusMean, muPlusMean_err), "");
    info->AddEntry((TObject*)nullptr,
                Form("<p_{T}>_{#mu^{-}} = %.2f #pm %.2f GeV",
                        muMinusMean, muMinusMean_err), "");
    info->Draw();

    //Since we scaled it.
    drawLatexText("#times 10^{2}", 0.05, 0.93, 0.042);
    //Selections and cuts
    drawLatexText(Form("p_{T} > %.0f GeV, |#eta| < %.1f", PTCUTVALUE, ETACUTVALUE), 0.7, 0.35, 0.03);
    drawLatexText(Form("|y| < %.1f", RAPIDITYCUTVALUE), 0.7, 0.25, 0.03);
    drawLatexText(Form("%.0f < M_{#mu#mu} < %.0f GeV", MINZ_MASS, MAXZ_MASS), 0.7, 0.2, 0.03);
    if(dataset.system == CollisionSystem::PbPb2024) drawLatexText(Form("%d < Cent. < %.0f%%", MIN_CENTRALITY, MAX_CENTRALITY / 2.0), 0.7, 0.3, 0.03);


    TLegend *leg = new TLegend(0.2, 0.7, 0.35, 0.85);
    basicLegendFormatting(leg);
    leg->AddEntry(histMuPlus, "#mu^{+}", "l");
    leg->AddEntry(histMuMinus, "#mu^{-}", "l");
    leg->Draw();

    pad1->Update();

    //Bottom pannel
    pad2->cd();
    
    TH1D* histRatio = new TH1D("histRatio", "histRatio", histMuPlus->GetNbinsX(), histMuPlus->GetXaxis()->GetXmin(), histMuPlus->GetXaxis()->GetXmax());
    histRatio->Divide(histMuPlus, histMuMinus, 1.0, 1.0, "B");
    
    basicHistFormatting(histRatio, true);
    histRatio->SetMarkerStyle(24);
    histRatio->SetMarkerSize(0.8);
    histRatio->SetMarkerColor(kBlack);
    histRatio->SetLineColor(kBlack);
    histRatio->GetYaxis()->SetTitle("R(p_{T})");
    histRatio->GetXaxis()->SetTitle("p_{T} [GeV]");
    histRatio->GetYaxis()->ChangeLabel(1, -1);//Hides the first y-axis label, which is 0.
    histRatio->Draw("P");

    //A horizontal line just to guide the eye.
    TLine *line = new TLine(histRatio->GetXaxis()->GetXmin(), 1.0, histRatio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineColor(kMagenta+2);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");

    pad2->Update();

    //Back to the canvas
    c->cd();
    drawLatexText("#bf{CMS}", 0.12, 0.95, 0.038);
    drawLatexText("#it{Work in Progress}", 0.2, 0.95, 0.03);
    //For now we are only working with PbPb2024 and ppRef2024 datasets. So we can just check the dataset.system to know which one we are plotting.
    if (dataset.system == CollisionSystem::PbPb2024) drawLatexText("PbPb 2024 (5.36 TeV)", 0.68, 0.95, 0.03);
    else drawLatexText("ppRef 2024 (5.36 TeV)", 0.68, 0.95, 0.03);
    //Save
    TString output = dataset.name + "_ChargeYields" + plot_extension;
    c->Update();
    c->SaveAs(output);

    delete histRatio;
    delete info;
    delete leg;
    delete line;
    delete c;
}

void FillChargeYields(const Dataset& dataset, TH1D*& histMuPlus, TH1D*& histMuMinus){

    // Load root file.
    std::string fullPath = dataset.basePath + dataset.filePattern;

    TChain *chain = new TChain(dataset.treeName.c_str());
    chain->Add(fullPath.c_str());

    std::cout << "> Number of files = " << chain->GetListOfFiles()->GetEntries() << "\n" << std::endl;

    std::cout << "> Opening files " << fullPath << "\n" << std::endl;

    std::cout << "> Running function " << __func__ << " on " << dataset.name << "\n" << std::endl;
    
    //Total number of events on Tree.
    Long64_t nEvents = chain->GetEntries();

    const int MAX_DIMUON = 1000;
    const int MAX_MUON   = 1000;

    //For now, writing ONLY branches that are relevant to the observable we want to measure.

    //Event-level variables
    Int_t Centrality;
    Float_t  zVtx;

    chain->SetBranchAddress("zVtx", &zVtx);
    if(dataset.system == CollisionSystem::PbPb2024) {
        chain->SetBranchAddress("Centrality", &Centrality);
    }

    //Dimuon-level variables
    Short_t Reco_Dimuon_size;

    Short_t Reco_Dimuon_sign[1000];
    Short_t Reco_Dimuon_muonPlusIndex[1000];
    Short_t Reco_Dimuon_muonMinusIndex[1000];

    ULong64_t Reco_Dimuon_trig[1000];

    std::vector<float>* Reco_Dimuon_pt = nullptr;
    std::vector<float>* Reco_Dimuon_eta = nullptr;
    std::vector<float>* Reco_Dimuon_rapidity = nullptr;
    std::vector<float>* Reco_Dimuon_phi = nullptr;
    std::vector<float>* Reco_Dimuon_invMass = nullptr;

    std::vector<float>* Reco_Dimuon_muonPtDiff = nullptr;
    std::vector<float>* Reco_Dimuon_muonPtRelDiff = nullptr;

    chain->SetBranchAddress("Reco_Dimuon_size", &Reco_Dimuon_size);

    chain->SetBranchAddress("Reco_Dimuon_sign", Reco_Dimuon_sign);

    chain->SetBranchAddress("Reco_Dimuon_pt", &Reco_Dimuon_pt);
    chain->SetBranchAddress("Reco_Dimuon_eta", &Reco_Dimuon_eta);
    chain->SetBranchAddress("Reco_Dimuon_rapidity", &Reco_Dimuon_rapidity);
    chain->SetBranchAddress("Reco_Dimuon_phi", &Reco_Dimuon_phi);
    chain->SetBranchAddress("Reco_Dimuon_invMass", &Reco_Dimuon_invMass);

    chain->SetBranchAddress("Reco_Dimuon_muonPtDiff", &Reco_Dimuon_muonPtDiff);
    chain->SetBranchAddress("Reco_Dimuon_muonPtRelDiff", &Reco_Dimuon_muonPtRelDiff);

    chain->SetBranchAddress("Reco_Dimuon_muonPlusIndex", Reco_Dimuon_muonPlusIndex);
    chain->SetBranchAddress("Reco_Dimuon_muonMinusIndex", Reco_Dimuon_muonMinusIndex);

    chain->SetBranchAddress("Reco_Dimuon_trig", Reco_Dimuon_trig);

    //Muon-level variables
    Short_t Reco_Muon_size;

    std::vector<float>* Reco_Muon_pt = nullptr;
    //std::vector<float>* Reco_Muon_ptErrTrk = nullptr;
    std::vector<float>* Reco_Muon_eta = nullptr;
    std::vector<float>* Reco_Muon_phi = nullptr;
    std::vector<float>* Reco_Muon_mass = nullptr;

    ULong64_t Reco_Muon_trig[1000];
    Bool_t Reco_Muon_isTightCutBased[1000];
    
    chain->SetBranchAddress("Reco_Muon_size", &Reco_Muon_size);

    chain->SetBranchAddress("Reco_Muon_pt", &Reco_Muon_pt);
    //chain->SetBranchAddress("Reco_Muon_ptErrTrk", &Reco_Muon_ptErrTrk);
    chain->SetBranchAddress("Reco_Muon_eta", &Reco_Muon_eta);
    chain->SetBranchAddress("Reco_Muon_phi", &Reco_Muon_phi);
    chain->SetBranchAddress("Reco_Muon_mass", &Reco_Muon_mass);

    chain->SetBranchAddress("Reco_Muon_trig", Reco_Muon_trig);
    chain->SetBranchAddress("Reco_Muon_isTightCutBased", Reco_Muon_isTightCutBased);
    //

    //Event-level cut values.
    const int minCentrality = MIN_CENTRALITY;
    const int maxCentrality = MAX_CENTRALITY; //Cent max = maxCentrality/2.
    const double maxZvtx = MAX_ZVTX;

    //Z-level cut values.
    const double minZ_Mass = MINZ_MASS; 
    const double maxZ_Mass = MAXZ_MASS; 
    const double RapidityCutValue = RAPIDITYCUTVALUE;

    //Muon-level cut values.
    const float EtaCutValue = ETACUTVALUE;
    const double ptCutValue = PTCUTVALUE;
    bool MuPlIsTight;
    bool MuMiIsTight;

    //Observables
    TLorentzVector Z;
    double ptplus, ptminus;
    double etaplus, etaminus;

    //New trigger selection 'L2SingleMu12'.
    ULong64_t triggerBit = 1ULL << 7;

    //Aug.26- For now the trigger matching selection is only applied to PbPb2024 dataset, since there's an issue with ppRef2024 trigger.

    for(Long64_t i = 0; i < nEvents; ++i){//Loop through all EVENTS in the CHAIN.

        chain->GetEntry(i); //Get event i.

        //Good event selection
        bool goodVertex = (std::abs(zVtx) < maxZvtx);
        bool goodCent = true;

        if (dataset.system == CollisionSystem::PbPb2024) {
            goodCent = (Centrality > minCentrality && Centrality < maxCentrality);
        }

        if (!goodVertex) continue;
        if (!goodCent) continue;

        for(Short_t j = 0; j < Reco_Dimuon_size; ++j){ //Loop through all reco dimuon candidates of event i.
            
            //Good Z selection
            bool goodMass = (Reco_Dimuon_invMass->at(j) > minZ_Mass && Reco_Dimuon_invMass->at(j) < maxZ_Mass);
            bool goodRapidity = (std::abs(Reco_Dimuon_rapidity->at(j)) < RapidityCutValue);
            bool goodCharge = (Reco_Dimuon_sign[j] == 0);
            
            bool isTriggerMatched = true;

                if (dataset.system == CollisionSystem::PbPb2024){
                    isTriggerMatched = (Reco_Dimuon_trig[j]) & triggerBit;
                }

            if (!goodMass) continue;
            if (!goodRapidity) continue;
            if (!goodCharge) continue;
            if (!isTriggerMatched) continue;

            Short_t muonPlusIndex = Reco_Dimuon_muonPlusIndex[j];
            Short_t muonMinusIndex = Reco_Dimuon_muonMinusIndex[j];

            //Good muon selection
            ptplus = Reco_Muon_pt->at(muonPlusIndex); //pT of antimuon.
            ptminus = Reco_Muon_pt->at(muonMinusIndex); //pT of corresponding muon.
            etaplus = Reco_Muon_eta->at(muonPlusIndex); //Pseudorapidity of antimuon.
            etaminus = Reco_Muon_eta->at(muonMinusIndex); //Pseudorapidity of corresponding muon.
            MuPlIsTight = Reco_Muon_isTightCutBased[muonPlusIndex];
            MuMiIsTight = Reco_Muon_isTightCutBased[muonMinusIndex];
            
            bool goodMuPl = (ptplus > ptCutValue)
                            && (std::abs(etaplus) < EtaCutValue)
                            && (MuPlIsTight);

            bool goodMuMi = (ptminus > ptCutValue)
                            && (std::abs(etaminus) < EtaCutValue)
                            && (MuMiIsTight);

            if (!goodMuPl || !goodMuMi) continue;
        
            //Fill histograms with yields of daughter muons in pT bins.
            histMuPlus->Fill(ptplus);
            histMuMinus->Fill(ptminus);
        }

    }//Exiting event-by-event loop.

    delete chain;
}

void PlotPtRelativeDiff(TH1D* histPbPb, TH1D* histppRef){

    //Plots the Pt Relative Difference for PbPb2024 and ppRef2024 datasets in the top pad and the difference in the bottom pad.
    TCanvas *c = new TCanvas("c", "Pt Relative Diff", 800, 800);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.30, 1, 1.0);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.00, 1, 0.30);
    basicCanvasFormatting(c, pad1, pad2);
    pad1->SetLeftMargin(0.13);
    pad2->SetLeftMargin(0.13);
    pad1->Draw();
    pad2->Draw();

    //Top pannel
    pad1->cd();
    basicHistFormatting(histPbPb, false);
    basicHistFormatting(histppRef, false);

    //Calculating histogram statistics
    //Mean of each distribution
    double PbPbRelDiff = histPbPb->GetMean()*1e3;
    double PbPbRelDiff_err = histPbPb->GetMeanError()*1e3;
    double ppRefRelDiff = histppRef->GetMean()*1e3;
    double ppRefRelDiff_err = histppRef->GetMeanError()*1e3;

    //We actually need to normalize it.
    histPbPb->Scale(1./histPbPb->Integral());
    histppRef->Scale(1./histppRef->Integral());

    //Histogram formatting
    histPbPb->SetMarkerStyle(20);
    histPbPb->SetMarkerSize(0.8);
    histPbPb->SetMarkerColor(kBlack);
    histPbPb->SetLineColor(kBlack);

    histppRef->SetMarkerStyle(24);
    histppRef->SetMarkerSize(0.8);
    histppRef->SetMarkerColor(kRed);
    histppRef->SetLineColor(kRed);

    histPbPb->GetYaxis()->CenterTitle(true);
    histPbPb->GetYaxis()->SetTitle("1/N dN/d(#Delta p_{T})");
    histPbPb->GetYaxis()->SetTitleOffset(1.2);

    histPbPb->Draw("P");
    histppRef->Draw("P SAME");

    //Write histogram statistics
    TLegend *info = new TLegend(0.62, 0.75, 0.84, 0.85);
    info->SetBorderSize(0);
    info->SetFillStyle(0);
    info->SetTextFont(42);
    info->SetTextSize(0.03);
    info->SetMargin(0.0);
    info->SetEntrySeparation(0.02);
    info->AddEntry((TObject*)nullptr,
                Form("#LT #Delta p_{T}#GT_{PbPb} = %.1f #pm %.1f #times 10^{%d} GeV",
                        PbPbRelDiff, PbPbRelDiff_err, 3), "");
    info->AddEntry((TObject*)nullptr,
                Form("#LT #Delta p_{T}#GT_{ppRef} = %.1f #pm %.1f #times 10^{%d} GeV",
                        ppRefRelDiff, ppRefRelDiff_err, 3), "");
    info->Draw();

    //Selections and cuts
    drawLatexText(Form("p_{T} > %.0f GeV, |#eta| < %.1f", PTCUTVALUE, ETACUTVALUE), 0.7, 0.35, 0.03);
    drawLatexText(Form("|y| < %.1f", RAPIDITYCUTVALUE), 0.7, 0.25, 0.03);
    drawLatexText(Form("%.0f < M_{#mu#mu} < %.0f GeV", MINZ_MASS, MAXZ_MASS), 0.7, 0.2, 0.03);
    drawLatexText(Form("%d < Cent. < %.0f%%", MIN_CENTRALITY, MAX_CENTRALITY / 2.0), 0.7, 0.3, 0.03);

    TLegend *leg = new TLegend(0.22, 0.65, 0.45, 0.81);
    basicLegendFormatting(leg);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.036);
    leg->AddEntry(histPbPb, "PbPb2024", "p");
    leg->AddEntry(histppRef, "ppRef2024", "p");
    leg->Draw();

    pad1->Update();

    //Bottom pannel
    pad2->cd();
    TH1D* histDiff = (TH1D*)histPbPb->Clone("histDiff");
    histDiff->Add(histppRef, -1.0);//Subtraction bin-by-bin.
    
    basicHistFormatting(histDiff, true);
    histDiff->SetMarkerStyle(24);
    histDiff->SetMarkerSize(0.8);
    histDiff->SetMarkerColor(kBlack);
    histDiff->SetLineColor(kBlack);
    histDiff->GetYaxis()->SetTitle("#Delta_{PbPb} - #Delta_{ppRef}");
    histDiff->GetXaxis()->SetTitle("#Delta p_{T}");
    histDiff->GetYaxis()->SetTitleOffset(0.6);
    histDiff->Draw("P");

    //A horizontal line just to guide the eye.
    TLine *line = new TLine(histDiff->GetXaxis()->GetXmin(), 0.0, histDiff->GetXaxis()->GetXmax(), 0.0);
    line->SetLineColor(kMagenta+2);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");

    pad2->Update();

    //Back to the canvas
    c->cd();
    drawLatexText("#bf{CMS}", 0.12, 0.95, 0.038);
    drawLatexText("#it{Work in Progress}", 0.2, 0.95, 0.03);
    drawLatexText("PbPb 2024, ppRef 2024 (5.36 TeV)", 0.52, 0.95, 0.03);

    //Save
    TString output = "PbPb_ppRef_RelativeDiff" + plot_extension;
    c->Update();
    c->SaveAs(output);

    delete histPbPb;
    delete histppRef;
    delete histDiff;
    delete line;
    delete c;
}

TH1D* MakePtRelativeDiff(const Dataset& dataset){

    // Load root file.
    std::string fullPath = dataset.basePath + dataset.filePattern;

    TChain *chain = new TChain(dataset.treeName.c_str());
    chain->Add(fullPath.c_str());

    std::cout << "> Number of files = " << chain->GetListOfFiles()->GetEntries() << "\n" << std::endl;

    std::cout << "> Opening files " << fullPath << "\n" << std::endl;

    std::cout << "> Running function " << __func__ << " on " << dataset.name << "\n" << std::endl;
    
    //Total number of events on Tree.
    Long64_t nEvents = chain->GetEntries();

    const int MAX_DIMUON = 1000;
    const int MAX_MUON   = 1000;

    //For now, writing ONLY branches that are relevant to the observable we want to measure.

    //Event-level variables
    Int_t Centrality;
    Float_t  zVtx;

    chain->SetBranchAddress("zVtx", &zVtx);
    if(dataset.system == CollisionSystem::PbPb2024) {
        chain->SetBranchAddress("Centrality", &Centrality);
    }

    //Dimuon-level variables
    Short_t Reco_Dimuon_size;

    Short_t Reco_Dimuon_sign[1000];
    Short_t Reco_Dimuon_muonPlusIndex[1000];
    Short_t Reco_Dimuon_muonMinusIndex[1000];

    ULong64_t Reco_Dimuon_trig[1000];

    std::vector<float>* Reco_Dimuon_pt = nullptr;
    std::vector<float>* Reco_Dimuon_eta = nullptr;
    std::vector<float>* Reco_Dimuon_rapidity = nullptr;
    std::vector<float>* Reco_Dimuon_phi = nullptr;
    std::vector<float>* Reco_Dimuon_invMass = nullptr;

    std::vector<float>* Reco_Dimuon_muonPtDiff = nullptr;
    std::vector<float>* Reco_Dimuon_muonPtRelDiff = nullptr;

    chain->SetBranchAddress("Reco_Dimuon_size", &Reco_Dimuon_size);

    chain->SetBranchAddress("Reco_Dimuon_sign", Reco_Dimuon_sign);

    chain->SetBranchAddress("Reco_Dimuon_pt", &Reco_Dimuon_pt);
    chain->SetBranchAddress("Reco_Dimuon_eta", &Reco_Dimuon_eta);
    chain->SetBranchAddress("Reco_Dimuon_rapidity", &Reco_Dimuon_rapidity);
    chain->SetBranchAddress("Reco_Dimuon_phi", &Reco_Dimuon_phi);
    chain->SetBranchAddress("Reco_Dimuon_invMass", &Reco_Dimuon_invMass);

    chain->SetBranchAddress("Reco_Dimuon_muonPtDiff", &Reco_Dimuon_muonPtDiff);
    chain->SetBranchAddress("Reco_Dimuon_muonPtRelDiff", &Reco_Dimuon_muonPtRelDiff);

    chain->SetBranchAddress("Reco_Dimuon_muonPlusIndex", Reco_Dimuon_muonPlusIndex);
    chain->SetBranchAddress("Reco_Dimuon_muonMinusIndex", Reco_Dimuon_muonMinusIndex);

    chain->SetBranchAddress("Reco_Dimuon_trig", Reco_Dimuon_trig);

    //Muon-level variables
    Short_t Reco_Muon_size;

    std::vector<float>* Reco_Muon_pt = nullptr;
    //std::vector<float>* Reco_Muon_ptErrTrk = nullptr;
    std::vector<float>* Reco_Muon_eta = nullptr;
    std::vector<float>* Reco_Muon_phi = nullptr;
    std::vector<float>* Reco_Muon_mass = nullptr;

    ULong64_t Reco_Muon_trig[1000];
    Bool_t Reco_Muon_isTightCutBased[1000];
    
    chain->SetBranchAddress("Reco_Muon_size", &Reco_Muon_size);

    chain->SetBranchAddress("Reco_Muon_pt", &Reco_Muon_pt);
    //chain->SetBranchAddress("Reco_Muon_ptErrTrk", &Reco_Muon_ptErrTrk);
    chain->SetBranchAddress("Reco_Muon_eta", &Reco_Muon_eta);
    chain->SetBranchAddress("Reco_Muon_phi", &Reco_Muon_phi);
    chain->SetBranchAddress("Reco_Muon_mass", &Reco_Muon_mass);

    chain->SetBranchAddress("Reco_Muon_trig", Reco_Muon_trig);
    chain->SetBranchAddress("Reco_Muon_isTightCutBased", Reco_Muon_isTightCutBased);
    //

    TString flag = dataset.name;
    TH1D* hist_muonPtRelDiff = new TH1D("hist_muonPtRelDiff"+flag, "", 20, -1., 1.);

    //Event-level cut values.
    const int minCentrality = MIN_CENTRALITY;
    const int maxCentrality = MAX_CENTRALITY; //Cent max = maxCentrality/2.
    const double maxZvtx = MAX_ZVTX;

    //Z-level cut values.
    const double minZ_Mass = MINZ_MASS; 
    const double maxZ_Mass = MAXZ_MASS; 
    const double RapidityCutValue = RAPIDITYCUTVALUE;

    //Muon-level cut values.
    const float EtaCutValue = ETACUTVALUE;
    const double ptCutValue = PTCUTVALUE;
    bool MuPlIsTight;
    bool MuMiIsTight;

    //Observables
    TLorentzVector Z;
    double ptplus, ptminus;
    double etaplus, etaminus;

    //New trigger selection 'L2SingleMu12'
    ULong64_t triggerBit = 1ULL << 7;

    for(Long64_t i = 0; i < nEvents; ++i){//Loop through all EVENTS in the CHAIN.

        chain->GetEntry(i); //Get event i.

        //Good event selection
        bool goodVertex = (std::abs(zVtx) < maxZvtx);
        bool goodCent = true;

        if (dataset.system == CollisionSystem::PbPb2024) {
            goodCent = (Centrality > minCentrality && Centrality < maxCentrality);
        }

        if (!goodVertex) continue;
        if (!goodCent) continue;

        for(Short_t j = 0; j < Reco_Dimuon_size; ++j){ //Loop through all reco dimuon candidates of event i.
            
            //Good Z selection
            bool goodMass = (Reco_Dimuon_invMass->at(j) > minZ_Mass && Reco_Dimuon_invMass->at(j) < maxZ_Mass);
            bool goodRapidity = (std::abs(Reco_Dimuon_rapidity->at(j)) < RapidityCutValue);
            bool goodCharge = (Reco_Dimuon_sign[j] == 0);
            
            bool isTriggerMatched = true;

                if (dataset.system == CollisionSystem::PbPb2024){
                    isTriggerMatched = (Reco_Dimuon_trig[j]) & triggerBit;
                }

            if (!goodMass) continue;
            if (!goodRapidity) continue;
            if (!goodCharge) continue;
            if (!isTriggerMatched) continue;

            Short_t muonPlusIndex = Reco_Dimuon_muonPlusIndex[j];
            Short_t muonMinusIndex = Reco_Dimuon_muonMinusIndex[j];

            //Good muon selection
            ptplus = Reco_Muon_pt->at(muonPlusIndex); //pT of antimuon.
            ptminus = Reco_Muon_pt->at(muonMinusIndex); //pT of corresponding muon.
            etaplus = Reco_Muon_eta->at(muonPlusIndex); //Pseudorapidity of antimuon.
            etaminus = Reco_Muon_eta->at(muonMinusIndex); //Pseudorapidity of corresponding muon.
            MuPlIsTight = Reco_Muon_isTightCutBased[muonPlusIndex];
            MuMiIsTight = Reco_Muon_isTightCutBased[muonMinusIndex];
            
            bool goodMuPl = (ptplus > ptCutValue)
                            && (std::abs(etaplus) < EtaCutValue)
                            && (MuPlIsTight);

            bool goodMuMi = (ptminus > ptCutValue)
                            && (std::abs(etaminus) < EtaCutValue)
                            && (MuMiIsTight);

            if (!goodMuPl || !goodMuMi) continue;
        
            //Fill histograms for muon pT distributions
            hist_muonPtRelDiff->Fill(Reco_Dimuon_muonPtRelDiff->at(j));
        }
    }//Exiting event-by-event loop.

    return hist_muonPtRelDiff;
}

TH1D* MakeChargeAsymmetryHist(const Dataset& dataset){

    // Load root file.
    std::string fullPath = dataset.basePath + dataset.filePattern;

    TChain *chain = new TChain(dataset.treeName.c_str());
    chain->Add(fullPath.c_str());

    std::cout << "> Number of files = " << chain->GetListOfFiles()->GetEntries() << "\n" << std::endl;

    std::cout << "> Opening files " << fullPath << "\n" << std::endl;

    std::cout << "> Running function " << __func__ << " on " << dataset.name << "\n" << std::endl;
    
    //Total number of events on Tree.
    Long64_t nEvents = chain->GetEntries();

    const int MAX_DIMUON = 1000;
    const int MAX_MUON   = 1000;

    //For now, writing ONLY branches that are relevant to the observable we want to measure.

    //Event-level variables
    Int_t Centrality;
    Float_t  zVtx;

    chain->SetBranchAddress("zVtx", &zVtx);
    if(dataset.system == CollisionSystem::PbPb2024) {
        chain->SetBranchAddress("Centrality", &Centrality);
    }

    //Dimuon-level variables
    Short_t Reco_Dimuon_size;

    Short_t Reco_Dimuon_sign[1000];
    Short_t Reco_Dimuon_muonPlusIndex[1000];
    Short_t Reco_Dimuon_muonMinusIndex[1000];

    ULong64_t Reco_Dimuon_trig[1000];

    std::vector<float>* Reco_Dimuon_pt = nullptr;
    std::vector<float>* Reco_Dimuon_eta = nullptr;
    std::vector<float>* Reco_Dimuon_rapidity = nullptr;
    std::vector<float>* Reco_Dimuon_phi = nullptr;
    std::vector<float>* Reco_Dimuon_invMass = nullptr;

    std::vector<float>* Reco_Dimuon_muonPtDiff = nullptr;
    std::vector<float>* Reco_Dimuon_muonPtRelDiff = nullptr;

    chain->SetBranchAddress("Reco_Dimuon_size", &Reco_Dimuon_size);

    chain->SetBranchAddress("Reco_Dimuon_sign", Reco_Dimuon_sign);

    chain->SetBranchAddress("Reco_Dimuon_pt", &Reco_Dimuon_pt);
    chain->SetBranchAddress("Reco_Dimuon_eta", &Reco_Dimuon_eta);
    chain->SetBranchAddress("Reco_Dimuon_rapidity", &Reco_Dimuon_rapidity);
    chain->SetBranchAddress("Reco_Dimuon_phi", &Reco_Dimuon_phi);
    chain->SetBranchAddress("Reco_Dimuon_invMass", &Reco_Dimuon_invMass);

    chain->SetBranchAddress("Reco_Dimuon_muonPtDiff", &Reco_Dimuon_muonPtDiff);
    chain->SetBranchAddress("Reco_Dimuon_muonPtRelDiff", &Reco_Dimuon_muonPtRelDiff);

    chain->SetBranchAddress("Reco_Dimuon_muonPlusIndex", Reco_Dimuon_muonPlusIndex);
    chain->SetBranchAddress("Reco_Dimuon_muonMinusIndex", Reco_Dimuon_muonMinusIndex);

    chain->SetBranchAddress("Reco_Dimuon_trig", Reco_Dimuon_trig);

    //Muon-level variables
    Short_t Reco_Muon_size;

    std::vector<float>* Reco_Muon_pt = nullptr;
    //std::vector<float>* Reco_Muon_ptErrTrk = nullptr;
    std::vector<float>* Reco_Muon_eta = nullptr;
    std::vector<float>* Reco_Muon_phi = nullptr;
    std::vector<float>* Reco_Muon_mass = nullptr;

    ULong64_t Reco_Muon_trig[1000];
    Bool_t Reco_Muon_isTightCutBased[1000];
    
    chain->SetBranchAddress("Reco_Muon_size", &Reco_Muon_size);

    chain->SetBranchAddress("Reco_Muon_pt", &Reco_Muon_pt);
    //chain->SetBranchAddress("Reco_Muon_ptErrTrk", &Reco_Muon_ptErrTrk);
    chain->SetBranchAddress("Reco_Muon_eta", &Reco_Muon_eta);
    chain->SetBranchAddress("Reco_Muon_phi", &Reco_Muon_phi);
    chain->SetBranchAddress("Reco_Muon_mass", &Reco_Muon_mass);

    chain->SetBranchAddress("Reco_Muon_trig", Reco_Muon_trig);
    chain->SetBranchAddress("Reco_Muon_isTightCutBased", Reco_Muon_isTightCutBased);
    //

    TH1D* histMuPl = new TH1D("histMuPl", "", 100, 0., 100.);
    TH1D* histMuMi = new TH1D("histMuMi", "", 100, 0., 100.);

    //Event-level cut values.
    const int minCentrality = MIN_CENTRALITY;
    const int maxCentrality = MAX_CENTRALITY; //Cent max = maxCentrality/2.
    const double maxZvtx = MAX_ZVTX;

    //Z-level cut values.
    const double minZ_Mass = MINZ_MASS; 
    const double maxZ_Mass = MAXZ_MASS; 
    const double RapidityCutValue = RAPIDITYCUTVALUE;

    //Muon-level cut values.
    const float EtaCutValue = ETACUTVALUE;
    const double ptCutValue = PTCUTVALUE;
    bool MuPlIsTight;
    bool MuMiIsTight;

    //Observables
    TLorentzVector Z;
    double ptplus, ptminus;
    double etaplus, etaminus;

    //New trigger selection 'L2SingleMu12'
    ULong64_t triggerBit = 1ULL << 7;

    for(Long64_t i = 0; i < nEvents; ++i){//Loop through all EVENTS in the CHAIN.

        chain->GetEntry(i); //Get event i.

        //Good event selection
        bool goodVertex = (std::abs(zVtx) < maxZvtx);
        bool goodCent = true;

        if (dataset.system == CollisionSystem::PbPb2024) {
            goodCent = (Centrality > minCentrality && Centrality < maxCentrality);
        }

        if (!goodVertex) continue;
        if (!goodCent) continue;

        for(Short_t j = 0; j < Reco_Dimuon_size; ++j){ //Loop through all reco dimuon candidates of event i.
            
            //Good Z selection
            bool goodMass = (Reco_Dimuon_invMass->at(j) > minZ_Mass && Reco_Dimuon_invMass->at(j) < maxZ_Mass);
            bool goodRapidity = (std::abs(Reco_Dimuon_rapidity->at(j)) < RapidityCutValue);
            bool goodCharge = (Reco_Dimuon_sign[j] == 0);
            
            bool isTriggerMatched = true;

                if (dataset.system == CollisionSystem::PbPb2024){
                    isTriggerMatched = (Reco_Dimuon_trig[j]) & triggerBit;
                }

            if (!goodMass) continue;
            if (!goodRapidity) continue;
            if (!goodCharge) continue;
            if (!isTriggerMatched) continue;

            Short_t muonPlusIndex = Reco_Dimuon_muonPlusIndex[j];
            Short_t muonMinusIndex = Reco_Dimuon_muonMinusIndex[j];

            //Good muon selection
            ptplus = Reco_Muon_pt->at(muonPlusIndex); //pT of antimuon.
            ptminus = Reco_Muon_pt->at(muonMinusIndex); //pT of corresponding muon.
            etaplus = Reco_Muon_eta->at(muonPlusIndex); //Pseudorapidity of antimuon.
            etaminus = Reco_Muon_eta->at(muonMinusIndex); //Pseudorapidity of corresponding muon.
            MuPlIsTight = Reco_Muon_isTightCutBased[muonPlusIndex];
            MuMiIsTight = Reco_Muon_isTightCutBased[muonMinusIndex];
            
            bool goodMuPl = (ptplus > ptCutValue)
                            && (std::abs(etaplus) < EtaCutValue)
                            && (MuPlIsTight);

            bool goodMuMi = (ptminus > ptCutValue)
                            && (std::abs(etaminus) < EtaCutValue)
                            && (MuMiIsTight);

            if (!goodMuPl || !goodMuMi) continue;
        
            //Fill histograms for muon pT distributions
            histMuPl->Fill(ptplus);
            histMuMi->Fill(ptminus);
        }
    }//Exiting event-by-event loop.

    //Create a new histogram for charge asymmetry
    TH1D* hist_ChargeAsymmetry = (TH1D*)histMuPl->Clone("hist_ChargeAsymmetry");

    for (int i = 1; i <= histMuPl->GetNbinsX(); i++) {
        
        double NPlus  = histMuPl->GetBinContent(i);
        double NMinus = histMuMi->GetBinContent(i);
        double Asymmetry = 0.0;
        
        if (NPlus + NMinus != 0) Asymmetry = (NPlus - NMinus) / (NPlus + NMinus);
        
        hist_ChargeAsymmetry->SetBinContent(i, Asymmetry);
    }
    
    delete histMuPl;
    delete histMuMi;
    delete chain;

    return hist_ChargeAsymmetry;
}

void PlotChargeAsymmetry(TH1D* histPbPb, TH1D* histppRef){

    //Plots charge asymmetry histograms of PbPb2024 and ppRef2024 datasets.
    TCanvas *c = new TCanvas("c", "Charge Asymmetry", 800, 800);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.30, 1, 1.0);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.00, 1, 0.30);
    basicCanvasFormatting(c, pad1, pad2);
    pad1->SetLeftMargin(0.13);
    pad2->SetLeftMargin(0.13);
    pad1->Draw();
    pad2->Draw();

    //Top pannel
    pad1->cd();
    basicHistFormatting(histPbPb, false);
    basicHistFormatting(histppRef, false);

    basicHistFormatting(histPbPb);
    basicHistFormatting(histppRef);

    //Calculating histogram statistics
    double AsymPbPbMean = histPbPb->GetMean();
    double AsymPbPbMean_err = histPbPb->GetMeanError();
    double AsymppRefMean = histppRef->GetMean();
    double AsymppRefMean_err = histppRef->GetMeanError();

    //Normalization
    //histPbPb->Scale(1./histPbPb->Integral());
    //histppRef->Scale(1./histppRef->Integral());

    //Histogram formatting
    histPbPb->SetMarkerStyle(24);
    histPbPb->SetMarkerSize(0.8);
    histPbPb->SetMarkerColor(kBlack);
    histPbPb->SetLineColor(kBlack);

    histppRef->SetMarkerStyle(24);
    histppRef->SetMarkerSize(0.8);
    histppRef->SetMarkerColor(kRed);
    histppRef->SetLineColor(kRed);

    histPbPb->GetYaxis()->CenterTitle(true);
    histPbPb->GetYaxis()->SetTitle("A(p_{T})");
    histPbPb->GetYaxis()->SetTitleOffset(1.2);

    histPbPb->Draw("P");
    histppRef->Draw("P SAME");

    //Write histogram statistics
    TLegend *info = new TLegend(0.2, 0.12, 0.28, 0.24);
    info->SetBorderSize(0);
    info->SetFillStyle(0);
    info->SetTextFont(42);
    info->SetTextSize(0.03);
    info->SetMargin(0.0);
    info->SetEntrySeparation(0.02);
    info->AddEntry((TObject*)nullptr,
                Form("#LT A_{PbPb} #GT = %.1f #pm %.1f GeV",
                        AsymPbPbMean, AsymPbPbMean_err), "");
    info->AddEntry((TObject*)nullptr,
                Form("#LT A_{ppRef} #GT = %.1f #pm %.1f GeV",
                        AsymppRefMean, AsymppRefMean_err), "");
    info->Draw();

    TLegend *leg = new TLegend(0.42, 0.71, 0.62, 0.85);
    basicLegendFormatting(leg);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.036);
    leg->AddEntry(histPbPb, "PbPb2024", "p");
    leg->AddEntry(histppRef, "ppRef2024", "p");
    leg->Draw();

    pad1->Update();

    //Bottom pannel
    pad2->cd();

    TH1D* histDiff = (TH1D*)histPbPb->Clone("histDiff");
    histDiff->Add(histppRef, -1.0);//Subtraction bin-by-bin.
    
    basicHistFormatting(histDiff, true);
    histDiff->SetMarkerStyle(24);
    histDiff->SetMarkerSize(0.8);
    histDiff->SetMarkerColor(kBlack);
    histDiff->SetLineColor(kBlack);
    histDiff->GetYaxis()->SetTitle("A_{PbPb} - A_{ppRef}");
    histDiff->GetXaxis()->SetTitle("p_{T} [GeV]");
    histDiff->GetYaxis()->SetTitleOffset(0.6);
    histDiff->Draw("P");

    //A horizontal line just to guide the eye.
    TLine *line = new TLine(histDiff->GetXaxis()->GetXmin(), 0.0, histDiff->GetXaxis()->GetXmax(), 0.0);
    line->SetLineColor(kMagenta+2);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");

    pad2->Update();

    //Back to the canvas
    c->cd();
    drawLatexText("#bf{CMS}", 0.12, 0.95, 0.038);
    drawLatexText("#it{Work in Progress}", 0.2, 0.95, 0.03);
    drawLatexText("PbPb 2024, ppRef 2024 (5.36 TeV)", 0.52, 0.95, 0.03);

    TString output = "ChargeAsymmetry_PbPb2024_vs_ppRef2024" + plot_extension;
    c->Update();
    c->SaveAs(output);

    delete histPbPb;
    delete histppRef;
    delete histDiff;
    delete line;
    delete c;
}

void basicCanvasFormatting(TCanvas* c, TPad* pad1, TPad* pad2){
    // Canvas
    c->SetFillColor(0);
    c->SetFrameFillColor(0);
    c->SetTickx(1);
    c->SetTicky(1);

    // Top pad
    pad1->SetLeftMargin(0.098);
    pad1->SetRightMargin(0.036);
    pad1->SetTopMargin(0.08);
    pad1->SetBottomMargin(0.02);

    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad1->SetFillColor(0);
    pad1->SetFrameFillColor(0);
    pad1->SetFrameLineWidth(1);

    // Bottom pad
    pad2->SetLeftMargin(0.098);
    pad2->SetRightMargin(0.036);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.30);

    pad2->SetTickx(1);
    pad2->SetTicky(0);
    pad2->SetFillColor(0);
    pad2->SetFrameFillColor(0);
    pad2->SetFrameLineWidth(1);
}

void basicHistFormatting(TH1D* hist, bool isRatio = false){

    hist->SetTitle("");
    hist->SetStats(0);

    hist->GetXaxis()->CenterTitle(false);
    hist->GetYaxis()->CenterTitle(false);

    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);

    if (!isRatio) {

        // Main plot
        hist->GetXaxis()->SetTitleSize(0.055);
        hist->GetYaxis()->SetTitleSize(0.046);

        hist->GetXaxis()->SetLabelSize(0.0);  // Hide x labels
        hist->GetXaxis()->SetTitleSize(0.0);  // Hide x title

        hist->GetYaxis()->SetLabelSize(0.036);

        hist->GetYaxis()->SetTitleOffset(0.82);

    } else {

        // Ratio plot
        hist->GetYaxis()->CenterTitle(true);
        hist->GetXaxis()->SetTitleSize(0.1);
        hist->GetYaxis()->SetTitleSize(0.098);

        hist->GetXaxis()->SetLabelSize(0.08);
        hist->GetYaxis()->SetLabelSize(0.08);

        hist->GetXaxis()->SetTitleOffset(1.1);
        hist->GetYaxis()->SetTitleOffset(0.4);

        hist->GetXaxis()->SetTickLength(0.03);
        hist->GetYaxis()->SetTickLength(0.025);
    }

    hist->SetMarkerStyle(20);
}

void basicLegendFormatting(TLegend* leg){
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.044);
    leg->SetMargin(0.2);
    leg->SetEntrySeparation(0.04);
}

void drawLatexText(TString latexText, double x, double y, double TextSize){
// We can add any pT or eta selection. If no text is passed to the function the CMS Header will be drawn.
TLatex latex;
latex.SetNDC();              // For normalized coordinates
latex.SetTextSize(TextSize);
latex.SetTextFont(42);       // Helvetica
latex.SetTextAlign(11);      // Left-top aligned.
latex.DrawLatex(x, y, latexText);
}

void PrintBinInfo(TH1D* hist)
{
    std::cout << "\n=== Histogram: " << hist->GetName() << " ===\n";

    for (int bin = 40; bin <= 43; ++bin) {

        double content = hist->GetBinContent(bin);
        double error   = hist->GetBinError(bin);

        std::cout << "Bin " << bin << "\n";
        std::cout << "  Content      = " << content << "\n";
        std::cout << "  ROOT error   = " << error << "\n";
        std::cout << "  sqrt(N)      = " << std::sqrt(content) << "\n";
        std::cout << "  Error^2      = " << error * error << "\n";
        std::cout << std::endl;
    }
}