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

The initial centrality bins would be:
Set 1:
- 0-10%
- 10-30%
- 30-50%
- 50-70%
- 70-100%

Or maybe another temptative list is:
Set 2:
- 0-10%
- 10-20%
- 20-30%
- 30-40%
- 40-100%

pT sectors:
- 20-30 GeV
- 30-40 GeV
- 40-50 GeV
- 50-60 GeV

Another thing to observe is that we don't need to analyze the whole pT spectrum.
The proposed effect is in the 30-50 GeV region, so we can focus on that region and maybe extend it a little bit to 20-60 GeV.

The observables that i was thinking about were:

5- Peak and mean of PbPb and ppRef dN/dpT distributions as function of centrality (i.e. per centrality bin).

The asymmetry A has the same binnage as the pT distributions, but i think it could be interesting
to vary the pT binning and see if the asymmetry is more pronounced in a specific pT region.
For example, we could check it per pT sector.
One possibility to study is to vary the bin width for A.

If there is a measurable asymmetry between the pT distributions, it should be only apparent in PbPb mid-central collisions.

That means A(ppRef) should be approximately 0 for all pT sectors, while A(PbPb) should have asymmetry in some pT sector,
thus producing a deviation ΔA = A(PbPb) - A(ppRef) non-zero in that pT sector.

The double ratio is also an option, even to plot on the same canvas, since it uses the same pT binning as the asymmetry A.

The other observable we have in hand is the event-by-event relative difference in pT between the two muons, namely

            ΔpT = (pT(mu+) - pT(mu-)) / (pT(mu+) + pT(mu-)).

This could be a valuable observable, but it needs to be treated with care.
A particularly natural observable would be like 

            <ΔpT(PbPb)> - <ΔpT(ppRef)>,

where the mean is calculated per centrality bin.

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
std::string BasePath = "/home/lucas/Documents/CMS/z_boson_analysis/"; //IFT
//std::string BasePath = "/home/lucasdoriac/z_boson_analysis/data/"; //Home

const double MAX_ZVTX = 15.0;

const double MINZ_MASS = 60.;
const double MAXZ_MASS = 120.;
const double RAPIDITYCUTVALUE = 2.4;

const float ETACUTVALUE = 2.4;
const double PTCUTVALUE = 20.;

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

//Centrality bins for PbPb2024 data
std::vector<std::pair<int, int>> CentralitySet = {
    {0, 10},
    {10, 20},
    {20, 30},
    {30, 100},
};

std::vector<std::pair<double, double>> pT_sectors = {
    {20., 30.},
    {30., 40.},
    {40., 50.},
    {50., 60.},
    {60., 70.}
};

//---Function declarations
std::array<double, 3> CalculatePeakAndMeanDiff(const Dataset& dataset, const std::pair<int, int>& centRange, TFile* outputFile = nullptr);
void PlotDeltaPeakAndMeanDiff();
void CalculateAsymmetryPt(TFile* outputFile, const std::vector<std::pair<int, int>>& centRanges, const std::vector<std::pair<double, double>>& pTRanges);
void PlotFinalAsymmetry(TFile* outputFile, const std::vector<std::pair<int, int>>& centRanges);

//---Histogram formatting
void basicCanvasFormatting(TCanvas *c, TPad *pad1, TPad *pad2);
void basicHistFormatting(TH1D *hist, bool isRatio = false);
void basicLegendFormatting(TLegend *leg);
void drawLatexText(TString latexText = "#bf{CMS}", double x = 0.15, double y = 0.93, double TextSize = 0.05);

//---Main function
void MakeAll_v2(){

    gROOT->SetBatch(kTRUE);

    TFile* outputFile = new TFile("SingleMuonChargedYields.root", "RECREATE");
    //datasets[1] = PbPb2024
    //datasets[2] = ppRef2024

    //Calculate the peak and mean difference for ppRef2024 dataset, which will be used as reference for PbPb2024 dataset.
    auto values = CalculatePeakAndMeanDiff(datasets[2], {0., 100.}, outputFile);
    double Peak_diff_ppRef = values[0];
    double Mean_diff_ppRef = values[1];
    double Mean_diff_err_ppRef = values[2];

    for(Short_t cBin = 0; cBin < CentralitySet.size(); ++cBin){

        std::cout << "\n -> Centrality bin: " << CentralitySet[cBin].first 
        << "-" << CentralitySet[cBin].second << "%" << std::endl;

        values = CalculatePeakAndMeanDiff(datasets[1], CentralitySet[cBin], outputFile);
        double Peak_diff_PbPb = values[0];
        double Mean_diff_PbPb = values[1];
        double Mean_diff_err_PbPb = values[2];

        DeltaPeak_pT.push_back(Peak_diff_PbPb - Peak_diff_ppRef);
        DeltaMean_pT.push_back(Mean_diff_PbPb - Mean_diff_ppRef);
        DeltaMean_err_pT.push_back(sqrt(pow(Mean_diff_err_PbPb, 2) + pow(Mean_diff_err_ppRef, 2)));
    }

    PlotDeltaPeakAndMeanDiff();
    /*std::cout << "DeltaPeak_pT size = " << DeltaPeak_pT.size() << std::endl;
    std::cout << "DeltaMean_pT size = " << DeltaMean_pT.size() << std::endl;

    for (size_t i = 0; i < DeltaPeak_pT.size(); ++i) {
        std::cout << i
                << "  DeltaPeak = " << DeltaPeak_pT[i]
                << "  DeltaMean = " << DeltaMean_pT[i]
                << std::endl;
    }*/

    CalculateAsymmetryPt(outputFile, CentralitySet, pT_sectors);
    PlotFinalAsymmetry(outputFile, CentralitySet);

    // Clean up
    outputFile->Close();
    DeltaPeak_pT.clear();
    DeltaMean_pT.clear();
    DeltaMean_err_pT.clear();
    delete outputFile;
}

void PlotFinalAsymmetry(TFile* outputFile, const std::vector<std::pair<int, int>>& centRanges){

    TGraphErrors* grFinalAsymmetry =
        new TGraphErrors(centRanges.size());

        grFinalAsymmetry->SetName("grFinalAsymmetry");
        grFinalAsymmetry->SetTitle(
            ";Centrality (%);#LT#Delta A#GT"
        );

    for (size_t cBin = 0; cBin < centRanges.size(); ++cBin) {

        TString histName = Form(
            "DeltaAsym_PbPb2024_%d_%d",
            centRanges[cBin].first,
            centRanges[cBin].second
        );

        TH1D* hDeltaAsymmetry =
            (TH1D*)outputFile->Get(histName);

        if (!hDeltaAsymmetry) {
            std::cerr << "Could not find histogram: "
                      << histName << std::endl;
            continue;
        }

        double fitMin = 20.0;
        double fitMax = 60.0;

        // Constant fit:
        // Delta A(pT) = constant
        TF1* constantFit = new TF1(
            Form("constantFit_%d_%d",
                 centRanges[cBin].first,
                 centRanges[cBin].second),
            "pol0",
            fitMin,
            fitMax
        );

        // Fit histogram
        hDeltaAsymmetry->Fit(constantFit, "RQ");

        // --------------------------------------------------
        // Extract <Delta A>
        // --------------------------------------------------

        double meanDeltaA = constantFit->GetParameter(0);
        double meanError  = constantFit->GetParError(0);

        // --------------------------------------------------
        // Centrality midpoint and width
        // --------------------------------------------------

        double centMin = centRanges[cBin].first;
        double centMax = centRanges[cBin].second;

        double centMid   = 0.5 * (centMin + centMax);
        double centWidth = 0.5 * (centMax - centMin);

        // --------------------------------------------------
        // Add point to TGraphErrors
        // --------------------------------------------------

        grFinalAsymmetry->SetPoint(
            cBin,
            centMid,
            meanDeltaA
        );

        grFinalAsymmetry->SetPointError(
            cBin,
            centWidth,
            meanError
        );

        // --------------------------------------------------
        // Print result
        // --------------------------------------------------

        std::cout << "Centrality "
                  << centMin << "-" << centMax << "% : "
                  << "<Delta A> = "
                  << meanDeltaA << " +/- "
                  << meanError
                  << std::endl;
                
    }

        // --------------------------------------------------
        // Draw final graph
        // --------------------------------------------------

        TCanvas* canvas = new TCanvas(
            "cFinalAsymmetry",
            "Final Delta A",
            800,
            600
        );
        canvas->SetLeftMargin(0.13);
        canvas->SetBottomMargin(0.15);
        canvas->SetRightMargin(0.05);
        canvas->SetTopMargin(0.05);
        canvas->SetTicks(1, 1);

        //Formatting
        grFinalAsymmetry->GetXaxis()->CenterTitle(false);
        grFinalAsymmetry->GetYaxis()->CenterTitle(false);
        grFinalAsymmetry->GetXaxis()->SetTitleOffset(.9);
        grFinalAsymmetry->GetYaxis()->SetTitleOffset(1.);
        grFinalAsymmetry->GetXaxis()->SetTitleFont(42);
        grFinalAsymmetry->GetYaxis()->SetTitleFont(42);
        grFinalAsymmetry->GetXaxis()->SetLabelFont(42);
        grFinalAsymmetry->GetYaxis()->SetLabelFont(42);
        grFinalAsymmetry->GetXaxis()->SetTitleSize(0.055);
        grFinalAsymmetry->GetYaxis()->SetTitleSize(0.055);
        grFinalAsymmetry->GetXaxis()->SetLabelSize(0.042);
        grFinalAsymmetry->GetYaxis()->SetLabelSize(0.042);
        grFinalAsymmetry->SetTitle("");

        grFinalAsymmetry->SetMarkerStyle(21);
        grFinalAsymmetry->SetMarkerSize(1.1);
        grFinalAsymmetry->SetLineWidth(2);
        grFinalAsymmetry->SetLineColor(kRed+1);
        grFinalAsymmetry->SetMarkerColor(kRed+1);

        grFinalAsymmetry->Draw("AP");

        drawLatexText("#bf{CMS}", 0.15, 0.93, 0.05);
        drawLatexText("Work in Progress", 0.3, 0.93, 0.05);
        drawLatexText("PbPb #sqrt{#it{s}_{NN}} = 5.36 TeV, 2024 Data", 0.55, 0.93, 0.05);

        // --------------------------------------------------
        // Horizontal zero line
        // --------------------------------------------------

        double xMin = centRanges.front().first;
        double xMax = centRanges.back().second;

        TLine* zeroLine = new TLine(
            xMin, 0.,
            xMax, 0.
        );

        zeroLine->SetLineStyle(2);
        zeroLine->SetLineWidth(2);
        zeroLine->Draw("same");

        // --------------------------------------------------
        // Save plot
        // --------------------------------------------------

        canvas->SaveAs(
            "FinalDeltaAsymmetry_vs_Centrality.png"
        );

        // --------------------------------------------------
        // Save graph to ROOT file
        // --------------------------------------------------

        outputFile->cd();
        grFinalAsymmetry->Write();
}


void CalculateAsymmetryPt(TFile* outputFile, const std::vector<std::pair<int, int>>& centRanges, const std::vector<std::pair<double, double>>& pTRanges){

    TH1D* histMuPlusRef = (TH1D*)outputFile->Get("histMuPlus_ppRef2024_Data_0_100");
    TH1D* histMuMinusRef = (TH1D*)outputFile->Get("histMuMinus_ppRef2024_Data_0_100");

    TH1D* hAsymmetryRef = new TH1D("ReferenceAsym", "Reference Asymmetry", histMuPlusRef->GetNbinsX(), 
                    histMuPlusRef->GetXaxis()->GetXmin(),
                    histMuPlusRef->GetXaxis()->GetXmax());

    //Get reference values
    for (int i = 1; i <= histMuPlusRef->GetNbinsX(); ++i) {

        double NplusRef  = histMuPlusRef->GetBinContent(i);
        double NminusRef = histMuMinusRef->GetBinContent(i);

        double errPlusRef  = histMuPlusRef->GetBinError(i);
        double errMinusRef = histMuMinusRef->GetBinError(i);

        double denominatorRef = NplusRef + NminusRef;

        if (denominatorRef > 0.) {

            double asymmetryRef = (NplusRef - NminusRef) / denominatorRef;
            double asymmetryErrorRef = 2./(denominatorRef*denominatorRef)*sqrt( pow(NminusRef*errPlusRef,2) + pow(NplusRef*errMinusRef,2));

            hAsymmetryRef->SetBinContent(i, asymmetryRef);
            hAsymmetryRef->SetBinError(i, asymmetryErrorRef);

        }

    }
    hAsymmetryRef->Write();

    //Get muon charge yield histograms from outputFile.
    for (size_t cBin = 0; cBin < centRanges.size(); ++cBin) {

        TString namePlus = Form(
            "histMuPlus_PbPb2024_Data_%d_%d",
            centRanges[cBin].first, centRanges[cBin].second
        );

        TString nameMinus = Form(
            "histMuMinus_PbPb2024_Data_%d_%d",
            centRanges[cBin].first, centRanges[cBin].second
        );

        TH1D* histMuPlus = (TH1D*)outputFile->Get(namePlus);
        TH1D* histMuMinus = (TH1D*)outputFile->Get(nameMinus);
        TH1D* hDeltaAsymmetry = new TH1D(
            Form("DeltaAsym_PbPb2024_%d_%d", centRanges[cBin].first, centRanges[cBin].second),
            Form("#Delta A, %d-%d%%;p_{T} (GeV/#it{c});#Delta A",
                centRanges[cBin].first, centRanges[cBin].second),
            histMuPlus->GetNbinsX(),
            histMuPlus->GetXaxis()->GetXmin(),
            histMuPlus->GetXaxis()->GetXmax()
        );

        //Now calculate for each centrality bin
        for (int i = 1; i <= histMuPlus->GetNbinsX(); ++i) {

            double Nplus  = histMuPlus->GetBinContent(i);
            double Nminus = histMuMinus->GetBinContent(i);

            double errPlus  = histMuPlus->GetBinError(i);
            double errMinus = histMuMinus->GetBinError(i);

            double denominator = Nplus + Nminus;
            double asymmetryRef = hAsymmetryRef->GetBinContent(i);
            double asymmetryErrorRef = hAsymmetryRef->GetBinError(i);

            if (denominator > 0.) {

                double asymmetry = (Nplus - Nminus) / denominator;
                double asymmetryError = 2./(denominator*denominator)*sqrt( pow(Nminus*errPlus,2) + pow(Nplus*errMinus,2));

                double DeltaA = asymmetry - asymmetryRef;
                double deltaAError = sqrt(
                                    asymmetryError * asymmetryError +
                                    asymmetryErrorRef * asymmetryErrorRef
                                );

                hDeltaAsymmetry->SetBinContent(i, DeltaA);
                hDeltaAsymmetry->SetBinError(i, deltaAError);
            }
        }

        hDeltaAsymmetry->Write();

    }//End of centrality bin loop.

}

void PlotDeltaPeakAndMeanDiff()
{
    std::vector<double> centX;
    std::vector<double> centXErr;

    for (const auto& bin : CentralitySet) {
        double center = 0.5 * (bin.first + bin.second);
        double halfWidth = 0.5 * (bin.second - bin.first);

        centX.push_back(center);
        centXErr.push_back(halfWidth);
    }

    /*TGraph* gDeltaPeak = new TGraph(
        centX.size(),
        centX.data(),
        DeltaPeak_pT.data()
    );*/

    TGraphErrors* gDeltaMean = new TGraphErrors(
        centX.size(),
        centX.data(),
        DeltaMean_pT.data(),
        centXErr.data(),
        DeltaMean_err_pT.data()
    );

    TCanvas* c1 = new TCanvas("c1","Delta pT",800,700);
    c1->SetLeftMargin(0.13);
    c1->SetRightMargin(0.04);
    c1->SetBottomMargin(0.13);
    c1->SetTopMargin(0.06);
    c1->SetTicks(1, 1);

    // Graph styles
    /*gDeltaPeak->SetMarkerStyle(20);
    gDeltaPeak->SetMarkerSize(1.1);
    gDeltaPeak->SetMarkerColor(kRed+1);
    gDeltaPeak->SetLineColor(kRed+1);
    gDeltaPeak->SetLineWidth(2);*/

    gDeltaMean->SetMarkerStyle(21);
    gDeltaMean->SetMarkerSize(1.1);
    gDeltaMean->SetMarkerColor(kBlue+1);
    gDeltaMean->SetLineColor(kBlue+1);
    gDeltaMean->SetLineWidth(2);

    gDeltaMean->SetTitle("");
    gDeltaMean->GetXaxis()->SetTitle("Centrality (%)");
    gDeltaMean->GetYaxis()->SetTitle("#LT #Delta p_{T} #GT [GeV]");

    //gDeltaMean->GetXaxis()->SetLimits(0., 100.);
    //gDeltaMean->GetYaxis()->SetRangeUser(-0.05, 0.05);

    gDeltaMean->GetXaxis()->SetTitleSize(0.045);
    gDeltaMean->GetYaxis()->SetTitleSize(0.045);
    gDeltaMean->GetXaxis()->SetLabelSize(0.040);
    gDeltaMean->GetYaxis()->SetLabelSize(0.040);
    gDeltaMean->GetXaxis()->SetTitleOffset(1.05);
    gDeltaMean->GetYaxis()->SetTitleOffset(1.15);

    gDeltaMean->Draw("AP");
    //gDeltaPeak->Draw("P SAME");

    //Zero line
    TLine* zeroLine = new TLine(0., 0.,100., 0.);
    zeroLine->SetLineStyle(2);
    zeroLine->SetLineColor(kGray+2);
    zeroLine->SetLineWidth(1);
    zeroLine->Draw("SAME");

    TLegend* leg = new TLegend(0.55, 0.72,0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.040);
    //leg->AddEntry(gDeltaPeak,"Peak difference","p");
    leg->AddEntry(gDeltaMean,"Mean difference","p");
    leg->Draw();

    c1->Update();
    c1->SaveAs("DeltaPeakAndMeanDiff.png");
}

std::array<double, 3> CalculatePeakAndMeanDiff(const Dataset& dataset, const std::pair<int, int>& centRange, TFile* outputFile){

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
    Float_t Reco_Dimuon_vtxProb[1000];

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
    chain->SetBranchAddress("Reco_Dimuon_vtxProb", Reco_Dimuon_vtxProb);


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

    TString namePlus  = Form("histMuPlus_%s_%d_%d", dataset.name.c_str(),centRange.first, centRange.second);
    TString nameMinus = Form("histMuMinus_%s_%d_%d", dataset.name.c_str(), centRange.first, centRange.second);
    TString titlePlus  = Form("#mu^{+}, %d-%d%%;p_{T} (GeV/#it{c});Counts",
                               centRange.first, centRange.second);
    TString titleMinus = Form("#mu^{-}, %d-%d%%;p_{T} (GeV/#it{c});Counts",
                                centRange.first, centRange.second);

    TH1D* histMuPlus = new TH1D(namePlus, titlePlus, 100, 0., 100.);
    TH1D* histMuMinus = new TH1D(nameMinus, titleMinus, 100, 0., 100.);

    //Event-level cut values.
    const int minCentrality = 2.*centRange.first;
    const int maxCentrality = 2.*centRange.second;
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
            bool goodVtxProb = (Reco_Dimuon_vtxProb[j] > 0.001);

            bool isTriggerMatched = true;

                if (dataset.system == CollisionSystem::PbPb2024){
                    isTriggerMatched = (Reco_Dimuon_trig[j]) & triggerBit;
                }

            if (!goodMass) continue;
            if (!goodRapidity) continue;
            if (!goodCharge) continue;
            if (!goodVtxProb) continue;
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

    histMuPlus->Write();
    histMuMinus->Write();

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

    double Mean_diff = muPlusMean - muMinusMean;
    double Mean_diff_err = sqrt(muPlusMean_err*muPlusMean_err + muMinusMean_err*muMinusMean_err);

    delete histMuPlus;
    delete histMuMinus;
    delete chain;

    return {Peak_diff, Mean_diff, Mean_diff_err};
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