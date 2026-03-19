/*
So, what i want to test here is if we are correctly choosing the dimuon pairs.

We can build a set of tests, including invMassSpectrum.

This is a test event-by-event like Cesar suggested.

USE ONLY SINGLE MUONS AND NOT QQ MUONS.
*/

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
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstring>
#include <vector>

// ### User Settings ###
double eta_LowLimit = -2.4;
double eta_HighLimit = 2.4;
double pt_LowLimit = 14.0; // GeV.
double pt_HighLimit = 100.0; // GeV. Later need to check if this limit is reasonable.
double minDimuonMass = 85.0;
double maxDimuonMass = 95.0; // GeV.
double maxCentrality = 180; // % = what_centrality/2.
double delta = 1e-5;
TString plot_extension = ".pdf";

// ##############################################################################
// ##############################################################################


// --- Structs ---
struct dataFile {
    std::string path;
    std::string name;
    std::string label;
};

const dataFile PbPb_5TeV_2024 = {
    "../data/Oniatree_PbPb2024_PromptReco.root",
    "Oniatree_PbPb2024_PromptReco.root",
    "PbPb_5TeV_2024",
};

/*
const dataFile ppRef_2024 = {
	"../../../OniaTree/Oniatree_ppRef2024_PromptReco.root"	
	"Oniatree_ppRef2024_PromptReco.root",
    "ppRef_2024",
}

const dataFile MC_PbPb_5TeV = {
	"../../../OniaTree/Oniatree_MC_PbPb_5TeV.root"	
	"MC_PbPb_5TeV.root",
    "MC_PbPb_5TeV",
}
*/

// --- ready-to-go functions
void printTreeContents(const std::string& PathToROOTFile, const char* dir = "hionia", const char* tree = "myTree");
void invMassSpectrum(const std::string& PathToROOTFile = "Oniatree_PbPb2024_PromptReco.root", int nBins = 40, int minMass = 70, int maxMass = 110);
void drawLatexText(const char *latexText = "#bf{CMS}", double x = 0.16, double y = 0.93, double TextSize = 0.05);
void basicCanvasFormatting(TCanvas *c);
void basicHistFormatting(TH1D *hist);
void basicLegendFormatting(TLegend *leg);
// --- still writing
void singleMuonTests(const std::string& PathToROOTFile = "Oniatree_PbPb2024_PromptReco.root");

void test(){

	gROOT->SetBatch(kTRUE);

	//printTreeContents(PbPb_5TeV_2024.path);
	//invMassSpectrum(PbPb_5TeV_2024.path);
	singleMuonTests(PbPb_5TeV_2024.path);
}

void singleMuonTests(const std::string& PathToROOTFile){
	// Load root file.
	TFile *rootFile = TFile::Open(PathToROOTFile.c_str(), "READ");
    rootFile->ls();

    // Load directory.
    TDirectoryFile *dir1 = dynamic_cast<TDirectoryFile*>(rootFile->Get("hionia"));

    // Load tree.
	TTree *dimuonTree = dynamic_cast<TTree*>(dir1->Get("myTree"));

	//Total number of events on Tree.
	Long64_t nEvents = dimuonTree->GetEntries();

	dimuonTree->Print();

	Short_t Reco_mu_size;
	dimuonTree->SetBranchAddress("Reco_mu_size", &Reco_mu_size);

	const int MAX_MU = 20;
	Short_t Reco_mu_charge[MAX_MU];
	dimuonTree->SetBranchAddress("Reco_mu_charge", Reco_mu_charge);

	std::vector<float>* Reco_mu_4mom_pt = nullptr;
	std::vector<float>* Reco_mu_4mom_eta = nullptr;
	std::vector<float>* Reco_mu_4mom_phi = nullptr;
	std::vector<float>* Reco_mu_4mom_m = nullptr;
	dimuonTree->SetBranchAddress("Reco_mu_4mom_pt", &Reco_mu_4mom_pt);
	dimuonTree->SetBranchAddress("Reco_mu_4mom_eta", &Reco_mu_4mom_eta);
	dimuonTree->SetBranchAddress("Reco_mu_4mom_phi", &Reco_mu_4mom_phi);
	dimuonTree->SetBranchAddress("Reco_mu_4mom_m", &Reco_mu_4mom_m);

	Short_t Reco_QQ_size;
	dimuonTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);

	const int MAX_QQ = 10;
	Short_t Reco_QQ_mupl_idx[MAX_QQ];
	Short_t Reco_QQ_mumi_idx[MAX_QQ];
	dimuonTree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
	dimuonTree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);

	std::vector<float>* Reco_QQ_4mom_pt = nullptr;
	std::vector<float>* Reco_QQ_4mom_eta = nullptr;
	std::vector<float>* Reco_QQ_4mom_phi = nullptr;
	std::vector<float>* Reco_QQ_4mom_m = nullptr;
	dimuonTree->SetBranchAddress("Reco_QQ_4mom_pt", &Reco_QQ_4mom_pt);
	dimuonTree->SetBranchAddress("Reco_QQ_4mom_eta", &Reco_QQ_4mom_eta);
	dimuonTree->SetBranchAddress("Reco_QQ_4mom_phi", &Reco_QQ_4mom_phi);
	dimuonTree->SetBranchAddress("Reco_QQ_4mom_m", &Reco_QQ_4mom_m);

	//TH2D *hist_ptplmi = new TH2D("hist_pt_plmi", "p_{T}(#mu^{+}) vs p_{T}(#mu^{-})", 100, 0, 100, 100, 0, 100);
	TH1D *hist_ptpl = new TH1D("hist_pt_pl", "p_{T}(#mu^{+})", 100, 0, 100);
	TH1D *hist_ptmi = new TH1D("hist_pt_mi", "p_{T}(#mu^{i})", 100, 0, 100);
	double ptplus, ptminus;

	for(Long64_t i = 0; i < nEvents; ++i){//Loop events. "event-by-event"...

		dimuonTree->GetEntry(i);//Get BRANCH values at the i-th event.
		for(Short_t j = 0; j < Reco_QQ_size; ++j){//If event has dimuon candidate, flagged by Reco_QQ_size > 0...

			if(Reco_QQ_4mom_m->at(j) > minDimuonMass && Reco_QQ_4mom_m->at(j) < maxDimuonMass){//Note that i am getting the mass not from Reco_mu_4mom_m but from QQ...
				
				ptplus = Reco_mu_4mom_pt->at(Reco_QQ_mupl_idx[j]); //pT of antimuon...
				ptminus = Reco_mu_4mom_pt->at(Reco_QQ_mumi_idx[j]); //pT of corresponding muon...

				hist_ptpl->Fill(ptplus);
				hist_ptmi->Fill(ptminus);
			}
			
		}

	}//Exiting event-by-event loop.

	//Calculate muon and antimuon yield...
	double mu_yield_pl, mu_yield_pl_err;
	double mu_yield_mi, mu_yield_mi_err;
	mu_yield_pl = hist_ptpl->IntegralAndError(1, 100, mu_yield_pl_err);
	mu_yield_mi = hist_ptmi->IntegralAndError(1, 100, mu_yield_mi_err);

	//Calculate mean of pT distributions (should be around ~ 40 GeV...)
	double mean_pl = hist_ptpl->GetMean();
	double mean_mi = hist_ptmi->GetMean();
	double mean_pl_err = hist_ptpl->GetMeanError();
	double mean_mi_err = hist_ptmi->GetMeanError();

	//Calculate standard deviation of pT distributions (should be around ?)
	double std_pl = hist_ptpl->GetStdDev();
	double std_mi = hist_ptmi->GetStdDev();
	double std_pl_err = hist_ptpl->GetStdDevError();
	double std_mi_err = hist_ptmi->GetStdDevError();

	//Draw lines at respective peaks.
	int binMaxPl = hist_ptpl->GetMaximumBin();
	int binMaxMi = hist_ptmi->GetMaximumBin();
	double xMaxPl = hist_ptpl->GetBinCenter(binMaxPl);
	double xMaxMi = hist_ptmi->GetBinCenter(binMaxMi);
	double ymaxPl = hist_ptpl->GetMaximum();
	double ymaxMi = hist_ptmi->GetMaximum();

	TLine *linePl = new TLine(xMaxPl, 0.0, xMaxPl, ymaxPl);
	TLine *lineMi = new TLine(xMaxMi, 0.0, xMaxMi, ymaxMi);

	linePl->SetLineColor(kRed);
	linePl->SetLineWidth(2);
	linePl->SetLineStyle(2);

	lineMi->SetLineColor(kBlue);
	lineMi->SetLineWidth(2);
	lineMi->SetLineStyle(2);

	TCanvas *myCanvas = new TCanvas("c_ptplmi", "", 700, 600);
	basicCanvasFormatting(myCanvas);
	basicHistFormatting(hist_ptpl);
	basicHistFormatting(hist_ptmi);

    //hist_ptpl->GetXaxis()->SetTitle("p_{T}^{#mu^{+}} (GeV)");
	//hist_ptpl->GetYaxis()->SetTitle("Number of events");
	hist_ptpl->SetLineColor(kRed + 1);
	hist_ptpl->SetLineWidth(2);
	hist_ptpl->SetFillStyle(3004);
	hist_ptpl->SetFillColorAlpha(kRed + 1, 0.35);
	hist_ptmi->SetLineColor(kBlue + 1);
	hist_ptmi->SetLineWidth(2);
	hist_ptmi->SetFillStyle(3005);
	hist_ptmi->SetFillColorAlpha(kBlue + 1, 0.35);
    hist_ptpl->GetXaxis()->SetTitle("p_{T} (GeV)");
	hist_ptpl->GetYaxis()->SetTitle("Number of Z^{0} candidates (GeV)^{-1}");
	//Zooming in for display only.
	hist_ptpl->GetXaxis()->SetRangeUser(10, 100);
	hist_ptpl->Draw("HIST");
	hist_ptmi->Draw("HIST SAME");

	linePl->Draw("SAME");
	lineMi->Draw("SAME");

	drawLatexText();
	drawLatexText("#it{Internal}", 0.25, 0.93, 0.042);
	drawLatexText("#it{partial} PbPb 2024 (5.36 TeV)", 0.55, 0.93, 0.042);
	//drawLatexText(Form("Z yield = %d #pm %d", TMath::Nint(z_yield), TMath::Nint(z_yield_err)), 0.18, 0.8, 0.04);
	//drawLatexText(Form("m_{#mu^{+}#mu^{-}} = (%.2f #pm %.2f) GeV", mean, mean_err), 0.18, 0.72, 0.04);
	//drawLatexText(Form("#sigma_{m_{#mu^{+}#mu^{-}}} = (%.2f #pm %.2f) GeV", std_dev, std_dev_err), 0.18, 0.64, 0.04);
	drawLatexText(Form("%d < m_{#mu^{+}#mu^{-}} < %d GeV", TMath::Nint(minDimuonMass), TMath::Nint(maxDimuonMass)), 0.54, 0.8, 0.042);
	//drawLatexText(Form("|#eta| < %.1f", eta_HighLimit), 0.23, 0.33, 0.04);

	TLegend* leg = new TLegend(0.54, 0.45, 0.87, 0.75);
	basicLegendFormatting(leg);
	leg->AddEntry(hist_ptpl, Form("#LT p_{T}(#mu^{+}) #GT = %.2f #pm %.2f", mean_pl, mean_pl_err),"f");
	leg->AddEntry(hist_ptmi, Form("#LT p_{T}(#mu^{-}) #GT = %.2f #pm %.2f", mean_mi, mean_mi_err),"f");
	leg->AddEntry(hist_ptpl, Form("#mu^{+} peak = %.1f", xMaxPl), "l");
	leg->AddEntry(hist_ptmi, Form("#mu^{-} peak = %.1f", xMaxMi), "l");
	leg->Draw();

	TString output = TString(__func__) + plot_extension;

	myCanvas->Update();
	myCanvas->SaveAs(output);

	rootFile->Close();
}

void invMassSpectrum(const std::string& PathToROOTFile, int nBins, int minMass, int maxMass){

	// Load root file.
	TFile *rootFile = TFile::Open(PathToROOTFile.c_str(), "READ");
    rootFile->ls();

    // Load directory.
    TDirectoryFile *dir1 = dynamic_cast<TDirectoryFile*>(rootFile->Get("hionia"));

    // Load tree.
	TTree *dimuonTree = dynamic_cast<TTree*>(dir1->Get("myTree"));

	//Total number of events on Tree.
	Long64_t nEvents = dimuonTree->GetEntries();

	Short_t Reco_QQ_size;
	dimuonTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);

	int Centrality;
	dimuonTree->SetBranchAddress("Centrality", &Centrality);

	const int MAX_QQ = 10;
	Short_t Reco_QQ_mupl_idx[MAX_QQ];
	Short_t Reco_QQ_mumi_idx[MAX_QQ];
	dimuonTree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
	dimuonTree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);
	
	std::vector<float>* Reco_QQ_4mom_pt = nullptr;
	std::vector<float>* Reco_QQ_4mom_eta = nullptr;
	std::vector<float>* Reco_QQ_4mom_phi = nullptr;
	std::vector<float>* Reco_QQ_4mom_m = nullptr;
	dimuonTree->SetBranchAddress("Reco_QQ_4mom_pt", &Reco_QQ_4mom_pt);
	dimuonTree->SetBranchAddress("Reco_QQ_4mom_eta", &Reco_QQ_4mom_eta);
	dimuonTree->SetBranchAddress("Reco_QQ_4mom_phi", &Reco_QQ_4mom_phi);
	dimuonTree->SetBranchAddress("Reco_QQ_4mom_m", &Reco_QQ_4mom_m);

	TH1D *h_invMass = new TH1D("h_invMass","Dimuon invariant mass", nBins, minMass, maxMass);
	h_invMass->SetDirectory(0);

	for(Long64_t i = 0; i < nEvents; ++i){

		dimuonTree->GetEntry(i);

		if (Reco_QQ_size > 0){//If the event has a dimuon candidate, runs the loop, if not, go to next event.
			for (size_t j = 0; j < Reco_QQ_4mom_m->size(); ++j){
				h_invMass->Fill(Reco_QQ_4mom_m->at(j));
			}
		}

	}//Exiting event-by-event loop...

	//Set specified mass range for Z candidates...
	h_invMass->GetXaxis()->SetRangeUser(minDimuonMass, maxDimuonMass);
	//Calculate Z boson yield...
	int minMass_bin = h_invMass->FindBin(minDimuonMass + delta);
	int maxMass_bin = h_invMass->FindBin(maxDimuonMass - delta);
	double z_yield, z_yield_err;
	z_yield = h_invMass->IntegralAndError(minMass_bin, maxMass_bin, z_yield_err);

	//Calculate mean of invMass distribution (should be around ~ 91 GeV...)
	double mean = h_invMass->GetMean();
	double mean_err = h_invMass->GetMeanError();

	//Calculate standard deviation of invMass distribution (should be around ~2GeV...)
	double std_dev = h_invMass->GetRMS();
	double std_dev_err = h_invMass->GetRMSError();

	TCanvas *myCanvas = new TCanvas("c_invMass", "", 700, 600);
	basicCanvasFormatting(myCanvas);
	basicHistFormatting(h_invMass);
	//basicLegendFormatting();

	//Zooming out for display only.
	h_invMass->GetXaxis()->SetRangeUser(minMass, maxMass);
	h_invMass->Scale(1e-3);
    h_invMass->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
    h_invMass->GetYaxis()->SetTitle("Number of events");
	h_invMass->SetLineColor(kOrange + 7);
	h_invMass->SetLineWidth(2);
	h_invMass->SetFillStyle(3005);
	h_invMass->SetFillColorAlpha(kOrange + 7, 0.35); //less bright. more control.

	h_invMass->Draw("E1");
	h_invMass->Draw("HIST SAME");

	drawLatexText();
	drawLatexText("#it{Internal}", 0.25, 0.93, 0.042);
	drawLatexText("#times 10^{3}", 0.08, 0.93, 0.042);
	drawLatexText("#it{partial} PbPb 2024 (5.36 TeV)", 0.55, 0.93, 0.042);
	drawLatexText(Form("Z yield = %d #pm %d", TMath::Nint(z_yield), TMath::Nint(z_yield_err)), 0.18, 0.8, 0.04);
	drawLatexText(Form("m_{#mu^{+}#mu^{-}} = (%.2f #pm %.2f) GeV", mean, mean_err), 0.18, 0.72, 0.04);
	drawLatexText(Form("#sigma_{m_{#mu^{+}#mu^{-}}} = (%.2f #pm %.2f) GeV", std_dev, std_dev_err), 0.18, 0.64, 0.04);
	drawLatexText(Form("%d < m_{#mu^{+}#mu^{-}} < %d GeV", TMath::Nint(minDimuonMass), TMath::Nint(maxDimuonMass)), 0.64, 0.8, 0.04);
	//drawLatexText(Form("|#eta| < %.1f", eta_HighLimit), 0.23, 0.33, 0.04);

	TString output = TString(__func__) + plot_extension;

	myCanvas->Update();
	myCanvas->SaveAs(output);

	rootFile->Close();
}

void drawLatexText(const char *latexText, double x, double y, double TextSize){
// We can add any pT or eta selection. If no text is passed to the function the CMS Header will be drawn.
TLatex latex;
latex.SetNDC();              // Use normalized coordinates
latex.SetTextSize(TextSize);
latex.SetTextFont(42);       // Helvetica
latex.SetTextAlign(11);      // Left-top aligned.
latex.DrawLatex(x, y, latexText);
}

void basicCanvasFormatting(TCanvas* c){
    c->SetLeftMargin(0.14);
    c->SetRightMargin(0.035);
    c->SetBottomMargin(0.14);
    c->SetTopMargin(0.08);
    c->SetTickx(1);
    c->SetTicky(1);
    c->SetFillColor(0);
    c->SetFrameFillColor(0);
    c->SetFrameLineWidth(2);
}

void basicHistFormatting(TH1D* hist){
	hist->GetXaxis()->CenterTitle(true);
    hist->GetYaxis()->CenterTitle(true);
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetYaxis()->SetTitleOffset(0.9);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetLabelSize(0.042);
    hist->GetYaxis()->SetLabelSize(0.042);
    hist->SetTitle("");
	hist->SetStats(0);
	hist->SetMarkerStyle(20);
	//hist->SetMarkerSize(0.9);
}

void basicLegendFormatting(TLegend *leg){
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.042);
    leg->SetTextFont(42);
    leg->SetMargin(0.2);
    leg->SetEntrySeparation(0.04);
}

void printTreeContents(const std::string& PathToROOTFile, const char* dir, const char* tree){
	TFile *rootFile = TFile::Open(PathToROOTFile.c_str(), "READ");
    rootFile->ls();
    TDirectoryFile *dir1 = dynamic_cast<TDirectoryFile*>(rootFile->Get(dir));
	TTree *dimuonTree = dynamic_cast<TTree*>(dir1->Get(tree));
	dimuonTree->Print(); // Prints content of TTree.
}