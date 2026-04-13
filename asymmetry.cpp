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
double minDimuonMass = 85.0; // GeV.
double maxDimuonMass = 95.0; // GeV.
double maxCentrality = 180; // cen_bin% = max_centrality/2. e.g. max_centrality = 180 -> 0-90% centrality events.
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
//void TH1D_HIST_muPLmuMI(const std::string& PathToROOTFile = "Oniatree_PbPb2024_PromptReco.root");
void printTreeContents(const std::string& PathToROOTFile, const char* dir = "hionia", const char* tree = "myTree");
void invMassSpectrum(const std::string& PathToROOTFile = "Oniatree_PbPb2024_PromptReco.root", int nBins = 40, int minMass = 70, int maxMass = 110);
void drawLatexText(const char *latexText = "#bf{CMS}", double x = 0.15, double y = 0.93, double TextSize = 0.05);
void basicCanvasFormatting(TCanvas *c);
void basicHistFormatting(TH1D *hist);
void basicLegendFormatting(TLegend *leg);
// --- still writing
void pt_asymmetry(const std::string& PathToROOTFile = "Oniatree_PbPb2024_PromptReco.root");
//void pt_ratioPAD(const std::string& PathToROOTFile = "Oniatree_PbPb2024_PromptReco.root");
//void pt_sideBand(const std::string& PathToROOTFile = "Oniatree_PbPb2024_PromptReco.root");

void asymmetry(){

	gROOT->SetBatch(kTRUE);

	//printTreeContents(PbPb_5TeV_2024.path);
	//invMassSpectrum(PbPb_5TeV_2024.path);
	//TH1D_HIST_muPLmuMI(PbPb_5TeV_2024.path);
	pt_asymmetry(PbPb_5TeV_2024.path);
	//pt_ratioPAD(PbPb_5TeV_2024.path);
	//pt_sideBand(PbPb_5TeV_2024.path);
}


void pt_asymmetry(const std::string& PathToROOTFile){
	// Load root file.
	TFile *rootFile = TFile::Open(PathToROOTFile.c_str(), "READ");
    rootFile->ls();

    // Load directory.
    TDirectoryFile *dir1 = dynamic_cast<TDirectoryFile*>(rootFile->Get("hionia"));

    // Load tree.
	TTree *dimuonTree = dynamic_cast<TTree*>(dir1->Get("myTree"));

	//Total number of events on Tree.
	Long64_t nEvents = dimuonTree->GetEntries();
	//dimuonTree->Print();

	//Event-level quantities
	float zVtx;
	dimuonTree->SetBranchAddress("zVtx", &zVtx);
	int Centrality;
	dimuonTree->SetBranchAddress("Centrality", &Centrality);

	//Reco muon quantities
	Short_t Reco_mu_size;
	dimuonTree->SetBranchAddress("Reco_mu_size", &Reco_mu_size);

	const int MAX_MU = 1000;
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

	//Reco Z candidates quantities
	Short_t Reco_QQ_size;
	dimuonTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);

	const int MAX_QQ = 100;
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

	TH1D *h_asym = new TH1D("h_asym", "(p_{T}(#mu^{+}) - p_{T}(#mu^{-}))/(p_{T}(#mu^{+}) + p_{T}(#mu^{-}))", 100, -1, 1);
	double ptplus, ptminus, etaplus, etaminus, pt_asym;

	for(Long64_t i = 0; i < nEvents; ++i){//Loop events. "event-by-event"...

		dimuonTree->GetEntry(i);//Get BRANCH values at the i-th event.
		bool goodVertex = (std::abs(zVtx) < 15.0);
		bool goodCent = (Centrality >= 60 && Centrality <= 140); //30 to 70%.

		if(goodVertex && goodCent){//Good event selection.

			for (Short_t j = 0; j < Reco_QQ_size; ++j){//If event has dimuon candidate, flagged by Reco_QQ_size > 0...

				if(Reco_QQ_4mom_m->at(j) >= 80.0 && Reco_QQ_4mom_m->at(j) <= 100.0){//Z mass window
					
						ptplus = Reco_mu_4mom_pt->at(Reco_QQ_mupl_idx[j]); //pT of antimuon...
						ptminus = Reco_mu_4mom_pt->at(Reco_QQ_mumi_idx[j]); //pT of corresponding muon...
						etaplus = Reco_mu_4mom_eta->at(Reco_QQ_mupl_idx[j]);
						etaminus = Reco_mu_4mom_eta->at(Reco_QQ_mumi_idx[j]);

						//Good muon selection
						bool goodMuonPlus = (ptplus >= 20.0) && (std::abs(etaplus) < 2.4);
						bool goodMuonMinus = (ptminus >= 20.0) && (std::abs(etaminus) < 2.4);

						if(goodMuonPlus && goodMuonMinus){
							pt_asym = (ptplus - ptminus)/(ptplus + ptminus);
							h_asym->Fill(pt_asym);
						}
				}
			}
		}

	}//Exiting event-by-event loop.

	//Calculate mean of asymmetry distribution (should be around ~ 0 GeV...)
	double mean = h_asym->GetMean();
	double mean_err = h_asym->GetMeanError();

	TCanvas *myCanvas = new TCanvas("c_ptplmi", "", 700, 600);
	basicCanvasFormatting(myCanvas);
	basicHistFormatting(h_asym);

    h_asym->GetXaxis()->SetTitle("p_{T}(#mu^{+}) - p_{T}(#mu^{+})/p_{T}(#mu^{+}) + p_{T}(#mu^{-})");
    h_asym->GetYaxis()->SetTitle("Dimuon candidates [GeV^{-1}]");
	h_asym->SetLineColor(kAzure);
	h_asym->SetLineWidth(1);
	h_asym->SetFillStyle(3001);
	h_asym->SetFillColorAlpha(kAzure, 0.35); //less bright. more control.
	h_asym->SetStats(1);
    
    h_asym->GetXaxis()->SetTitleSize(0.04);
    h_asym->GetXaxis()->SetTitleOffset(1.1);
	h_asym->Draw("HIST");

	// Draw vertical line at x = 0
	TLine *line = new TLine(0, h_asym->GetMinimum(), 0, h_asym->GetMaximum());
	line->SetLineColor(kRed);
	line->SetLineStyle(2); // dashed
	line->SetLineWidth(1);
	line->Draw("SAME");

	drawLatexText();
	drawLatexText("#it{Internal}", 0.25, 0.93, 0.042);
	drawLatexText("PbPb 2024 (5.36 TeV)", 0.65, 0.93, 0.042);

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

	h_invMass->GetXaxis()->SetRangeUser(minMass, maxMass);
	h_invMass->Scale(1e-3);
    h_invMass->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} [GeV]");
    h_invMass->GetYaxis()->SetTitle("Dimuon candidates [GeV^{-1}]");
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
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.035);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    c->SetTickx(1);
    c->SetTicky(1);
    c->SetFillColor(0);
    c->SetFrameFillColor(0);
    c->SetFrameLineWidth(2);
}

void basicHistFormatting(TH1D* hist){
	hist->GetXaxis()->CenterTitle(false);
    hist->GetYaxis()->CenterTitle(false);
    hist->GetXaxis()->SetTitleOffset(.9);
    hist->GetYaxis()->SetTitleOffset(1.);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetTitleSize(0.055);
    hist->GetYaxis()->SetTitleSize(0.055);
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