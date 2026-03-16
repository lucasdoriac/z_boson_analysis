#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <iostream>
#include <string>
#include <vector>


TString plot_extension = ".pdf";
TString output_name = "out";
TString base_output_path = "./";

// --- Structs ---
struct dataFile {
    std::string path;
    std::string name;
    std::string label;
    std::string sufix;
};

// --- struct to 5TeV dataset
const dataFile rootFile_5TeV = {
    "./Oniatree_PbPb2024_PromptReco.root",
    "Oniatree_PbPb2024_PromptReco.root",
    "5TeV_data",
    "_5TeV"
};

void invMassSpectrum(){

	TFile *rootFile = TFile::Open(rootFile_5TeV.path.c_str(), "READ");
	if (!rootFile || rootFile->IsZombie()) {
        std::cerr << "Error in " << __func__ << ": Could not open "
          << rootFile_5TeV.name << " file." << std::endl;
        return;
    }
    rootFile->ls();

    TDirectoryFile *dir1 = dynamic_cast<TDirectoryFile*>(rootFile->Get("hionia"));
    if(!dir1){
		std::cerr << "Error in " << __func__ << ": Could not open dir1." 
		  << std::endl;
		return;
    }

    // Selecting tree.
	TTree *dimuonTree = dynamic_cast<TTree*>(dir1->Get("myTree"));
	if(!dimuonTree){
		std::cerr << "Error in " << __func__ << ": Could not load Tree." 
		  << std::endl;
		return;
	}
	//dimuonTree->Print(); // Prints content of TTree.
	//dimuonTree->Show(0); // Prints content of n-th event.
	Long64_t nEvents = dimuonTree->GetEntries(); 
	//std::cout << "Total number of events on the TTree -> " << nEvents << std::endl;// Prints the number of total events on the TTree.

	TH1D *h_invMass = new TH1D("h_invMass","Dimuon invariant mass", 20, 80, 100);
	h_invMass->SetDirectory(0);
	
	std::vector<float>* Reco_QQ_4mom_pt  = nullptr;
	std::vector<float>* Reco_QQ_4mom_eta = nullptr;
	std::vector<float>* Reco_QQ_4mom_phi = nullptr;
	std::vector<float>* Reco_QQ_4mom_m = nullptr;
	dimuonTree->SetBranchAddress("Reco_QQ_4mom_pt",  &Reco_QQ_4mom_pt);
	dimuonTree->SetBranchAddress("Reco_QQ_4mom_eta", &Reco_QQ_4mom_eta);
	dimuonTree->SetBranchAddress("Reco_QQ_4mom_phi", &Reco_QQ_4mom_phi);
	dimuonTree->SetBranchAddress("Reco_QQ_4mom_m", &Reco_QQ_4mom_m);

	for(Long64_t i = 0; i < nEvents; ++i){

		dimuonTree->GetEntry(i);
		if (Reco_QQ_4mom_m){		
			
			for (size_t j = 0; j < Reco_QQ_4mom_m->size(); ++j){ //If the event has a dimuon candidate, runs the loop, if not, go to next event.
				double mass = Reco_QQ_4mom_m->at(j);
				h_invMass->Fill(mass);
			}
		}
	
	}

	/*
	int bin_min = FindBin(80);
	int bin_max = FindBin(110);
	double z_yield_err = 0;
	double z_yield = IntegralAndError(bin_min, bin_max, z_yield_err);
	// Calculate mean of mass distribution.

	// Calculate standard deviation.
	*/

	double delta = 1e-5;
	int bin_min = h_invMass->FindBin(80 + delta);
	int bin_max = h_invMass->FindBin(100 - delta);

	// Calculate Z yield.
	double z_yield_err = 0;
	double z_yield = h_invMass->IntegralAndError(bin_min, bin_max, z_yield_err);

	TCanvas *c = new TCanvas("c","Invariant Mass",700,600);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.035);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    c->SetTickx(1);
    c->SetTicky(1);
    c->SetFillColor(0);
    c->SetFrameFillColor(0);
    c->SetFrameLineWidth(1);

    h_invMass->GetXaxis()->CenterTitle(true);
    h_invMass->GetYaxis()->CenterTitle(true);
    h_invMass->GetXaxis()->SetTitleOffset(1.);
    h_invMass->GetYaxis()->SetTitleOffset(1.);
    h_invMass->GetXaxis()->SetTitleFont(42);
    h_invMass->GetYaxis()->SetTitleFont(42);
    h_invMass->GetXaxis()->SetLabelFont(42);
    h_invMass->GetYaxis()->SetLabelFont(42);
    h_invMass->GetXaxis()->SetTitleSize(0.05);
    h_invMass->GetYaxis()->SetTitleSize(0.044);
    h_invMass->GetXaxis()->SetLabelSize(0.036);
    h_invMass->GetYaxis()->SetLabelSize(0.036);

    h_invMass->SetTitle("");
    h_invMass->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    h_invMass->GetYaxis()->SetTitle("Number of Events");

	h_invMass->SetMarkerStyle(20);
	h_invMass->SetMarkerSize(0.9);
	h_invMass->SetStats(0);
	h_invMass->Draw("E");

	/*TLegend *leg = new TLegend(0.60,0.68,0.81,0.85);
	leg->AddEntry(h_invMass,"dimuon mass","lep");
	leg->AddEntry()	
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.045);
    leg->SetTextFont(42);
    leg->SetMargin(0.2);
    leg->SetEntrySeparation(0.04);
    leg->Draw();*/

	TLatex latex;
	latex.SetNDC(); // For normalized coordinates.
	latex.SetTextSize(0.04);
	latex.SetTextFont(42);
	latex.SetTextAlign(11); // left-top aligned.
	TString cmsText = "#bf{CMS} #it{Internal}";
	latex.DrawLatex(0.12, 0.93, cmsText);

	TLatex latex2;
	latex2.SetNDC(); // For normalized coordinates.
	latex2.SetTextSize(0.04);
	latex2.SetTextFont(42);
	latex2.SetTextAlign(31); // right-top aligned.
	TString lumiText = "partial PbPb sample (5.36 TeV)";
	latex2.DrawLatex(0.96, 0.93, lumiText);

	TLatex latex3;
	latex.SetNDC(); // For normalized coordinates.
	latex.SetTextSize(0.04);
	latex.SetTextFont(42);
	latex.SetTextAlign(11); // left-top aligned.

	latex.DrawLatex(0.2, 0.8, Form("Z yield = %d #pm %d", TMath::Nint(z_yield), TMath::Nint(z_yield_err)));

	c->Update();
    c->SaveAs("./output.pdf");

	rootFile->Close();
}