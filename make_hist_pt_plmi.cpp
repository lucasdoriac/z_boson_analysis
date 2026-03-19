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

double minMass = 80.0;
double maxMass = 100.0;
double delta = 1e-5;
TString plot_extension = ".png";
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
    "../../../OniaTree/Oniatree_PbPb2024_PromptReco.root",
    "Oniatree_PbPb2024_PromptReco.root",
    "5TeV_data",
    "_5TeV"
};

void make_hist_pt_plmi(){

	gROOT->SetBatch(kTRUE);

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
	//dimuonTree->Show(26); // Prints content of n-th event.
	Long64_t nEvents = dimuonTree->GetEntries(); 
	
	TH2D *hist_pt_plmi = new TH2D("hist_pt_plmi", "p_{T}(#mu^{+}) vs p_{T}(#mu^{-})", 100, 0, 200, 100, 0, 200);
	//hist_pt_plmi->SetDirectory(0);

	Short_t Reco_QQ_size;
	dimuonTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);

	Short_t Reco_mu_size;
	dimuonTree->SetBranchAddress("Reco_mu_size", &Reco_mu_size);

	const int MAX_MU = 100;
	const int MAX_QQ = 10;

	Short_t Reco_QQ_mupl_idx[MAX_QQ];
	Short_t Reco_QQ_mumi_idx[MAX_QQ];

	dimuonTree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
	dimuonTree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);

	Short_t Reco_mu_charge[MAX_MU];
	dimuonTree->SetBranchAddress("Reco_mu_charge", Reco_mu_charge);

	std::vector<float>* Reco_mu_4mom_pt  = nullptr;
	std::vector<float>* Reco_mu_4mom_m = nullptr;	
	std::vector<float>* Reco_QQ_4mom_pt  = nullptr;
	std::vector<float>* Reco_QQ_4mom_m = nullptr;
	
	dimuonTree->SetBranchAddress("Reco_mu_4mom_pt", &Reco_mu_4mom_pt);
	dimuonTree->SetBranchAddress("Reco_mu_4mom_m", &Reco_mu_4mom_m);
	dimuonTree->SetBranchAddress("Reco_QQ_4mom_pt", &Reco_QQ_4mom_pt);
	dimuonTree->SetBranchAddress("Reco_QQ_4mom_m", &Reco_QQ_4mom_m);

	for(Long64_t i = 0; i < nEvents; ++i){// Loop over events...
		
		dimuonTree->GetEntry(i);
		
		if(Reco_QQ_size > 0){//If event has dimuon candidate...
			
			for(Short_t j = 0; j < Reco_QQ_size; ++j){

				// The corresponding dimuon pair is indexed with idx variables.
				double ptPlus = Reco_mu_4mom_pt->at(Reco_QQ_mupl_idx[j]);
				double ptMinus = Reco_mu_4mom_pt->at(Reco_QQ_mumi_idx[j]);
				hist_pt_plmi->Fill(ptPlus, ptMinus);
			}
		}
	}

    TCanvas *c = new TCanvas("c", "p_{T}(#mu^{+}) vs p_{T}(#mu^{-}", 900, 700);
	c->SetLeftMargin(0.11);
    c->SetRightMargin(0.1);
    c->SetBottomMargin(0.11);
    c->SetTopMargin(0.08);
    c->SetTickx(1);
    c->SetTicky(1);
    c->SetFillColor(0);
    c->SetFrameFillColor(0);
    c->SetFrameLineWidth(2);

	hist_pt_plmi->GetXaxis()->CenterTitle(true);
	hist_pt_plmi->GetYaxis()->CenterTitle(true);
	hist_pt_plmi->GetZaxis()->CenterTitle(true);

	hist_pt_plmi->GetXaxis()->SetTitleOffset(1.);
	hist_pt_plmi->GetYaxis()->SetTitleOffset(1.);
	hist_pt_plmi->GetZaxis()->SetTitleOffset(1.);

	hist_pt_plmi->GetXaxis()->SetTitleFont(42);
	hist_pt_plmi->GetYaxis()->SetTitleFont(42);
	hist_pt_plmi->GetZaxis()->SetTitleFont(42);

	hist_pt_plmi->GetXaxis()->SetLabelFont(42);
	hist_pt_plmi->GetYaxis()->SetLabelFont(42);
	hist_pt_plmi->GetZaxis()->SetLabelFont(42);

	hist_pt_plmi->GetXaxis()->SetTitleSize(0.05);
	hist_pt_plmi->GetYaxis()->SetTitleSize(0.05);
	hist_pt_plmi->GetZaxis()->SetTitleSize(0.05);

	hist_pt_plmi->GetXaxis()->SetLabelSize(0.038);
	hist_pt_plmi->GetYaxis()->SetLabelSize(0.038);
	hist_pt_plmi->GetZaxis()->SetLabelSize(0.035);

	hist_pt_plmi->GetXaxis()->SetTitle("p_{T}(#mu^{+}) [GeV]");
	hist_pt_plmi->GetYaxis()->SetTitle("p_{T}(#mu^{-}) [GeV]");
	hist_pt_plmi->GetZaxis()->SetTitle("Entries");
	hist_pt_plmi->GetXaxis()->SetRangeUser(10., 200.);
	hist_pt_plmi->GetYaxis()->SetRangeUser(10., 200.);
	hist_pt_plmi->SetStats(0);

	hist_pt_plmi->Draw("COLZ");
	c->Update();
	c->SaveAs("output.png");

	rootFile->Close();
}