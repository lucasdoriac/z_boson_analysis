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

void ptPL_vs_ptMI(){

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
	
	TH2D *hist_pt_plmi = new TH2D("hist_pt_plmi", "p_{T}(#mu^{+}) vs p_{T}(#mu^{-})", 100, 0, 100, 100, 0, 100);
	//hist_pt_plmi->SetDirectory(0);

	Short_t Reco_QQ_size;
	dimuonTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);

	const int MAX_QQ = 10;
	Short_t Reco_QQ_mupl_idx[MAX_QQ];
	Short_t Reco_QQ_mumi_idx[MAX_QQ];
	dimuonTree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
	dimuonTree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);

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
			
			for(Short_t j = 0; j < Reco_QQ_size; ++j){//Maybe event has more than one candidate...

				if(Reco_QQ_4mom_m->at(j) < 95.0 && Reco_QQ_4mom_m->at(j) > 85.0){//If j-th candidate has mass in [85,95] GeV...
					
					// The corresponding dimuon pair is indexed with idx variables.
					double ptPlus = Reco_mu_4mom_pt->at(Reco_QQ_mupl_idx[j]);
					double ptMinus = Reco_mu_4mom_pt->at(Reco_QQ_mumi_idx[j]);
					hist_pt_plmi->Fill(ptPlus, ptMinus);
				}
			}
		}
	}
    
	//Now we create two TH1 histograms for each pT spectrum.
	TH1D *hist_ptPl = (TH1D*)hist_pt_plmi->ProjectionX();
	TH1D *hist_ptMi = (TH1D*)hist_pt_plmi->ProjectionY();

    TCanvas *c = new TCanvas("c", "p_{T}(#mu^{+}) vs p_{T}(#mu^{-}", 900, 700);
	c->SetLeftMargin(0.11);
    c->SetRightMargin(0.04);
    c->SetBottomMargin(0.11);
    c->SetTopMargin(0.08);
    c->SetTickx(1);
    c->SetTicky(1);
    c->SetFillColor(0);
    c->SetFrameFillColor(0);
    c->SetFrameLineWidth(2);

    hist_ptPl->GetXaxis()->CenterTitle(true);
	hist_ptPl->GetYaxis()->CenterTitle(true);

	hist_ptPl->GetXaxis()->SetTitleOffset(1.);
	hist_ptPl->GetYaxis()->SetTitleOffset(1.);

	hist_ptPl->GetXaxis()->SetTitleFont(42);
	hist_ptPl->GetYaxis()->SetTitleFont(42);

	hist_ptPl->GetXaxis()->SetLabelFont(42);
	hist_ptPl->GetYaxis()->SetLabelFont(42);

	hist_ptPl->GetXaxis()->SetTitleSize(0.05);
	hist_ptPl->GetYaxis()->SetTitleSize(0.05);

	hist_ptPl->GetXaxis()->SetLabelSize(0.038);
	hist_ptPl->GetYaxis()->SetLabelSize(0.038);

	hist_ptMi->GetXaxis()->CenterTitle(true);
	hist_ptMi->GetYaxis()->CenterTitle(true);

	hist_ptMi->GetXaxis()->SetTitleOffset(1.);
	hist_ptMi->GetYaxis()->SetTitleOffset(1.);

	hist_ptMi->GetXaxis()->SetTitleFont(42);
	hist_ptMi->GetYaxis()->SetTitleFont(42);

	hist_ptMi->GetXaxis()->SetLabelFont(42);
	hist_ptMi->GetYaxis()->SetLabelFont(42);

	hist_ptMi->GetXaxis()->SetTitleSize(0.05);
	hist_ptMi->GetYaxis()->SetTitleSize(0.05);

	hist_ptMi->GetXaxis()->SetLabelSize(0.038);
	hist_ptMi->GetYaxis()->SetLabelSize(0.038);

	hist_ptPl->SetStats(0);
	hist_ptMi->SetStats(0);

	hist_ptPl->SetLineColor(kRed);
	hist_ptMi->SetLineColor(kBlue);

	hist_ptPl->SetLineWidth(2);
	hist_ptMi->SetLineWidth(2);

	hist_ptPl->GetXaxis()->SetTitle("p_{T} [GeV]");
	hist_ptPl->GetYaxis()->SetTitle("Reco Z^{0} candidates");
	//hist_ptPl->GetXaxis()->SetRangeUser(20., 60.);
	//hist_ptPl->GetYaxis()->SetRangeUser(100., 850.);

	int binMaxPl = hist_ptPl->GetMaximumBin();
	int binMaxMi = hist_ptMi->GetMaximumBin();
	double xMaxPl = hist_ptPl->GetBinCenter(binMaxPl);
	double xMaxMi = hist_ptMi->GetBinCenter(binMaxMi);
	double ymaxPl = hist_ptPl->GetMaximum();
	double ymaxMi = hist_ptMi->GetMaximum();
	double errPl  = hist_ptPl->GetBinError(binMaxPl);
	double errMi  = hist_ptMi->GetBinError(binMaxMi);

	TLine *linePl = new TLine(xMaxPl, 0.0, xMaxPl, ymaxPl);
	TLine *lineMi = new TLine(xMaxMi, 0.0, xMaxMi, ymaxMi);

	linePl->SetLineColor(kRed);
	linePl->SetLineWidth(2);
	linePl->SetLineStyle(2);

	lineMi->SetLineColor(kBlue);
	lineMi->SetLineWidth(2);
	lineMi->SetLineStyle(2);

	hist_ptPl->Draw("HIST E");
	linePl->Draw("SAME");
	hist_ptMi->Draw("HIST E SAME");
	lineMi->Draw("SAME");

	TLegend *leg1 = new TLegend(0.6, 0.7, 0.85, 0.85);
	leg1->SetBorderSize(0);
	leg1->SetFillStyle(0);
	leg1->SetTextFont(42);
	leg1->SetTextSize(0.04);
	leg1->AddEntry(hist_ptPl, Form("#mu^{+} peak = %.1f #pm %.1f", xMaxPl, errPl), "l");
	leg1->AddEntry(hist_ptMi, Form("#mu^{-} peak = %.1f #pm %.1f", xMaxMi, errMi), "l");
	leg1->Draw();

	TLatex latex;
	latex.SetNDC();              // Use normalized coordinates
	latex.SetTextSize(0.04);
	latex.SetTextFont(42);       // Helvetica
	latex.SetTextAlign(11);      // Left-top aligned.	
	latex.DrawLatex(0.12, 0.93, "#bf{CMS} #it{Internal}");

	TLatex latex2;
	latex2.SetNDC();              // Use normalized coordinates
	latex2.SetTextSize(0.04);
	latex2.SetTextFont(42);       // Helvetica
	latex2.SetTextAlign(11);      // Left-top aligned.	
	latex2.DrawLatex(0.62, 0.55, "85 < m_{#mu^{+}#mu^{-}} < 95 GeV");

	TLatex latex3;
	latex3.SetNDC();              // Use normalized coordinates
	latex3.SetTextSize(0.04);
	latex3.SetTextFont(42);       // Helvetica
	latex3.SetTextAlign(11);      // Left-top aligned.	
	latex3.DrawLatex(0.7, 0.93, "Partial PbPb 2024");

	c->Update();
	c->SaveAs("output3.png");

	rootFile->Close();
}