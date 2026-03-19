/*
Section for comments.

17-03-2026
- Initial goal is to write a macro that plots mass, pT, eta, phi and centrality distributions for PbPb 2024 CMS data sample.
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
double eta_LowLimit = -1.5;
double eta_HighLimit = 1.5;
double pt_LowLimit = 20.0; // GeV.
double pt_HighLimit = 1000.0; // GeV. Later need to check if this limit is reasonable.
TString plot_extension = ".png";
TString output_name = "out";
TString base_output_path = "./";

// ##############################################################################
// ##############################################################################


// --- Structs ---
struct dataFile {
    std::string path;
    std::string name;
    std::string label;
    std::string sufix;
};

const dataFile PbPb_5TeV_2024 = {
    "../../../OniaTree/Oniatree_PbPb2024_PromptReco.root",
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

void invMassSpectrum(const std::string& PathToROOTFile = "Oniatree_PbPb2024_PromptReco.root");
void drawLatexText(const char *latexText = "#bf{CMS} #it{Internal}", double x = 0.12, double y = 0.93);
void basicHistFormatting(TH1D *hist);

void controlPlots(){

	invMassSpectrum(PbPb_5TeV_2024.path);
}

void invMassSpectrum(const std::string& PathToROOTFile){

	TFile *rootFile = TFile::Open(PathToROOTFile.c_str(), "READ");
	if (!rootFile || rootFile->IsZombie()) {
        std::cerr << "Error in " << __func__ << ": Could not open file." << std::endl;
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

	dimuonTree->Print(); // Prints content of TTree.
	dimuonTree->Show(0); // Prints content of n-th event.
	Long64_t nEvents = dimuonTree->GetEntries(); 
	//std::cout << "Total number of events on the TTree -> " << nEvents << std::endl;// Prints the number of total events on the TTree.

	//TH1D *h_invMass = new TH1D("h_invMass","Dimuon invariant mass", 20, 80, 100);
	//h_invMass->SetDirectory(0);
	
	std::vector<float>* Reco_QQ_4mom_pt  = nullptr;
	std::vector<float>* Reco_QQ_4mom_eta = nullptr;
	std::vector<float>* Reco_QQ_4mom_phi = nullptr;
	std::vector<float>* Reco_QQ_4mom_m = nullptr;
	dimuonTree->SetBranchAddress("Reco_QQ_4mom_pt",  &Reco_QQ_4mom_pt);
	dimuonTree->SetBranchAddress("Reco_QQ_4mom_eta", &Reco_QQ_4mom_eta);
	dimuonTree->SetBranchAddress("Reco_QQ_4mom_phi", &Reco_QQ_4mom_phi);
	dimuonTree->SetBranchAddress("Reco_QQ_4mom_m", &Reco_QQ_4mom_m);

	//dimuonTree->GetEntry(0);
	//std::cout << "Print single-event Reco_QQ_4momm_m content-> " << Reco_QQ_4mom_m->at(0) << std::endl;// Prints the number of total events on the TTree.

	for(Long64_t i = 0; i < nEvents; ++i){

		dimuonTree->GetEntry(i);
		if (Reco_QQ_4mom_m && Reco_QQ_4mom_m->size() > 0){ //If the event has a dimuon candidate, runs the loop, if not, go to next event.
			
			for (size_t j = 0; j < Reco_QQ_4mom_m->size(); ++j) {
    			std::cout << "Event " << i << " Reco_QQ_4mom_m content-> " << Reco_QQ_4mom_m->at(j) << std::endl;
    			std::cout << "Event " << i << " Reco_QQ_4mom_pt content-> ";
    		for (auto pt : *Reco_QQ_4mom_pt) {
        		std::cout << pt << " ";
    		}
    		std::cout << std::endl;
    		std::cout << "Event " << i << " Reco_QQ_4mom_pt size-> " << Reco_QQ_4mom_pt->size() << std::endl;
			}

		}
	}

/*double get_n_events(const std::string& filename, double EHFmin, double EHFmax){
    std::string this_function = __func__;
    TFile *f = TFile::Open(filename.c_str(), "READ");
    if (!f || f->IsZombie()) {
    	macro_log("From " + this_function + ": error while opening the file " + filename);
    	return -1;
    }

    TDirectory *dir = (TDirectory*)f->Get("QA_histograms");
    if (!dir) {
    	macro_log("From " + this_function + ": error while opening 'QA_histograms' directory.");
        f->Close();
        return -1;
    }

    TH1D *hf_hist = (TH1D*)dir->Get("hfSumEtPb");
    if (!hf_hist) {
    	macro_log("From " + this_function + ": histogram 'hfSumEtPb' not found!");
        f->Close();
        return -1;
    }
    hf_hist->SetDirectory(0);
    f->Close();

    int bin_min = hf_hist->GetXaxis()->FindBin(EHFmin + delta);
    int bin_max = hf_hist->GetXaxis()->FindBin(EHFmax - delta);

    // If one wants to check interval of integration over transversal E_HF energy.
    double lowEdge   = hf_hist->GetXaxis()->GetBinLowEdge(bin_min);
    double highEdge  = hf_hist->GetXaxis()->GetBinUpEdge(bin_min);
    double lowEdge_   = hf_hist->GetXaxis()->GetBinLowEdge(bin_max);
    double highEdge_  = hf_hist->GetXaxis()->GetBinUpEdge(bin_max);
    printf("\n-> From get_n_events: Integrating 'hfSumEtPb' from bin [%.1f,%.1f] GeV to bin [%.1f,%.1f] GeV \n\n\n", lowEdge, highEdge, lowEdge_, highEdge_);
    //

    // Calculates the number of events n_events detected by the Pb side HF, for the specified centrality class.
    double n_events = hf_hist->Integral(bin_min, bin_max);

    return n_events;
}*/
}

void drawLatexText(const char *latexText, double x, double y){

// Add any pT or eta selection. If no text is passed to the function the CMS Header will be drawn.
TLatex latex;
latex.SetNDC();              // Use normalized coordinates
latex.SetTextSize(0.04);
latex.SetTextFont(42);       // Helvetica
latex.SetTextAlign(11);      // Left-top aligned.
latex.DrawLatex(x, y, latexText);
}