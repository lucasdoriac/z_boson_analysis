/*
Code to analyze MC datasets available at cern box.

Step 1:
1. Create a TChain, add the files, and print the total number of entries.
2. Read a single branch, e.g. Reco_QQ_size and print its contents.
3. Make a simple histogram.
4. Read and plot a vector branch.

Step 2:
1. Inspect, through simple histograms, the reconstructed objects.
	1. Make histograms of simple distributions such as Reco_QQ_size, dimuon pT, dimuon pseudorapidity, dimuon phi angle and dimuon charge.
2. Inspect, through simple histograms, the single muon distributions.
	1. Make histograms of simple distributions such as Reco_mu_size, muon pT, muon pseudorapidity, muon phi angle and muon charge.
3. Reconstruct the invariant mass spectrum out of muon variables and compare with result from OniaTree (Reco_QQ_4mom_m).
4. Same as Step 3, but now applying muon-level and pair-level cuts in order to improve quality.
	1. pT > 20 GeV
	2. fabs(eta) < 2.4
	3. quality criteria
	4. trigger matching
5. NOT DONE HERE!!! Apply specific pair-level cuts in order to optmize the final observable.
	1. opposite charge
	2. 80 < M_pair < 100 GeV
	3. rapidity selection
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
#include <TChain.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>

// ### User Settings ###
double eta_LowLimit = -2.4;
double eta_HighLimit = 2.4;
double pt_LowLimit = 20.0; // GeV.
double pt_HighLimit = 100.0; // GeV. Later need to check if this limit is reasonable.
double minDimuonMass = 80.0; // GeV.
double maxDimuonMass = 100.0; // GeV.
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

const dataFile MC_PbPb2024 = {
    "/home/lucasdoriac/z_boson_analysis/MC/PbPb2024/DYto2Mu_MLL-50_TuneCP5_5p36TeV_powheg-pythia8/PowhegEmbedded_March9/260309_143939/0000/",
    "",
    "MC_PbPb_2024",
};

//---simple histogram functions for dimuon candidates
void test_chain(); //Histogram of number of reco dimuon candidates.
void dimuon_pt(); //Histogram of reco dimuon candidates pT distribution.
void dimuon_pseudorap(); //Histogram of reco dimuon candidates eta distribution.
void dimuon_phi(); //Histogram of reco dimuon candidates phi angle distribution.
void dimuon_mass(); //Histogram of reco dimuon candidates mass distribution.
void dimuon_rapidity();

//---simple histogram functions for single reco muons.
void reco_mu_size();
void muon_pt();
void muon_pseudorap();
void muon_phi();
void muon_mass();
void muon_charge();

//---function to reconstruct (calculate) inv mass spectrum from mu variables and compare with Reco_QQ_4mom_m.
void mu_inv_mass();
void mu_inv_mass_WithSelection();
void check_triggerBranch();

//---main
void monte_carlo_analysis(){

	//Part 4:
	//mu_inv_mass_WithSelection();
	//check_triggerBranch();
	//muon_charge();

	//Part 3:
	//mu_inv_mass();

	// Part 2:
	//muon_phi();
	//muon_mass();
	//muon_pseudorap();
	//reco_mu_size();
	//muon_pt();
	
	// Part 1:
	//dimuon_mass();
	//dimuon_phi();
	//dimuon_pseudorap();
	//dimuon_pt();
	dimuon_rapidity();
	//test_chain();

}

void muon_charge(){

	TChain *chain = new TChain("hionia/myTree");
	chain->Add((MC_PbPb2024.path + "Oniatree_PowhegZtoMuMu_PbPb2024_*.root").c_str());

	int nFiles = chain->GetListOfFiles()->GetEntries();
	std::cout << "Number of files = " << nFiles << std::endl;

		if (nFiles != 20) //Just a check to make sure all 20 files were loaded.
	{
	   	std::cout << "ERROR: Wrong number of files were added to the chain!" << std::endl;
    	return;
	}

	Long64_t nEntries = chain->GetEntries(); //Get number of EVENTS.

	const int MAX_MU = 1000;
    Short_t Reco_mu_charge[MAX_MU]; //Charge of reco muons.
    Short_t Reco_mu_size; //Number of reconstructed muons in the event.

	chain->SetBranchAddress("Reco_mu_size", &Reco_mu_size);
    chain->SetBranchAddress("Reco_mu_charge", Reco_mu_charge);

 	TH1D *h = new TH1D("h", "Reco muon charge distribution;charge (e);n of reco muons", 10, -5., +5.);

    for(Long64_t i = 0; i < nEntries; ++i){

    	chain->GetEntry(i);

    	for(Short_t j = 0; j < Reco_mu_size; ++j){

    		if (Reco_mu_charge[j] != -1 && Reco_mu_charge[j] != 1) {
			    std::cout << "Unexpected charge = "
			              << Reco_mu_charge[j]
			              << " in event " << i
			              << std::endl;
			}

    		h->Fill(Reco_mu_charge[j]);
    	}

		if (i % (nEntries / 100) == 0){ //Just to check the progress.

			        int percent = static_cast<int>(100.0 * i / nEntries+0.5);

		    	    std::cout << "\r"
		        	          << percent
		            	      << "% complete..."
		                	  << std::flush;
		}

    }

    std::cout << "Negative muons: "
          << h->GetBinContent(h->FindBin(-1))
          << std::endl;

	std::cout << "Positive muons: "
          << h->GetBinContent(h->FindBin(+1))
          << std::endl;

    TCanvas *c = new TCanvas();
	h->Draw();
	c->SaveAs("../../../mnt/c/Users/lucas/OneDrive/Documentos/charge of reco single muons.png");

    delete chain;
    delete c;
    delete h;


}

void mu_inv_mass(){

	TChain *chain = new TChain("hionia/myTree");
	chain->Add((MC_PbPb2024.path + "Oniatree_PowhegZtoMuMu_PbPb2024_*.root").c_str());

	int nFiles = chain->GetListOfFiles()->GetEntries();
	std::cout << "Number of files = " << nFiles << std::endl;

	if (nFiles != 20) //Just a check to make sure all 20 files were loaded.
	{
	   	std::cout << "ERROR: Wrong number of files were added to the chain!" << std::endl;
    	return;
	}

	Long64_t nEntries = chain->GetEntries(); //Get number of EVENTS.

	std::vector<float> *Reco_mu_4mom_m = nullptr; //Mass of reco muons.
	std::vector<float> *Reco_mu_4mom_pt = nullptr; //pT of reco muons.
	std::vector<float> *Reco_mu_4mom_eta = nullptr; //eta of reco muons.
	std::vector<float> *Reco_mu_4mom_phi = nullptr; //phi of reco muons.

 	chain->SetBranchAddress("Reco_mu_4mom_m", &Reco_mu_4mom_m);
 	chain->SetBranchAddress("Reco_mu_4mom_pt", &Reco_mu_4mom_pt);
 	chain->SetBranchAddress("Reco_mu_4mom_eta", &Reco_mu_4mom_eta);
 	chain->SetBranchAddress("Reco_mu_4mom_phi", &Reco_mu_4mom_phi);

	const int MAX_MU = 1000;
	Short_t Reco_mu_charge[MAX_MU];
	chain->SetBranchAddress("Reco_mu_charge", Reco_mu_charge);
 	
 	std::vector<float> *Reco_QQ_4mom_m = nullptr; //Mass of reco Z boson candidates.
	std::vector<float> *Reco_QQ_4mom_pt = nullptr; //pT of reco Z boson candidates.
	std::vector<float> *Reco_QQ_4mom_eta = nullptr; //eta of reco Z boson candidates.
	std::vector<float> *Reco_QQ_4mom_phi = nullptr; //phi of reco Z boson candidates.

 	chain->SetBranchAddress("Reco_QQ_4mom_m", &Reco_QQ_4mom_m);
 	chain->SetBranchAddress("Reco_QQ_4mom_pt", &Reco_QQ_4mom_pt);
 	chain->SetBranchAddress("Reco_QQ_4mom_eta", &Reco_QQ_4mom_eta);
 	chain->SetBranchAddress("Reco_QQ_4mom_phi", &Reco_QQ_4mom_phi);
 	
 	TH1D *h_mu = new TH1D("h_mu", "", 100, 20., 120.);
 	TH1D *h_QQ = new TH1D("h_QQ", "", 100, 20., 120.);

 	for(Long64_t i = 0; i < nEntries; i++){

 		chain->GetEntry(i); //Get EVENT i.
 		Short_t n_of_muons = Reco_mu_4mom_pt->size(); //Number of reco muons on the EVENT i;

 		if (n_of_muons < 2) continue; //If there are less than 2 muons in the event, skip event...

 		for(size_t j = 0; j < n_of_muons; j++){

 			for(size_t k = j+1; k < n_of_muons; k++){ //Double for loop because we have to check every opposite sign pair inside EVENT i.

 				if(Reco_mu_charge[j] * Reco_mu_charge[k] < 0){ //If they have opposite sign...

 					//TLorentzVector object, one for each muon.
					TLorentzVector mu1;
					TLorentzVector mu2;

					mu1.SetPtEtaPhiM(
					Reco_mu_4mom_pt->at(j),
					Reco_mu_4mom_eta->at(j),
					Reco_mu_4mom_phi->at(j),
					Reco_mu_4mom_m->at(j)
					);

					mu2.SetPtEtaPhiM(
				    Reco_mu_4mom_pt->at(k),
				    Reco_mu_4mom_eta->at(k),
				    Reco_mu_4mom_phi->at(k),
				    Reco_mu_4mom_m->at(k)
					);

					TLorentzVector Z = mu1 + mu2;
					double mass = Z.M();
 					
 					h_mu->Fill(mass);
 				}
 			}
 		}

 		for(size_t j = 0; j < Reco_QQ_4mom_m->size(); j++){

 			h_QQ->Fill(Reco_QQ_4mom_m->at(j));
 		}

 		if (i % (nEntries / 100) == 0){ //Just to check the progress.

	        int percent = static_cast<int>(100.0 * i / nEntries+0.5);

    	    std::cout << "\r"
        	          << percent
            	      << "% complete..."
                	  << std::flush;
	    }
 	}

 	std::cout << "\r100% complete!" << std::endl;
	
	TCanvas *c = new TCanvas();
	h_mu->GetXaxis()->SetTitle("Mass (GeV)");
	h_mu->GetYaxis()->SetTitle("n of dimuon candidates");
	h_mu->SetLineColor(kOrange + 7);
	h_mu->Draw();
	h_QQ->Draw("same");
	c->SaveAs("../../../mnt/c/Users/lucas/OneDrive/Documentos/TLorentzVector test.png");

    delete chain;
    delete c;
    delete h_mu;
    delete h_QQ;	
}

void muon_phi(){

	TChain *chain = new TChain("hionia/myTree");
	chain->Add((MC_PbPb2024.path + "Oniatree_PowhegZtoMuMu_PbPb2024_*.root").c_str());

	int nFiles = chain->GetListOfFiles()->GetEntries();
	std::cout << "Number of files = " << nFiles << std::endl;
	
	if (nFiles != 20)
	{
	   	std::cout << "ERROR: Wrong number of files were added to the chain!" << std::endl;
    	return;
	}

	std::vector<float> *Reco_mu_4mom_phi = nullptr;

	Long64_t nEntries = chain->GetEntries();
 	chain->SetBranchAddress("Reco_mu_4mom_phi", &Reco_mu_4mom_phi);
 	
 	TH1D *h = new TH1D("h", "Reco muon phi distribution;#phi (rad);n of reco muons", 64, -3.2, 3.2);

 	for(Long64_t i = 0; i < nEntries; i++){

 		chain->GetEntry(i); //Get event i.

 		for(size_t j = 0; j < Reco_mu_4mom_phi->size(); j++){

 			h->Fill(Reco_mu_4mom_phi->at(j)); //One fill per reco muon...
 		}

 		if (i % (nEntries / 100) == 0){ //Just to check the progress.

	        int percent = static_cast<int>(100.0 * i / nEntries+0.5);

    	    std::cout << "\r"
        	          << percent
            	      << "% complete..."
                	  << std::flush;
	    }
 	}

 	std::cout << "\r100% complete!" << std::endl;
	
	TCanvas *c = new TCanvas();
	//h->GetXaxis()->SetTitle("pT of dimuon candidates");
	//h->GetYaxis()->SetTitle("Events");
	h->Draw();
	c->SaveAs("../../../mnt/c/Users/lucas/OneDrive/Documentos/phi of reco single muons.png");

    delete chain;
    delete c;
    delete h;
}


void muon_mass(){

	TChain *chain = new TChain("hionia/myTree");
	chain->Add((MC_PbPb2024.path + "Oniatree_PowhegZtoMuMu_PbPb2024_*.root").c_str());

	int nFiles = chain->GetListOfFiles()->GetEntries();
	std::cout << "Number of files = " << nFiles << std::endl;
	
	if (nFiles != 20)
	{
	   	std::cout << "ERROR: Wrong number of files were added to the chain!" << std::endl;
    	return;
	}

	std::vector<float> *Reco_mu_4mom_m = nullptr;

	Long64_t nEntries = chain->GetEntries();
 	chain->SetBranchAddress("Reco_mu_4mom_m", &Reco_mu_4mom_m);
 	
 	TH1D *h = new TH1D("h", "Reco muon mass distribution;m (GeV);n of reco muons", 100, 0., 1.);

 	for(Long64_t i = 0; i < nEntries; i++){

 		chain->GetEntry(i); //Get event i.

 		for(size_t j = 0; j < Reco_mu_4mom_m->size(); j++){

 			h->Fill(Reco_mu_4mom_m->at(j)); //One fill per reco muon...
 		}

 		if (i % (nEntries / 100) == 0){ //Just to check the progress.

	        int percent = static_cast<int>(100.0 * i / nEntries+0.5);

    	    std::cout << "\r"
        	          << percent
            	      << "% complete..."
                	  << std::flush;
	    }
 	}

 	std::cout << "\r100% complete!" << std::endl;
	
	TCanvas *c = new TCanvas();
	//h->GetXaxis()->SetTitle("pT of dimuon candidates");
	//h->GetYaxis()->SetTitle("Events");
	h->Draw();
	c->SaveAs("../../../mnt/c/Users/lucas/OneDrive/Documentos/mass of reco single muons.png");

    delete chain;
    delete c;
    delete h;	
}

void muon_pseudorap(){

	TChain *chain = new TChain("hionia/myTree");
	chain->Add((MC_PbPb2024.path + "Oniatree_PowhegZtoMuMu_PbPb2024_*.root").c_str());

	int nFiles = chain->GetListOfFiles()->GetEntries();
	std::cout << "Number of files = " << nFiles << std::endl;
	
	if (nFiles != 20)
	{
	   	std::cout << "ERROR: Wrong number of files were added to the chain!" << std::endl;
    	return;
	}

	std::vector<float> *Reco_mu_4mom_eta = nullptr;

	Long64_t nEntries = chain->GetEntries();
 	chain->SetBranchAddress("Reco_mu_4mom_eta", &Reco_mu_4mom_eta);
	if (chain->SetBranchAddress("Reco_mu_4mom_eta",&Reco_mu_4mom_eta) < 0)
	{
    	std::cout << "Couldn't find branch Reco_mu_4mom_eta\n";
    	return;
	}
 	
 	TH1D *h = new TH1D("h", "Reco muon #eta distribution;#eta;n of reco muons", 100, -5., 5.);

 	for(Long64_t i = 0; i < nEntries; i++){

 		chain->GetEntry(i); //Get event i.

 		for(size_t j = 0; j < Reco_mu_4mom_eta->size(); j++){

 			h->Fill(Reco_mu_4mom_eta->at(j)); //One fill per reco muon...
 		}

 		if (i % (nEntries / 100) == 0){ //Just to check the progress.

	        int percent = static_cast<int>(100.0 * i / nEntries+0.5);

    	    std::cout << "\r"
        	          << percent
            	      << "% complete..."
                	  << std::flush;
	    }
 	}

 	std::cout << "\r100% complete!" << std::endl;
	
	TCanvas *c = new TCanvas();
	//h->GetXaxis()->SetTitle("pT of dimuon candidates");
	//h->GetYaxis()->SetTitle("Events");
	h->Draw();
	c->SaveAs("../../../mnt/c/Users/lucas/OneDrive/Documentos/pseudorapidity of reco single muons.png");

    delete chain;
    delete c;
    delete h;
}

void reco_mu_size(){

	TChain *chain = new TChain("hionia/myTree");
	chain->Add((MC_PbPb2024.path + "Oniatree_PowhegZtoMuMu_PbPb2024_*.root").c_str());

	int nFiles = chain->GetListOfFiles()->GetEntries();
	std::cout << "Number of files = " << nFiles << std::endl;
	
	if (nFiles != 20)
	{
	   	std::cout << "ERROR: Wrong number of files were added to the chain!" << std::endl;
    	return;
	}

	Short_t Reco_mu_size;

	Long64_t nEntries = chain->GetEntries();
 	chain->SetBranchAddress("Reco_mu_size", &Reco_mu_size);
 	
 	TH1D *h = new TH1D("h", "Distribution of reco muons;number of reco muons (GeV);events", 50, 0, 50);

 	for(Long64_t i = 0; i < nEntries; i++){

 		chain->GetEntry(i); //Get event i.

 		h->Fill(Reco_mu_size); //One fill per reco muon...

 		if (i % (nEntries / 100) == 0){ //Just to check the progress.

	        int percent = static_cast<int>(100.0 * i / nEntries+0.5);

    	    std::cout << "\r"
        	          << percent
            	      << "% complete..."
                	  << std::flush;
	    }
 	}

 	std::cout << "\r100% complete!" << std::endl;
	
	TCanvas *c = new TCanvas();
	//h->GetXaxis()->SetTitle("pT of dimuon candidates");
	//h->GetYaxis()->SetTitle("Events");
	h->Draw();
	c->SaveAs("../../../mnt/c/Users/lucas/OneDrive/Documentos/number of reco single muons.png");

    delete chain;
    delete c;
    delete h;
}

void muon_pt(){

	TChain *chain = new TChain("hionia/myTree");
	chain->Add((MC_PbPb2024.path + "Oniatree_PowhegZtoMuMu_PbPb2024_*.root").c_str());

	int nFiles = chain->GetListOfFiles()->GetEntries();
	std::cout << "Number of files = " << nFiles << std::endl;
	
	if (nFiles != 20)
	{
	   	std::cout << "ERROR: Wrong number of files were added to the chain!" << std::endl;
    	return;
	}

	std::vector<float> *Reco_mu_4mom_pt = nullptr;

	Long64_t nEntries = chain->GetEntries();
 	chain->SetBranchAddress("Reco_mu_4mom_pt", &Reco_mu_4mom_pt);
	if (chain->SetBranchAddress("Reco_mu_4mom_pt",&Reco_mu_4mom_pt) < 0)
	{
    	std::cout << "Couldn't find branch Reco_mu_4mom_pt\n";
    	return;
	}
 	
 	TH1D *h = new TH1D("h", "Reco muon pT distribution;p_{T} (GeV);n of reco muons", 100, 0., 200.);

 	for(Long64_t i = 0; i < nEntries; i++){

 		chain->GetEntry(i); //Get event i.

 		for(size_t j = 0; j < Reco_mu_4mom_pt->size(); j++){

 			h->Fill(Reco_mu_4mom_pt->at(j)); //One fill per reco muon...
 		}

 		if (i % (nEntries / 100) == 0){ //Just to check the progress.

	        int percent = static_cast<int>(100.0 * i / nEntries+0.5);

    	    std::cout << "\r"
        	          << percent
            	      << "% complete..."
                	  << std::flush;
	    }
 	}

 	std::cout << "\r100% complete!" << std::endl;
	
	TCanvas *c = new TCanvas();
	//h->GetXaxis()->SetTitle("pT of dimuon candidates");
	//h->GetYaxis()->SetTitle("Events");
	c->SetLogy();
	h->Draw();
	c->SaveAs("../../../mnt/c/Users/lucas/OneDrive/Documentos/pt of reco single muons.png");

    delete chain;
    delete c;
    delete h;
}

void dimuon_mass(){

	TChain *chain = new TChain("hionia/myTree");
	chain->Add((MC_PbPb2024.path + "Oniatree_PowhegZtoMuMu_PbPb2024_*.root").c_str());

	int nFiles = chain->GetListOfFiles()->GetEntries();
	std::cout << "Number of files = " << nFiles << std::endl;
	
	if (nFiles != 20)
	{
	   	std::cout << "ERROR: Wrong number of files were added to the chain!" << std::endl;
    	return;
	}

	std::vector<float> *Reco_QQ_4mom_m = nullptr;

	Long64_t nEntries = chain->GetEntries();
 	chain->SetBranchAddress("Reco_QQ_4mom_m", &Reco_QQ_4mom_m);
	if (chain->SetBranchAddress("Reco_QQ_4mom_m",&Reco_QQ_4mom_m) < 0)
	{
    	std::cout << "Couldn't find branch Reco_QQ_4mom_m\n";
    	return;
	}
 	
 	TH1D *h = new TH1D("h", "Dimuon mass distribution;mass (GeV);Candidates", 180, 30., 120.);

 	for(Long64_t i = 0; i < nEntries; i++){

 		chain->GetEntry(i); //Get event i.

 		for(size_t j = 0; j < Reco_QQ_4mom_m->size(); j++){

 			h->Fill(Reco_QQ_4mom_m->at(j)); //One fill per dimuon candidate...
 		}

 		if (i % (nEntries / 100) == 0){ //Just to check the progress.

	        int percent = static_cast<int>(100.0 * i / nEntries+0.5);

    	    std::cout << "\r"
        	          << percent
            	      << "% complete..."
                	  << std::flush;
	    }
 	}

 	std::cout << "\r100% complete!" << std::endl;
	
	TCanvas *c = new TCanvas();
	//h->GetXaxis()->SetTitle("pT of dimuon candidates");
	//h->GetYaxis()->SetTitle("Events");
	h->Draw();
	c->SaveAs("../../../mnt/c/Users/lucas/OneDrive/Documentos/mass of reco dimuon candidates.png");

    delete chain;
    delete c;
    delete h;
}

void dimuon_phi(){

	TChain *chain = new TChain("hionia/myTree");
	chain->Add((MC_PbPb2024.path + "Oniatree_PowhegZtoMuMu_PbPb2024_*.root").c_str());

	int nFiles = chain->GetListOfFiles()->GetEntries();
	std::cout << "Number of files = " << nFiles << std::endl;
	
	if (nFiles != 20)
	{
	   	std::cout << "ERROR: Wrong number of files were added to the chain!" << std::endl;
    	return;
	}

	std::vector<float> *Reco_QQ_4mom_phi = nullptr;

	Long64_t nEntries = chain->GetEntries();
 	chain->SetBranchAddress("Reco_QQ_4mom_phi", &Reco_QQ_4mom_phi);
	if (chain->SetBranchAddress("Reco_QQ_4mom_phi",&Reco_QQ_4mom_phi) < 0)
	{
    	std::cout << "Couldn't find branch Reco_QQ_4mom_phi\n";
    	return;
	}
 	
 	TH1D *h = new TH1D("h", "Dimuon azimuthal angle;#phi (rad);Candidates", 64, -3.2, 3.2);

 	for(Long64_t i = 0; i < nEntries; i++){

 		chain->GetEntry(i); //Get event i.

 		for(size_t j = 0; j < Reco_QQ_4mom_phi->size(); j++){

 			h->Fill(Reco_QQ_4mom_phi->at(j)); //One fill per dimuon candidate...
 		}

 		if (i % (nEntries / 100) == 0){ //Just to check the progress.

	        int percent = static_cast<int>(100.0 * i / nEntries+0.5);

    	    std::cout << "\r"
        	          << percent
            	      << "% complete..."
                	  << std::flush;
	    }
 	}

 	std::cout << "\r100% complete!" << std::endl;
	
	TCanvas *c = new TCanvas();
	//h->GetXaxis()->SetTitle("pT of dimuon candidates");
	//h->GetYaxis()->SetTitle("Events");
	h->Draw();
	c->SaveAs("../../../mnt/c/Users/lucas/OneDrive/Documentos/phi of reco dimuon candidates.png");

    delete chain;
    delete c;
    delete h;
}

void dimuon_pseudorap(){

	TChain *chain = new TChain("hionia/myTree");
	chain->Add((MC_PbPb2024.path + "Oniatree_PowhegZtoMuMu_PbPb2024_*.root").c_str());

	int nFiles = chain->GetListOfFiles()->GetEntries();
	std::cout << "Number of files = " << nFiles << std::endl;
	
	if (nFiles != 20)
	{
	   	std::cout << "ERROR: Wrong number of files were added to the chain!" << std::endl;
    	return;
	}

	std::vector<float> *Reco_QQ_4mom_eta = nullptr;

	Long64_t nEntries = chain->GetEntries();
 	chain->SetBranchAddress("Reco_QQ_4mom_eta", &Reco_QQ_4mom_eta);
	if (chain->SetBranchAddress("Reco_QQ_4mom_eta",&Reco_QQ_4mom_eta) < 0)
	{
    	std::cout << "Couldn't find branch Reco_QQ_4mom_eta\n";
    	return;
	}
 	
 	TH1D *h = new TH1D("h","Dimuon pseudorapidity;#eta;Candidates", 100, -5., 5.);

 	for(Long64_t i = 0; i < nEntries; i++){

 		chain->GetEntry(i); //Get event i.

 		for(size_t j = 0; j < Reco_QQ_4mom_eta->size(); j++){

 			h->Fill(Reco_QQ_4mom_eta->at(j)); //One fill per dimuon candidate...
 		}

 		if (i % (nEntries / 100) == 0){ //Just to check the progress.

	        int percent = static_cast<int>(100.0 * i / nEntries+0.5);

    	    std::cout << "\r"
        	          << percent
            	      << "% complete..."
                	  << std::flush;
	    }
 	}

 	std::cout << "\r100% complete!" << std::endl;
	
	TCanvas *c = new TCanvas();
	//h->GetXaxis()->SetTitle("pT of dimuon candidates");
	//h->GetYaxis()->SetTitle("Events");
	h->Draw();
	c->SaveAs("../../../mnt/c/Users/lucas/OneDrive/Documentos/eta of reco dimuon candidates.png");

    delete chain;
    delete c;
    delete h;
}

void dimuon_pt(){

	TChain *chain = new TChain("hionia/myTree");
	chain->Add((MC_PbPb2024.path + "Oniatree_PowhegZtoMuMu_PbPb2024_*.root").c_str());
	
	int nFiles = chain->GetListOfFiles()->GetEntries();
	std::cout << "Number of files = " << nFiles << std::endl;
	
	if (nFiles != 20)
	{
	   	std::cout << "ERROR: Wrong number of files were added to the chain!" << std::endl;
    	return;
	}

	Short_t Reco_QQ_size;
	std::vector<float> *Reco_QQ_4mom_pt = nullptr;

	Long64_t nEntries = chain->GetEntries();
 	chain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
 	chain->SetBranchAddress("Reco_QQ_4mom_pt", &Reco_QQ_4mom_pt);

 	TH1D *h = new TH1D("h","Dimuon p_{T};p_{T} (GeV);Candidates", 100, 0, 100);

 	for(Long64_t i = 0; i < nEntries; i++){

 		chain->GetEntry(i); //Get event i.

 		for(Short_t j = 0; j < Reco_QQ_4mom_pt->size(); j++){

 			h->Fill(Reco_QQ_4mom_pt->at(j)); //One fill per dimuon candidate...
 		}

 		if (i % (nEntries / 100) == 0){ //Just to check the progress.

	        int percent = static_cast<int>(100.0 * i / nEntries+0.5);

    	    std::cout << "\r"
        	          << percent
            	      << "% complete..."
                	  << std::flush;
	    }
 	}

 	std::cout << "\r100% complete!" << std::endl;
	
	TCanvas *c = new TCanvas();
	//h->GetXaxis()->SetTitle("pT of dimuon candidates");
	//h->GetYaxis()->SetTitle("Events");
	h->Draw();
	c->SetLogy();
	c->Update();
	c->SaveAs("../../../mnt/c/Users/lucas/OneDrive/Documentos/pT of reco dimuon candidates.png");

    delete chain;
    delete c;
    delete h;
}

void dimuon_rapidity(){

		TChain *chain = new TChain("hionia/myTree");
	chain->Add((MC_PbPb2024.path + "Oniatree_PowhegZtoMuMu_PbPb2024_*.root").c_str());

	int nFiles = chain->GetListOfFiles()->GetEntries();
	std::cout << "Number of files = " << nFiles << std::endl;

	if (nFiles != 20) //Just a check to make sure all 20 files were loaded.
	{
	   	std::cout << "ERROR: Wrong number of files were added to the chain!" << std::endl;
    	return;
	}

	Long64_t nEntries = chain->GetEntries(); //Get number of EVENTS.

 	Short_t Reco_QQ_size;
 	std::vector<float> *Reco_QQ_4mom_m = nullptr; //Mass of reco Z boson candidates.
	std::vector<float> *Reco_QQ_4mom_pt = nullptr; //pT of reco Z boson candidates.
	std::vector<float> *Reco_QQ_4mom_eta = nullptr; //eta of reco Z boson candidates.
	std::vector<float> *Reco_QQ_4mom_phi = nullptr; //phi of reco Z boson candidates.

 	chain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
 	chain->SetBranchAddress("Reco_QQ_4mom_m", &Reco_QQ_4mom_m);
 	chain->SetBranchAddress("Reco_QQ_4mom_pt", &Reco_QQ_4mom_pt);
 	chain->SetBranchAddress("Reco_QQ_4mom_eta", &Reco_QQ_4mom_eta);
 	chain->SetBranchAddress("Reco_QQ_4mom_phi", &Reco_QQ_4mom_phi);
 	
 	TH1D *h = new TH1D("h_mu", "", 100, -5., 5.);

 	for(Long64_t i = 0; i < nEntries; i++){

 		chain->GetEntry(i); //Get EVENT i.

 		for(Short_t j = 0; j < Reco_QQ_size; j++){

			TLorentzVector Z;
			//Set 4-momentum vector for j-th dimuon candidate.
        	Z.SetPtEtaPhiM(
                Reco_QQ_4mom_pt->at(j),
                Reco_QQ_4mom_eta->at(j),
                Reco_QQ_4mom_phi->at(j),
                Reco_QQ_4mom_m->at(j)
        	);

			h->Fill(Z.Rapidity());
		}

 		if (i % (nEntries / 100) == 0){ //Just to check the progress.

	        int percent = static_cast<int>(100.0 * i / nEntries+0.5);

    	    std::cout << "\r"
        	          << percent
            	      << "% complete..."
                	  << std::flush;
	    }
 	}

 	std::cout << "\r100% complete!" << std::endl;
	
	TCanvas *c = new TCanvas();
	h->Draw();
	c->SaveAs("../../../mnt/c/Users/lucas/OneDrive/Documentos/rapidity distribution of reco dimuon candidates.png");

    delete chain;
    delete c;
    delete h;
}

void test_chain(){

	TChain *chain = new TChain("hionia/myTree");

	chain->Add((MC_PbPb2024.path + "Oniatree_PowhegZtoMuMu_PbPb2024_*.root").c_str());

	int nFiles = chain->GetListOfFiles()->GetEntries();

	std::cout << "Number of files = " << nFiles << std::endl;
	if (nFiles != 0)
	{
	   	std::cout << "ERROR: Wrong number of files were added to the chain!" << std::endl;
    	return;
	}

    std::cout << "Total entries = " << chain->GetEntries() << std::endl;
    //

    Short_t Reco_QQ_size;
	std::vector<float> *Reco_QQ_4mom_m = nullptr;

    Short_t Reco_mu_charge;          // if it's a scalar
	std::vector<float> *Reco_mu_4mom_pt = nullptr;

	Long64_t nEntries = chain->GetEntries();
 	chain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
 	chain->SetBranchAddress("Reco_QQ_4mom_pt", &Reco_mu_4mom_pt);

 	TH1D *h = new TH1D("h","QQ candidates", 20, 0, 20);

	for (Long64_t i = 0; i < nEntries; i++){
    
    	chain->GetEntry(i);

	    h->Fill(Reco_QQ_size);

	    if (i % (nEntries / 100) == 0){ //Just to check the progress.

	        int percent = static_cast<int>(100.0 * i / nEntries+0.5);

    	    std::cout << "\r"
        	          << percent
            	      << "% complete..."
                	  << std::flush;
	    	}
	
	}

	std::cout << "\r100% complete!" << std::endl;
	
	TCanvas *c = new TCanvas();
	h->GetXaxis()->SetTitle("Number of dimuon candidates");
	h->GetYaxis()->SetTitle("Events");
	h->Draw();
	c->SaveAs("../../../mnt/c/Users/lucas/OneDrive/Documentos/number of reco dimuon candidates.png");

    //chain->GetListOfFiles()->Print();
	//chain->GetListOfBranches()->Print(); // 83 branches. Verify consistency with data trees later. Should be equal.

    delete chain;
    delete c;
    delete h;
}