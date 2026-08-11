/*I want to test a better way to read the trees without writing 200 lines of code in the functions everytime.

Is the best way to use a library, a class or a struct?

The branches between each collisionsystem and sampletype can be different.

I want something like:

1. Get chain
2. Declare variables.
3. Read tree.

Suggestion to consider in the future:

include/

    Dataset.h

    EventBranches.h
    MuonBranches.h
    DimuonBranches.h

    TreeReader.h

src/

    TreeReader.cpp

analysis/

    plotPtSpectrum.cpp
    invariantMass.cpp
    asymmetry.cpp
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

//---Class
class TreeReader {	
	
	//Testing with event branches first.
	
	public:
	    Float_t zVtx;
	    Short_t nPV;
	    Int_t Centrality;
	    Short_t NpixelTracks;

	    void connectBranches(TChain* chain){
	    	chain->SetBranchAddress("zVtx", &zVtx);
		    chain->SetBranchAddress("nPV", &nPV);
		    chain->SetBranchAddress("Centrality", &Centrality);
		    chain->SetBranchAddress("NpixelTracks", &NpixelTracks);
	    }
		
//	private:

};

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
        "/home/lucasdoriac/z_boson_analysis/Data/PbPb2023/"
    },

    {
        "PbPb2024_Data",
        SampleType::Data,
        CollisionSystem::PbPb2024,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_PbPb2024Data.root",
        "/home/lucasdoriac/z_boson_analysis/Data/PbPb2024/"
    },

    {
        "ppRef2024_Data",
        SampleType::Data,
        CollisionSystem::ppRef2024,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_ppRef2024.root",
        "/home/lucasdoriac/z_boson_analysis/Data/ppRef2024/"
    },

    {
        "PbPb2024_MC",
        SampleType::MC,
        CollisionSystem::PbPb2024,
        "hionia/myTree",
        "Oniatree_PowhegZtoMuMu_PbPb2024_*.root",
        "/home/lucasdoriac/z_boson_analysis/MC/PbPb2024/DYto2Mu_MLL-50_TuneCP5_5p36TeV_powheg-pythia8/PowhegEmbedded_March9/260309_143939/0000/"
    },

    {
        "ppRef2024_MC",
        SampleType::MC,
        CollisionSystem::ppRef2024,
        "hionia/myTree",
        "Oniatree_PowhegZtoMuMu_ppRef2024_*.root",
        "/home/lucasdoriac/z_boson_analysis/MC/ppRef2024/DYToMuMu_M-50_TuneCP5_5p36TeV_powheg-pythia8/Powheg_ppRefPileup_March20/260320_125046/0000/"
    }
};

//---Miscellaneous
void printTreeContents(const Dataset &dataset);

//---main()
void test_treeReaderClass(){

	gROOT->SetBatch(kTRUE);

	printTreeContents(datasets[0]);

}

//---Miscellaneous
void printTreeContents(const Dataset &dataset){

    std::string fullPath = dataset.basePath + dataset.filePattern;

    TChain *chain = new TChain(dataset.treeName.c_str());
    chain->Add(fullPath.c_str());

    chain->Print();
    std::cout << "Number of files = " << chain->GetListOfFiles()->GetEntries() << std::endl;

    Long64_t nEvents = chain->GetEntries();
    std::cout << "Total number of events = " << nEvents << std::endl;

    std::cout << "Dataset: " << dataset.name << std::endl;
    std::cout << "Tree: " << dataset.treeName << std::endl;
    std::cout << "Path: " << fullPath << std::endl;

    std::cout << "Testing new class: \n" << std::endl;

    TreeReader reader;
    reader.connectBranches(chain);

    for(int i = 0; i < 10; ++i){
    	chain->GetEntry(i);
	    std::cout << "Event " << i << " | zVtx = " << reader.zVtx << std::endl;
    }

}