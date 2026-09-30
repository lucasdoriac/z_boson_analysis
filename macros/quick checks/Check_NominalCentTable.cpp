
//---Libraries
#include <TROOT.h>
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
#include <TVector2.h>
#include <algorithm>


//Location of datasets
std::string BasePath = "/home/lucas/Documents/CMS/z_boson_analysis/"; //IFT
//std::string BasePath = "/home/lucasdoriac/z_boson_analysis/data/"; //Home


// ##############################################################################
// ##############################################################################

//---Enumerates
enum class SampleType {
    Data,
    MC
};

enum class CollisionSystem {
    ppRef,
    PbPb
};

//---Structs
struct Dataset {
    std::string name;
    SampleType type;
    CollisionSystem system;
    int year;
    std::string treeName;
    std::string filePattern;
    std::string basePath;

    bool hasCentrality;//Or maybe is AA
    bool applyTrigger;
    ULong64_t triggerBit;
};

Dataset datasets[] = {
    
    {
        "ppRef2024_Data",
        SampleType::Data,
        CollisionSystem::ppRef,
        2024,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_ppRef2024.root",
        BasePath + "Data/ppRef2024/",
        false,
        false,
        0ULL
    },
    
    {
        "PbPb2023_Data",
        SampleType::Data,
        CollisionSystem::PbPb,
        2023,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_PbPb2023.root",
        BasePath + "Data/PbPb2023/",
        true,
        true,
        1ULL << 6 //'HLT_HIL2SingleMu7_v'
    },

    {
        "PbPb2024_Data",
        SampleType::Data,
        CollisionSystem::PbPb,
        2024,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_PbPb2024Data.root",
        BasePath + "Data/PbPb2024/",
        true,
        true,
        1ULL << 7 //'HLT_HIL2SingleMu12_v'
    },

    {
        "PbPb2025_Data",
        SampleType::Data,
        CollisionSystem::PbPb,
        2025,
        "hionia/DimuonTree",
        "HighPtMuon_PbPb2025Data.root",
        BasePath + "Data/PbPb2025/",
        true,
        true,
        1ULL << 7 //'HLT_HIL2SingleMu12_v'
    },

    {
        "PbPb2026_Data",
        SampleType::Data,
        CollisionSystem::PbPb,
        2026,
        "hionia/DimuonTree",
        "HighPtMuon_PbPb2026Data.root",
        BasePath + "Data/PbPb2026/",
        true,
        true,
        1ULL << 7 //'HLT_HIL2SingleMu12_v'
    },

    {
        "PbPb2024_MC",
        SampleType::MC,
        CollisionSystem::PbPb,
        2024,
        "hionia/myTree",
        "Oniatree_PowhegZtoMuMu_PbPb2024_*.root",
        BasePath + "MC/PbPb2024/DYto2Mu_MLL-50_TuneCP5_5p36TeV_powheg-pythia8/PowhegEmbedded_March9/260309_143939/0000/",
        false,
        false,
        0ULL
    },

    {
        "ppRef2024_MC",
        SampleType::MC,
        CollisionSystem::ppRef,
        2024,
        "hionia/myTree",
        "Oniatree_PowhegZtoMuMu_ppRef2024_*.root",
        BasePath + "MC/ppRef2024/DYToMuMu_M-50_TuneCP5_5p36TeV_powheg-pythia8/Powheg_ppRefPileup_March20/260320_125046/0000/",
        false,
        false,
        0ULL
    }
};


// Bin edges: paste the table for your year and variation from the sections above.
const Int_t nBins = 200;
const Double_t binTable[nBins+1] = {0, 11.5689, 12.4389, 13.3363, 14.2201, 15.1293, 16.0823, 16.9587, 17.9343, 18.8936, 19.8677,
20.8638, 21.9125, 22.9434, 24.0678, 25.2431, 26.3926, 27.6413, 28.946, 30.2756, 31.6543, 33.1128, 34.6179, 36.1481, 37.7914, 
39.5531, 41.3402, 43.179, 45.1796, 47.2191, 49.3825, 51.6477, 53.9232, 56.3972, 58.8486, 61.4979, 64.2856, 67.1077, 69.8402, 
72.9557, 76.1467, 79.3769, 82.7643, 86.4625, 90.1211, 93.9367, 98.1249, 102.035, 106.015, 110.114, 114.383, 118.792, 123.387, 
128.15, 133.094, 138.236, 143.543, 149.056, 154.767, 160.674, 166.76, 173.102, 179.617, 186.345, 193.344, 200.559, 208.029, 
215.75, 223.692, 231.898, 240.317, 249.005, 257.953, 267.175, 276.639, 286.349, 296.401, 306.774, 317.492, 328.484, 339.789, 
351.358, 363.302, 375.515, 388.086, 400.922, 414.203, 427.638, 441.548, 455.909, 470.475, 485.474, 500.927, 516.731, 532.781, 
549.247, 566.125, 583.356, 600.881, 618.868, 637.213, 655.996, 675.17, 694.91, 714.993, 735.468, 756.41, 777.771, 799.49, 821.753, 
844.528, 867.693, 891.25, 915.28, 939.779, 964.756, 990.071, 1016.04, 1042.28, 1068.94, 1096.18, 1124.05, 1152.41, 1181.53, 
1211.12, 1241.12, 1271.41, 1302.34, 1333.98, 1366.26, 1398.99, 1432.18, 1466.15, 1500.5, 1535.59, 1571.27, 1607.6, 1644.25, 
1681.44, 1719.25, 1757.98, 1797.16, 1837.19, 1877.66, 1918.74, 1960.51, 2003.43, 2046.46, 2090.34, 2134.96, 2180.46, 2226.49, 
2273.22, 2320.56, 2368.93, 2418.37, 2468.43, 2519.45, 2571.34, 2624, 2677.87, 2732.18, 2787.62, 2844.07, 2901.51, 2959.73, 
3019.17, 3079.22, 3140.48, 3203.46, 3266.95, 3331.62, 3397.3, 3464.65, 3532.86, 3602.4, 3673.34, 3745.72, 3819.65, 3894.02, 
3970.59, 4048.6, 4128.86, 4210.48, 4294.01, 4379.37, 4466.51, 4555.7, 4647.1, 4741.5, 4838.52, 4937.76, 5039.49, 5144.82, 
5253.13, 5364.45, 5480.54, 5601.13, 5733.06, 5892.36, 9397.83};


void myFunction(const Dataset& dataset);
Int_t getHiBinFromhiHF(const Double_t hiHF);

void Check_NominalCentTable(){

    myFunction(datasets[3]);//PbPb2025_Data
}

void myFunction(const Dataset& dataset){

    //Load root file.
    std::string fullPath = dataset.basePath + dataset.filePattern;

    TChain *chain = new TChain(dataset.treeName.c_str());
    chain->Add(fullPath.c_str());

    std::cout << "\n> Number of files added to TChain = " << chain->GetListOfFiles()->GetEntries() << "\n" << std::endl;
    std::cout << "> Opening files " << fullPath << "\n" << std::endl;
    
    //Total number of events on Tree.
    Long64_t nEvents = chain->GetEntries();

    std::cout << "> Total number of events on tree = " << nEvents << "\n" << std::endl;

    Float_t SumET_HF;
    chain->SetBranchAddress("SumET_HF", &SumET_HF);

    TH1D* h1D_centrality = new TH1D("h1D_centrality",
        "Centrality from SumET_HF (PbPb2025);Centrality;Entries",
        200, 0, 200);

        
    for(Long64_t i = 0; i < nEvents; ++i){//Loop through all EVENTS in the CHAIN.

        chain->GetEntry(i); //Get event i.

        //Fill centrality according to the CentTable.
        Int_t centBin = getHiBinFromhiHF(SumET_HF);
        h1D_centrality->Fill(centBin);

    }//Exiting event-by-event loop.

    TCanvas *c = new TCanvas("c", "c", 800, 600);
    //c->SetLogy();
    h1D_centrality->Draw();
    c->SaveAs("CentralityDistribution.png");

}    

Int_t getHiBinFromhiHF(const Double_t hiHF){

    Int_t binPos = -1;

    for (int i = 0; i < nBins; ++i) {
        if (hiHF >= binTable[i] && hiHF < binTable[i+1]) {
            binPos = i;
            break;
        }
    }

    binPos = nBins - 1 - binPos;
    return (Int_t)(200*((Double_t)binPos)/((Double_t)nBins));
}
