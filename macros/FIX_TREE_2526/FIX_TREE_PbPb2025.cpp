/*
Let's fix the PbPb2025 data TTree.
The centrality branch from the TTree is wrong. So i'll create a new branch called centrality_fixed and fill it with the correct values.
The centrality_fixed values will be calculated from the SumET_HF using the 2025 calibration table.
Tables available at: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideHeavyIonCentrality#Nominal_data_AN2
*/
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TBranch.h>
#include <iostream>
#include <string>

//Location of root file.
std::string PathToROOTFile = "/home/lucas/Documents/CMS/z_boson_analysis/Data/PbPb2025/HighPtMuon_PbPb2025Data.root"; //IFT

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


Int_t getHiBinFromhiHF(const Double_t hiHF);
void FixTREE();


void FIX_TREE_PbPb2025(){

    FixTREE();
}


void FixTREE(){

    TFile* rootFile = new TFile(PathToROOTFile.c_str(), "READ");
    TFile* outputFile = new TFile("Fixed_Tree_PbPb2025.root", "RECREATE");

    TTree* OldTree = (TTree*)rootFile->Get("hionia/DimuonTree");

    //Copy all branches except the centrality branch to a new tree.
    OldTree->SetBranchStatus("*", 1);
    OldTree->SetBranchStatus("Centrality", 0);

    TDirectory* dir = outputFile->mkdir("hionia");
    dir->cd();

    //Clone the old tree, keeping only the active branches.
    TTree* newTree = OldTree->CloneTree(-1, "fast");

    //Get the SumET_HF branch from the old tree.
    Float_t SumET_HF;
    OldTree->SetBranchAddress("SumET_HF", &SumET_HF);

    //Create branch for Centrality.
    Int_t Centrality = -1;
    TBranch* newBranch = newTree->Branch("Centrality", &Centrality, "Centrality/I");

    for(Long64_t i = 0; i < OldTree->GetEntries(); ++i){//Loop through all EVENTS in the TREE.

        OldTree->GetEntry(i); //Get event i.

        //Fill centrality according to the CentTable.
        Centrality = getHiBinFromhiHF(SumET_HF);//Returns a value between 0 and 199.
        newBranch->Fill();
    }//Exiting event-by-event loop.

    std::cout
    << "Old tree: " << OldTree->GetEntries() << '\n'
    << "New tree: " << newTree->GetEntries() << '\n'
    << "New Centrality: " << newBranch->GetEntries() << '\n';

    newTree->Write();

    outputFile->Close();
    rootFile->Close();

    delete outputFile;
    delete rootFile;

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
