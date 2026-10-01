/*
Fix the PbPb2026 data TTree. This is a COPY of the macro that fixes the PbPb2025 data TTree.
The centrality_fixed values will be calculated from the SumET_HF using the 2026 calibration table.
Tables available at: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideHeavyIonCentrality#Nominal_data_AN2
*/
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TBranch.h>
#include <iostream>
#include <string>

//Location of root file.
std::string PathToROOTFile = "/home/lucas/Documents/CMS/z_boson_analysis/Data/PbPb2026/HighPtMuon_PbPb2026Data.root"; //IFT

// Bin edges: paste the table for your year and variation from the sections above.
const Int_t nBins = 200;
const Double_t binTable[nBins+1] = {0, 12.767, 13.751, 14.7043, 15.6582, 16.6397, 17.6616, 18.6705, 19.6969, 20.7669, 21.8597,
    22.9919, 24.1573, 25.3552, 26.6099, 27.9199, 29.2841, 30.6558, 32.077, 33.5808, 35.1456, 36.782, 38.4701, 40.232, 42.0559,
    43.9605, 45.981, 48.0874, 50.2921, 52.5807, 54.9396, 57.4183, 60.0068, 62.7031, 65.5132, 68.4721, 71.5323, 74.7414, 78.024,
    81.4764, 84.965, 88.6944, 92.5398, 96.5507, 100.612, 104.479, 108.484, 112.636, 116.948, 121.432, 126.058, 130.906, 
    135.923, 141.112, 146.553, 152.153, 157.942, 163.953, 170.222, 176.713, 183.398, 190.347, 197.517, 204.875, 212.537, 
    220.463, 228.629, 237.056, 245.746, 254.745, 264.039, 273.629, 283.464, 293.673, 304.174, 314.967, 326.061, 337.498, 
    349.252, 361.295, 373.734, 386.497, 399.606, 413.105, 426.885, 441.092, 455.655, 470.513, 485.821, 501.512, 517.572, 534.12, 
    551.108, 568.368, 586.024, 604.069, 622.682, 641.697, 661.059, 680.939, 701.218, 721.863, 743.086, 764.793, 786.859, 809.547, 
    832.446, 855.851, 879.8, 904.355, 929.408, 954.93, 980.838, 1007.34, 1034.43, 1061.95, 1090.07, 1118.64, 1147.77, 1177.47, 
    1207.69, 1238.49, 1269.81, 1301.88, 1334.66, 1367.75, 1401.56, 1435.87, 1470.49, 1505.69, 1541.86, 1578.75, 1616.2, 1654.09, 
    1692.54, 1732.07, 1772.06, 1812.68, 1853.88, 1895.75, 1938.2, 1981.31, 2025.06, 2069.85, 2115.4, 2161.83, 2208.65, 2256.63, 
    2305.06, 2354.58, 2404.85, 2455.93, 2507.86, 2560.24, 2613.65, 2668.22, 2723.59, 2779.64, 2837.08, 2895.21, 2954.08, 3014.03, 
    3075.11, 3137.47, 3201.18, 3265.85, 3331.38, 3397.75, 3465.58, 3534.59, 3604.56, 3676.17, 3749.19, 3822.85, 3898.35, 3974.95, 
    4053.36, 4133.43, 4214.92, 4298.2, 4383.14, 4469.64, 4558.53, 4648.9, 4741.61, 4836.09, 4932.74, 5031.83, 5132.76, 5236.97, 
    5343.17, 5452.6, 5566.15, 5681.9, 5802.18, 5926.03, 6054.17, 6188.47, 6333.47, 6510.25, 12150.3};


Int_t getHiBinFromhiHF(const Double_t hiHF);
void FixTREE();


void FIX_TREE_PbPb2026(){

    FixTREE();
}


void FixTREE(){

    TFile* rootFile = new TFile(PathToROOTFile.c_str(), "READ");
    TFile* outputFile = new TFile("Fixed_Tree_PbPb2026.root", "RECREATE");

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
