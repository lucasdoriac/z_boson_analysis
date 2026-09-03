/*
I want to test a new way to declare and initialize the branch variables.
*/

#include "../headers/PbPbBranches.h"
#include "../headers/TheDatasets.h"
#include "../headers/basicFormatting.h"
#include "../headers/Miscellaneous.h"

using namespace std;

//---Macro settings
std::string plot_extension = ".png"; // ".png" for regular development and ".pdf" for final quality plots

//Cut values
const int MIN_CENTRALITY = 0.;
const int MAX_CENTRALITY = 180;
const double MAX_ZVTX = 15.0;

const double MINZ_MASS = 80.;
const double MAXZ_MASS = 100.;
const double RAPIDITYCUTVALUE = 2.4;

const float ETACUTVALUE = 2.4;
const double PTCUTVALUE = 10.;

//Map:
/*
datasets[0] = PbPb2023_Data
datasets[1] = PbPb2024_Data
datasets[2] = ppRef2024_Data
datasets[3] = PbPb2024_MC
datasets[4] = ppRef2024_MC
*/

void ApplyGoodSelection(const Dataset &dataset);

void testBranches(){

    ApplyGoodSelection(datasets[1]); //PbPb2024_Data
}

void ApplyGoodSelection(const Dataset &dataset){

    //Load root file.
    std::string fullPath = dataset.path + dataset.filePattern;

    TChain *chain = new TChain(dataset.treeName.c_str());
    chain->Add(fullPath.c_str());

    std::cout << "> Number of files = " << chain->GetListOfFiles()->GetEntries() << "\n" << std::endl;

    std::cout << "> Opening files " << fullPath << "\n" << std::endl;

    std::cout << "> Running function " << __func__ << " on " << dataset.name << "\n" << std::endl;
    
    //Total number of events on Tree.
    Long64_t nEvents = chain->GetEntries();    

    //Initialize the PbPb branches.
    SetPbPbBranches PbPb2024;
    PbPb2024.SetBranches(chain);

    cout << "> Number of events in the tree = " << nEvents << endl;

    //Event-level cut values.
    const int minCentrality = MIN_CENTRALITY;
    const int maxCentrality = MAX_CENTRALITY; //Cent max = maxCentrality/2.
    const double maxZvtx = MAX_ZVTX;

    //Z-level cut values.
    const double minZ_Mass = MINZ_MASS; 
    const double maxZ_Mass = MAXZ_MASS; 
    const double RapidityCutValue = RAPIDITYCUTVALUE;

    //Muon-level cut values.
    const float EtaCutValue = ETACUTVALUE;
    const double ptCutValue = PTCUTVALUE;
    bool MuPlIsTight;
    bool MuMiIsTight;

    //Observables
    TLorentzVector Z;
    double ptplus, ptminus;
    double etaplus, etaminus;

    //Trigger selection 'L2SingleMu12'.
    ULong64_t triggerBit = 1ULL << 7;

    int goodDimuonCount = 0;

    for(Long64_t iEvent = 0; iEvent < nEvents; iEvent++){

        chain->GetEntry(iEvent);

        //Good event selection
        bool goodCent = (PbPb2024.Centrality > minCentrality && PbPb2024.Centrality < maxCentrality);
        bool goodVertex = (std::abs(PbPb2024.zVtx) < maxZvtx);
        //Still missing 2 HF towers of 4 GeV of energy in the event.
        //Shapes of clusters compatibility.

        if (!goodCent) continue;
        if (!goodVertex) continue;

        for(Short_t j = 0; j < PbPb2024.Reco_Dimuon_size; j++){

            //Good Z selection
            bool goodMass = (PbPb2024.Reco_Dimuon_invMass->at(j) > minZ_Mass && PbPb2024.Reco_Dimuon_invMass->at(j) < maxZ_Mass);
            bool goodRapidity = (std::abs(PbPb2024.Reco_Dimuon_rapidity->at(j)) < RapidityCutValue);
            bool goodCharge = (PbPb2024.Reco_Dimuon_sign[j] == 0);
            bool isTriggerMatched = (PbPb2024.Reco_Dimuon_trig[j]) & triggerBit;
            
            if (!goodMass) continue;
            if (!goodRapidity) continue;
            if (!goodCharge) continue;
            if (!isTriggerMatched) continue;

            //Good muon selection
            Short_t muonPlusIndex = PbPb2024.Reco_Dimuon_muonPlusIndex[j];//Index of the muon+ in the Reco_Muon arrays.
            Short_t muonMinusIndex = PbPb2024.Reco_Dimuon_muonMinusIndex[j];//Index of the muon- in the Reco_Muon arrays.

            ptplus = PbPb2024.Reco_Muon_pt->at(muonPlusIndex); //pT of antimuon.
            ptminus = PbPb2024.Reco_Muon_pt->at(muonMinusIndex); //pT of corresponding muon.
            etaplus = PbPb2024.Reco_Muon_eta->at(muonPlusIndex); //Pseudorapidity of antimuon.
            etaminus = PbPb2024.Reco_Muon_eta->at(muonMinusIndex); //Pseudorapidity of corresponding muon.
            MuPlIsTight = PbPb2024.Reco_Muon_isTightCutBased[muonPlusIndex];
            MuMiIsTight = PbPb2024.Reco_Muon_isTightCutBased[muonMinusIndex];

            bool goodMuPl = (ptplus > ptCutValue)
                            && (std::abs(etaplus) < EtaCutValue)
                            && (MuPlIsTight);
            
            bool goodMuMi = (ptminus > ptCutValue)
                            && (std::abs(etaminus) < EtaCutValue)
                            && (MuMiIsTight);
            
            if (!goodMuPl || !goodMuMi) continue;

            goodDimuonCount++;
        }

        //Just to check the progress.
        if (iEvent % (nEvents / 100) == 0){

            int percent = static_cast<int>(100.0 * iEvent / nEvents+0.5);

            std::cout << "\r"
                      << percent
                      << "% complete..."
                      << std::flush;
        }//

    }//End of event loop.

    std::cout << "> Number of good dimuon candidates = " << goodDimuonCount << std::endl;

}