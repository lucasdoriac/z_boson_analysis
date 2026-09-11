/*
1- First candidate plot let's do the usual N(\mu+)(pT) and N(\mu-)(pT).
This is just the distribution of dN/dpT for the individual muon daughters and it's the distribution that
the theory paper uses to predict the shift in the peak of THESE distributions in the 30-50 GeV region.
On the bottom panel we will plot R(pT) = N(\mu+)/N(\mu-). We can do this for both PbPb and ppRef.

2- We would like to also look at the double ratio of the R(pT) distributions, i.e. R(PbPb)/R(ppRef).
For this we can plot the usual N(\mu+)(pT) and N(\mu-)(pT), as we did in the first candidate plot, but now for PbPb and ppRef in the same canvas. 
This would give four distributions in the top pannel, and in the bottom panel we can plot the double ratio R(PbPb)/R(ppRef).

Of course, we could also do that in the first candidate plot, since we are calculating R(pT) = N(\mu+)/N(\mu-) for each dataset,
i.e. we are already calculating R(PbPb) and R(ppRef). However, it is interesting at a first moment to plot both datasets separately,
and then later plot them in the same canvas.

3- Another observable option is the event-by-event relative difference in pT between the two muons. This is defined as:
    ΔpT = (pT(mu+) - pT(mu-)) / (pT(mu+) + pT(mu-)).
If we call that variable ΔpT, then we can calculate ΔpT(PbPb) and ΔpT(ppRef), event-by-event.
In an ideal scenario, one might expect ΔpT(PbPb) to be modified by the presence of the B field in the QGP,
with maximum effects on the 40-70% centrality bin (and maybe in the low pT region?).

The sweet observable would be the mean of ΔpT(PbPb), <ΔpT(PbPb)> per centrality bin,
as well as their difference with respect to the ppRef dataset, i.e. <ΔpT(PbPb)> - <ΔpT(ppRef)>.

4- The last simple observable we can try to check is the charge asymmetry of the muons in each dataset. This is defined as:
    A(X) = (N(mu+) - N(mu-)) / (N(mu+) + N(mu-)).
Notice that this is a YIELD asymmetry, NOT a pT asymmetry. It is the difference in the number of muons of each charge.
(Normalized by the total number of muons?)
We can calculate A(PbPb) and A(ppRef) and plot them in the same canvas.
In the bottom panel we can plot the difference in the charge asymmetry between PbPb and ppRef, i.e. ΔA = A(PbPb) - A(ppRef).

-------------------------------------
After doing the four candidate plots, now we move to another plot that gives a centrality-dependent observable.

The initial centrality bins would be:
Set 1:
- 0-10%
- 10-30%
- 30-50%
- 50-70%
- 70-100%

Or maybe another temptative list is:
Set 2:
- 0-10%
- 10-20%
- 20-30%
- 30-40%
- 40-100%

pT sectors:
- 20-30 GeV
- 30-40 GeV
- 40-50 GeV
- 50-60 GeV

Another thing to observe is that we don't need to analyze the whole pT spectrum.
The proposed effect is in the 30-50 GeV region, so we can focus on that region and maybe extend it a little bit to 20-60 GeV.

The observables that i was thinking about were:

5- Peak and mean of PbPb and ppRef dN/dpT distributions as function of centrality (i.e. per centrality bin).

The asymmetry A has the same binnage as the pT distributions, but i think it could be interesting
to vary the pT binning and see if the asymmetry is more pronounced in a specific pT region.
For example, we could check it per pT sector.
One possibility to study is to vary the bin width for A.

If there is a measurable asymmetry between the pT distributions, it should be only apparent in PbPb mid-central collisions.

That means A(ppRef) should be approximately 0 for all pT sectors, while A(PbPb) should have asymmetry in some pT sector,
thus producing a deviation ΔA = A(PbPb) - A(ppRef) non-zero in that pT sector.

The double ratio is also an option, even to plot on the same canvas, since it uses the same pT binning as the asymmetry A.

The other observable we have in hand is the event-by-event relative difference in pT between the two muons, namely

            ΔpT = (pT(mu+) - pT(mu-)) / (pT(mu+) + pT(mu-)).

This could be a valuable observable, but it needs to be treated with care.
A particularly natural observable would be like 

            <ΔpT(PbPb)> - <ΔpT(ppRef)>,

where the mean is calculated per centrality bin.

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
#include <tuple>
#include "../headers/basicFormatting.h"


//---Macro settings
std::string plot_extension = ".png"; // ".png" for regular development and ".pdf" for final quality plots
std::string BasePath = "/home/lucas/Documents/CMS/z_boson_analysis/"; //IFT
//std::string BasePath = "/home/lucasdoriac/z_boson_analysis/data/"; //Home


//Good selection threshold values
const double MAX_ZVTX = 15.0;

const double MINZ_MASS = 60.;
const double MAXZ_MASS = 120.;
const double RAPIDITYCUTVALUE = 2.4;

const float ETACUTVALUE = 2.4;
const double PTCUTVALUE = 20.;

// ##############################################################################
// ##############################################################################

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
        BasePath + "Data/PbPb2023/"
    },

    {
        "PbPb2024_Data",
        SampleType::Data,
        CollisionSystem::PbPb2024,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_PbPb2024Data.root",
        BasePath + "Data/PbPb2024/"
    },

    {
        "ppRef2024_Data",
        SampleType::Data,
        CollisionSystem::ppRef2024,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_ppRef2024.root",
        BasePath + "Data/ppRef2024/"
    },

    {
        "PbPb2024_MC",
        SampleType::MC,
        CollisionSystem::PbPb2024,
        "hionia/myTree",
        "Oniatree_PowhegZtoMuMu_PbPb2024_*.root",
        BasePath + "MC/PbPb2024/DYto2Mu_MLL-50_TuneCP5_5p36TeV_powheg-pythia8/PowhegEmbedded_March9/260309_143939/0000/"
    },

    {
        "ppRef2024_MC",
        SampleType::MC,
        CollisionSystem::ppRef2024,
        "hionia/myTree",
        "Oniatree_PowhegZtoMuMu_ppRef2024_*.root",
        BasePath + "MC/ppRef2024/DYToMuMu_M-50_TuneCP5_5p36TeV_powheg-pythia8/Powheg_ppRefPileup_March20/260320_125046/0000/"
    }
};

//Centrality bins for PbPb2024 data
std::vector<std::pair<int, int>> CentralitySet = {
    {0, 10},
    {10, 20},
    {20, 30},
    {30, 100}
};

//Vector with ppRef values. Zcount, Mean, MeanError, Variance. Only diff to the result vector is that on the latter we have cent bin string as first element.
std::vector<std::tuple<int, double, double, double>> ppRefResults; 

//Vector of 5-tuples with PbPb centrality bin values.
std::vector<std::tuple<std::string, int, double, double, double>> PbPbResults;

//Vector of 5-tuples with final results.
std::vector<std::tuple<std::string, int, double, double, double>> Results;


//---Function declarations
void CalculatePtRelativeDiff(const Dataset& dataset, std::pair<int, int> centralityBins);
void PlotPtRelativeDiff();

//---Main function
void PtRelDiffFromTree(){

    ppRefResults.clear();
    PbPbResults.clear();
    Results.clear();

    // The centrality argument is irrelevant for ppRef dataset. 
    CalculatePtRelativeDiff(datasets[2], {0, 100});

    for(const auto& centralityBin : CentralitySet){
        std::cout << "\nCalculating PtRelDiff_PbPb for centrality bin: " << centralityBin.first << "-" << centralityBin.second << "%\n";
        CalculatePtRelativeDiff(datasets[1], centralityBin); //PbPb2024 dataset
    }

    PlotPtRelativeDiff();

    //Print results to file
    std::ofstream outFile("PtRelDiffResults.dat");
    outFile << "# CentralityBin Z-count Mean MeanError Variance\n";
    outFile << std::setprecision(17);

    for (const auto& result : Results) {
        outFile << std::get<0>(result) << '\t'
                << std::get<1>(result) << '\t'
                << std::get<2>(result) << '\t'
                << std::get<3>(result) << '\t'
                << std::get<4>(result) << '\n';
    }

    outFile.close();
}

void PlotPtRelativeDiff(){

    //Before the actual plot we have to take the difference of PbPb values from ppRef values.
    for (size_t i = 0; i < PbPbResults.size(); ++i) { //Index of mean = 2 and meanerror = 3 in PbPb vector. i-1 for ppref vector.

        const auto& PbPbResult = PbPbResults[i];
        const auto& ppRefResult = ppRefResults[i];

        std::string centralityBinStr = std::get<0>(PbPbResult);
        int Zcount = std::get<1>(PbPbResult);
        double meanDiff = std::get<2>(PbPbResult) - std::get<1>(ppRefResult); //Mean difference
        double meanErrorDiff = std::sqrt(std::pow(std::get<3>(PbPbResult), 2) + std::pow(std::get<2>(ppRefResult), 2)); //Mean error difference
        double varianceDiff = std::get<4>(PbPbResult) - std::get<3>(ppRefResult); //Variance difference

        //Make 5-tuple
        std::tuple<std::string, int, double, double, double> resultTuple(centralityBinStr, Zcount, meanDiff, meanErrorDiff, varianceDiff);
        Results.push_back(resultTuple);
    }


    //Values for TGraphErrors are just mean diff and mean error diff so
    std::vector<double> xValues, yValues, xErrors, yErrors;
    int idx = 0;
    for (const auto& result : Results) {
        std::string centralityBinStr = std::get<0>(result);
        double meanDiff = std::get<2>(result);
        double meanErrorDiff = std::get<3>(result);
        float x_centrality = (CentralitySet[idx].second + CentralitySet[idx].first) / 2.0; // Midpoint of the centrality bin

        xValues.push_back(x_centrality);
        yValues.push_back(meanDiff);
        xErrors.push_back(0.5 * (CentralitySet[idx].second - CentralitySet[idx].first)); // Half-width of the centrality bin
        yErrors.push_back(meanErrorDiff);
        ++idx;
    }

    //TGraphErrors
    TGraphErrors* graph = new TGraphErrors(xValues.size(), xValues.data(), yValues.data(), xErrors.data(), yErrors.data());
    basicGraphFormatting(graph);

    //Plot
    TCanvas* canvas = new TCanvas("canvas", "Pt Relative Difference vs Centrality", 800, 600);
    graph->Draw("AP");

    canvas->SaveAs("output.png");

}

void CalculatePtRelativeDiff(const Dataset& dataset, std::pair<int, int> centralityBins){

    // Load root file.
    std::string fullPath = dataset.basePath + dataset.filePattern;

    TChain *chain = new TChain(dataset.treeName.c_str());
    chain->Add(fullPath.c_str());

    std::cout << "> Number of files = " << chain->GetListOfFiles()->GetEntries() << "\n" << std::endl;

    std::cout << "> Opening files " << fullPath << "\n" << std::endl;

    std::cout << "> Running function " << __func__ << " on " << dataset.name << "\n" << std::endl;
    
    //Total number of events on Tree.
    Long64_t nEvents = chain->GetEntries();

    const int MAX_DIMUON = 1000;
    const int MAX_MUON   = 1000;

    //For now, writing ONLY branches that are relevant to the observable we want to measure.

    //Event-level variables
    Int_t Centrality;
    Float_t  zVtx;

    chain->SetBranchAddress("zVtx", &zVtx);

    if(dataset.system == CollisionSystem::PbPb2024) {//Centrality is only defined for PbPb2024 dataset.
        chain->SetBranchAddress("Centrality", &Centrality);
    }

    //Dimuon-level variables
    Short_t Reco_Dimuon_size;

    Short_t Reco_Dimuon_sign[MAX_DIMUON];
    Short_t Reco_Dimuon_muonPlusIndex[MAX_DIMUON];
    Short_t Reco_Dimuon_muonMinusIndex[MAX_DIMUON];

    ULong64_t Reco_Dimuon_trig[MAX_DIMUON];
    Float_t Reco_Dimuon_vtxProb[MAX_DIMUON];

    std::vector<float>* Reco_Dimuon_pt = nullptr;
    std::vector<float>* Reco_Dimuon_eta = nullptr;
    std::vector<float>* Reco_Dimuon_rapidity = nullptr;
    std::vector<float>* Reco_Dimuon_phi = nullptr;
    std::vector<float>* Reco_Dimuon_invMass = nullptr;

    std::vector<float>* Reco_Dimuon_muonPtDiff = nullptr;
    std::vector<float>* Reco_Dimuon_muonPtRelDiff = nullptr;

    chain->SetBranchAddress("Reco_Dimuon_size", &Reco_Dimuon_size);

    chain->SetBranchAddress("Reco_Dimuon_sign", Reco_Dimuon_sign);
    chain->SetBranchAddress("Reco_Dimuon_muonPlusIndex", Reco_Dimuon_muonPlusIndex);
    chain->SetBranchAddress("Reco_Dimuon_muonMinusIndex", Reco_Dimuon_muonMinusIndex);

    chain->SetBranchAddress("Reco_Dimuon_trig", Reco_Dimuon_trig);
    chain->SetBranchAddress("Reco_Dimuon_vtxProb", Reco_Dimuon_vtxProb);

    chain->SetBranchAddress("Reco_Dimuon_pt", &Reco_Dimuon_pt);
    chain->SetBranchAddress("Reco_Dimuon_eta", &Reco_Dimuon_eta);
    chain->SetBranchAddress("Reco_Dimuon_rapidity", &Reco_Dimuon_rapidity);
    chain->SetBranchAddress("Reco_Dimuon_phi", &Reco_Dimuon_phi);
    chain->SetBranchAddress("Reco_Dimuon_invMass", &Reco_Dimuon_invMass);

    chain->SetBranchAddress("Reco_Dimuon_muonPtDiff", &Reco_Dimuon_muonPtDiff);
    chain->SetBranchAddress("Reco_Dimuon_muonPtRelDiff", &Reco_Dimuon_muonPtRelDiff);

    //Muon-level variables
    Short_t Reco_Muon_size;

    std::vector<float>* Reco_Muon_pt = nullptr;
    //std::vector<float>* Reco_Muon_ptErrTrk = nullptr; //I think we need to study this branch further.
    std::vector<float>* Reco_Muon_eta = nullptr;
    std::vector<float>* Reco_Muon_phi = nullptr;
    std::vector<float>* Reco_Muon_mass = nullptr;

    ULong64_t Reco_Muon_trig[MAX_MUON];
    Bool_t Reco_Muon_isTightCutBased[MAX_MUON];
    
    chain->SetBranchAddress("Reco_Muon_size", &Reco_Muon_size);

    chain->SetBranchAddress("Reco_Muon_pt", &Reco_Muon_pt);
    //chain->SetBranchAddress("Reco_Muon_ptErrTrk", &Reco_Muon_ptErrTrk);
    chain->SetBranchAddress("Reco_Muon_eta", &Reco_Muon_eta);
    chain->SetBranchAddress("Reco_Muon_phi", &Reco_Muon_phi);
    chain->SetBranchAddress("Reco_Muon_mass", &Reco_Muon_mass);

    chain->SetBranchAddress("Reco_Muon_trig", Reco_Muon_trig);
    chain->SetBranchAddress("Reco_Muon_isTightCutBased", Reco_Muon_isTightCutBased);


    //Event-level cut values.
    const float minCentrality = 2.*centralityBins.first;
    const float maxCentrality = 2.*centralityBins.second;
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

    //Trigger selection: 'L2SingleMu12'.
    ULong64_t triggerBit = 1ULL << 7;

    //Helpers to calculate mean of PtRelDiff event-by-event.
    Long64_t n = 0;
    double mean = 0;
    double M2 = 0;

    for(Long64_t i = 0; i < nEvents; ++i){//Loop through all EVENTS in the CHAIN.

        chain->GetEntry(i); //Get event i.

        //Good event selection
        bool goodVertex = (std::abs(zVtx) < maxZvtx);
        bool goodCent = true;

        if (dataset.system == CollisionSystem::PbPb2024) {
            goodCent = (Centrality >= minCentrality && Centrality < maxCentrality);
        }

        if (!goodVertex) continue;
        if (!goodCent) continue;

        for(Short_t j = 0; j < Reco_Dimuon_size; ++j){ //Loop through all reco dimuon candidates of event i.
            
            //Good Z selection
            bool goodMass = (Reco_Dimuon_invMass->at(j) > minZ_Mass && Reco_Dimuon_invMass->at(j) < maxZ_Mass);
            bool goodRapidity = (std::abs(Reco_Dimuon_rapidity->at(j)) < RapidityCutValue);
            bool goodCharge = (Reco_Dimuon_sign[j] == 0);
            bool goodVtxProb = (Reco_Dimuon_vtxProb[j] > 0.001); //Vertex probability cut of .1% for dimuon candidates.
            bool isTriggerMatched = true;

                if (dataset.system == CollisionSystem::PbPb2024){
                    isTriggerMatched = (Reco_Dimuon_trig[j] & triggerBit);//**At least** one of the daughter muons must be matched to the trigger.
                }

            if (!goodMass) continue;
            if (!goodRapidity) continue;
            if (!goodCharge) continue;
            if (!goodVtxProb) continue;
            if (!isTriggerMatched) continue;

            //Good muon selection
            Short_t muonPlusIndex = Reco_Dimuon_muonPlusIndex[j]; //Index of antimuon in the reco muon arrays.
            Short_t muonMinusIndex = Reco_Dimuon_muonMinusIndex[j]; //Index of corresponding muon in the reco muon arrays.

            ptplus = Reco_Muon_pt->at(muonPlusIndex); //pT of antimuon.
            ptminus = Reco_Muon_pt->at(muonMinusIndex); //pT of corresponding muon.
            etaplus = Reco_Muon_eta->at(muonPlusIndex); //Pseudorapidity of antimuon.
            etaminus = Reco_Muon_eta->at(muonMinusIndex); //Pseudorapidity of corresponding muon.
            MuPlIsTight = Reco_Muon_isTightCutBased[muonPlusIndex];
            MuMiIsTight = Reco_Muon_isTightCutBased[muonMinusIndex];
            
            bool goodMuPl = (ptplus > ptCutValue)
                            && (std::abs(etaplus) < EtaCutValue)
                            && (MuPlIsTight);

            bool goodMuMi = (ptminus > ptCutValue)
                            && (std::abs(etaminus) < EtaCutValue)
                            && (MuMiIsTight);

            if (!goodMuPl || !goodMuMi) continue;
            //End of good selection for dimuon candidate j of event i.

            //Cesar pointed out that a measure of PtRelDiff from a histogram may introduce a bias because of the binning.
            //Thus, we will calculate the PtRelDiff event-by-event, and calculate its mean without binning it.
            //We just need to study exactly how to best calculate the mean or another representative variable of the distribution.
            //Calculate mean, standard deviation and skewness.

            ++n;
            double x = Reco_Dimuon_muonPtRelDiff->at(j);
            double Delta = x - mean;
            mean += Delta / n;
            M2 += Delta * (x - mean);

        }//End of dimuon candidate loop.


        //Track progress of event loop.
        Long64_t progressStep = std::max<Long64_t>(1, nEvents / 100);
        if (i % progressStep == 0){

            int percent = static_cast<int>(100.0 * i / nEvents+0.5);

            std::cout << "\r"
                      << percent
                      << "% complete..."
                      << std::flush;
        }

    }//Exiting event-by-event loop.

    double variance = M2 / (n - 1);
    double meanError = std::sqrt(variance / n);

    if(dataset.system == CollisionSystem::PbPb2024){
        //Centrality to string
        std::string centralityBinStr = std::to_string(centralityBins.first) + "-" + std::to_string(centralityBins.second);

        //Make 4-tuple
        std::tuple<std::string, int, double, double, double> PbPbresultTuple(centralityBinStr, n, mean, meanError, variance);
        PbPbResults.push_back(PbPbresultTuple);
    }

    else {
        //Store ppRef results on its own vector.
        std::tuple<int, double, double, double> ppRefResultTuple(n, mean, meanError, variance);
        ppRefResults.push_back(ppRefResultTuple);
    }

}