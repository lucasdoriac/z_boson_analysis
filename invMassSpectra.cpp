/*
The goal is to build n invMassSpectrum plots, where n = N of good selections.
Apply each selection individually until all selections are applied.
Keep track of Z yield at every applied basic selection.

For example: goodCent -> Plot invMassSpectrum & Z yield = N -> goodCent + goodVertex -> Plot invMassSpectrum & Z yield = N-a.

13/08 - Começando macro.
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

//---Macro settings

TString plot_extension = ".pdf";

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
        "/home/lucas/Documents/Analysis/Data/PbPb2023/"
    },

    {
        "PbPb2024_Data",
        SampleType::Data,
        CollisionSystem::PbPb2024,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_PbPb2024Data.root",
        "/home/lucas/Documents/Analysis/Data/PbPb2024/"
    },

    {
        "ppRef2024_Data",
        SampleType::Data,
        CollisionSystem::ppRef2024,
        "hionia/DimuonTree",
        "HighPtMuons_HLTL2SingleMu_ppRef2024.root",
        "/home/lucas/Documents/Analysis/Data/ppRef2024/"
    },

    {
        "PbPb2024_MC",
        SampleType::MC,
        CollisionSystem::PbPb2024,
        "hionia/myTree",
        "Oniatree_PowhegZtoMuMu_PbPb2024_*.root",
        "/home/lucas/Documents/Analysis/MC/PbPb2024/DYto2Mu_MLL-50_TuneCP5_5p36TeV_powheg-pythia8/PowhegEmbedded_March9/260309_143939/0000/"
    },

    {
        "ppRef2024_MC",
        SampleType::MC,
        CollisionSystem::ppRef2024,
        "hionia/myTree",
        "Oniatree_PowhegZtoMuMu_ppRef2024_*.root",
        "/home/lucas/Documents/Analysis/MC/ppRef2024/DYToMuMu_M-50_TuneCP5_5p36TeV_powheg-pythia8/Powheg_ppRefPileup_March20/260320_125046/0000/"
    }
};


//---Observables

//---Histogram formatting
void basicCanvasFormatting(TCanvas *c);
void basicHistFormatting(TH1D *hist);
void basicLegendFormatting(TLegend *leg);
void drawLatexText(TString latexText = "#bf{CMS}", double x = 0.15, double y = 0.93, double TextSize = 0.05);

//---Miscellaneous
void printTreeContents(const Dataset &dataset);

//---Still writing !!DONT USE!!
void Make_invMassSpectrum(const Dataset &dataset);

//---main()
void invMassSpectra(){

	gROOT->SetBatch(kTRUE);

    //Make_invMassSpectrum(datasets[1]);
}

//---Observables


//---Histogram formatting

void basicCanvasFormatting(TCanvas* c){
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.035);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    c->SetTickx(1);
    c->SetTicky(1);
    c->SetFillColor(0);
    c->SetFrameFillColor(0);
    c->SetFrameLineWidth(2);
}

void basicHistFormatting(TH1D* hist){
	hist->GetXaxis()->CenterTitle(false);
    hist->GetYaxis()->CenterTitle(false);
    hist->GetXaxis()->SetTitleOffset(.9);
    hist->GetYaxis()->SetTitleOffset(1.);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetTitleSize(0.055);
    hist->GetYaxis()->SetTitleSize(0.055);
    hist->GetXaxis()->SetLabelSize(0.042);
    hist->GetYaxis()->SetLabelSize(0.042);
    hist->SetTitle("");
	hist->SetStats(0);
	hist->SetMarkerStyle(20);
	//hist->SetMarkerSize(0.9);
}

void basicLegendFormatting(TLegend *leg){
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.042);
    leg->SetTextFont(42);
    leg->SetMargin(0.2);
    leg->SetEntrySeparation(0.04);
}

void drawLatexText(TString latexText, double x, double y, double TextSize){
// We can add any pT or eta selection. If no text is passed to the function the CMS Header will be drawn.
TLatex latex;
latex.SetNDC();              // For normalized coordinates
latex.SetTextSize(TextSize);
latex.SetTextFont(42);       // Helvetica
latex.SetTextAlign(11);      // Left-top aligned.
latex.DrawLatex(x, y, latexText);
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
}

//---Still writing !!DONT USE!!
void Make_invMassSpectrum(const Dataset &dataset){

    // Load root file.
    std::string fullPath = dataset.basePath + dataset.filePattern;

    TChain *chain = new TChain(dataset.treeName.c_str());
    chain->Add(fullPath.c_str());

    std::cout << "> Number of files = " << chain->GetListOfFiles()->GetEntries() << "\n" << std::endl;

    std::cout << "> Opening files " << fullPath << "\n" << std::endl;

    std::cout << "> Running function " << __func__ << " on " << dataset.name << "\n" << std::endl;
    
    //Total number of events on Tree.
    Long64_t nEvents = chain->GetEntries();

    //Event-level variables
    UInt_t eventNb;
    UInt_t runNb;
    UInt_t LS;

    Float_t zVtx;
    Short_t nPV;
    Int_t Centrality;
    Short_t NpixelTracks;

    Int_t trigPrescale[8];
    ULong64_t HLTriggers;

    Float_t SumET_HF;
    Float_t SumET_HFplus;
    Float_t SumET_HFminus;
    Float_t SumET_HFplusEta4;
    Float_t SumET_HFminusEta4;
    Float_t SumET_ET;
    Float_t SumET_EE;
    Float_t SumET_EB;

    Int_t nEP;

    Float_t rpAng[100];
    Float_t rpSin[100];
    Float_t rpCos[100];

    Float_t rpAng_origin[100];
    Float_t rpSin_origin[100];
    Float_t rpCos_origin[100];

    chain->SetBranchAddress("eventNb", &eventNb);
    chain->SetBranchAddress("runNb", &runNb);
    chain->SetBranchAddress("LS", &LS);

    chain->SetBranchAddress("zVtx", &zVtx);
    chain->SetBranchAddress("nPV", &nPV);
    chain->SetBranchAddress("Centrality", &Centrality);
    chain->SetBranchAddress("NpixelTracks", &NpixelTracks);

    chain->SetBranchAddress("trigPrescale", trigPrescale);
    chain->SetBranchAddress("HLTriggers", &HLTriggers);

    chain->SetBranchAddress("SumET_HF", &SumET_HF);
    chain->SetBranchAddress("SumET_HFplus", &SumET_HFplus);
    chain->SetBranchAddress("SumET_HFminus", &SumET_HFminus);
    chain->SetBranchAddress("SumET_HFplusEta4", &SumET_HFplusEta4);
    chain->SetBranchAddress("SumET_HFminusEta4", &SumET_HFminusEta4);
    chain->SetBranchAddress("SumET_ET", &SumET_ET);
    chain->SetBranchAddress("SumET_EE", &SumET_EE);
    chain->SetBranchAddress("SumET_EB", &SumET_EB);

    chain->SetBranchAddress("nEP", &nEP);

    chain->SetBranchAddress("rpAng", rpAng);
    chain->SetBranchAddress("rpSin", rpSin);
    chain->SetBranchAddress("rpCos", rpCos);

    chain->SetBranchAddress("rpAng_origin", rpAng_origin);
    chain->SetBranchAddress("rpSin_origin", rpSin_origin);
    chain->SetBranchAddress("rpCos_origin", rpCos_origin);
    
    //Reconstructed dimuons
    Short_t Reco_Dimuon_size;

    Short_t Reco_Dimuon_sign[1000];
    Short_t Reco_Dimuon_muonPlusIndex[1000];
    Short_t Reco_Dimuon_muonMinusIndex[1000];

    ULong64_t Reco_Dimuon_trig[1000];

    Float_t Reco_Dimuon_ctau[1000];
    Float_t Reco_Dimuon_ctauErr[1000];
    Float_t Reco_Dimuon_cosAlpha[1000];

    Float_t Reco_Dimuon_ctau3D[1000];
    Float_t Reco_Dimuon_ctauErr3D[1000];
    Float_t Reco_Dimuon_cosAlpha3D[1000];

    Float_t Reco_Dimuon_vtxProb[1000];
    Float_t Reco_Dimuon_dca[1000];
    Float_t Reco_Dimuon_MassErr[1000];

    std::vector<float>* Reco_Dimuon_pt = nullptr;
    std::vector<float>* Reco_Dimuon_eta = nullptr;
    std::vector<float>* Reco_Dimuon_rapidity = nullptr;
    std::vector<float>* Reco_Dimuon_phi = nullptr;
    std::vector<float>* Reco_Dimuon_invMass = nullptr;

    std::vector<float>* Reco_Dimuon_muonPtDiff = nullptr;
    std::vector<float>* Reco_Dimuon_muonPtRelDiff = nullptr;

    std::vector<float>* Reco_Dimuon_vtx_xpos = nullptr;
    std::vector<float>* Reco_Dimuon_vtx_ypos = nullptr;
    std::vector<float>* Reco_Dimuon_vtx_zpos = nullptr;

    chain->SetBranchAddress("Reco_Dimuon_size", &Reco_Dimuon_size);

    chain->SetBranchAddress("Reco_Dimuon_sign", Reco_Dimuon_sign);

    chain->SetBranchAddress("Reco_Dimuon_pt", &Reco_Dimuon_pt);
    chain->SetBranchAddress("Reco_Dimuon_eta", &Reco_Dimuon_eta);
    chain->SetBranchAddress("Reco_Dimuon_rapidity", &Reco_Dimuon_rapidity);
    chain->SetBranchAddress("Reco_Dimuon_phi", &Reco_Dimuon_phi);
    chain->SetBranchAddress("Reco_Dimuon_invMass", &Reco_Dimuon_invMass);

    chain->SetBranchAddress("Reco_Dimuon_muonPtDiff", &Reco_Dimuon_muonPtDiff);
    chain->SetBranchAddress("Reco_Dimuon_muonPtRelDiff", &Reco_Dimuon_muonPtRelDiff);

    chain->SetBranchAddress("Reco_Dimuon_muonPlusIndex", Reco_Dimuon_muonPlusIndex);
    chain->SetBranchAddress("Reco_Dimuon_muonMinusIndex", Reco_Dimuon_muonMinusIndex);

    chain->SetBranchAddress("Reco_Dimuon_trig", Reco_Dimuon_trig);

    chain->SetBranchAddress("Reco_Dimuon_ctau", Reco_Dimuon_ctau);
    chain->SetBranchAddress("Reco_Dimuon_ctauErr", Reco_Dimuon_ctauErr);
    chain->SetBranchAddress("Reco_Dimuon_cosAlpha", Reco_Dimuon_cosAlpha);

    chain->SetBranchAddress("Reco_Dimuon_ctau3D", Reco_Dimuon_ctau3D);
    chain->SetBranchAddress("Reco_Dimuon_ctauErr3D", Reco_Dimuon_ctauErr3D);
    chain->SetBranchAddress("Reco_Dimuon_cosAlpha3D", Reco_Dimuon_cosAlpha3D);

    chain->SetBranchAddress("Reco_Dimuon_vtxProb", Reco_Dimuon_vtxProb);
    chain->SetBranchAddress("Reco_Dimuon_dca", Reco_Dimuon_dca);
    //chain->SetBranchAddress("Reco_Dimuon_MassErr", Reco_Dimuon_MassErr);

    chain->SetBranchAddress("Reco_Dimuon_vtx_xpos", &Reco_Dimuon_vtx_xpos);
    chain->SetBranchAddress("Reco_Dimuon_vtx_ypos", &Reco_Dimuon_vtx_ypos);
    chain->SetBranchAddress("Reco_Dimuon_vtx_zpos", &Reco_Dimuon_vtx_zpos);

    //Reconstructed muons
    Short_t Reco_Muon_size;

    std::vector<float>* Reco_Muon_pt = nullptr;
    std::vector<float>* Reco_Muon_ptErrTrk = nullptr;
    std::vector<float>* Reco_Muon_eta = nullptr;
    std::vector<float>* Reco_Muon_phi = nullptr;
    std::vector<float>* Reco_Muon_mass = nullptr;

    std::vector<float>* Reco_Muon_etaL1 = nullptr;
    std::vector<float>* Reco_Muon_phiL1 = nullptr;

    ULong64_t Reco_Muon_trig[1000];

    Bool_t Reco_Muon_isPF[1000];
    Bool_t Reco_Muon_isTracker[1000];
    Bool_t Reco_Muon_isGlobal[1000];

    Bool_t Reco_Muon_isSoftCutBased[1000];
    Bool_t Reco_Muon_isHybridSoft[1000];
    Bool_t Reco_Muon_isLooseCutBased[1000];
    Bool_t Reco_Muon_isMediumCutBased[1000];
    Bool_t Reco_Muon_isTightCutBased[1000];

    Float_t Reco_Muon_softMVAValue[1000];
    Float_t Reco_Muon_muonMVAValue[1000];

    std::vector<bool>* Reco_Muon_passesPFIsoLoose = nullptr;
    std::vector<bool>* Reco_Muon_passesPFIsoMedium = nullptr;
    std::vector<bool>* Reco_Muon_passesPFIsoTight = nullptr;
    std::vector<bool>* Reco_Muon_passesPFIsoVeryTight = nullptr;

    std::vector<float>* Reco_Muon_isoTrackSumPt = nullptr;

    std::vector<bool>* Reco_Muon_passesMultiIsoMedium = nullptr;

    std::vector<float>* Reco_Muon_HIMVAIso = nullptr;

    std::vector<bool>* Reco_Muon_HIMVAIsoWP80 = nullptr;
    std::vector<bool>* Reco_Muon_HIMVAIsoWP85 = nullptr;
    std::vector<bool>* Reco_Muon_HIMVAIsoWP90 = nullptr;
    std::vector<bool>* Reco_Muon_HIMVAIsoWP95 = nullptr;

    chain->SetBranchAddress("Reco_Muon_size", &Reco_Muon_size);

    chain->SetBranchAddress("Reco_Muon_pt", &Reco_Muon_pt);
    chain->SetBranchAddress("Reco_Muon_ptErrTrk", &Reco_Muon_ptErrTrk);
    chain->SetBranchAddress("Reco_Muon_eta", &Reco_Muon_eta);
    chain->SetBranchAddress("Reco_Muon_phi", &Reco_Muon_phi);
    chain->SetBranchAddress("Reco_Muon_mass", &Reco_Muon_mass);

    chain->SetBranchAddress("Reco_Muon_etaL1", &Reco_Muon_etaL1);
    chain->SetBranchAddress("Reco_Muon_phiL1", &Reco_Muon_phiL1);

    chain->SetBranchAddress("Reco_Muon_trig", Reco_Muon_trig);

    chain->SetBranchAddress("Reco_Muon_isPF", Reco_Muon_isPF);
    chain->SetBranchAddress("Reco_Muon_isTracker", Reco_Muon_isTracker);
    chain->SetBranchAddress("Reco_Muon_isGlobal", Reco_Muon_isGlobal);

    chain->SetBranchAddress("Reco_Muon_isSoftCutBased", Reco_Muon_isSoftCutBased);
    chain->SetBranchAddress("Reco_Muon_isHybridSoft", Reco_Muon_isHybridSoft);
    chain->SetBranchAddress("Reco_Muon_isLooseCutBased", Reco_Muon_isLooseCutBased);
    chain->SetBranchAddress("Reco_Muon_isMediumCutBased", Reco_Muon_isMediumCutBased);
    chain->SetBranchAddress("Reco_Muon_isTightCutBased", Reco_Muon_isTightCutBased);

    chain->SetBranchAddress("Reco_Muon_softMVAValue", Reco_Muon_softMVAValue);
    chain->SetBranchAddress("Reco_Muon_muonMVAValue", Reco_Muon_muonMVAValue);

    chain->SetBranchAddress("Reco_Muon_passesPFIsoLoose", &Reco_Muon_passesPFIsoLoose);
    chain->SetBranchAddress("Reco_Muon_passesPFIsoMedium", &Reco_Muon_passesPFIsoMedium);
    chain->SetBranchAddress("Reco_Muon_passesPFIsoTight", &Reco_Muon_passesPFIsoTight);
    chain->SetBranchAddress("Reco_Muon_passesPFIsoVeryTight", &Reco_Muon_passesPFIsoVeryTight);

    chain->SetBranchAddress("Reco_Muon_isoTrackSumPt", &Reco_Muon_isoTrackSumPt);

    chain->SetBranchAddress("Reco_Muon_passesMultiIsoMedium", &Reco_Muon_passesMultiIsoMedium);

    chain->SetBranchAddress("Reco_Muon_HIMVAIso", &Reco_Muon_HIMVAIso);

    chain->SetBranchAddress("Reco_Muon_HIMVAIsoWP80", &Reco_Muon_HIMVAIsoWP80);
    chain->SetBranchAddress("Reco_Muon_HIMVAIsoWP85", &Reco_Muon_HIMVAIsoWP85);
    chain->SetBranchAddress("Reco_Muon_HIMVAIsoWP90", &Reco_Muon_HIMVAIsoWP90);
    chain->SetBranchAddress("Reco_Muon_HIMVAIsoWP95", &Reco_Muon_HIMVAIsoWP95);

    //Histograms
    TH1D *h_invMass = new TH1D("hist_invMass", "Z0 invMassSpectrum", 40, 70, 110);
    
    //Event-level cut values.
    const int minCentrality = 0.0;
    const int maxCentrality = 180; //Cent max = maxCentrality/2.
    const double maxZvtx = 15.0;

    //Z-level cut values.
    const double minZ_Mass = 80.0; 
    const double maxZ_Mass = 100.0; 
    const double RapidityCutValue = 2.4;

    //Muon-level cut values.
    const float EtaCutValue = 2.4;
    const double ptCutValue = 10.0;
    bool MuPlIsTight;
    bool MuMiIsTight;

    //Observables
    TLorentzVector Z;
    double ptplus, ptminus;
    double etaplus, etaminus;
    //

    for(Long64_t i = 0; i < nEvents; ++i){ //Loop through all EVENTS in the CHAIN.

        chain->GetEntry(i); //Get event i.

        //Good event selection
        //bool goodCent = (Centrality > minCentrality && Centrality < maxCentrality);
        //bool goodVertex = (std::abs(zVtx) < maxZvtx);
        //Still missing 2 HF towers of 4 GeV of energy in the event.
        //Shapes of clusters compatibility.

        //if (!goodCent) continue;
        //if (!goodVertex) continue;

        for(Short_t j = 0; j < Reco_Dimuon_size; ++j){ //Loop through all reco dimuon candidates of event i.
            
            //Good Z selection
            //bool goodMass = (Reco_Dimuon_invMass->at(j) > minZ_Mass && Reco_Dimuon_invMass->at(j) < maxZ_Mass);
            //bool goodRapidity = (std::abs(Reco_Dimuon_rapidity->at(j)) < RapidityCutValue);

            //if (!goodMass) continue;
            //if (!goodRapidity) continue;

            //Good muon selection
            //ptplus = Reco_Muon_pt->at(Reco_Dimuon_muonPlusIndex[j]); //pT of antimuon.
            //ptminus = Reco_Muon_pt->at(Reco_Dimuon_muonMinusIndex[j]); //pT of corresponding muon.
            //etaplus = Reco_Muon_eta->at(Reco_Dimuon_muonPlusIndex[j]); //Pseudorapidity of antimuon.
            //etaminus = Reco_Muon_eta->at(Reco_Dimuon_muonMinusIndex[j]); //Pseudorapidity of corresponding muon.
            //MuPlIsTight = Reco_Muon_isTightCutBased[Reco_Dimuon_muonPlusIndex[j]];
            //MuMiIsTight = Reco_Muon_isTightCutBased[Reco_Dimuon_muonMinusIndex[j]];
            
            //Still missing trigger matching selection. "At least one muon matched to L1 and HLT trigger".
            //L1 matching, HLT matching...
            //Reco_mu_trig[Reco_mu_size]/I.
            //Muon ID selection (isGlobal, isTracker, isTight...)
            //Track quality such as Number of Track Hits NTrkHits and Chi^2/ndf.
            //Vertex quality with dxy, dz.
            //Muon isolation. Maybe not needed...?
            //One thing we need to verify is the possibility of the same muon being counted in different Z^0 candidates.
            //e.g. 
            //     Muon A + Muon B → candidate 1
            //     Muon A + Muon C → candidate 2

            //bool goodMuPl = (ptplus > ptCutValue) 
            //                && (std::abs(etaplus) < EtaCutValue) 
            //                && (MuPlIsTight);

            //bool goodMuMi = (ptminus > ptCutValue) 
            //                && (std::abs(etaminus) < EtaCutValue) 
            //                && (MuMiIsTight);

            //if (!goodMuPl || !goodMuMi) continue;

            h_invMass->Fill(Reco_Dimuon_invMass->at(j));
        }

        //Just to check the progress.
        if (i % (nEvents / 100) == 0){

            int percent = static_cast<int>(100.0 * i / nEvents+0.5);

            std::cout << "\r"
                      << percent
                      << "% complete..."
                      << std::flush;
        }//
    
    }//Exiting event-by-event loop.

    //End of histogram filling.
    std::cout << "\r100% complete!" << std::endl;

    //Starting histogram statistics calculation.

    TCanvas *myCanvas = new TCanvas("c_ptplmi", "", 700, 600);
    basicCanvasFormatting(myCanvas);
    basicHistFormatting(hist_ptpl);
    basicHistFormatting(hist_ptmi);

    hist_ptpl->SetLineColor(kRed + 1);
    hist_ptpl->SetLineWidth(2);
    hist_ptmi->SetLineColor(kBlue + 1);
    hist_ptmi->SetLineWidth(2);
    hist_ptpl->GetXaxis()->SetTitle("p_{T} [GeV]");
    hist_ptpl->GetYaxis()->SetTitle("Number of Z^{0} candidates [GeV^{-1}]");

    //Zooming in for display only.
    hist_ptpl->GetXaxis()->SetRangeUser(10, 100);
    hist_ptpl->GetYaxis()->SetRangeUser(0, hist_ptpl->GetMaximum() + 60);
    hist_ptpl->Draw("HIST");
    hist_ptmi->Draw("HIST SAME");

    linePl->Draw("SAME");
    lineMi->Draw("SAME");

    drawLatexText();
    drawLatexText("#it{Internal}", 0.24, 0.93, 0.042);
    drawLatexText("PbPb 2024 (5.36 TeV)", 0.65, 0.93, 0.042);        
    drawLatexText("Cent. 30-60%", 0.65, 0.35, 0.044);

    TLegend* leg = new TLegend(0.5, 0.60, 0.72, 0.85);
    basicLegendFormatting(leg);
    leg->AddEntry(hist_ptpl, Form("#LT p_{T}(#mu^{+}) #GT = (%.2f #pm %.2f) GeV", mean_pl, mean_pl_err),"f");
    leg->AddEntry(hist_ptmi, Form("#LT p_{T}(#mu^{-}) #GT = (%.2f #pm %.2f) GeV", mean_mi, mean_mi_err),"f");
    leg->AddEntry(hist_ptpl, Form("p_{T}(#mu^{+}) peak = (%.1f #pm %.1f) GeV", xMaxPl, hist_ptpl->GetBinWidth(binMaxPl)/2.), "l");
    leg->AddEntry(hist_ptmi, Form("p_{T}(#mu^{-}) peak = (%.1f #pm %.1f) GeV", xMaxMi, hist_ptmi->GetBinWidth(binMaxMi)/2.), "l");
    leg->SetTextSize(0.038);
    leg->SetEntrySeparation(0.038);
    leg->Draw();

    TString output = TString(dataset.name) + TString(__func__) + plot_extension;

    myCanvas->Update();
    myCanvas->SaveAs(output);

    delete chain;
    delete myCanvas;
}

//---