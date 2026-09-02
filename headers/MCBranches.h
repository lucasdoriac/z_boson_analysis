//MC branches
#ifndef MCBranches_h
#define MCBranches_h

#include <vector>

const int MAX_DIMUON = 1000;
const int MAX_MUON   = 1000;

//Event-level variables
UInt_t eventNb;
Float_t zVtx;
Short_t nPV;

Int_t Centrality;
Short_t NpixelTracks;

Short_t Npix;
Short_t Ntracks;

Int_t trigPrescale[11];
ULong64_t HLTriggers;

Float_t SumET_HF;
Float_t SumET_HFplus;
Float_t SumET_HFminus;
Float_t SumET_HFplusEta4;
Float_t SumET_HFminusEta4;
Float_t SumET_ET;
Float_t SumET_EE;
Float_t SumET_EB;

chain->SetBranchAddress("eventNb", &eventNb);
chain->SetBranchAddress("zVtx", &zVtx);
chain->SetBranchAddress("nPV", &nPV);
chain->SetBranchAddress("Centrality", &Centrality);
chain->SetBranchAddress("Npix", &Npix);
chain->SetBranchAddress("NpixelTracks", &NpixelTracks);
chain->SetBranchAddress("Ntracks", &Ntracks);
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

//Dimuon-level variables
Short_t Reco_QQ_size;

Short_t Reco_QQ_sign[MAX_QQ];
Short_t Reco_QQ_mupl_idx[MAX_QQ];
Short_t Reco_QQ_mumi_idx[MAX_QQ];

ULong64_t Reco_QQ_trig[MAX_QQ];

Float_t Reco_QQ_ctau[MAX_QQ];
Float_t Reco_QQ_ctauErr[MAX_QQ];
Float_t Reco_QQ_cosAlpha[MAX_QQ];

Float_t Reco_QQ_ctau3D[MAX_QQ];
Float_t Reco_QQ_ctauErr3D[MAX_QQ];
Float_t Reco_QQ_cosAlpha3D[MAX_QQ];

std::vector<float>* Reco_QQ_4mom_pt  = nullptr;
std::vector<float>* Reco_QQ_4mom_eta = nullptr;
std::vector<float>* Reco_QQ_4mom_y   = nullptr;
std::vector<float>* Reco_QQ_4mom_phi = nullptr;
std::vector<float>* Reco_QQ_4mom_m   = nullptr;

Float_t Reco_QQ_VtxProb[MAX_QQ];
Float_t Reco_QQ_MassErr[MAX_QQ];

Short_t Reco_QQ_whichGen[MAX_QQ];

TClonesArray* Reco_QQ_vtx = nullptr;

chain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
chain->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign);

chain->SetBranchAddress("Reco_QQ_4mom_pt", &Reco_QQ_4mom_pt);
chain->SetBranchAddress("Reco_QQ_4mom_eta", &Reco_QQ_4mom_eta);
chain->SetBranchAddress("Reco_QQ_4mom_y", &Reco_QQ_4mom_y);
chain->SetBranchAddress("Reco_QQ_4mom_phi", &Reco_QQ_4mom_phi);
chain->SetBranchAddress("Reco_QQ_4mom_m", &Reco_QQ_4mom_m);

chain->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
chain->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);
chain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig);

chain->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau);
chain->SetBranchAddress("Reco_QQ_ctauErr", Reco_QQ_ctauErr);
chain->SetBranchAddress("Reco_QQ_cosAlpha", Reco_QQ_cosAlpha);
chain->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D);
chain->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D);
chain->SetBranchAddress("Reco_QQ_cosAlpha3D", Reco_QQ_cosAlpha3D);

chain->SetBranchAddress("Reco_QQ_whichGen", Reco_QQ_whichGen);

chain->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb);
chain->SetBranchAddress("Reco_QQ_MassErr", Reco_QQ_MassErr);

chain->SetBranchAddress("Reco_QQ_vtx", &Reco_QQ_vtx);

//Muon-level variables
Short_t Reco_mu_size;

ULong64_t Reco_mu_trig[MAX_MU];

Short_t Reco_mu_whichGen[MAX_MU];
Short_t Reco_mu_charge[MAX_MU];

std::vector<float>* Reco_mu_4mom_pt  = nullptr;
std::vector<float>* Reco_mu_4mom_eta = nullptr;
std::vector<float>* Reco_mu_4mom_phi = nullptr;
std::vector<float>* Reco_mu_4mom_m   = nullptr;

std::vector<float>* Reco_mu_L1_4mom_pt  = nullptr;
std::vector<float>* Reco_mu_L1_4mom_eta = nullptr;
std::vector<float>* Reco_mu_L1_4mom_phi = nullptr;
std::vector<float>* Reco_mu_L1_4mom_m   = nullptr;

Bool_t Reco_mu_isPF[MAX_MU];
Bool_t Reco_mu_isTracker[MAX_MU];
Bool_t Reco_mu_isGlobal[MAX_MU];

Float_t Reco_mu_softMvaRun3Value[MAX_MU];

Bool_t Reco_mu_isMediumCutBased[MAX_MU];
Bool_t Reco_mu_isTightCutBased[MAX_MU];

Int_t Reco_mu_nPixValHits[MAX_MU];
Int_t Reco_mu_nMuValHits[MAX_MU];
Int_t Reco_mu_nTrkHits[MAX_MU];

Float_t Reco_mu_normChi2_inner[MAX_MU];

Int_t Reco_mu_nPixWMea[MAX_MU];
Int_t Reco_mu_nTrkWMea[MAX_MU];
Int_t Reco_mu_StationsMatched[MAX_MU];

Float_t Reco_mu_dxy[MAX_MU];
Float_t Reco_mu_dxyErr[MAX_MU];
Float_t Reco_mu_dz[MAX_MU];
Float_t Reco_mu_dzErr[MAX_MU];

chain->SetBranchAddress("Reco_mu_size", &Reco_mu_size);
chain->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen);
chain->SetBranchAddress("Reco_mu_charge", Reco_mu_charge);

chain->SetBranchAddress("Reco_mu_4mom_pt", &Reco_mu_4mom_pt);
chain->SetBranchAddress("Reco_mu_4mom_eta", &Reco_mu_4mom_eta);
chain->SetBranchAddress("Reco_mu_4mom_phi", &Reco_mu_4mom_phi);
chain->SetBranchAddress("Reco_mu_4mom_m", &Reco_mu_4mom_m);

chain->SetBranchAddress("Reco_mu_L1_4mom_pt", &Reco_mu_L1_4mom_pt);
chain->SetBranchAddress("Reco_mu_L1_4mom_eta", &Reco_mu_L1_4mom_eta);
chain->SetBranchAddress("Reco_mu_L1_4mom_phi", &Reco_mu_L1_4mom_phi);
chain->SetBranchAddress("Reco_mu_L1_4mom_m", &Reco_mu_L1_4mom_m);

chain->SetBranchAddress("Reco_mu_trig", Reco_mu_trig);

chain->SetBranchAddress("Reco_mu_isPF", Reco_mu_isPF);
chain->SetBranchAddress("Reco_mu_isTracker", Reco_mu_isTracker);
chain->SetBranchAddress("Reco_mu_isGlobal", Reco_mu_isGlobal);

chain->SetBranchAddress("Reco_mu_softMvaRun3Value", Reco_mu_softMvaRun3Value);
chain->SetBranchAddress("Reco_mu_isMediumCutBased", Reco_mu_isMediumCutBased);
chain->SetBranchAddress("Reco_mu_isTightCutBased", Reco_mu_isTightCutBased);

chain->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits);
chain->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits);
chain->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits);

chain->SetBranchAddress("Reco_mu_normChi2_inner", Reco_mu_normChi2_inner);

chain->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea);
chain->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea);
chain->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched);

chain->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy);
chain->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr);
chain->SetBranchAddress("Reco_mu_dz", Reco_mu_dz);
chain->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr);

//Generator Dimuon
Float_t Gen_weight;
Float_t Gen_pthat;

Short_t Gen_QQ_size;
Short_t Gen_QQ_mupl_idx[MAX_QQ];
Short_t Gen_QQ_mumi_idx[MAX_QQ];

std::vector<float>* Gen_QQ_4mom_pt  = nullptr;
std::vector<float>* Gen_QQ_4mom_eta = nullptr;
std::vector<float>* Gen_QQ_4mom_y   = nullptr;
std::vector<float>* Gen_QQ_4mom_phi = nullptr;
std::vector<float>* Gen_QQ_4mom_m   = nullptr;

Short_t Gen_QQ_whichRec[MAX_QQ];

chain->SetBranchAddress("Gen_weight", &Gen_weight);
chain->SetBranchAddress("Gen_pthat", &Gen_pthat);
chain->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size);

chain->SetBranchAddress("Gen_QQ_4mom_pt", &Gen_QQ_4mom_pt);
chain->SetBranchAddress("Gen_QQ_4mom_eta", &Gen_QQ_4mom_eta);
chain->SetBranchAddress("Gen_QQ_4mom_y", &Gen_QQ_4mom_y);
chain->SetBranchAddress("Gen_QQ_4mom_phi", &Gen_QQ_4mom_phi);
chain->SetBranchAddress("Gen_QQ_4mom_m", &Gen_QQ_4mom_m);

chain->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx);
chain->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx);
chain->SetBranchAddress("Gen_QQ_whichRec", Gen_QQ_whichRec);

//Generator muons
Short_t Gen_mu_size;
Short_t Gen_mu_charge[MAX_MU];

std::vector<float>* Gen_mu_4mom_pt  = nullptr;
std::vector<float>* Gen_mu_4mom_eta = nullptr;
std::vector<float>* Gen_mu_4mom_phi = nullptr;
std::vector<float>* Gen_mu_4mom_m   = nullptr;

Short_t Gen_mu_whichRec[MAX_MU];

chain->SetBranchAddress("Gen_mu_size", &Gen_mu_size);
chain->SetBranchAddress("Gen_mu_charge", Gen_mu_charge);

chain->SetBranchAddress("Gen_mu_4mom_pt", &Gen_mu_4mom_pt);
chain->SetBranchAddress("Gen_mu_4mom_eta", &Gen_mu_4mom_eta);
chain->SetBranchAddress("Gen_mu_4mom_phi", &Gen_mu_4mom_phi);
chain->SetBranchAddress("Gen_mu_4mom_m", &Gen_mu_4mom_m);

chain->SetBranchAddress("Gen_mu_whichRec", Gen_mu_whichRec);

#endif // ifndef MCBranches_h