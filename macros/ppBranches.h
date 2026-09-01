//ppRef branches.
#ifndef ppBranches_h
#define ppBranches_h

#include <vector>

const int MAX_DIMUON = 1000;
const int MAX_MUON   = 1000;

//Event-level variables
UInt_t eventNb;
UInt_t runNb;
UInt_t LS;

Float_t zVtx;
Short_t nPV;

Int_t trigPrescale[9];
ULong64_t HLTriggers;

//Int_t nEP;

//Float_t rpAng[100];
//Float_t rpSin[100];
//Float_t rpCos[100];

//Float_t rpAng_origin[100];
//Float_t rpSin_origin[100];
//Float_t rpCos_origin[100];

chain->SetBranchAddress("eventNb", &eventNb);
chain->SetBranchAddress("runNb", &runNb);
chain->SetBranchAddress("LS", &LS);

chain->SetBranchAddress("zVtx", &zVtx);
chain->SetBranchAddress("nPV", &nPV);
//chain->SetBranchAddress("NpixelTracks", &NpixelTracks);

chain->SetBranchAddress("trigPrescale", trigPrescale);
chain->SetBranchAddress("HLTriggers", &HLTriggers);

//chain->SetBranchAddress("nEP", &nEP);

//chain->SetBranchAddress("rpAng", rpAng);
//chain->SetBranchAddress("rpSin", rpSin);
//chain->SetBranchAddress("rpCos", rpCos);

//chain->SetBranchAddress("rpAng_origin", rpAng_origin);
//chain->SetBranchAddress("rpSin_origin", rpSin_origin);
//chain->SetBranchAddress("rpCos_origin", rpCos_origin);

//Dimuon-level variables
Short_t Reco_Dimuon_size;

Short_t Reco_Dimuon_sign[MAX_DIMUON];
Short_t Reco_Dimuon_muonPlusIndex[MAX_DIMUON];
Short_t Reco_Dimuon_muonMinusIndex[MAX_DIMUON];

ULong64_t Reco_Dimuon_trig[MAX_DIMUON];

Float_t Reco_Dimuon_ctau[MAX_DIMUON];
Float_t Reco_Dimuon_ctauErr[MAX_DIMUON];
Float_t Reco_Dimuon_cosAlpha[MAX_DIMUON];

Float_t Reco_Dimuon_ctau3D[MAX_DIMUON];
Float_t Reco_Dimuon_ctauErr3D[MAX_DIMUON];
Float_t Reco_Dimuon_cosAlpha3D[MAX_DIMUON];

Float_t Reco_Dimuon_vtxProb[MAX_DIMUON];
Float_t Reco_Dimuon_dca[MAX_DIMUON];
Float_t Reco_Dimuon_invMassErr[MAX_DIMUON];

std::vector<float>* Reco_Dimuon_pt       = nullptr;
std::vector<float>* Reco_Dimuon_eta      = nullptr;
std::vector<float>* Reco_Dimuon_rapidity = nullptr;
std::vector<float>* Reco_Dimuon_phi      = nullptr;
std::vector<float>* Reco_Dimuon_invMass  = nullptr;

std::vector<float>* Reco_Dimuon_muonPtDiff    = nullptr;
std::vector<float>* Reco_Dimuon_muonPtRelDiff = nullptr;

std::vector<float>* Reco_Dimuon_vtx_xpos = nullptr;
std::vector<float>* Reco_Dimuon_vtx_ypos = nullptr;
std::vector<float>* Reco_Dimuon_vtx_zpos = nullptr;

chain->SetBranchAddress("Reco_Dimuon_size", &Reco_Dimuon_size);

chain->SetBranchAddress("Reco_Dimuon_sign", Reco_Dimuon_sign);
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

chain->SetBranchAddress("Reco_Dimuon_pt", &Reco_Dimuon_pt);
chain->SetBranchAddress("Reco_Dimuon_eta", &Reco_Dimuon_eta);
chain->SetBranchAddress("Reco_Dimuon_rapidity", &Reco_Dimuon_rapidity);
chain->SetBranchAddress("Reco_Dimuon_phi", &Reco_Dimuon_phi);
chain->SetBranchAddress("Reco_Dimuon_invMass", &Reco_Dimuon_invMass);

chain->SetBranchAddress("Reco_Dimuon_muonPtDiff", &Reco_Dimuon_muonPtDiff);
chain->SetBranchAddress("Reco_Dimuon_muonPtRelDiff", &Reco_Dimuon_muonPtRelDiff);

chain->SetBranchAddress("Reco_Dimuon_vtx_xpos", &Reco_Dimuon_vtx_xpos);
chain->SetBranchAddress("Reco_Dimuon_vtx_ypos", &Reco_Dimuon_vtx_ypos);
chain->SetBranchAddress("Reco_Dimuon_vtx_zpos", &Reco_Dimuon_vtx_zpos);

//Muon-level variables
Short_t Reco_Muon_size;

std::vector<float>* Reco_Muon_pt       = nullptr;
std::vector<float>* Reco_Muon_ptErrTrk = nullptr;
std::vector<float>* Reco_Muon_eta      = nullptr;
std::vector<float>* Reco_Muon_phi      = nullptr;
std::vector<float>* Reco_Muon_mass     = nullptr;

std::vector<float>* Reco_Muon_etaL1 = nullptr;
std::vector<float>* Reco_Muon_phiL1 = nullptr;

ULong64_t Reco_Muon_trig[MAX_MUON];

Bool_t Reco_Muon_isPF[MAX_MUON];
Bool_t Reco_Muon_isTracker[MAX_MUON];
Bool_t Reco_Muon_isGlobal[MAX_MUON];

Bool_t Reco_Muon_isSoftCutBased[MAX_MUON];
Bool_t Reco_Muon_isHybridSoft[MAX_MUON];
Bool_t Reco_Muon_isLooseCutBased[MAX_MUON];
Bool_t Reco_Muon_isMediumCutBased[MAX_MUON];
Bool_t Reco_Muon_isTightCutBased[MAX_MUON];

Float_t Reco_Muon_softMVAValue[MAX_MUON];
Float_t Reco_Muon_muonMVAValue[MAX_MUON];

std::vector<bool>* Reco_Muon_passesPFIsoLoose     = nullptr;
std::vector<bool>* Reco_Muon_passesPFIsoMedium    = nullptr;
std::vector<bool>* Reco_Muon_passesPFIsoTight     = nullptr;
std::vector<bool>* Reco_Muon_passesPFIsoVeryTight = nullptr;

std::vector<float>* Reco_Muon_isoTrackSumPt = nullptr;

std::vector<bool>* Reco_Muon_passesMultiIsoMedium = nullptr;

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

#endif // ifndef ppBranches_h