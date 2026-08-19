void print_pT_range(){

	TFile *rootFile = TFile::Open("Oniatree_PbPb2024_PromptReco.root", "READ");
    TDirectoryFile *dir1 = dynamic_cast<TDirectoryFile*>(rootFile->Get("hionia"));
	TTree *dimuonTree = dynamic_cast<TTree*>(dir1->Get("myTree"));

	std::vector<float>* Reco_mu_4mom_pt  = nullptr;
	std::vector<float>* Reco_mu_4mom_m = nullptr;	
	std::vector<float>* Reco_mu_4mom_eta = nullptr;	
	dimuonTree->SetBranchAddress("Reco_mu_4mom_pt", &Reco_mu_4mom_pt);
	dimuonTree->SetBranchAddress("Reco_mu_4mom_eta", &Reco_mu_4mom_eta);
	dimuonTree->SetBranchAddress("Reco_mu_4mom_m", &Reco_mu_4mom_m);

	float min_pt  = std::numeric_limits<float>::max();
    float max_pt  = std::numeric_limits<float>::lowest();

    float min_eta = std::numeric_limits<float>::max();
    float max_eta = std::numeric_limits<float>::lowest();

    float min_m   = std::numeric_limits<float>::max();
    float max_m   = std::numeric_limits<float>::lowest();

	for(Long64_t i = 0; i < dimuonTree->GetEntries(); ++i){

		dimuonTree->GetEntry(i);

		// Loop over all muons in this event
        for(size_t j = 0; j < Reco_mu_4mom_pt->size(); ++j){

            float pt  = Reco_mu_4mom_pt->at(j);
            float eta = Reco_mu_4mom_eta->at(j);
            float m   = Reco_mu_4mom_m->at(j);

            if(pt  < min_pt)  min_pt  = pt;
            if(pt  > max_pt)  max_pt  = pt;

            if(eta < min_eta) min_eta = eta;
            if(eta > max_eta) max_eta = eta;

            if(m   < min_m)   min_m   = m;
            if(m   > max_m)   max_m   = m;
        }
	}

	// Print results
    std::cout << "pT range:  [" << min_pt  << ", " << max_pt  << "]\n";
    std::cout << "eta range: [" << min_eta << ", " << max_eta << "]\n";
    std::cout << "mass range:[" << min_m   << ", " << max_m   << "]\n";
}