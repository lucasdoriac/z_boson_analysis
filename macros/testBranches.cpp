/*
I want to test a new way to declare and initialize the branch variables.
*/
#include "PbPbBranches.h"
#include "TheDatasets.h"

using namespace std;

//---Macro settings
std::string plot_extension = ".png"; // ".png" for regular development and ".pdf" for final quality plots

void testBranches(){

    Dataset dataset = datasets[1]; //PbPb2024_Data

    //Load root file.
    std::string fullPath = dataset.basePath + dataset.filePattern;

    TChain *chain = new TChain(dataset.treeName.c_str());
    chain->Add(fullPath.c_str());

    std::cout << "> Number of files = " << chain->GetListOfFiles()->GetEntries() << "\n" << std::endl;

    std::cout << "> Opening files " << fullPath << "\n" << std::endl;

    std::cout << "> Running function " << __func__ << " on " << dataset.name << "\n" << std::endl;
    
    //Total number of events on Tree.
    Long64_t nEvents = chain->GetEntries();    

    SetPbPbBranches PbPb2024;
    PbPb2024.SetBranches(chain);

    cout << "> Number of events in the tree = " << nEvents << endl;
    cout << "> Test succesfully completed." << endl;
}