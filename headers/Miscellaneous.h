#ifndef MISCELLANEOUS_H
#define MISCELLANEOUS_H

#include <iostream>
#include <TChain.h>
#include "TheDatasets.h"

//---Miscellaneous useful functions
void printTreeContents(const Dataset &dataset){

    std::string fullPath = dataset.path + dataset.filePattern;

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

#endif //MISCELLANEOUS_H