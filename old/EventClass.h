#ifndef EVENTCLASS_H
#define EVENTCLASS_H

#include <iostream>
#include <TChain.h>

// class definition here
class EventClass{

public:

    Float_t zVtx;

	EventClass(TChain *chain){
	    ConnectBranches(chain);
	}

private:
	void ConnectBranches(TChain *chain){
    	chain->SetBranchAddress("zVtx", &zVtx);

    	std::cout << "Branch zVtx connected succesfully." << std::endl;
	}
};


class MyClass{

public:

	std::string name;
	float number;

	//Constructor definition.
	MyClass(std::string n, float x){
		name = n;
		number = x;			
		std::cout << "Hello im " << name << " and im returning the value " << squaredNumber() << std::endl;
	}

	float squaredNumber(){
		return pow(number,2);
	}
//	private:
};

#endif // EVENTCLASS_H