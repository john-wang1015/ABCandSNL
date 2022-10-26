// The code to simulate the VCBM originates from Examining the efficacy of localised gemcitabine therapy 
// for the treatment of pancreatic cancer using a hybrid agent-based model
#include "Model.h"
#include <iostream>
#include <Vector>

extern "C" {
Pancreas* SeedAndGrowToStartVolume(double p0, double psc, int dmax, int gage, int page, double startVolume)
{
	Params* parameters = new Params(p0, psc, dmax, gage, page);
	vector<Cell*> empty;
	Pancreas* pancreas = new Pancreas(empty, parameters);
	// start with just one infected cancer cell nearest to (0, 0)
	pancreas->CreateInitialTumour();

	// pre-observation phase - run until tumour reaches start volume
	double volume = 0;
	int days = 0;
	while (volume < startVolume && days < 200)
		volume = pancreas->SimulateOneDay(days++);

	// who disposes parameters???
	return pancreas;
}


Pancreas* CreateNewParticle(double p0, double psc, int dmax, int gage, int page, Pancreas* pancreas)
{
	return pancreas->CreateNewParticle(new Params(p0, psc, dmax, gage, page));
}

void UpdateParticle(double p0, double psc, int dmax, int gage, int page, Pancreas* pancreas)
{
    pancreas->UpdateParameters(new Params(p0, psc, dmax, gage, page));
}

void TumourGrowthData(double* input, double p0, double psc, int dmax, int gage, int page, double startVolume, int days){
	Pancreas* pancreas = SeedAndGrowToStartVolume(p0, psc, dmax, gage, page,startVolume);
	
	for (int i = 0; i < days; i++)
	{
		pancreas->SimulateOneDay(1);
		input[i]=pancreas->TumourVolume();
		//output.push_back(pancreas->SimulateOneDay(1));
	}
	
}


void test(){
	std::cout << "Test";
}
}