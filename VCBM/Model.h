#pragma once

#include "Pancreas.hpp"

extern "C" {
Pancreas* SeedAndGrowToStartVolume(double p0, double psc, int dmax, int gage, int page, double startVolume);
Pancreas* CreateNewParticle(double p0, double psc, int dmax, int gage, int page, Pancreas* pancreas);
void UpdateParticle(double p0, double psc, int dmax, int gage, int page, Pancreas* pancreas);
void TumourGrowthData(double* input, double p0, double psc, int dmax, int gage, int page, double startVolume, int days);
}