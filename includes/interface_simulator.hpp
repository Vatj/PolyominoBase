#include "interface_model.hpp"
#include <iostream>
#include <fstream>




std::vector<simulation_params::population_size_type> RouletteWheelSelection(std::vector<double>& fitnesses);

void EvolvePopulation(std::string run_details); 
void RandomStrings();
void SetRuntimeConfigurations(int argc, char* argv[]);
