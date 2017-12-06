#include "interface_model.hpp"
#include <iostream>
#include <fstream>

namespace simulation_params
{
  typedef uint32_t population_size_type;
  extern population_size_type population_size;
  extern uint32_t generation_limit;
  extern uint8_t n_tiles,phenotype_builds;

}


std::vector<simulation_params::population_size_type> RouletteWheelSelection(std::vector<double>& fitnesses);

void EvolvePopulation(std::string run_details); 
RandomStrings()
void SetRuntimeConfigurations(int argc, char* argv[]);
