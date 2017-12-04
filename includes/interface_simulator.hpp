#include "interface_model.hpp"
#include <iostream>
#include <fstream>

namespace simulation_params
{
  typedef uint16_t population_size_type;
  extern population_size_type population_size;
  extern uint32_t generation_limit;
  extern uint8_t n_tiles,phenotype_builds;

}


std::vector<uint16_t> RouletteWheelSelection(std::vector<double>& fitnesses);
void EvolvePopulation();
void EvolvePopulation_Off();

void SetRuntimeConfigurations(int argc, char* argv[]);
