#include "interface_model.hpp"

namespace simulation_params
{
  typedef uint16_t population_size_type;
  extern population_size_type population_size;

}


std::vector<uint16_t> RouletteWheelSelection(std::vector<double>& fitnesses);
void EvolvePopulation();
