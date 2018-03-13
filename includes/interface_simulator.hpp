#include "interface_model.hpp"
#include <functional>


struct PopulationGenotype {
  std::vector<interface_type> genotype;
  Phenotype_ID pid;
  //std::vector<uint8_t> interacting_interfaces;
  PopulationGenotype(void) : genotype(simulation_params::n_tiles*4), pid{0,0} {};
};



std::vector<uint8_t> SequenceDifference(const std::vector<interface_type>& parent, const std::vector<interface_type>& child);
uint8_t ConjugateInterface(std::vector<interface_type>& genotype,uint8_t mutation_site);



void RandomStrings();

/* Fitness selection */
std::vector<uint16_t> RouletteWheelSelection(std::vector<double>& fitnesses);

/* Main evolution runners */
void EvolvePopulation(std::string run_details); 
void EvolutionRunner();
void SetRuntimeConfigurations(int argc, char* argv[]);
