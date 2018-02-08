#include "interface_model.hpp"
#include <fstream>
#include <functional>



std::vector<uint16_t> RouletteWheelSelection(std::vector<double>& fitnesses);

void EvolvePopulation(std::string run_details); 
void RandomStrings();
void SetRuntimeConfigurations(int argc, char* argv[]);
std::vector<uint8_t> SequenceDifference(const std::vector<interface_type>& parent, const std::vector<interface_type>& child);
