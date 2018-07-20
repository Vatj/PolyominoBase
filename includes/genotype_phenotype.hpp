#include "stochastic_model.hpp"
#include <iostream>

namespace simulation_params
{
  extern uint16_t phenotype_builds, preprocess_builds;
}

std::vector<Phenotype_ID> GetSetPIDs(Genotype genotype, PhenotypeTable* pt_it);
std::map<Phenotype_ID, uint8_t> GetPIDCounter(Genotype genotype, PhenotypeTable* pt_it);

void PreProcessSampled(std::vector<Genotype> genomes, Set_to_Genome& set_to_genome, PhenotypeTable* pt);
void FilterExhaustive(std::vector<Genotype> genomes, PhenotypeTable* pt);
