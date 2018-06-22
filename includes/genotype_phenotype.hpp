#include "stochastic_model.hpp"
#include <iostream>

std::vector<Phenotype_ID> GetPhenotypeIDs(Genotype& genotype, uint8_t k_builds, PhenotypeTable* pt_it);
std::vector<Phenotype_ID> GetSetPIDs(Genotype genotype, PhenotypeTable* pt_it);

void PreProcessSampled(std::vector<Genotype> genomes, Set_to_Genome& set_to_genome, PhenotypeTable* pt);
void FilterExhaustive(std::vector<Genotype> genomes, PhenotypeTable* pt);
