#include "stochastic_model.hpp"
#include <iostream>

/*External wrappers for python integration */
extern "C"
{
  void GetPhenotypesIDs(const char* file_path_c,const char* file_name_c, uint8_t n_genes, uint8_t colours);
  // void PreProcessGenotypesTopology(const char* file_path_c, uint8_t n_genes, uint8_t colours);
  // void PreProcessGenotypes(const char* file_path_c, uint8_t n_genes, uint8_t colours);
  // void PreProcessWrite(const char* file_path_c, const char* file_name_c, uint8_t n_genes, uint8_t colours);
}

std::vector<Phenotype_ID> GetPhenotypeIDs(Genotype& genotype, uint8_t k_builds, StochasticPhenotypeTable* pt_it);
std::vector<Phenotype_ID> GetSetPIDs(Genotype& genotype, uint8_t k_builds, StochasticPhenotypeTable* pt_it);

void PreProcessSampled(std::vector<Genotype> genomes, Set_to_Genome& set_to_genome, StochasticPhenotypeTable* pt);
