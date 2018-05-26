#include "stochastic_model.hpp"
#include <iostream>

/*External wrappers for python integration */
extern "C"
{
  void GetPhenotypesIDs(const char* file_path_c,const char* file_name_c, uint8_t n_genes, uint8_t colours);
}

std::vector<Phenotype_ID> GetPhenotypeIDs(Genotype& genotype, uint8_t k_builds, StochasticPhenotypeTable* pt_it);
std::vector<Phenotype_ID> GetSetPIDs(Genotype& genotype, uint8_t k_builds, StochasticPhenotypeTable* pt_it);
