#include "stochastic_model.hpp"
// #include "interface_model.hpp"
#include <iostream>

/*External wrappers for python integration */
// extern "C"
// {
//   void GetPhenotypesIDs(const char* file_path_c,const char* file_name_c, uint8_t n_genes, uint8_t colours);
// }

std::vector<Phenotype_ID> GetPhenotypeIDs(Genotype& genotype, uint8_t k_builds, PhenotypeTable* pt_it);
std::vector<Phenotype_ID> GetSetPIDs(Genotype& genotype, PhenotypeTable* pt_it);

void PreProcessSampled(std::vector<Genotype> genomes, Set_to_Genome& set_to_genome, PhenotypeTable* pt);
