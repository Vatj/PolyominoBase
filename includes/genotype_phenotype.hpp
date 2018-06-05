#include "stochastic_model.hpp"
#include <iostream>

/*External wrappers for python integration */
extern "C"
{
  void GetPhenotypesIDs(const char* file_path_c,const char* file_name_c, uint8_t n_genes, uint8_t colours);
  void PreProcessGenotypesTopology(const char* file_path_c, uint8_t n_genes, uint8_t colours);
  void PreProcessGenotypes(const char* file_path_c, uint8_t n_genes, uint8_t colours);
}

std::vector<Phenotype_ID> GetPhenotypeIDs(Genotype& genotype, uint8_t k_builds, StochasticPhenotypeTable* pt_it);
std::vector<Phenotype_ID> GetSetPIDs(Genotype& genotype, uint8_t k_builds, StochasticPhenotypeTable* pt_it);

void LoadExistingTable(std::ifstream& fin, StochasticPhenotypeTable* pt_it);
void PreProcessSampled(std::ifstream& genome_in, Set_to_Genome& set_to_genome, StochasticPhenotypeTable* pt, std::ofstream& pIDset_out);
void PrintPreProcessFile(std::ofstream& fout, Set_to_Genome& set_to_genome);
void PrintSetTable(std::ofstream& fout, Set_to_Genome& set_to_genome);
