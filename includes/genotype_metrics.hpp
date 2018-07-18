#include "genotype_generate.hpp"
#include <iostream>

namespace simulation_params
{
  extern uint8_t n_genes, colours, metric_colours;
  extern uint32_t n_jiggle;
  extern std::mt19937 RNG_Engine;
  extern bool dup_aware;
}

namespace io_params
{
  extern std::string set_metric_file, genome_metric_file;
}

void GP_MapSampler(std::vector<Set_Metrics>& metrics, Set_to_Genome& set_to_genome, PhenotypeTable* pt);
void GP_MapSimple(std::vector<Set_Metrics>& metrics, Set_to_Genome& set_to_genome, PhenotypeTable* pt);

std::vector<Genotype> genotype_neighbourhood(const Genotype& genome);
void JiggleGenotype(Genotype& genotype);
void JiggleDuplicateGenotype(Genotype& genotype);

/*Neutral size calculations*/
uint64_t NeutralSize(Genotype genotype,uint32_t N_neutral_colours,uint32_t N_possible_interacting_colours);
uint64_t combination_with_repetiton(uint8_t space_size , uint8_t sample_size);
uint64_t nChoosek(uint8_t n, uint8_t k);
