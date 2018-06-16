#include "GP_generators.hpp"
#include <sstream>
#include <iterator>

extern "C"
{
void GP_MapSampler(const char* file_path_c,uint8_t n_genes, uint8_t rcolours,uint8_t colours,bool file_of_genotypes);
}

/*GP map calculations*/
std::vector<Phenotype_ID> GetPhenotypeIDs(Genotype& genotype, uint8_t k_builds, PhenotypeTable* pt_it);
std::vector<Genotype> genotype_neighbourhood(const Genotype& genome, uint8_t ngenes, uint8_t colours);
void JiggleGenotype(Genotype& genotype, uint8_t max_colour);


/*Neutral size calculations*/
uint64_t NeutralSize(Genotype genotype,uint32_t N_neutral_colours,uint32_t N_possible_interacting_colours);
uint64_t combination_with_repetiton(uint8_t space_size , uint8_t sample_size);
uint64_t nChoosek(uint8_t n, uint8_t k);
