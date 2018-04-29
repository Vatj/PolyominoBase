#include "genotype_phenotype.hpp"
#include <iostream>

/*External wrappers for python integration */
extern "C"
{
  void GP_MapSampler(const char* file_path_c,uint8_t n_genes, uint8_t rcolours,uint8_t colours,bool file_of_genotypes);
}

/*GP map calculations*/
std::vector<Phenotype_ID> GetPhenotypeIDs(Genotype& genotype, uint8_t k_builds, StochasticPhenotypeTable* pt_it);
std::vector<Genotype> genotype_neighbourhood(const Genotype& genome, uint8_t ngenes, uint8_t colours);
void JiggleGenotype(Genotype& genotype, uint8_t max_colour);

/*Neutral size calculations*/
uint64_t NeutralSize(Genotype genotype,uint32_t N_neutral_colours,uint32_t N_possible_interacting_colours);
uint64_t combination_with_repetiton(uint8_t space_size , uint8_t sample_size);
uint64_t nChoosek(uint8_t n, uint8_t k);

struct Genotype_Metrics
{
  uint8_t n_genes;
  double neutral_size;

  Genotype ref_genotype;
  std::vector <Phenotype_ID> ref_pIDs;

  Phenotype_ID pID_InfiniteLoop{0, 255}, pID_NonDeterministic{0, 0};

  std::vector <double> robustness, new_evolvability;
  std::vector <double> death_InfiniteLoop, death_NonDeterministic, death;

  std::vector < std::set <Phenotype_ID> > diversity;

  Genotype_Metrics(uint8_t ngenes);

  void set_reference(Genotype& ref_genotype, std::vector <Phenotype_ID>& ref_pIDs);

  void clear();

  void analyse_pIDs(std::vector <Phenotype_ID>& pIDs);

  void save_to_file(std::ofstream& fout);
};

typedef struct Genotype_Metrics Genotype_Metrics;

struct Phenotype_Metrics
{
  uint8_t n_genes, colours;
  uint32_t N_JIGGLE=100;

  std::vector <Phenotype_ID> pIDs;

  std::vector <std::vector <double>> robustnesses;
  std::vector <std::vector <double>> new_evolvabilities;
  std::vector <std::vector <double>> deaths_InfiniteLoop;
  std::vector <std::vector <double>> deaths_NonDeterministic;
  std::vector <std::vector <double>> deaths;

  std::vector<uint64_t> neutral_weightings;

  std::vector <std::set <Phenotype_ID>> diversity;

  Phenotype_Metrics(uint8_t n_genes, uint8_t colours);

  void add_genotype_metrics(Genotype_Metrics& gmetrics);

  void save_to_file(std::ofstream& fout);

  void clear();

};

typedef struct Phenotype_Metrics Phenotype_Metrics;
