#include "genotype_phenotype.hpp"
#include <iostream>

/*External wrappers for python integration */
extern "C"
{
  void GP_MapSampler(const char* file_path_c,uint8_t n_genes, uint8_t rcolours,uint8_t colours);
  void GP_MapSampler_new(const char* file_path_c, uint8_t n_genes, uint8_t rcolours, uint8_t colours);
}

/*GP map calculations*/
// std::vector<Phenotype_ID> GetPhenotypeIDs(Genotype& genotype, uint8_t k_builds, StochasticPhenotypeTable* pt_it);
std::vector<Genotype> genotype_neighbourhood(const Genotype& genome, uint8_t ngenes, uint8_t colours);
void JiggleGenotype(Genotype& genotype, uint8_t max_colour);

/*Neutral size calculations*/
uint64_t NeutralSize(Genotype genotype,uint32_t N_neutral_colours,uint32_t N_possible_interacting_colours);
uint64_t combination_with_repetiton(uint8_t space_size , uint8_t sample_size);
uint64_t nChoosek(uint8_t n, uint8_t k);

struct Shape_Metrics
{
  Phenotype_ID pID;
  uint32_t robustness;

  Shape_Metrics(Phenotype_ID pID);

  void robust_pID(std::vector <Phenotype_ID> pIDs);
};

typedef struct Shape_Metrics Shape_Metrics;

struct Genotype_Metrics
{
  uint8_t n_genes, colours;
  uint64_t neutral_weight;

  Genotype ref_genotype;
  std::vector <Phenotype_ID> ref_pIDs;

  Phenotype_ID death_pID = {0, 0}, loop_pID = {255, 0};

  double number_of_neighbours;
  double strict_robustness = 0, intersection_robustness = 0, union_evolvability = 0;
  double death = 0, loop = 0;

  std::vector <Shape_Metrics> shapes;
  std::set <Phenotype_ID> diversity;

  Genotype_Metrics(uint8_t ngenes, uint8_t colours);

  void set_reference(Genotype& genotype, std::vector <Phenotype_ID> pIDs);

  void clear();

  void analyse_pIDs(std::vector <Phenotype_ID>& pIDs);

  void save_to_file(std::ofstream& fout);
};

typedef struct Genotype_Metrics Genotype_Metrics;

struct Set_Metrics
{
  uint8_t n_genes, colours;

  std::vector <Phenotype_ID> ref_pIDs;

  std::vector <double> strict_robustnesses, union_evolvabilities;
  std::vector <double> intersection_robustnesses, deaths, loops;
  std::vector <uint64_t> neutral_weightings;
  std::vector <Genotype_Metrics> genome_metrics;

  std::set <Phenotype_ID> diversity;

  Set_Metrics(uint8_t n_genes, uint8_t colours);

  void add_genotype_metrics(Genotype_Metrics& gmetrics);

  void save_to_file(std::ofstream& set_out, std::ofstream& genome_out);

  void clear();
};

typedef struct Set_Metrics Set_Metrics;
