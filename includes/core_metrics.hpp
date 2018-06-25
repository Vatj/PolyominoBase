#include "core_phenotype.hpp"


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

  Genotype ref_genotype, original;
  std::vector <Phenotype_ID> ref_pIDs;
  std::map <Phenotype_ID, uint8_t> pID_counter;

  Phenotype_ID rare_pID = {0, 0}, unbound_pID = {255, 0};

  double number_of_neighbours;
  double strict_robustness = 0, intersection_robustness = 0, union_evolvability = 0;
  double rare = 0, unbound = 0, neutral_weight=0;

  std::vector <Shape_Metrics> shapes;
  std::set <Phenotype_ID> diversity;

  Genotype_Metrics(uint8_t ngenes, uint8_t colours);

  void set_reference(Genotype& genotype, std::vector <Phenotype_ID> pIDs, double neutral);

  void clear();

  void analyse_pIDs(std::vector <Phenotype_ID>& pIDs);

  void save_to_file(std::ofstream& fout);
};

typedef struct Genotype_Metrics Genotype_Metrics;

struct Set_Metrics
{
  uint8_t n_genes, colours;
  uint32_t analysed;

  std::vector <Phenotype_ID> ref_pIDs;
  std::vector <Genotype> originals;
  std::map <Genotype, std::vector<Phenotype_ID>> misclassified;

  std::vector <double> strict_robustnesses, union_evolvabilities;
  std::vector <double> intersection_robustnesses, rares, unbounds;
  std::vector <double> neutral_weightings;
  std::vector <Genotype_Metrics> genome_metrics;
  std::vector <uint32_t> diversity_tracker;

  std::set <Phenotype_ID> diversity;

  Set_Metrics(uint8_t n_genes, uint8_t colours);

  void add_genotype_metrics(Genotype_Metrics& gmetrics);

  void save_to_file(std::ofstream& set_out, std::ofstream& genome_out);

  void clear();
};

typedef struct Set_Metrics Set_Metrics;
