#include "genotype_duplicate.hpp"
#include <iostream>

namespace model_params
{
  extern uint8_t n_genes, colours, metric_colours;
}

namespace simulation_params
{
  extern uint8_t phenotype_builds;
  extern uint32_t n_samples, n_jiggle;
  extern double UND_threshold;
  extern bool allow_duplicates;
}


void PrintConfigFile(std::string config_file);
void PrintGenomeFile(std::string genome_file, std::vector<Genotype>& genomes);
void PrintPreProcessFile(std::string preprocess_file, Set_to_Genome& set_to_genome);
void PrintSetTable(std::string set_file, Set_to_Genome& set_to_genome);
void PrintMetrics(std::string set_metric_file, std::string genome_metric_file, std::vector<Set_Metrics> metrics);

void LoadGenomeFile(std::string genome_file, std::vector<Genotype>& genomes);
void LoadPreProcessFile(std::string preprocess_file, Set_to_Genome& set_to_genome);
void LoadPhenotypeTable(std::string phenotype_file, PhenotypeTable* pt_it);
