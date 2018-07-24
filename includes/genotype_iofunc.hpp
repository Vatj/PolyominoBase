#include "genotype_duplicate.hpp"
#include <iostream>

namespace simulation_params
{
  extern uint16_t n_genes, colours, metric_colours;
  extern uint16_t phenotype_builds;
  extern uint32_t n_samples, n_jiggle;
  extern double UND_threshold;
  extern bool allow_duplicates, iso, duplicate_exhaustive;
}

std::string full_filename(std::string file_name, bool extra, bool iso);
void infer_file_details();
void all_files_to_full_names();

// void PrintConfigFile(std::string config_file);
void PrintGenomeFile(std::string genome_file, std::vector<Genotype>& genomes);
void PrintPreProcessFile(std::string preprocess_file, Set_to_Genome& set_to_genome);
void PrintSetTable(std::string set_file, Set_to_Genome& set_to_genome);
void header_metric_files(std::ofstream& set_metric_out, std::ofstream& genome_metric_out);

void LoadGenomeFile(std::string genome_file, std::vector<Genotype>& genomes);
void LoadPreProcessFile(std::string preprocess_file, Set_to_Genome& set_to_genome);
void LoadPhenotypeTable(std::string phenotype_file, PhenotypeTable* pt_it);
