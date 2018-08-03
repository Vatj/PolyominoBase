#include "genotype_metrics.hpp"

namespace io_params
{
  extern std::string file_path, config_file;
  extern std::string out_genome_file, in_genome_file, duplicate_file;
  extern std::string out_phenotype_file, in_phenotype_file;
  extern std::string set_file, preprocess_file;
  extern std::string set_metric_file, genome_metric_file;
}

int main(int argc, char *argv[]);

void JustExhaustive();
void QuickRandom();
void QuickFromFile();
void ExhaustiveMetricsPrintAll();
// void DuplicateJiggle();
void DuplicateExhaustive();
void QuickMapFromFile();
void MakePhenotypeTableFromFile();
