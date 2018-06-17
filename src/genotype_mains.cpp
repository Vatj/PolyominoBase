#include "genotype_mains.hpp"

namespace model_params
{
  uint8_t n_genes=3, colours=7, metric_colours=9;
}

namespace simulation_params
{
  uint8_t phenotype_builds=50;
  uint32_t n_samples = 10, n_jiggle = 10;
  std::mt19937 RNG_Engine(std::random_device{}());
  double UND_threshold=0.25;
  bool allow_duplicates = true;
}

namespace io_params
{
  std::string file_path = "/rscratch/vatj2/public_html/Polyominoes/data/gpmap/V4/exhaustive/threshold25p/";
  std::string file_details = "_N" + std::to_string(model_params::n_genes) + "_C" + std::to_string(model_params::colours);
  std::string extra="_Cx" + std::to_string(model_params::metric_colours), ending=".txt";

  std::string genome_file = file_path + "SampledGenotypes" + file_details + ending;
  std::string phenotype_file = file_path + "PhenotypeTable" + file_details + ending;
  std::string set_file = file_path + "SetTable" + file_details + ending;
  std::string preprocess_file = file_path + "PreProcessGenotypes" + file_details + ending;
  std::string set_metric_file = file_path + "SetMetrics" + file_details + extra + ending;
  std::string genome_metric_file = file_path + "GenomeMetrics" + file_details + extra + ending;
  std::string duplicate_file = file_path + "DuplicateGenotypes" + file_details + ending;

  std::string config_file = file_path + "Config" + file_details + extra + ending;
}

int main()
{
  std::cout << "Global path : " << io_params::file_path  << "\n";

  // ExhaustiveFast();
  ExhaustiveMetricsAndPrint();
  // Quick();
  // Duplicate();
  // FullFast();
  // QuickFromIso();

  std::cout << "Back to sleep!" << std::endl;

  PrintConfigFile(io_params::config_file);

  return 0;
}

// void Duplicate()
// {
//   PhenotypeTable pt;
//   std::vector<Genotype> genomes, duplicates;
//   Set_to_Genome set_to_genome;
//   std::vector<Set_Metrics> metrics;
//
//   LoadPhenotypeTable(io_params::phenotype_file, &pt);
//   LoadGenomeFile(io_params::genome_file, genomes);
//
//   duplicates = GenomesDuplication(genomes);
//
//   PrintGenomeFile(io_params::duplicate_file, duplicates);
//   // PreProcessSampled(duplicates, set_to_genome, &pt);
//   // PrintPreProcessFile(io_params::preprocess_file, set_to_genome);
//   // PrintSetTable(io_params::set_file, set_to_genome);
//
//   // GP_MapSampler(metrics, set_to_genome, &pt);
//   // PrintMetrics(io_params::set_metric_file, io_params::genome_metric_file, metrics);
//
//   // pt.PrintTable(io_params::phenotype_file);
//
//   std::cout << "Back to sleep!" << std::endl;
// }


void ExhaustiveFast()
{
  PhenotypeTable pt;
  std::vector<Genotype> genomes;
  Set_to_Genome set_to_genome;

  genomes = ExhaustiveMinimalGenotypes(&pt);
  // genomes = SampleMinimalGenotypes(n_genes, colours, N_SAMPLES, allow_duplicates, &pt);
  PrintGenomeFile(io_params::genome_file, genomes);
  pt.PrintTable(io_params::phenotype_file);

  // PreProcessSampled(genomes, set_to_genome, &pt);
  // PrintPreProcessFile(io_params::preprocess_file, set_to_genome);
  // PrintSetTable(io_params::set_file, set_to_genome);
}


// void FullFast()
// {
//   PhenotypeTable pt;
//   std::vector<Genotype> genomes;
//   Set_to_Genome set_to_genome;
//   std::vector<Set_Metrics> metrics;
//
//   genomes = ExhaustiveFullGenotypes2(colours, &pt);
//   // genomes = SampleMinimalGenotypes(n_genes, colours, N_SAMPLES, allow_duplicates, &pt);
//   // PrintGenomeFile(io_params::genome_file, genomes);
//   // pt.PrintTable(io_params::phenotype_file);
//
//   PreProcessSampled(genomes, set_to_genome, &pt);
//   // PrintPreProcessFile(io_params::preprocess_file, set_to_genome);
//   // PrintSetTable(io_params::set_file, set_to_genome);
//
//   GP_MapSampler(metrics, set_to_genome, &pt);
//   PrintMetrics(io_params::set_metric_file, io_params::genome_metric_file, metrics);
//
// }

void ExhaustiveMetricsAndPrint()
{
  PhenotypeTable pt;
  std::vector<Genotype> genomes;
  Set_to_Genome set_to_genome;
  std::vector<Set_Metrics> metrics;

  genomes = ExhaustiveMinimalGenotypes(&pt);
  // genomes = SampleMinimalGenotypes(model_params::n_genes, colours, N_SAMPLES, allow_duplicates, &pt);
  // LoadGenomeFile(genome_file, genomes);
  PrintGenomeFile(io_params::genome_file, genomes);
  pt.PrintTable(io_params::phenotype_file);

  PreProcessSampled(genomes, set_to_genome, &pt);
  PrintPreProcessFile(io_params::preprocess_file, set_to_genome);
  PrintSetTable(io_params::set_file, set_to_genome);

  GP_MapSampler(metrics, set_to_genome, &pt);
  PrintMetrics(io_params::set_metric_file, io_params::genome_metric_file, metrics);
}

void Quick()
{
  PhenotypeTable pt;
  std::vector<Genotype> genomes;
  Set_to_Genome set_to_genome;
  std::vector<Set_Metrics> metrics;

  genomes = SampleMinimalGenotypes(&pt);
  PrintGenomeFile(io_params::genome_file, genomes);
  PreProcessSampled(genomes, set_to_genome, &pt);
  GP_MapSampler(metrics, set_to_genome, &pt);
  PrintMetrics(io_params::set_metric_file, io_params::genome_metric_file, metrics);
}

void QuickFromIso()
{
  PhenotypeTable pt;
  std::vector<Genotype> genomes;
  Set_to_Genome set_to_genome;
  std::vector<Set_Metrics> metrics;

  LoadGenomeFile(io_params::genome_file, genomes);
  PreProcessSampled(genomes, set_to_genome, &pt);
  GP_MapSampler(metrics, set_to_genome, &pt);
  PrintMetrics(io_params::set_metric_file, io_params::genome_metric_file, metrics);
}


// int main()
// {
//   PhenotypeTable pt;
//   std::vector<Genotype> genomes;
//   Set_to_Genome set_to_genome;
//   std::vector<Set_Metrics> metrics;
//
//   bool print_genome = true, print_preprocess = true;
//   bool print_set_table = true, print_phenotype_table = true;
//   bool load_genome = false, load_preprocess = false, load_phenotype_table = false;
//
//   if(load_phenotype_table)
//     LoadPhenotypeTable(io_params::phenotype_file, &pt);
//
//   if(load_genome)
//     LoadGenomeFile(io_params::genome_file, genomes);
//   else
//     genomes = SampleMinimalGenotypes(n_genes, colours, N_SAMPLES, allow_duplicates, &pt);
//
//   if(print_genome)
//     PrintGenomeFile(io_params::genome_file, genomes);
//
//   // std::vector<Genotype> dups = GenomesDuplication(genomes, n_genes);
//   // PrintGenomeFile(io_params::duplicate_file, dups);
//
//   if(print_phenotype_table)
//     pt.PrintTable(io_params::phenotype_file);
//
//   if(load_preprocess)
//     LoadPreProcessFile(io_params::preprocess_file, set_to_genome);
//   else
//     PreProcessSampled(genomes, set_to_genome, &pt);
//
//   if(print_preprocess)
//     PrintPreProcessFile(io_params::preprocess_file, set_to_genome);
//
//   if(print_set_table)
//     PrintSetTable(io_params::set_file, set_to_genome);
//
//   GP_MapSampler(metrics, set_to_genome, &pt);
//   PrintMetrics(io_params::set_metric_file, io_params::genome_metric_file, metrics);
//
//   std::cout << "Back to sleep!" << std::endl;
//
//   return 0;
// }
