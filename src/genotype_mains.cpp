#include "genotype_mains.hpp"

namespace simulation_params
{
  uint8_t n_genes=3, colours=7, metric_colours=9;
  uint8_t phenotype_builds=40;
  uint32_t n_samples = 10, n_jiggle = 201;
  std::mt19937 RNG_Engine(std::random_device{}());
  double UND_threshold=0.25;
  bool allow_duplicates = true, STERIC_FORBIDDEN = false;
}

namespace io_params
{
  std::string file_path = "/rscratch/vatj2/public_html/Polyominoes/data/gpmap/V6/experiment/";
  std::string threshold = "_T" + std::to_string((int) ceil(100 * simulation_params::UND_threshold));
  std::string builds = "_B" + std::to_string(simulation_params::phenotype_builds);
  std::string file_details = "_N" + std::to_string(simulation_params::n_genes) + "_C" + std::to_string(simulation_params::colours) + threshold + builds;
  std::string extra="_Cx" + std::to_string(simulation_params::metric_colours) + "_J" + std::to_string(simulation_params::n_jiggle);
  std::string ending=".txt", iso_ending="_Iso.txt";
  // std::string file_path3 = "/rscratch/vatj2/public_html/Polyominoes/data/gpmap/V6/reproducibility/";

  std::string genome_file = file_path + "SampledGenotypes" + file_details + ending;
  std::string phenotype_file = file_path + "PhenotypeTable" + ending;
  std::string set_file = file_path + "SetTable" + file_details + ending;
  std::string preprocess_file = file_path + "PreProcessGenotypes" + file_details + ending;
  std::string set_metric_file = file_path + "SetMetrics" + file_details + extra + ending;
  std::string genome_metric_file = file_path + "GenomeMetrics" + file_details + extra + ending;

  // std::string file_path2 = "/rscratch/vatj2/public_html/Polyominoes/data/gpmap/V6/duplication/";
  std::string file_details2 = "_N" + std::to_string(simulation_params::n_genes + 1) + "_C" + std::to_string(simulation_params::colours) + threshold + builds;
  // std::string duplicate_file = file_path + "DuplicateGenotypes" + file_details2 + ending;
  std::string duplicate_file = file_path + "JiggleDuplicateGenotypes" + file_details2 + extra + ending;
  std::string dup_set_metric_file = file_path + "JiggleDuplicateSetMetrics" + file_details2 + extra + ending;
  std::string dup_genome_metric_file = file_path + "JiggleDuplicateGenomeMetrics" + file_details2 + extra + ending;

  std::string jiggle_file = file_path + "JiggleGenotypes" + file_details + extra + ending;
  std::string jiggle_set_metric_file = file_path + "JiggleSetMetrics" + file_details + extra + ending;
  std::string jiggle_genome_metric_file = file_path + "JiggleGenomeMetrics" + file_details + extra + ending;
}

int main()
{
  std::cout << "Global path : " << io_params::file_path  << "\n";

  // JustExhaustive();
  // ExhaustiveMetricsPrintAll();
  // QuickFromFile();
  // QuickRandom();
  DuplicateJiggle();

  std::cout << "Back to sleep!" << std::endl;

  // PrintConfigFile(io_params::config_file);

  return 0;
}

void JustExhaustive()
{
  PhenotypeTable pt;
  std::vector<Genotype> genomes;
  Set_to_Genome set_to_genome;

  // genomes = ExhaustiveMinimalGenotypesIL(&pt);
  genomes = ExhaustiveMinimalGenotypesFiltered(&pt);
  // genomes = SampleMinimalGenotypes(n_genes, colours, N_SAMPLES, allow_duplicates, &pt);
  PrintGenomeFile(io_params::genome_file, genomes);
  pt.PrintTable(io_params::phenotype_file);

  uint64_t nphen = 0;
  for(auto iter = std::begin(pt.known_phenotypes); iter != std::end(pt.known_phenotypes); iter++)
    nphen += (iter->second).size();
  std::cout << "The phenotype table has " <<+ nphen << " entries \n";
}

void ExhaustiveMetricsPrintAll()
{
  PhenotypeTable pt;
  std::vector<Genotype> genomes;
  Set_to_Genome set_to_genome;
  std::vector<Set_Metrics> metrics;

  // genomes = ExhaustiveMinimalGenotypesIL(&pt);
  genomes = ExhaustiveMinimalGenotypesFiltered(&pt);
  // genomes = SampleMinimalGenotypes(&pt);
  // LoadGenomeFile(genome_file, genomes);
  // FilterExhaustive(genomes, &pt);
  PrintGenomeFile(io_params::genome_file, genomes);
  pt.PrintTable(io_params::phenotype_file);

  PreProcessSampled(genomes, set_to_genome, &pt);
  PrintPreProcessFile(io_params::preprocess_file, set_to_genome);
  PrintSetTable(io_params::set_file, set_to_genome);

  GP_MapSampler(metrics, set_to_genome, &pt);
  PrintMetrics(io_params::set_metric_file, io_params::genome_metric_file, metrics);
}

void QuickRandom()
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

void QuickFromFile()
{
  PhenotypeTable pt;
  std::vector<Genotype> genomes;
  Set_to_Genome set_to_genome;
  std::vector<Set_Metrics> metrics;

  LoadPhenotypeTable(io_params::phenotype_file, &pt);
  LoadGenomeFile(io_params::genome_file, genomes);
  PreProcessSampled(genomes, set_to_genome, &pt);
  GP_MapSampler(metrics, set_to_genome, &pt);
  PrintMetrics(io_params::set_metric_file, io_params::genome_metric_file, metrics);
}

void DuplicateJiggle()
{
  PhenotypeTable pt, pt2;
  std::vector<Genotype> genomes, genomes_jiggle, duplicates;
  Set_to_Genome set_to_genome;
  std::vector<Set_Metrics> metrics;

  LoadGenomeFile(io_params::genome_file, genomes);
  GenomesJiggleDuplication(genomes, genomes_jiggle, duplicates);
  genomes.clear();

  LoadPhenotypeTable(io_params::phenotype_file, &pt);
  PreProcessSampled(genomes_jiggle, set_to_genome, &pt);
  GP_MapSimple(metrics, set_to_genome, &pt);
  PrintMetrics(io_params::jiggle_set_metric_file, io_params::jiggle_genome_metric_file, metrics);

  simulation_params::n_genes++;
  genomes_jiggle.clear(), set_to_genome.clear(), metrics.clear();

  LoadPhenotypeTable(io_params::phenotype_file, &pt2);
  PreProcessSampled(duplicates, set_to_genome, &pt2);
  GP_MapSimple(metrics, set_to_genome, &pt2);
  PrintMetrics(io_params::dup_set_metric_file, io_params::dup_genome_metric_file, metrics);
}
