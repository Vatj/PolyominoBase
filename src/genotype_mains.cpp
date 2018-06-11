#include "genotype_mains.hpp"

int main()
{
  // ExhaustiveFast();
  // ExhaustiveMetricsAndPrint();
  // Quick();
  // Duplicate();
  // FullFast();
  QuickFromIso();

  return 0;
}

void Duplicate()
{
  uint8_t n_genes = 3, colours = 7, metric_colours = 9;
  StochasticPhenotypeTable pt;
  std::vector<Genotype> genomes, duplicates;
  Set_to_Genome set_to_genome;
  std::vector<Set_Metrics> metrics;

  std::string file_path("/rscratch/vatj2/Polyominoes/data/gpmap/V4/exhaustive/");
  std::string file_details = "_N" + std::to_string(n_genes) + "_C" + std::to_string(colours);
  std::string extra="_Cx" + std::to_string(metric_colours), ending=".txt";

  std::string genome_file = file_path + "SampledGenotypes" + file_details + ending;
  std::string phenotype_file = file_path + "PhenotypeTable" + file_details + ending;

  LoadPhenotypeTable(phenotype_file, &pt);
  LoadGenomeFile(genome_file, genomes);

  duplicates = GenomesDuplication(genomes, n_genes);

  n_genes++;
  file_path = "/rscratch/vatj2/Polyominoes/data/gpmap/V4/duplication2/";
  file_details = "_N" + std::to_string(n_genes) + "_C" + std::to_string(colours);
  phenotype_file = file_path + "PhenotypeTable" + file_details + ending;

  std::string set_file = file_path + "SetTable" + file_details + ending;
  std::string preprocess_file = file_path + "PreProcessGenotypes" + file_details + ending;
  std::string duplicate_file = file_path + "DuplicateGenotypes" + file_details + ending;
  std::string set_metric_file = file_path + "SetMetrics" + file_details + extra + ending;
  std::string genome_metric_file = file_path + "GenomeMetrics" + file_details + extra + ending;

  PrintGenomeFile(duplicate_file, duplicates);
  PreProcessSampled(duplicates, set_to_genome, &pt);
  PrintPreProcessFile(preprocess_file, set_to_genome);
  PrintSetTable(set_file, set_to_genome);

  GP_MapSampler(n_genes, metric_colours, metrics, set_to_genome, &pt);
  PrintMetrics(set_metric_file, genome_metric_file, metrics);

  pt.PrintTable(phenotype_file);

  std::cout << "Back to sleep!" << std::endl;
}


void ExhaustiveFast()
{
  uint8_t n_genes = 4, colours = 7;
  // uint32_t N_SAMPLES = 300;
  StochasticPhenotypeTable pt;
  // bool allow_duplicates = true;
  std::vector<Genotype> genomes;
  Set_to_Genome set_to_genome;

  std::string file_path("/rscratch/vatj2/Polyominoes/data/gpmap/V4/exhaustive2/");
  std::string file_details = "_N" + std::to_string(n_genes) + "_C" + std::to_string(colours), ending=".txt";

  std::string genome_file = file_path + "SampledGenotypes" + file_details + ending;
  std::string phenotype_file = file_path + "PhenotypeTable" + file_details + ending;
  std::string set_file = file_path + "SetTable" + file_details + ending;
  std::string preprocess_file = file_path + "PreProcessGenotypes" + file_details + ending;

  std::cout << "Global path : " << file_path  << "\n";

  genomes = ExhaustiveMinimalGenotypes(n_genes, colours, &pt);
  // genomes = SampleMinimalGenotypes(n_genes, colours, N_SAMPLES, allow_duplicates, &pt);
  PrintGenomeFile(genome_file, genomes);
  pt.PrintTable(phenotype_file);

  PreProcessSampled(genomes, set_to_genome, &pt);
  PrintPreProcessFile(preprocess_file, set_to_genome);
  PrintSetTable(set_file, set_to_genome);

  std::cout << "Back to sleep!" << std::endl;
}


void FullFast()
{
  uint8_t n_genes = 2, colours = 3, metric_colours=3;
  // uint32_t N_SAMPLES = 300;
  StochasticPhenotypeTable pt;
  // bool allow_duplicates = true;
  std::vector<Genotype> genomes;
  Set_to_Genome set_to_genome;
  std::vector<Set_Metrics> metrics;

  std::string file_path("/rscratch/vatj2/public_html/Polyominoes/data/gpmap/V4/full/");
  std::string file_details = "_N" + std::to_string(n_genes) + "_C" + std::to_string(colours);
  std::string extra="_Cx" + std::to_string(metric_colours), ending=".txt";

  std::string genome_file = file_path + "SampledGenotypes" + file_details + ending;
  std::string phenotype_file = file_path + "PhenotypeTable" + file_details + ending;
  std::string set_file = file_path + "SetTable" + file_details + ending;
  std::string preprocess_file = file_path + "PreProcessGenotypes" + file_details + ending;
  std::string set_metric_file = file_path + "SetMetrics" + file_details + extra + ending;
  std::string genome_metric_file = file_path + "GenomeMetrics" + file_details + extra + ending;

  std::cout << "Global path : " << file_path  << "\n";

  genomes = ExhaustiveFullGenotypes2(colours, &pt);
  // genomes = SampleMinimalGenotypes(n_genes, colours, N_SAMPLES, allow_duplicates, &pt);
  // PrintGenomeFile(genome_file, genomes);
  // pt.PrintTable(phenotype_file);

  PreProcessSampled(genomes, set_to_genome, &pt);
  // PrintPreProcessFile(preprocess_file, set_to_genome);
  // PrintSetTable(set_file, set_to_genome);

  GP_MapSampler(n_genes, metric_colours, metrics, set_to_genome, &pt);
  PrintMetrics(set_metric_file, genome_metric_file, metrics);

  std::cout << "Back to sleep!" << std::endl;
}

void ExhaustiveMetricsAndPrint()
{
  uint8_t n_genes = 3, colours = 7, metric_colours = 9;
  StochasticPhenotypeTable pt;
  std::vector<Genotype> genomes;
  Set_to_Genome set_to_genome;
  std::vector<Set_Metrics> metrics;

  std::string file_path("/rscratch/vatj2/Polyominoes/data/gpmap/V4/exhaustive/");
  std::string file_details = "_N" + std::to_string(n_genes) + "_C" + std::to_string(colours);
  std::string extra="_Cx" + std::to_string(metric_colours), ending=".txt";

  std::string genome_file = file_path + "SampledGenotypes" + file_details + ending;
  std::string phenotype_file = file_path + "PhenotypeTable" + file_details + ending;
  std::string set_file = file_path + "SetTable" + file_details + ending;
  std::string preprocess_file = file_path + "PreProcessGenotypes" + file_details + ending;
  std::string set_metric_file = file_path + "SetMetrics" + file_details + extra + ending;
  std::string genome_metric_file = file_path + "GenomeMetrics" + file_details + extra + ending;

  std::cout << "Global path : " << file_path  << "\n";

  genomes = ExhaustiveMinimalGenotypes(n_genes, colours, &pt);
  PrintGenomeFile(genome_file, genomes);
  pt.PrintTable(phenotype_file);

  PreProcessSampled(genomes, set_to_genome, &pt);
  PrintPreProcessFile(preprocess_file, set_to_genome);
  PrintSetTable(set_file, set_to_genome);

  GP_MapSampler(n_genes, metric_colours, metrics, set_to_genome, &pt);
  PrintMetrics(set_metric_file, genome_metric_file, metrics);

  std::cout << "Back to sleep!" << std::endl;
}

void Quick()
{
  uint8_t n_genes = 3, colours = 7, metric_colours = 9;
  uint32_t N_SAMPLES = 20;
  StochasticPhenotypeTable pt;
  bool allow_duplicates = true;
  std::vector<Genotype> genomes;
  Set_to_Genome set_to_genome;
  std::vector<Set_Metrics> metrics;

  std::string file_path("/rscratch/vatj2/public_html/Polyominoes/data/gpmap/V4/experiment/");
  std::string file_details = "_N" + std::to_string(n_genes) + "_C" + std::to_string(colours);
  std::string extra="_Cx" + std::to_string(metric_colours), ending=".txt";

  std::string genome_file = file_path + "SampledGenotypes" + file_details + ending;
  std::string set_metric_file = file_path + "SetMetrics" + file_details + extra + ending;
  std::string genome_metric_file = file_path + "GenomeMetrics" + file_details + extra + ending;

  std::cout << "Global path : " << file_path  << "\n";

  genomes = SampleMinimalGenotypes(n_genes, colours, N_SAMPLES, allow_duplicates, &pt);
  PrintGenomeFile(genome_file, genomes);
  PreProcessSampled(genomes, set_to_genome, &pt);
  GP_MapSampler(n_genes, metric_colours, metrics, set_to_genome, &pt);
  PrintMetrics(set_metric_file, genome_metric_file, metrics);

  std::cout << "Back to sleep!" << std::endl;
}

void QuickFromIso()
{
  uint8_t n_genes = 3, colours = 7, metric_colours = 9;
  StochasticPhenotypeTable pt;
  std::vector<Genotype> genomes;
  Set_to_Genome set_to_genome;
  std::vector<Set_Metrics> metrics;

  std::string file_path("/rscratch/vatj2/public_html/Polyominoes/data/gpmap/V4/exhaustive/");
  std::string file_details = "_N" + std::to_string(n_genes) + "_C" + std::to_string(colours);
  std::string extra="_Cx" + std::to_string(metric_colours), ending=".txt";

  std::string genome_file = file_path + "SampledGenotypes" + file_details + "_Iso" + ending;
  std::string set_metric_file = file_path + "SetMetrics" + file_details + extra + "_Iso" + ending;
  std::string genome_metric_file = file_path + "GenomeMetrics" + file_details + extra + "_Iso" + ending;

  std::cout << "Global path : " << file_path  << "\n";

  LoadGenomeFile(genome_file, genomes);
  PreProcessSampled(genomes, set_to_genome, &pt);
  GP_MapSampler(n_genes, metric_colours, metrics, set_to_genome, &pt);
  PrintMetrics(set_metric_file, genome_metric_file, metrics);

  std::cout << "Back to sleep!" << std::endl;
}



// int main()
// {
//   uint8_t n_genes = 3, colours = 7, metric_colours = 9;
//   uint32_t N_SAMPLES = 10000;
//   StochasticPhenotypeTable pt;
//   bool allow_duplicates = true;
//   std::vector<Genotype> genomes;
//   Set_to_Genome set_to_genome;
//   std::vector<Set_Metrics> metrics;
//
//   bool print_genome = true, print_preprocess = true;
//   bool print_set_table = true, print_phenotype_table = true;
//   bool load_genome = false, load_preprocess = false, load_phenotype_table = false;
//
//   std::string file_path("/rscratch/vatj2/Polyominoes/data/gpmap/V4/experiment/");
//   std::string file_details = "_N" + std::to_string(n_genes) + "_C" + std::to_string(colours);
//   std::string extra="_Cx" + std::to_string(metric_colours), ending=".txt";
//
//   std::string genome_file = file_path + "SampledGenotypes" + file_details + ending;
//   std::string phenotype_file = file_path + "PhenotypeTable" + file_details + ending;
//   std::string set_file = file_path + "SetTable" + file_details + ending;
//   std::string preprocess_file = file_path + "PreProcessGenotypes" + file_details + ending;
//   std::string duplicate_file = file_path + "DuplicateGenotypes" + file_details + ending;
//   std::string set_metric_file = file_path + "SetMetrics" + file_details + extra + ending;
//   std::string genome_metric_file = file_path + "GenomeMetrics" + file_details + extra + ending;
//
//   if(load_phenotype_table)
//     LoadPhenotypeTable(phenotype_file, &pt);
//
//   if(load_genome)
//     LoadGenomeFile(genome_file, genomes);
//   else
//     genomes = SampleMinimalGenotypes(n_genes, colours, N_SAMPLES, allow_duplicates, &pt);
//
//   if(print_genome)
//     PrintGenomeFile(genome_file, genomes);
//
//   // std::vector<Genotype> dups = GenomesDuplication(genomes, n_genes);
//   // PrintGenomeFile(duplicate_file, dups);
//
//   if(print_phenotype_table)
//     pt.PrintTable(phenotype_file);
//
//   if(load_preprocess)
//     LoadPreProcessFile(preprocess_file, set_to_genome);
//   else
//     PreProcessSampled(genomes, set_to_genome, &pt);
//
//   if(print_preprocess)
//     PrintPreProcessFile(preprocess_file, set_to_genome);
//
//   if(print_set_table)
//     PrintSetTable(set_file, set_to_genome);
//
//   GP_MapSampler(n_genes, metric_colours, metrics, set_to_genome, &pt);
//   PrintMetrics(set_metric_file, genome_metric_file, metrics);
//
//   std::cout << "Back to sleep!" << std::endl;
//
//   return 0;
// }
