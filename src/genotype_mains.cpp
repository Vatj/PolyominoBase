#include "genotype_mains.hpp"
// #include <boost/lambda/lambda.hpp>
// #include <boost/regex.hpp>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <iterator>


namespace simulation_params
{
  uint16_t n_genes, colours, metric_colours;
  uint16_t phenotype_builds, preprocess_builds;
  uint32_t n_samples, n_jiggle;
  double UND_threshold;
  bool allow_duplicates, STERIC_FORBIDDEN, iso, dup_aware;
  bool exhaustive, metrics, duplicate_exhaustive;

  std::mt19937 RNG_Engine(std::random_device{}());
}

int main (int argc, char *argv[])
{

  try {

    po::options_description config("Configuration");
    config.add_options()
      ("help", "produce help message")
      ("config",
        po::value<std::string>(&io_params::config_file)->default_value("configuration.cfg"),
        "Specify a configuration file")
      ("n_genes,n",
        po::value<uint16_t>(&simulation_params::n_genes)->default_value(3),
        "Number of gene in genome (modulo 4)")
      ("generate_colours,c",
        po::value<uint16_t>(&simulation_params::colours)->default_value(7),
        "Allowed colours in generated genomes")
      ("metric_colours,m",
        po::value<uint16_t>(&simulation_params::metric_colours)->default_value(9),
        "Allowed colours to jiggle genomes")
      ("builds,b",
        po::value<uint16_t>(&simulation_params::phenotype_builds)->default_value(40),
        "Number of random assembly trial per gene")
      ("threshold,t",
        po::value<double>(&simulation_params::UND_threshold)->default_value(0.25),
        "Threshold percentage to save a shape")
      ("n_jiggle,j",
        po::value<uint32_t>(&simulation_params::n_jiggle)->default_value(3),
        "Number of analysed genome per representant")
      ("n_samples",
        po::value<uint32_t>(&simulation_params::n_samples)->default_value(10),
        "Number of minimal represant to save. Only apply to random sampling")
      ("dup_aware",
        po::value<bool>(&simulation_params::dup_aware)->default_value(false),
        "Duplicate genes in a genome are jiggled in the same way")
      ("iso",
        po::value<bool>(&simulation_params::iso)->default_value(false),
        "Run for isomorphism trimmed files")
    ;

    po::options_description execution("Execution");
    execution.add_options()
      ("exhaustive",
        po::value<bool>(&simulation_params::exhaustive)->default_value(false),
        "Generate all minimal genomes and save pseudo-deterministic")
      ("metrics",
        po::value<bool>(&simulation_params::metrics)->default_value(false),
        "Compute the metric distributions by jiggling minimal genomes")
      ("duplicate_exhaustive",
        po::value<bool>(&simulation_params::duplicate_exhaustive)->default_value(false),
        "Generate all minimal genomes and save pseudo-deterministic genomes with a dupicate gene")
    ;

    po::options_description io_options("IO");
    io_options.add_options()
      ("file_path",
        po::value<std::string>(&io_params::file_path),
        "Main file path")
      ("in_genome_file",
        po::value<std::string>(&io_params::in_genome_file),
        "file containing the input genomes")
      ("out_genome_file",
        po::value<std::string>(&io_params::out_genome_file),
        "file containing the output genomes")
      ("duplicate_genome_file",
        po::value<std::string>(&io_params::duplicate_file),
        "file containing the output genomes with a duplicated gene")
      ("in_phenotype_file",
        po::value<std::string>(&io_params::in_phenotype_file),
        "file containing the input phenotype table")
      ("out_phenotype_file",
        po::value<std::string>(&io_params::out_phenotype_file),
        "file containing the output phenotype table")
      ("genome_metric_file",
        po::value<std::string>(&io_params::genome_metric_file),
        "file containing the output genome metrics")
      ("set_metric_file",
        po::value<std::string>(&io_params::set_metric_file),
        "file containing the output set metrics")
      ("set_file",
        po::value<std::string>(&io_params::set_file),
        "file containing allowed set found when preprocessing")
      ("preprocess_file",
        po::value<std::string>(&io_params::preprocess_file),
        "file containing the output map between pID set and genomes")
    ;

    po::options_description hidden("Hidden");
    hidden.add_options()
      ("preprocess_builds",
        po::value<uint16_t>(&simulation_params::preprocess_builds)->default_value(250),
        "Same as builds but only for preprocessing function")
      ("allow_duplicates",
        po::value<bool>(&simulation_params::allow_duplicates)->default_value(false),
        "Gene duplication is allowed when generating minimal genomes")
      ("steric_forbidden",
        po::value<bool>(&simulation_params::STERIC_FORBIDDEN)->default_value(false),
        "Steric constraint is relaxed in the assembly process")
    ;

    po::options_description cmdline_options;
    cmdline_options.add(config).add(execution).add(io_options).add(hidden);

    po::options_description config_file_options;
    config_file_options.add(config).add(execution).add(io_options).add(hidden);

    po::options_description visible("Allowed options");
    visible.add(config).add(execution).add(io_options);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
    po::notify(vm);

    std::ifstream ifs(io_params::config_file.c_str());
    if (!ifs)
    {
        std::cout << "can not open config file: " << io_params::config_file << "\n";
        return 0;
    }
    else
    {
        po::store(parse_config_file(ifs, config_file_options), vm);
        po::notify(vm);
    }

    if (vm.count("help")) {
          std::cout << visible << "\n";
          return 0;
    }

    all_files_to_full_names();

    if(vm["exhaustive"].as<bool>())
      JustExhaustive();
    else if(vm["metrics"].as<bool>())
      QuickFromFile();
    else if(vm["duplicate_exhaustive"].as<bool>())
      DuplicateExhaustive();
    else
      std::cout << "Why waking me up?" << std::endl;
    }
    catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!\n";
    }

  std::cout << "Global path : " << io_params::file_path  << "\n";


  // JustExhaustive();
  // // ExhaustiveMetricsPrintAll();
  // // QuickFromFile();
  // // QuickRandom();
  // // DuplicateExhaustive();
  //
  std::cout << "Back to sleep!" << std::endl;

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
  PrintGenomeFile(io_params::out_genome_file, genomes);
  pt.PrintTable(io_params::out_phenotype_file);
}

void ExhaustiveMetricsPrintAll()
{
  PhenotypeTable pt;
  std::vector<Genotype> genomes;
  Set_to_Genome set_to_genome;

  // genomes = ExhaustiveMinimalGenotypesIL(&pt);
  genomes = ExhaustiveMinimalGenotypesFiltered(&pt);
  // genomes = SampleMinimalGenotypes(&pt);
  // LoadGenomeFile(simulation_params::in_genome_file, genomes);
  // FilterExhaustive(genomes, &pt);
  PrintGenomeFile(io_params::out_genome_file, genomes);
  pt.PrintTable(io_params::out_phenotype_file);

  PreProcessSampled(genomes, set_to_genome, &pt);
  PrintPreProcessFile(io_params::preprocess_file, set_to_genome);
  PrintSetTable(io_params::set_file, set_to_genome);

  GP_MapSampler(set_to_genome, &pt);
}

void QuickRandom()
{
  PhenotypeTable pt;
  std::vector<Genotype> genomes;
  Set_to_Genome set_to_genome;

  genomes = SampleMinimalGenotypes(&pt);
  PrintGenomeFile(io_params::out_genome_file, genomes);
  PreProcessSampled(genomes, set_to_genome, &pt);
  GP_MapSampler(set_to_genome, &pt);
}

void QuickFromFile()
{
  PhenotypeTable pt;
  std::vector<Genotype> genomes;
  Set_to_Genome set_to_genome;

  LoadPhenotypeTable(io_params::in_phenotype_file, &pt);
  LoadGenomeFile(io_params::in_genome_file, genomes);
  PreProcessSampled(genomes, set_to_genome, &pt);
  GP_MapSampler(set_to_genome, &pt);
}


void DuplicateExhaustive()
{
  PhenotypeTable pt;
  std::vector<Genotype> genomes, duplicates;

  duplicates = ExhaustiveMinimalGenotypesFilteredDuplicate(genomes, &pt);
  PrintGenomeFile(io_params::out_genome_file, genomes);
  PrintGenomeFile(io_params::duplicate_file, duplicates);
}
