#include "genotype_iofunc.hpp"
#include <sstream>
#include <iterator>
#include <string>
#include <utility>

namespace io_params
{
  std::string file_path, file_details, extra, config_file;
  std::string out_genome_file, in_genome_file;
  std::string duplicate_file;
  std::string out_phenotype_file, in_phenotype_file;
  std::string set_file, preprocess_file;
  std::string set_metric_file, genome_metric_file;
  std::string ending=".txt", iso_ending="_Iso.txt";
}

void infer_file_details()
{
  using namespace simulation_params;

  std::string file_details = "_N" + std::to_string(n_genes);
  file_details += "_C" + std::to_string(colours);
  file_details += "_T" + std::to_string((int) ceil(100 * UND_threshold));
  file_details += "_B" + std::to_string(phenotype_builds);
  io_params::file_details = file_details;

  std::string extra = "_Cx" + std::to_string(metric_colours);
  extra += "_J" + std::to_string(n_jiggle);
  io_params::extra = extra;
}

void all_files_to_full_names()
{
  using namespace io_params;
  infer_file_details();

  in_phenotype_file = file_path + in_phenotype_file + ending;

  in_genome_file = full_filename(in_genome_file, false, simulation_params::iso);
  out_phenotype_file = full_filename(out_phenotype_file, false, false);
  genome_metric_file = full_filename(genome_metric_file, true, simulation_params::iso);
  set_metric_file = full_filename(set_metric_file, true, simulation_params::iso);
  set_file = full_filename(set_file, false, simulation_params::iso);
  preprocess_file = full_filename(preprocess_file, false, simulation_params::iso);
  duplicate_file = full_filename(duplicate_file, false, simulation_params::iso);

  if(simulation_params::duplicate_exhaustive)
  {
    simulation_params::n_genes -= 1;
    infer_file_details();
    out_genome_file = full_filename(out_genome_file, false, false);
    simulation_params::n_genes += 1;
  } else
    out_genome_file = full_filename(out_genome_file, false, false);
}

std::string full_filename(std::string file_name, bool extra, bool iso)
{
  std::string full = io_params::file_path + file_name + io_params::file_details;
  if(extra)
    full += io_params::extra;
  if(iso)
    full += io_params::iso_ending;
  else
    full += io_params::ending;

  return full;
}


// void PrintConfigFile(std::string config_file)
// {
//   std::ofstream fout(config_file);
//
//   fout << "n_genes : " <<+ simulation_params::n_genes << "\n";
//   fout << "colours : " <<+ simulation_params::colours << "\n";
//   fout << "metric_colours : " <<+ simulation_params::metric_colours << "\n";
//   fout << "n_samples : " <<+ simulation_params::n_samples << "\n";
//   fout << "Threshold : " <<+ simulation_params::UND_threshold << "\n";
//   fout << "Phenotype builds : " <<+ simulation_params::phenotype_builds << "\n";
//   fout << "n_jiggle : " <<+ simulation_params::n_jiggle << "\n";
// }

void PrintGenomeFile(std::string genome_file, std::vector<Genotype>& genomes)
{
  std::cout << "Printing " <<+ genomes.size() << " genomes to file : " << genome_file << "\n";
  std::ofstream fout(genome_file);
  for(auto genome: genomes)
  {
    for(auto base: genome)
      fout <<+ base << " ";
    fout << "\n";
  }
}

void PrintPreProcessFile(std::string preprocess_file, Set_to_Genome& set_to_genome)
{
  std::cout << "Printing preprocessed genomes to file : " << preprocess_file << "\n";
  std::ofstream fout(preprocess_file);

  for(Set_to_Genome::iterator iter = std::begin(set_to_genome); iter != std::end(set_to_genome); iter++)
  {
    fout << "x ";
    for (auto pID: iter->first)
      fout <<+ pID.first << " " <<+ pID.second << " ";
    fout << "\n";

    for(auto genome: iter->second)
    {
      for(auto index: genome)
        fout <<+ index << " ";
      fout << "\n";
    }
  }
}

void PrintPreProcessFile2(std::string preprocess_file, Set_to_Genome& set_to_genome)
{
  std::cout << "Printing preprocessed genomes to file : " << preprocess_file << "\n";
  std::ofstream fout(preprocess_file);

  for(Set_to_Genome::iterator iter = std::begin(set_to_genome); iter != std::end(set_to_genome); iter++)
  {
    for(auto genome: iter->second)
    {
      for(auto index: genome)
        fout <<+ index << " ";
      fout << "{";
      for (auto pID: iter->first)
        fout <<+ pID.first << " " <<+ pID.second << " ";
      fout << "}\n";
    }
  }
}

void PrintSetTable(std::string set_file, Set_to_Genome& set_to_genome)
{
  std::cout << "Printing set table to file : " << set_file << "\n";
  std::ofstream fout(set_file);

  for(Set_to_Genome::iterator iter = std::begin(set_to_genome); iter != std::end(set_to_genome); iter++)
  {
    fout << "{";
    for (auto pID: iter->first)
      fout <<+ "(" <<+ pID.first << "," <<+ pID.second << "),";
    fout.seekp((long) fout.tellp() - 1);
    fout << "} ";

    fout << "[";
    for (auto original: iter->second)
    {
      fout << "(";
      for (auto face: original)
        fout <<+ face << ",";
      fout.seekp((long) fout.tellp() - 1);
      fout << "),";
    }
    fout.seekp((long) fout.tellp() - 1);
    fout << "]\n";
  }
}

void header_metric_files(std::ofstream& set_metric_out, std::ofstream& genome_metric_out)
{
  // Logging
  std::cout << "Print metrics to files : \n";
  std::cout << io_params::set_metric_file << "\n" << io_params::genome_metric_file << "\n";

  // Header for the metric files
  set_metric_out << "srobustness irobustness evolvability";
  set_metric_out << " robust_evolvability complex_evolvability";
  set_metric_out << " complex_diversity rare unbound";
  set_metric_out << " analysed misclassified neutral_size";
  set_metric_out << " diversity diversity_tracker originals misclassified_details pIDs\n";

  genome_metric_out << "genome original srobustness irobustness";
  genome_metric_out << " evolvability complex_evolvability robust_evolvability";
  genome_metric_out << " rare unbound complex_diversity diversity";
  genome_metric_out << " neutral_weight frequencies pIDs\n";
}



void LoadGenomeFile(std::string genome_file, std::vector<Genotype>& genomes)
{
  std::string str;
  Genotype genotype;
  std::ifstream genome_in(genome_file);

  std::cout << "Loading genome from file : " << genome_file << "\n";

  while (std::getline(genome_in, str))
  {
    std::istringstream is( str );
    genotype.assign( std::istream_iterator<int>( is ), std::istream_iterator<int>() );
    genomes.emplace_back(genotype);
  }
}


void LoadPreProcessFile(std::string preprocess_file, Set_to_Genome& set_to_genome)
{
  std::string str;
  Genotype genotype;
  std::vector<int> pre_pIDs;
  std::vector<Phenotype_ID> pIDs;
  std::ifstream fin(preprocess_file);

  std::cout << "Loading preprocess file : " << preprocess_file << "\n";

  while (std::getline(fin, str))
  {
    if(str.compare(0, 1, "x"))
    {
      std::istringstream is(str);
      pre_pIDs.assign(std::istream_iterator<int>(is), std::istream_iterator<int>());

      for(uint8_t index=0; index < pre_pIDs.size() - 1; index+=2)
        pIDs.emplace_back(std::make_pair(pre_pIDs[index], pre_pIDs[index + 1]));
      continue;
    }

    std::istringstream is( str );
    genotype.assign( std::istream_iterator<int>( is ), std::istream_iterator<int>() );
    set_to_genome[pIDs].emplace_back(genotype);
  }
}

void LoadPhenotypeTable(std::string phenotype_file, PhenotypeTable* pt_it)
{
  std::ifstream pheno_in(phenotype_file);
  std::string str;

  std::cout << "Loading phenotype table : " << phenotype_file << "\n";

  while (std::getline(pheno_in, str))
  {
    std::stringstream iss(str);
    int number;
    std::vector<uint8_t> phenotype_line;
    while (iss>>number)
      phenotype_line.push_back(static_cast<uint8_t>(number));
    Phenotype phen;
    phen.dx=phenotype_line[2];
    phen.dy=phenotype_line[3];
    phen.tiling=std::vector<uint8_t>(phenotype_line.begin()+4, phenotype_line.end());
    pt_it->known_phenotypes[phenotype_line[0]].emplace_back(phen);
  }
}
