#include "genotype_iofunc.hpp"
#include <python2.7/Python.h>
#include <sstream>
#include <iterator>
#include <string>
#include <utility>

// void IsoCall(std::string file_name)
// {
//   Py_Initialize();
//   PyObject* graph_string = PyString_FromString((char *) "./scripts/graph_methods.py");
//   PyObject* graph_script = PyImport_Import(graph_string);
//
//   PyObject* trim_topo = PyObject_GetAttrString(graph_script, (char*) "Trim_Topologies");
//   PyObject* args = PyString_FromString(file_name.c_str());
//
//   PyObject_CallObject(trim_topo, args);
//   Py_Finalize();
// }

void PrintConfigFile(std::string config_file)
{
  std::ofstream fout(config_file);

  fout << "n_genes : " <<+ simulation_params::n_genes << "\n";
  fout << "colours : " <<+ simulation_params::colours << "\n";
  fout << "metric_colours : " <<+ simulation_params::metric_colours << "\n";
  fout << "n_samples : " <<+ simulation_params::n_samples << "\n";
  fout << "Threshold : " <<+ simulation_params::UND_threshold << "\n";
  fout << "Phenotype builds : " <<+ simulation_params::phenotype_builds << "\n";
  fout << "n_jiggle : " <<+ simulation_params::n_jiggle << "\n";
}

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
    fout << "{";
    for (auto pID: ref_pIDs)
      fout <<+ "(" <<+ pID.first << "," <<+ pID.second << "),";
    fout.seekp((long) set_out.tellp() - 1);
    fout << "} ";

    fout << "[";
    for (auto original: originals)
    {
      set_out << "(";
      for (auto face: original)
        set_out <<+ face << ",";
      set_out.seekp((long) set_out.tellp() - 1);
      set_out << "),";
    }
    set_out.seekp((long) set_out.tellp() - 1);
    set_out << "]\n";
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
    fout << "} "<<+ (iter->second).size() << std::endl;
  }
}

void PrintMetrics(std::string set_metric_file, std::string genome_metric_file, std::vector<Set_Metrics> metrics)
{
  // Create new files and open them for writing (erase previous data)
  std::ofstream set_metric_out(set_metric_file);
  std::ofstream genome_metric_out(genome_metric_file);

  // Logging
  std::cout << "Print metrics to files : \n";
  std::cout << set_metric_file << "\n" << genome_metric_file << "\n";

  // Header for the metric files
  set_metric_out << "srobustness irobustness meta_evolvability evolvability";
  set_metric_out << " rare unbound analysed misclassified neutral_size";
  set_metric_out << " diversity diversity_tracker originals pIDs\n";

  genome_metric_out << "genome original srobustness irobustness";
  genome_metric_out << " meta_evolvability evolvability rare unbound diversity";
  genome_metric_out << " neutral_weight frequencies pIDs\n";

  for (auto metric: metrics)
    metric.save_to_file(set_metric_out, genome_metric_out);
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
  "Defunct Code for now."
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
