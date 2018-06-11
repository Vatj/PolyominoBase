#include "genotype_iofunc.hpp"
#include <sstream>
#include <iterator>
#include <string>
#include <utility>


void PrintGenomeFile(std::string genome_file, std::vector<Genotype>& genomes)
{
  std::cout << "Printing genomes to file : " << genome_file << "\n";
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
    // fout << "{";
    // for (auto pID: iter->first)
    //   fout <<+ "(" <<+ pID.first << "," <<+ pID.second << "),";
    // fout.seekp((long) fout.tellp() - 1);
    // fout << "}" << "\n";
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
  std::ofstream set_metric_out(set_metric_file);
  std::ofstream genome_metric_out(genome_metric_file);

  std::cout << "Print metrics to files : \n" << set_metric_file << "\n" << genome_metric_file << "\n";

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

void LoadPhenotypeTable(std::string phenotype_file, StochasticPhenotypeTable* pt_it)
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
