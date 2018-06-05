#include "genotype_phenotype.hpp"
#include <sstream>
#include <iterator>
#include <functional>
#include <set>

// std::mt19937 RNG_Engine1(std::random_device{}());

std::vector<Phenotype_ID> GetPhenotypeIDs(Genotype& genotype, uint8_t k_builds, StochasticPhenotypeTable* pt_it)
{
  std::vector<Phenotype_ID> pIDs;
  Clean_Genome(genotype,0,false);

  const uint8_t n_genes=genotype.size()/4;
  std::map<uint8_t,uint8_t> dups=DuplicateGenes(genotype);
  std::map<uint8_t,Phenotype_ID> seed_m;

  for(uint8_t seed=0;seed<genotype.size()/4;++seed)
    seed_m[seed]=Stochastic::Analyse_Genotype_Outcome(genotype,k_builds,pt_it,seed);
  uint8_t index=0;

  Phenotype_ID phen_id;

  for(uint8_t gene=0;gene<n_genes; ++gene) {
    pIDs.emplace_back(dups.count(gene) ? seed_m[dups[gene]] : phen_id=seed_m[index++]);
  }
  std::sort(pIDs.begin(),pIDs.end());
  return pIDs;
}

std::vector<Phenotype_ID> GetSetPIDs(Genotype& genotype, uint8_t k_builds, StochasticPhenotypeTable* pt_it)
{
  Clean_Genome(genotype,0,false);
  double UND_frac = 0.1;

  std::vector<Phenotype_ID> pIDs = Stochastic::AssemblePlasticGenotype(genotype, k_builds, pt_it, UND_frac);

  std::sort(pIDs.begin(),pIDs.end());
  return pIDs;
}

void GetPhenotypesIDs(const char* file_path_c, const char* file_name_c, uint8_t n_genes, uint8_t colours)
{
  std::string file_path(file_path_c),file_name(file_name_c),str,details="_N"+std::to_string(n_genes)+"_C"+std::to_string(colours)+".txt";
  std::ifstream file_in(file_path+file_name);
  std::ofstream gfout(file_path+"Genotype_Codes"+details, std::ios_base::out);
  std::ofstream pfout(file_path+"Phenotype_Table"+details, std::ios_base::out);

  StochasticPhenotypeTable pt;
  uint8_t k_builds=10;
  Genotype genotype(n_genes*4);

  while (std::getline(file_in, str))
  {
    std::istringstream is( str );
    genotype.assign( std::istream_iterator<int>( is ), std::istream_iterator<int>() );

    for(auto phen_id : GetPhenotypeIDs(genotype, k_builds, &pt))
      gfout<<+phen_id.first<<" "<<+phen_id.second<<" ";
    gfout<<"\n";
  }
  pt.PrintTable(pfout);
}

void PreProcessGenotypesTopology(const char* file_path_c, uint8_t n_genes, uint8_t colours)
{
  std::string str,file_path(file_path_c),file_details="_N"+std::to_string(n_genes)+"_C"+std::to_string(colours);
  std::ifstream fin(file_path+"SampledGenotypes"+file_details+"_Iso.txt");
  std::ofstream fout(file_path+"SampledGenotypes"+file_details+"_Processed.txt");

  StochasticPhenotypeTable pt;
  uint8_t k_builds=5;
  std::map<std::vector<Phenotype_ID>, std::vector<Genotype>> phen_sets;
  Genotype genotype;

  while (std::getline(fin, str))
  {
    std::istringstream is( str );
    genotype.assign( std::istream_iterator<int>( is ), std::istream_iterator<int>() );

    const Genotype g_write(genotype);
    std::vector<Phenotype_ID> pIDs = GetSetPIDs(genotype, k_builds, &pt);
    phen_sets[pIDs].emplace_back(g_write);
  }

  for(auto kv : phen_sets)
  {
    for(Genotype& gen : kv.second)
    {
      for(uint8_t base : gen)
      {
	       fout<<+base<<" ";
      }
      fout<<"\n";
    }

    fout << "x ";
    for(auto pID: kv.first)
      fout <<+ pID.first << " " <<+ pID.second << " ";
    fout << "\n";

  }
}

void PreProcessGenotypes(const char* file_path_c, uint8_t n_genes, uint8_t colours)
{
  std::string str,file_path(file_path_c),file_details="_N"+std::to_string(n_genes)+"_C"+std::to_string(colours);
  std::ifstream fin(file_path+"SampledGenotypes"+file_details+".txt");
  std::ofstream fout(file_path+"SampledGenotypes"+file_details+"_Processed.txt");

  StochasticPhenotypeTable pt;
  uint8_t k_builds=5;
  std::map<std::vector<Phenotype_ID>, std::vector<Genotype>> phen_sets;
  Genotype genotype;

  while (std::getline(fin, str))
  {
    std::istringstream is( str );
    genotype.assign( std::istream_iterator<int>( is ), std::istream_iterator<int>() );

    const Genotype g_write(genotype);
    std::vector<Phenotype_ID> pIDs = GetSetPIDs(genotype, k_builds, &pt);
    phen_sets[pIDs].emplace_back(g_write);
  }

  for(auto kv : phen_sets)
  {
    for(Genotype& gen : kv.second)
    {
      for(uint8_t base : gen)
      {
	       fout<<+base<<" ";
      }
      fout<<"\n";
    }
    fout << "x ";
    for(auto pID: kv.first)
      fout <<+ pID.first << " " <<+ pID.second << " ";
    fout << "\n";

  }
}

void PreProcessSampled(std::ifstream& genome_in, Set_to_Genome& set_to_genome, StochasticPhenotypeTable* pt, std::ofstream& pIDset_out)
{
  std::string str;
  uint8_t k_builds=100;
  Genotype genotype;

  while (std::getline(genome_in, str))
  {
    std::istringstream is( str );
    genotype.assign( std::istream_iterator<int>( is ), std::istream_iterator<int>() );

    const Genotype g_write(genotype);
    std::vector<Phenotype_ID> pIDs = GetSetPIDs(genotype, k_builds, pt);
    set_to_genome[pIDs].emplace_back(genotype);
  }

  PrintSetTable(pIDset_out, set_to_genome);
}

void PrintPreProcessFile(std::ofstream& fout, Set_to_Genome& set_to_genome)
{
  for(Set_to_Genome::iterator iter = std::begin(set_to_genome); iter != std::end(set_to_genome); iter++)
  {
    fout << "x ";
    fout << "{";
    for (auto pID: iter->first)
      fout <<+ "(" <<+ pID.first << "," <<+ pID.second << "),";
    fout.seekp((long) fout.tellp() - 1);
    fout << "}" << "\n";

    for(auto genome: iter->second)
    {
      for(auto index: genome)
        fout <<+ index << " ";
      fout << "\n";
    }
  }
}

void PrintSetTable(std::ofstream& fout, Set_to_Genome& set_to_genome)
{
  for(Set_to_Genome::iterator iter = std::begin(set_to_genome); iter != std::end(set_to_genome); iter++)
  {
    fout << "{";
    for (auto pID: iter->first)
      fout <<+ "(" <<+ pID.first << "," <<+ pID.second << "),";
    fout.seekp((long) fout.tellp() - 1);
    fout << "} "<<+ (iter->second).size() << std::endl;
  }
}


void LoadExistingTable(std::ifstream& fin, StochasticPhenotypeTable* pt_it)
{
  std::string str;
  while (std::getline(fin, str)) {
    std::stringstream iss(str);
    int number;
    std::vector<uint8_t> phenotype_line;
    while (iss>>number)
      phenotype_line.push_back(static_cast<uint8_t>(number));
    Phenotype phen;
    phen.dx=phenotype_line[2];
    phen.dy=phenotype_line[3];
    phen.tiling=std::vector<uint8_t>(phenotype_line.begin()+4,phenotype_line.end());
    pt_it->known_phenotypes[phenotype_line[0]].emplace_back(phen);
  }
}
