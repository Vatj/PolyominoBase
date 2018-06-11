#include "genotype_phenotype.hpp"
#include <sstream>
#include <iterator>
#include <functional>
#include <set>

// deprecated functions

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

void GetPhenotypesIDs(const char* file_path_c, const char* file_name_c, uint8_t n_genes, uint8_t colours)
{
  std::string file_path(file_path_c),file_name(file_name_c),str,details="_N"+std::to_string(n_genes)+"_C"+std::to_string(colours)+".txt";
  std::ifstream file_in(file_path+file_name);
  std::ofstream gfout(file_path+"Genotype_Codes"+details, std::ios_base::out);
  // std::ofstream pfout(file_path+"Phenotype_Table"+details, std::ios_base::out);

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

  pt.PrintTable(file_path + "Phenotype_Table" + details);
}

// New pID set and PreProcessing functions

std::vector<Phenotype_ID> GetSetPIDs(Genotype& genotype, uint8_t k_builds, StochasticPhenotypeTable* pt_it)
{
  Clean_Genome(genotype, 0, false);
  double UND_frac = 0.1;

  std::vector<Phenotype_ID> pIDs = Stochastic::AssemblePlasticGenotype(genotype, k_builds, pt_it, UND_frac);

  std::sort(pIDs.begin(),pIDs.end());
  return pIDs;
}

void PreProcessSampled(std::vector<Genotype> genomes, Set_to_Genome& set_to_genome, StochasticPhenotypeTable* pt)
{
  uint8_t k_builds=100;
  Genotype genotype;
  std::vector<Phenotype_ID> pIDs;

  std::cout << "PreProcessing " <<+ genomes.size() << " genomes\n";

  #pragma omp parallel for schedule(dynamic) firstprivate(pIDs, genotype)
  for(uint64_t index=0; index < genomes.size(); index++)
  {
    genotype = genomes[index];
    pIDs = GetSetPIDs(genotype, k_builds, pt);

    if(index % 100 ==0)
      std::cout << "Currently preprocessing genome : " <<+ index << " out of " <<+ genomes.size() << "\n";

    #pragma omp critical
      set_to_genome[pIDs].emplace_back(genotype);
  }
}
