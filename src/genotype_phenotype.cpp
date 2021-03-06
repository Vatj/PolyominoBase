#include "genotype_phenotype.hpp"
#include <sstream>
#include <iterator>
#include <functional>
#include <set>


// New pID set and PreProcessing functions

std::vector<Phenotype_ID> GetSetPIDs(Genotype genotype, PhenotypeTable* pt_it)
{
  Clean_Genome(genotype, false);
  std::vector<Phenotype_ID> pIDs = Stochastic::AssemblePlasticGenotype(genotype, pt_it);
  std::sort(pIDs.begin(), pIDs.end());
  return pIDs;
}

std::map<Phenotype_ID, uint8_t> GetPIDCounter(Genotype genotype, PhenotypeTable* pt_it)
{
  Clean_Genome(genotype, false);
  std::map<Phenotype_ID, uint8_t> pID_counter = Stochastic::AssemblePlasticGenotypeFrequency(genotype, pt_it);
  return pID_counter;
}


void PreProcessSampled(std::vector<Genotype> genomes, Set_to_Genome& set_to_genome, PhenotypeTable* pt)
{
  Genotype genotype;
  std::vector<Phenotype_ID> pIDs;
  uint8_t save_config_builds = simulation_params::phenotype_builds;
  simulation_params::phenotype_builds = simulation_params::preprocess_builds;

  std::cout << "PreProcessing " <<+ genomes.size() << " genomes, building ";
  std::cout <<+ simulation_params::preprocess_builds << "th times\n";

  #pragma omp parallel for schedule(dynamic) firstprivate(pIDs, genotype)
  for(uint64_t index=0; index < genomes.size(); index++)
  {
    genotype = genomes[index];
    pIDs = GetSetPIDs(genotype, pt);

    if(index % 100 ==0)
      std::cout << "Currently preprocessing genome : " <<+ index << " out of " <<+ genomes.size() << "\n";

    #pragma omp critical
      set_to_genome[pIDs].emplace_back(genotype);
  }
  simulation_params::phenotype_builds = save_config_builds;
  std::cout << "Set back build parameters to config value : " <<+ simulation_params::phenotype_builds << std::endl;
}

void FilterExhaustive(std::vector<Genotype> genomes, PhenotypeTable* pt)
{
  Genotype genotype;
  std::vector<Genotype> new_genomes;
  std::vector<Phenotype_ID> pIDs;
  Phenotype_ID rare_pID = {0, 0}, unbound_pID = {255, 0};

  std::cout << "Threshold is : " << (ceil(simulation_params::phenotype_builds * simulation_params::UND_threshold));
  std::cout << " out of " <<+ simulation_params::phenotype_builds << " builds \n";

  #pragma omp parallel for schedule(dynamic) firstprivate(pIDs, genotype)
  for(uint64_t index=0; index < genomes.size(); index++)
  {
    genotype = genomes[index];
    pIDs = GetSetPIDs(genotype, pt);
    if(pIDs.front() == rare_pID || pIDs.back() == unbound_pID)
      continue;
    else
      new_genomes.emplace_back(genotype);

    if(index % 100 ==0)
      std::cout << "Currently filtering genome : " <<+ index << " out of " <<+ genomes.size() << "\n";
  }
  genomes = new_genomes;
}
