#include "core_phenotype.hpp"

#include <iostream>
#include <functional>
#include <tuple>
#include <numeric>
#include <map>
#include <random>
#include <climits>
#include <cmath>







namespace simulation_params
{
  extern uint8_t n_tiles,phenotype_builds;
  extern uint16_t population_size;
  extern uint32_t generation_limit,independent_trials,run_offset;

  extern bool fitness_selection,random_initilisation;
}

namespace model_params
{
  extern double temperature,mu_prob,misbinding_rate,fitness_factor,unbound_factor,UND_threshold;
  extern const uint8_t interface_size;

  extern std::binomial_distribution<uint8_t> b_dist;
  extern std::uniform_real_distribution<double> real_dist;

}

typedef uint8_t interface_type;

/* SPATIAL */
Phenotype SpatialGrid(std::vector<int8_t>& placed_tiles);


/* PRINTING */
void PrintShape(Phenotype phen);

void MinimalTilingRepresentation(std::vector<uint8_t>& tiling);
uint8_t PhenotypeSymmetryFactor(std::vector<uint8_t>& original_shape, uint8_t dx, uint8_t dy);
void DistributionStatistics(std::vector<double>& intf, double& mean, double& variance);
void InterfaceStrengths(std::vector<interface_type>& interfaces, std::vector<uint32_t>& strengths);


namespace interface_model
{
  
  
  struct InterfacePhenotypeTable;

  extern std::random_device rd;
  extern std::mt19937 RNG_Engine;
  interface_type reverse_bits(interface_type v);
  uint8_t ArbitraryPopcount(interface_type face1);
  uint8_t SammingDistance(interface_type face1,interface_type face2);
  void MutateInterfaces(std::vector<interface_type>& binary_genome);

  /* ASSEMBLY */
  double ProteinAssemblyOutcome(std::vector<interface_type> binary_genome, InterfacePhenotypeTable* pt,Phenotype_ID& pid);
  std::vector<int8_t> AssembleProtein(const std::vector<interface_type>& binary_genome);
  void PerimeterGrowth(int8_t x,int8_t y,int8_t theta,int8_t direction, int8_t tile_type,std::vector<int8_t>& growing_perimeter,std::vector<int8_t>& placed_tiles);


  

  struct InterfacePhenotypeTable : PhenotypeTable {

    std::unordered_map<uint8_t,std::vector<Phenotype> > undiscovered_phenotypes;
    std::unordered_map<uint8_t,std::vector<double> > phenotype_fitnesses{{0,{0}}};
    std::unordered_map<uint8_t,std::vector<uint16_t> > new_phenotype_xfer;
    std::vector<uint16_t> undiscovered_phenotype_counts;

    uint16_t PhenotypeCheck(Phenotype& phen) {
      uint8_t phenotype_size=std::count_if(phen.tiling.begin(),phen.tiling.end(),[](const int c){return c != 0;});
      for(uint16_t phenotype_index=0; phenotype_index != known_phenotypes[phenotype_size].size();++phenotype_index) {
        if(ComparePolyominoes(phen,known_phenotypes[phenotype_size][phenotype_index])) 
	  return phenotype_index;
      }
      uint8_t new_phenotype_index=0;
      for(Phenotype phen_p : undiscovered_phenotypes[phenotype_size]) {
        if(ComparePolyominoes(phen,phen_p)) {
          if(++undiscovered_phenotype_counts[new_phenotype_index]>=ceil(model_params::UND_threshold*simulation_params::phenotype_builds)) {
            new_phenotype_xfer[phenotype_size].emplace_back(phenotype_fitnesses[phenotype_size].size()+new_phenotype_index+simulation_params::phenotype_builds);
            known_phenotypes[phenotype_size].push_back(phen);
            //std::gamma_distribution<double> fitness_dist(sqrt(static_cast<double>(phenotype_size)),1);
            std::gamma_distribution<double> fitness_dist(pow(static_cast<double>(phenotype_size),2),1);
            phenotype_fitnesses[phenotype_size].emplace_back(fitness_dist(RNG_Engine));
            new_phenotype_xfer[phenotype_size].emplace_back(phenotype_fitnesses[phenotype_size].size()-1);
            
            return phenotype_fitnesses[phenotype_size].size()-1;
          }
          else
            return phenotype_fitnesses[phenotype_size].size()+new_phenotype_index+simulation_params::phenotype_builds;
        }
        ++new_phenotype_index;
      }
      
      MinimalTilingRepresentation(phen.tiling);
      undiscovered_phenotypes[phenotype_size].emplace_back(phen);
      undiscovered_phenotype_counts.emplace_back(1);
      return phenotype_fitnesses[phenotype_size].size()+new_phenotype_index+simulation_params::phenotype_builds;
    }
    
    std::map<Phenotype_ID,uint8_t> PhenotypeFrequencies(std::vector<Phenotype_ID >& Phenotype_IDs) {
      /* Replace previously undiscovered phenotype IDs with new one */
      for(std::unordered_map<uint8_t,std::vector<uint16_t> >::iterator xfer_iter=new_phenotype_xfer.begin();xfer_iter!=new_phenotype_xfer.end();++xfer_iter) 
        for(std::vector<uint16_t>::iterator rep_iter=xfer_iter->second.begin();rep_iter!=xfer_iter->second.end();rep_iter+=2) 
          std::replace(Phenotype_IDs.begin(),Phenotype_IDs.end(),std::make_pair(xfer_iter->first,*(rep_iter)),std::make_pair(xfer_iter->first,*(rep_iter+1)));
      /* Count each ID frequency */
      std::map<Phenotype_ID, uint8_t> ID_counter;
      for(std::vector<Phenotype_ID >::const_iterator ID_iter = Phenotype_IDs.begin(); ID_iter!=Phenotype_IDs.end(); ++ID_iter)
        if(ID_iter->second < phenotype_fitnesses[ID_iter->first].size())
          ++ID_counter[std::make_pair(ID_iter->first,ID_iter->second)];
      
      undiscovered_phenotypes.clear();
      undiscovered_phenotype_counts.clear();
      new_phenotype_xfer.clear();
      return ID_counter;

    }
    double GenotypeFitness(std::map<Phenotype_ID,uint8_t> ID_counter) {
      double fitness=0;
      /* Add fitness contribution from each phenotype */
      for(std::map<Phenotype_ID,uint8_t >::const_iterator p_it =ID_counter.begin();p_it!=ID_counter.end();++p_it)
        if(p_it->second >= ceil(model_params::UND_threshold*simulation_params::phenotype_builds))
	    fitness+=phenotype_fitnesses[p_it->first.first][p_it->first.second] * std::pow(static_cast<double>(p_it->second)/simulation_params::phenotype_builds,model_params::fitness_factor);     
      return fitness;
    }

    void ReassignFitness() {
      for(std::unordered_map<uint8_t,std::vector<double> >::iterator fit_iter=phenotype_fitnesses.begin();fit_iter!=phenotype_fitnesses.end();++fit_iter) {
	if(fit_iter->first) {
	  std::gamma_distribution<double> fitness_dist(sqrt(static_cast<double>(fit_iter->first)),1);
	  for(double& fitness : fit_iter->second)
	    fitness=fitness_dist(RNG_Engine);
	}
      }
    }

  };
}





