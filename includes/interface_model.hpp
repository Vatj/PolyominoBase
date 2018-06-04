#include "core_phenotype.hpp"

#include <functional>
#include <random>
#include <climits>
#include <cmath>

#include <set>
#include <array>

#include <iostream>


typedef uint64_t interface_type;
typedef std::vector<interface_type> BGenotype;
typedef std::pair<uint8_t,uint8_t> interaction_pair;

namespace simulation_params
{
  extern uint8_t n_tiles,phenotype_builds;
  extern uint16_t population_size;
  extern uint32_t generation_limit,independent_trials,run_offset;

  extern bool random_initilisation;
}

namespace model_params
{
  constexpr uint8_t interface_size=CHAR_BIT*sizeof(interface_type);
  
  extern double temperature,binding_threshold,mu_prob,fitness_factor,UND_threshold,interface_threshold;
  extern bool fixed_seed;
  
  extern std::binomial_distribution<uint8_t> b_dist;
  extern std::uniform_real_distribution<double> real_dist;
  extern std::array<double,model_params::interface_size+1> binding_probabilities;
}

std::array<double,model_params::interface_size+1> GenBindingProbsLUP();



/* SPATIAL */
Phenotype SpatialGrid(std::vector<int8_t>& placed_tiles);


/* PRINTING */
void PrintShape(Phenotype phen);

uint8_t PhenotypeSymmetryFactor(std::vector<uint8_t>& original_shape, uint8_t dx, uint8_t dy);
void DistributionStatistics(std::vector<double>& intf, double& mean, double& variance);
void InterfaceStrengths(BGenotype& interfaces, std::vector<uint32_t>& strengths);


namespace interface_model
{  
  struct InterfacePhenotypeTable;

  extern std::mt19937 RNG_Engine;
  constexpr uint8_t GAUGE=4;
  interface_type reverse_bits(interface_type v);
  uint8_t ArbitraryPopcount(interface_type face1);
  uint8_t SammingDistance(interface_type face1,interface_type face2);
  void MutateInterfaces(BGenotype& binary_genome);

  /* ASSEMBLY */
  double ProteinAssemblyOutcome(BGenotype binary_genome, InterfacePhenotypeTable* pt,Phenotype_ID& pid,std::vector<interaction_pair>& pid_interactions);
  std::vector<int8_t> AssembleProtein(const BGenotype& binary_genome,std::set<interaction_pair>& interacting_indices);
  void PerimeterGrowth(int8_t x,int8_t y,int8_t theta,int8_t direction, int8_t tile_type,std::vector<int8_t>& growing_perimeter,std::vector<int8_t>& placed_tiles);

  std::vector<int8_t> AssembleProteinNew(const BGenotype& binary_genome,std::set<interaction_pair>& interacting_indices);
  void ExtendPerimeter(const BGenotype& binary_genome,uint8_t tile_detail, int8_t x,int8_t y, std::vector<int8_t>& placed_tiles,std::vector<int8_t>& potential_sites,std::vector<double>& binding_strengths,std::vector<interaction_pair>& interaction_pairs);
  


  

  struct InterfacePhenotypeTable : PhenotypeTable {

    std::unordered_map<uint8_t,std::vector<Phenotype> > undiscovered_phenotypes;
    std::unordered_map<uint8_t,std::vector<double> > phenotype_fitnesses{{0,{0}}};
    std::unordered_map<uint8_t,std::vector<uint16_t> > new_phenotype_xfer;
    std::vector<uint16_t> undiscovered_phenotype_counts;

    uint16_t PhenotypeCheck(Phenotype& phen) {
      uint8_t phenotype_size=std::count_if(phen.tiling.begin(),phen.tiling.end(),[](const int c){return c != 0;});
      for(uint16_t phenotype_index=0; phenotype_index != known_phenotypes[phenotype_size].size();++phenotype_index) {
        if(ComparePolyominoes(phen,known_phenotypes[phenotype_size][phenotype_index],GAUGE)) 
	  return phenotype_index;
      }
      uint8_t new_phenotype_index=0;
      for(Phenotype phen_p : undiscovered_phenotypes[phenotype_size]) {
        if(ComparePolyominoes(phen,phen_p,GAUGE)) {
          if(++undiscovered_phenotype_counts[new_phenotype_index]>=ceil(model_params::UND_threshold*simulation_params::phenotype_builds)) {
            new_phenotype_xfer[phenotype_size].emplace_back(phenotype_fitnesses[phenotype_size].size()+new_phenotype_index+simulation_params::phenotype_builds);
            known_phenotypes[phenotype_size].push_back(phen);
            //std::gamma_distribution<double> fitness_dist(sqrt(static_cast<double>(phenotype_size)),1);
            std::gamma_distribution<double> fitness_dist(sqrt(phenotype_size),1);
            phenotype_fitnesses[phenotype_size].emplace_back(fitness_dist(RNG_Engine));
            new_phenotype_xfer[phenotype_size].emplace_back(phenotype_fitnesses[phenotype_size].size()-1);
            
            return phenotype_fitnesses[phenotype_size].size()-1;
          }
          else
            return phenotype_fitnesses[phenotype_size].size()+new_phenotype_index+simulation_params::phenotype_builds;
        }
        ++new_phenotype_index;
      }
      
      MinimizePhenRep(phen.tiling,GAUGE);
      undiscovered_phenotypes[phenotype_size].emplace_back(phen);
      undiscovered_phenotype_counts.emplace_back(1);
      return phenotype_fitnesses[phenotype_size].size()+new_phenotype_index+simulation_params::phenotype_builds;
    }

    /* Replace previously undiscovered phenotype IDs with new phenotype ID */
    void RelabelPhenotypes(std::vector<Phenotype_ID >& pids,std::map<Phenotype_ID, std::map<interaction_pair, uint8_t> >& p_ints) { 
      for(std::unordered_map<uint8_t,std::vector<uint16_t> >::iterator x_iter=new_phenotype_xfer.begin();x_iter!=new_phenotype_xfer.end();++x_iter) 
        for(std::vector<uint16_t>::iterator r_iter=x_iter->second.begin();r_iter!=x_iter->second.end();r_iter+=2) {
          std::replace(pids.begin(),pids.end(),std::make_pair(x_iter->first,*(r_iter)),std::make_pair(x_iter->first,*(r_iter+1)));
	  for(auto& imap :p_ints[std::make_pair(x_iter->first,*(r_iter))])
	    p_ints[std::make_pair(x_iter->first,*(r_iter+1))][imap.first]+=imap.second;
	}
      undiscovered_phenotypes.clear();
      undiscovered_phenotype_counts.clear();
      new_phenotype_xfer.clear();

    }
    
    /* Count each ID frequency */
    std::map<Phenotype_ID,uint8_t> PhenotypeFrequencies(std::vector<Phenotype_ID >& pids) {
      std::map<Phenotype_ID, uint8_t> ID_counter;
      for(std::vector<Phenotype_ID >::const_iterator ID_iter = pids.begin(); ID_iter!=pids.end(); ++ID_iter)
        if(ID_iter->second < phenotype_fitnesses[ID_iter->first].size())
          ++ID_counter[std::make_pair(ID_iter->first,ID_iter->second)];
      return ID_counter;
    }

    /* Add fitness contribution from each phenotype */
    double GenotypeFitness(std::map<Phenotype_ID,uint8_t> ID_counter) {
      double fitness=0; 
      for(std::map<Phenotype_ID,uint8_t >::const_iterator p_it =ID_counter.begin();p_it!=ID_counter.end();++p_it)
        if(p_it->second >= ceil(model_params::UND_threshold*simulation_params::phenotype_builds))
	    fitness+=phenotype_fitnesses[p_it->first.first][p_it->first.second] * std::pow(static_cast<double>(p_it->second)/simulation_params::phenotype_builds,model_params::fitness_factor);     
      return fitness;
    }

    void ReassignFitness() {
      for(std::unordered_map<uint8_t,std::vector<double> >::iterator fit_iter=phenotype_fitnesses.begin();fit_iter!=phenotype_fitnesses.end();++fit_iter) {
	if(fit_iter->first) {
	  std::gamma_distribution<double> fitness_dist(sqrt(fit_iter->first),1);
	  for(double& fitness : fit_iter->second)
	    fitness=fitness_dist(RNG_Engine);
	}
      }
    }

  };
}





