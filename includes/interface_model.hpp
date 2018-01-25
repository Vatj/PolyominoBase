#include <algorithm>
#include <vector>
#include <iostream>
#include <utility>
#include <tuple>
#include <array>
#include <numeric>
#include <unordered_map>
#include <random>
#include <climits>
#include <cstdint>
#include <cmath>

namespace simulation_params
{
  typedef uint8_t population_size_type;
  extern population_size_type population_size;
  extern uint32_t generation_limit,independent_trials,run_offset;
  extern uint8_t n_tiles,phenotype_builds;
  extern bool fitness_selection;
}

namespace model_params
{
  extern double temperature,mu_prob,misbinding_rate,fitness_factor,unbound_factor,UND_threshold;
  extern const uint8_t interface_size;

  extern std::binomial_distribution<uint8_t> b_dist;
  extern std::uniform_real_distribution<double> real_dist;

}
/* SPATIAL */
std::vector<uint8_t> SpatialGrid(std::vector<int8_t>& placed_tiles, uint8_t& dx,uint8_t& dy);
bool ComparePolyominoes(std::vector<uint8_t>& phenotype,uint8_t& dx,uint8_t& dy, std::vector<uint8_t>& phenotype_prime,uint8_t dx_prime,uint8_t dy_prime);

/* ROTATIONS */
void ClockwiseRotation(std::vector<uint8_t>& phenotype,uint8_t& dx,uint8_t& dy);
std::vector<uint8_t> ClockwisePiRotation(std::vector<uint8_t>& phenotype);

/* PRINTING */
void PrintShape(std::vector<uint8_t>& spatial_information,uint8_t dx,uint8_t dy);

namespace interface_model
{
  typedef uint32_t interface_type;
  struct PhenotypeTable;

  extern std::random_device rd;
  extern std::mt19937 RNG_Engine;
  inline interface_type reverse_bits(interface_type v);
  inline uint8_t ArbitraryPopcount(interface_type face1);
  inline double SammingDistance(interface_type face1,interface_type face2);
  void MutateInterfaces(std::vector<interface_type>& binary_genome);
  //double SymmetryFactor(interface_type face1);

  /* ASSEMBLY */
  double ProteinAssemblyOutcome(std::vector<interface_type> binary_genome, PhenotypeTable* pt);
  std::vector<int8_t> AssembleProtein(const std::vector<interface_type>& binary_genome);
  void PerimeterGrowth(int8_t x,int8_t y,int8_t theta,int8_t direction, int8_t tile_type,std::vector<int8_t>& growing_perimeter,std::vector<int8_t>& placed_tiles);

  

  struct PhenotypeTable {
    //uint32_t n_phenotypes;
    std::unordered_map<uint8_t,std::vector<uint8_t> > known_phenotypes;
    std::unordered_map<uint8_t,std::vector<uint8_t> > undiscovered_phenotypes;
    std::unordered_map<uint8_t,std::vector<double> > phenotype_fitnesses{{0,{0}}};
    std::unordered_map<uint8_t,std::vector<uint32_t> > new_phenotype_xfer;
    std::vector<uint32_t> undiscovered_phenotype_counts;
    
    //std::vector<uint32_t> phenotype_size_count_accepted(model_params::unbound_factor*simulation_params::n_tiles*simulation_params::n_tiles);
    
    //PhenotypeTable(void) : n_phenotypes(0) {};

    uint32_t PhenotypeCheck(std::vector<uint8_t>& phenotype, uint8_t dx, uint8_t dy) {
      if(phenotype.empty())
        return 0;
      uint32_t phenotype_ID=0;
      uint8_t phenotype_size=std::accumulate(phenotype.begin(),phenotype.end(),0);
      std::vector<uint8_t> temp_phenotype;
      for(std::vector<uint8_t>::iterator phen_iter = known_phenotypes[phenotype_size].begin();phen_iter!=known_phenotypes[phenotype_size].end();) {
        temp_phenotype.assign(phen_iter+2,phen_iter+2+*phen_iter* *(phen_iter+1));
        if(ComparePolyominoes(phenotype,dx,dy,temp_phenotype,*phen_iter,*(phen_iter+1))) {
	  return phenotype_ID;
	}
        phen_iter+=*phen_iter* *(phen_iter+1)+2;
        ++phenotype_ID;
      }
      uint8_t new_phenotype_index=0;
      for(std::vector<uint8_t>::iterator phen_iter = undiscovered_phenotypes[phenotype_size].begin();phen_iter!=undiscovered_phenotypes[phenotype_size].end();) {
        temp_phenotype.assign(phen_iter+2,phen_iter+2+*phen_iter* *(phen_iter+1));
        if(ComparePolyominoes(phenotype,dx,dy,temp_phenotype,*phen_iter,*(phen_iter+1))) {
          if(++undiscovered_phenotype_counts[new_phenotype_index]>=model_params::UND_threshold*simulation_params::phenotype_builds) {
	    new_phenotype_xfer[phenotype_size].emplace_back(phenotype_fitnesses[phenotype_size].size()+new_phenotype_index+simulation_params::phenotype_builds);
	    known_phenotypes[phenotype_size].insert(known_phenotypes[phenotype_size].end(),phen_iter,phen_iter+2+*phen_iter* *(phen_iter+1));
	    std::gamma_distribution<double> fitness_dist(sqrt(static_cast<double>(phenotype_size)),1);
	    phenotype_fitnesses[phenotype_size].emplace_back(fitness_dist(RNG_Engine));
	    new_phenotype_xfer[phenotype_size].emplace_back(phenotype_fitnesses[phenotype_size].size()-1);
	    return phenotype_fitnesses[phenotype_size].size()-1;
	  }
	  else {
	    return phenotype_fitnesses[phenotype_size].size()+new_phenotype_index+simulation_params::phenotype_builds;
	  }
	  
	}
        phen_iter+=*phen_iter* *(phen_iter+1)+2;
        ++new_phenotype_index;
      }
      
      undiscovered_phenotypes[phenotype_size].emplace_back(dx);
      undiscovered_phenotypes[phenotype_size].emplace_back(dy);
      undiscovered_phenotypes[phenotype_size].insert(undiscovered_phenotypes[phenotype_size].end(),phenotype.begin(),phenotype.end());
      undiscovered_phenotype_counts.emplace_back(1);

      return phenotype_fitnesses[phenotype_size].size()+new_phenotype_index+simulation_params::phenotype_builds;
    }
    
    double GenotypeFitness(std::vector<std::pair<uint8_t, uint32_t> >& phenotype_IDs) {
      for(std::unordered_map<uint8_t,std::vector<uint32_t> >::iterator xfer_iter=new_phenotype_xfer.begin();xfer_iter!=new_phenotype_xfer.end();++xfer_iter) 
	for(std::vector<uint32_t>::iterator rep_iter=xfer_iter->second.begin();rep_iter!=xfer_iter->second.end();rep_iter+=2) 
	std::replace(phenotype_IDs.begin(),phenotype_IDs.end(),std::make_pair(xfer_iter->first,*(rep_iter)),std::make_pair(xfer_iter->first,*(rep_iter+1)));
	  
      std::unordered_map<uint8_t, std::unordered_map<uint32_t,uint8_t> > ID_counter;
      for(std::vector<std::pair<uint8_t,uint32_t> >::const_iterator ID_iter = phenotype_IDs.begin(); ID_iter!=phenotype_IDs.end(); ++ID_iter)
        ++ID_counter[ID_iter->first][ID_iter->second];
      double fitness=0;
      
      for(std::unordered_map<uint8_t, std::unordered_map<uint32_t,uint8_t> >::const_iterator size_iter =ID_counter.begin();size_iter!=ID_counter.end();++size_iter) 
        for(std::unordered_map<uint32_t,uint8_t>::const_iterator f_iter =size_iter->second.begin();f_iter!=size_iter->second.end();++f_iter)
	  if(f_iter->first < phenotype_fitnesses[size_iter->first].size() && f_iter->second >= model_params::UND_threshold*simulation_params::phenotype_builds)
	    fitness+=phenotype_fitnesses[size_iter->first][f_iter->first] * std::pow(static_cast<double>(f_iter->second)/phenotype_IDs.size(),model_params::fitness_factor);
           
      undiscovered_phenotypes.clear();
      undiscovered_phenotype_counts.clear();
      new_phenotype_xfer.clear();
      
      return fitness;
    }

    void ReassignFitness() {
      for(std::unordered_map<uint8_t,std::vector<double> >::iterator fit_iter=phenotype_fitnesses.begin();fit_iter!=phenotype_fitnesses.end();++fit_iter) {
	if(fit_iter->first) {
	  std::gamma_distribution<double> fitness_dist(sqrt(static_cast<double>(fit_iter->first)),1);
	  for(double& fitness :fit_iter->second)
	    fitness=fitness_dist(RNG_Engine);
	}
      }
    }

  };
}

uint8_t PhenotypeSymmetryFactor(std::vector<uint8_t>& original_shape, uint8_t dx, uint8_t dy);
void DistributionStatistics(std::vector<double>& intf, double& mean, double& variance);
std::vector<uint16_t> InterfaceStrengths(std::vector<interface_model::interface_type>& interfaces);


