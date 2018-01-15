#include <algorithm>
#include <vector>
#include <iostream>
#include <utility>
#include <tuple>
#include <array>
#include <numeric>
#include <unordered_map>
#include <random>


		

namespace simulation_params
{
  typedef uint32_t population_size_type;
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
  extern std::mt19937_64 RNG_Engine;
  inline interface_type reverse_bits(interface_type v);
  inline uint8_t ArbitraryPopcount(interface_type face1);
  inline double SammingDistance(interface_type face1,interface_type face2);
  void MutateInterfaces(std::vector<interface_type>& binary_genome);
  //double SymmetryFactor(interface_type face1);

  /* ASSEMBLY */
  double ProteinAssemblyOutcome(std::vector<interface_type> binary_genome,uint8_t N_repeats,PhenotypeTable* pt);
  std::vector<int8_t> AssembleProtein(const std::vector<interface_type>& binary_genome);
  void PerimeterGrowth(int8_t x,int8_t y,int8_t theta,int8_t direction, int8_t tile_type,std::vector<int8_t>& growing_perimeter,std::vector<int8_t>& placed_tiles);

  

  struct PhenotypeTable {
    uint32_t n_phenotypes;
    std::unordered_map<uint8_t,std::vector<uint8_t> > known_phenotypes;
    std::unordered_map<uint8_t,std::vector<uint8_t> > undiscovered_phenotypes;
    std::unordered_map<uint8_t,std::vector<double> > phenotype_fitnesses{{0,{0}}};
    std::unordered_map<uint8_t,std::vector<uint32_t> > undiscovered_phenotype_xfer;
    std::vector<uint32_t> undiscovered_phenotype_counts;
    
    //std::vector<uint32_t> phenotype_size_count_accepted(model_params::unbound_factor*simulation_params::n_tiles*simulation_params::n_tiles);
    
    PhenotypeTable(void) : n_phenotypes(0) {};

    uint32_t PhenotypeCheck(std::vector<uint8_t>& phenotype, uint8_t dx, uint8_t dy) {
      if(phenotype.empty())
        return 0;
      uint32_t phenotype_ID=0;
      uint8_t phenotype_size=std::accumulate(phenotype.begin(),phenotype.end(),0);
      std::vector<uint8_t> temp_phenotype;
      for(std::vector<uint8_t>::iterator phen_iter = known_phenotypes[phenotype_size].begin();phen_iter!=known_phenotypes[phenotype_size].end();) {
        temp_phenotype.assign(phen_iter+2,phen_iter+2+*phen_iter* *(phen_iter+1));
        if(ComparePolyominoes(phenotype,dx,dy,temp_phenotype,*phen_iter,*(phen_iter+1)))
          return phenotype_ID;
        phen_iter+=*phen_iter* *(phen_iter+1)+2;
        ++phenotype_ID;
      }
      uint8_t new_phenotype_index=0;
      for(std::vector<uint8_t>::iterator phen_iter = undiscovered_phenotypes[phenotype_size].begin();phen_iter!=undiscovered_phenotypes[phenotype_size].end();) {
        temp_phenotype.assign(phen_iter+2,phen_iter+2+*phen_iter* *(phen_iter+1));
        if(ComparePolyominoes(phenotype,dx,dy,temp_phenotype,*phen_iter,*(phen_iter+1))) {
          if(++undiscovered_phenotype_counts[new_phenotype_index]>=model_params::UND_threshold*simulation_params::phenotype_builds) {
	    known_phenotypes[phenotype_size].insert(known_phenotypes[phenotype_size].end(),phen_iter,phen_iter+2+*phen_iter* *(phen_iter+1));
	    phenotype_fitnesses[phenotype_size].emplace_back(model_params::real_dist(RNG_Engine));
	    undiscovered_phenotype_xfer[phenotype_size].emplace_back(new_phenotype_index);
	    undiscovered_phenotype_xfer[phenotype_size].emplace_back(known_phenotypes[phenotype_size].size()-1);
	    return known_phenotypes[phenotype_size].size()-1+simulation_params::phenotype_builds;
	  }
	  else {
	    return known_phenotypes[phenotype_size].size()+new_phenotype_index+simulation_params::phenotype_builds;
	  }
	  
	}
        phen_iter+=*phen_iter* *(phen_iter+1)+2;
        ++new_phenotype_index;
      }
      
      undiscovered_phenotypes[phenotype_size].emplace_back(dx);
      undiscovered_phenotypes[phenotype_size].emplace_back(dy);
      undiscovered_phenotypes[phenotype_size].insert(undiscovered_phenotypes[phenotype_size].end(),phenotype.begin(),phenotype.end());
      undiscovered_phenotype_counts.emplace_back(1);

      /*
      known_phenotypes[phenotype_size].emplace_back(dx);
      known_phenotypes[phenotype_size].emplace_back(dy);
      known_phenotypes[phenotype_size].insert(known_phenotypes[phenotype_size].end(),phenotype.begin(),phenotype.end());
      phenotype_fitnesses[phenotype_size].emplace_back(model_params::real_dist(RNG_Engine));
      ++n_phenotypes;
      return phenotype_ID;
      */
      //++phenotype_size_count_active[phenotype_size];
      ++n_phenotypes;
      return known_phenotypes[phenotype_size].size()+new_phenotype_index-1+simulation_params::phenotype_builds;
    }
    
    double GenotypeFitness(std::vector<std::pair<uint8_t, uint32_t> >& phenotype_IDs) {
      for(std::unordered_map<uint8_t,std::vector<uint32_t> >::iterator xfer_iter=undiscovered_phenotype_xfer.begin();xfer_iter!=undiscovered_phenotype_xfer.end();++xfer_iter) {
	for(std::vector<uint32_t>::iterator replacing_iter=xfer_iter->second.begin();replacing_iter!=xfer_iter->second.end();replacing_iter+=2) {
	std::replace(phenotype_IDs.begin(),phenotype_IDs.end(),std::make_pair(xfer_iter->first,*(replacing_iter)),std::make_pair(xfer_iter->first,*(replacing_iter+1)));
	  }
      }
      

      std::unordered_map<uint8_t, std::unordered_map<uint32_t,uint8_t> > ID_counter;
      for(std::vector<std::pair<uint8_t,uint32_t> >::const_iterator ID_iter = phenotype_IDs.begin(); ID_iter!=phenotype_IDs.end(); ++ID_iter)
        ++ID_counter[ID_iter->first][ID_iter->second];
      double fitness=0;
      
      for(std::unordered_map<uint8_t, std::unordered_map<uint32_t,uint8_t> >::iterator size_iter =ID_counter.begin();size_iter!=ID_counter.end();++size_iter) 
        for(std::unordered_map<uint32_t,uint8_t>::iterator frequency_iter =size_iter->second.begin();frequency_iter!=size_iter->second.end();++frequency_iter)
	  if(frequency_iter->first < known_phenotypes[size_iter->first].size()) 
	    fitness+=phenotype_fitnesses[size_iter->first][frequency_iter->first] * std::pow(static_cast<double>(frequency_iter->second)/phenotype_IDs.size(),model_params::fitness_factor);


      //temporary_phenotype_indices.clear();
      undiscovered_phenotypes.clear();
      undiscovered_phenotype_counts.clear();
      undiscovered_phenotype_xfer.clear();
      
      return fitness;
    }
      

  };
}
uint8_t PhenotypeSymmetryFactor(std::vector<uint8_t>& original_shape, uint8_t dx, uint8_t dy);
void DistributionStatistics(std::vector<double>& intf, double& mean, double& variance);

std::vector<uint16_t> InterfaceStrengths(std::vector<interface_model::interface_type>& interfaces);


