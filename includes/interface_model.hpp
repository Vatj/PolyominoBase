#include "xorshift.hpp"


#include <algorithm>
#include <vector>
#include <iostream>
#include <utility>
#include <tuple>
#include <array>
#include <numeric>
#include <unordered_map>


		


namespace model_params
{
  extern double temperature,mu_prob,misbinding_rate,fitness_factor,unbound_factor;
  extern const uint8_t interface_size;

  extern std::binomial_distribution<uint8_t> b_dist;
  extern std::uniform_real_distribution<double> real_dist;

}

namespace interface_model
{
  typedef uint16_t interface_type;
  struct PhenotypeTable;

  
  extern xorshift RNG_Engine;
  inline interface_type reverse_bits(interface_type v);
  inline uint8_t ArbitraryPopcount(interface_type face1);
  double SammingDistance(interface_type face1,interface_type face2);
  void MutateInterfaces(std::vector<interface_type>& binary_genome);
  double SymmetryFactor(interface_type face1);

  /* ASSEMBLY */
  double ProteinAssemblyOutcome(std::vector<interface_type> binary_genome,uint8_t N_repeats,PhenotypeTable* pt);
  std::vector<int8_t> AssembleProtein(const std::vector<interface_type>& binary_genome);
  void PerimeterGrowth(int8_t x,int8_t y,int8_t theta,int8_t direction, int8_t tile_type,std::vector<int8_t>& growing_perimeter,std::vector<int8_t>& placed_tiles);

  /* SPATIAL */
  std::vector<uint8_t> SpatialGrid(std::vector<int8_t>& placed_tiles, uint8_t& dx,uint8_t& dy);
  bool ComparePolyominoes(std::vector<uint8_t>& phenotype,uint8_t& dx,uint8_t& dy, std::vector<uint8_t>& phenotype_prime,uint8_t dx_prime,uint8_t dy_prime);

  /* ROTATIONS */
  void ClockwiseRotation(std::vector<uint8_t>& phenotype,uint8_t& dx,uint8_t& dy);
  std::vector<uint8_t> ClockwisePiRotation(std::vector<uint8_t>& phenotype);

  /* PRINTING */
  void PrintShape(std::vector<uint8_t>& spatial_information,uint8_t dx,uint8_t dy);

  struct PhenotypeTable {
    uint32_t n_phenotypes;
    std::unordered_map<uint8_t,std::vector<uint8_t> > known_phenotypes;
    std::unordered_map<uint8_t,std::vector<double> > phenotype_fitnesses{{0,{0}}};

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
      known_phenotypes[phenotype_size].emplace_back(dx);
      known_phenotypes[phenotype_size].emplace_back(dy);
      
      known_phenotypes[phenotype_size].insert(known_phenotypes[phenotype_size].end(),phenotype.begin(),phenotype.end());
      phenotype_fitnesses[phenotype_size].emplace_back(model_params::real_dist(RNG_Engine));
      ++n_phenotypes;
      return phenotype_ID;
    }
    
    double GenotypeFitness(std::vector<std::pair<uint8_t, uint32_t> >& phenotype_IDs) {
      std::unordered_map<uint8_t, std::unordered_map<uint32_t,uint8_t> > ID_counter;
      for(std::vector<std::pair<uint8_t,uint32_t> >::const_iterator ID_iter = phenotype_IDs.begin(); ID_iter!=phenotype_IDs.end(); ++ID_iter)
        ++ID_counter[ID_iter->first][ID_iter->second];
      double fitness=0;
      for(std::unordered_map<uint8_t, std::unordered_map<uint32_t,uint8_t> >::iterator size_iter =ID_counter.begin();size_iter!=ID_counter.end();++size_iter) 
        for(std::unordered_map<uint32_t,uint8_t>::iterator frequency_iter =size_iter->second.begin();frequency_iter!=size_iter->second.end();++frequency_iter) 
          fitness+=phenotype_fitnesses[size_iter->first][frequency_iter->first] * std::pow(static_cast<double>(frequency_iter->second)/phenotype_IDs.size(),model_params::fitness_factor);
      return fitness;
    }
      

  };
}

void DistributionStatistics(std::vector<double>& intf, double& mean, double& variance);

std::vector<uint16_t> InterfaceStrengths(std::vector<interface_model::interface_type>& interfaces);
