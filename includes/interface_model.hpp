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
  extern uint8_t interface_size;

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
    uint16_t n_phenotypes;
    std::vector<uint8_t> known_phenotypes;
    std::vector<double> phenotype_fitnesses;

    PhenotypeTable(void) : n_phenotypes(0) {known_phenotypes.reserve(1000);phenotype_fitnesses.reserve(1000); phenotype_fitnesses.emplace_back(0);};


    uint16_t PhenotypeCheck(std::vector<uint8_t>& phenotype, uint8_t dx, uint8_t dy) {
      int phenotype_ID=1;
      std::vector<uint8_t> temp_phenotype;
      for(std::vector<uint8_t>::iterator phen_iter = known_phenotypes.begin();phen_iter!=known_phenotypes.end();) {
        temp_phenotype.assign(phen_iter+2,phen_iter+2+*phen_iter* *(phen_iter+1));
        if(ComparePolyominoes(phenotype,dx,dy,temp_phenotype,*phen_iter,*(phen_iter+1)))
          return phenotype_ID;
        phen_iter+=*phen_iter* *(phen_iter+1)+2;
        ++phenotype_ID;
      }      
      known_phenotypes.emplace_back(dx);
      known_phenotypes.emplace_back(dy);
      
      known_phenotypes.insert(known_phenotypes.end(),phenotype.begin(),phenotype.end());
      phenotype_fitnesses.emplace_back(model_params::real_dist(RNG_Engine));
      ++n_phenotypes;
      return phenotype_ID;
    }
    double GenotypeFitness(std::vector<uint8_t>& phenotype_IDs) {
      //std::cout<<"IDs :";
      std::unordered_map<uint8_t,uint16_t> ID_counter;
      for(std::vector<uint8_t>::const_iterator ID_iter = phenotype_IDs.begin(); ID_iter!=phenotype_IDs.end(); ++ID_iter)
        ++ID_counter[*ID_iter];
      
      if(ID_counter.size()==1) 
        return phenotype_fitnesses[ID_counter.begin()->first];
      
      double fitness=0;//,fitness_normalisation=0;
      /*
      std::cout<<"IDs ";
      for(std::unordered_map<uint8_t,uint16_t>::iterator frequency_iter =ID_counter.begin();frequency_iter!=ID_counter.end();++frequency_iter) {
        std::cout<<"ID "<<+frequency_iter->first<<" count "<<+frequency_iter->second<<" fitness "<<phenotype_fitnesses[frequency_iter->first]<<std::endl;
 
      }
      std::cout<<"\n";
      */
      for(std::unordered_map<uint8_t,uint16_t>::iterator frequency_iter =ID_counter.begin();frequency_iter!=ID_counter.end();++frequency_iter) {
        fitness+=phenotype_fitnesses[frequency_iter->first] * std::pow(static_cast<double>(frequency_iter->second)/phenotype_IDs.size(),model_params::fitness_factor);
        
        //std::cout<<"F "<<frequency_iter->second<<" "<<phenotype_IDs.size()<<" "<<model_params::fitness_factor<<" = "<<std::pow(static_cast<double>(frequency_iter->second)/phenotype_IDs.size(),model_params::fitness_factor)<<std::endl;
	//fitness+=phenotype_fitnesses[frequency_iter->first] * std::log(1-static_cast<double>(frequency_iter->second)/phenotype_IDs.size());//frequency_iter->second /phenotype_IDs.size();
        //fitness_normalisation+=std::log(1-static_cast<double>(frequency_iter->second) /phenotype_IDs.size());
        //std::cout<<"Norm "<<std::log(1-static_cast<double>(frequency_iter->second) /phenotype_IDs.size())<<std::endl;
	//std::cout<<"ID "<<+frequency_iter->first<<" count "<<+frequency_iter->second<<std::endl;
      }    

	
      //for(std::vector<uint8_t>::const_iterator ID = phenotype_IDs.begin(); ID != phenotype_IDs.end(); ++ID) {
      //fitness+=phenotype_fitnesses[*ID];
      //	std::cout<<+*ID<<" ";
      //}
      //std::cout<<"\n";
      //std::cout<<fitness_normalisation<<", "<<fitness<<std::endl;
      return fitness;//fitness_normalisation;
    }
      

  };
}

void DistributionStatistics(std::vector<double>& intf, double& mean, double& variance);

std::vector<double> InterfaceStrengths(std::vector<interface_model::interface_type>& interfaces);
