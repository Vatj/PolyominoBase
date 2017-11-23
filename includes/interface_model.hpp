#include "xorshift.hpp"


#include <algorithm>
#include <vector>
#include <iostream>
#include <utility>
#include <tuple>
#include <array>
#include <numeric>

		


namespace params
{
  extern double temperature;
  extern uint8_t interface_size;
  //extern std::array<uint8_t, 16> interface_indices;
  //extern std::array<uint8_t, 4> faces;
  extern std::binomial_distribution<uint8_t> b_dist;
  extern std::uniform_real_distribution<double> real_dist;

}

namespace interface_model
{
  typedef uint16_t interface_type;
  struct PhenotypeTable;

  
  extern xorshift RNG_Engine;
  inline interface_type reverse_bits(interface_type v);
  uint8_t SammingDistance(uint16_t face1,uint16_t face2);
  void MutateInterfaces(std::vector<interface_type>& binary_genome);

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

    PhenotypeTable(void) : n_phenotypes(0) {known_phenotypes.reserve(1000);phenotype_fitnesses.reserve(1000);};

    uint16_t PhenotypeCheck(std::vector<uint8_t>& phenotype, uint8_t dx, uint8_t dy) {
      int phenotype_ID=0;
      std::vector<uint8_t> temp_phenotype;
      for(std::vector<uint8_t>::const_iterator phen_iter = known_phenotypes.begin();phen_iter!=known_phenotypes.end();) {
        temp_phenotype.assign(phen_iter+2,phen_iter+2+*phen_iter* *(phen_iter+1));
        if(ComparePolyominoes(phenotype,dx,dy,temp_phenotype,*phen_iter,*(phen_iter+1)))
          return phenotype_ID;
        phen_iter+=*phen_iter* *(phen_iter+1)+2;
        ++phenotype_ID;
      }      
      known_phenotypes.emplace_back(dx);
      known_phenotypes.emplace_back(dy);
      known_phenotypes.insert(known_phenotypes.end(),phenotype.begin(),phenotype.end());
      phenotype_fitnesses.emplace_back(params::real_dist(RNG_Engine));
      ++n_phenotypes;
      return phenotype_ID;
    }

  };
}
