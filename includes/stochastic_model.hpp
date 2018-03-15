#include <graph_analysis.hpp>
#include <core_phenotype.hpp>


struct StochasticPhenotypeTable : PhenotypeTable {
  Phenotype_ID GetPhenotypeID(Phenotype phen) {
     uint8_t phen_size=std::count_if(phen.tiling.begin(),phen.tiling.end(),[](const uint8_t s) {return s>0;});
    for(size_t index=0;index<known_phenotypes[phen_size].size();++index) {
      Phenotype known=known_phenotypes[phen_size][index];
      if(ComparePolyominoes(phen,known))
        return std::make_pair(phen_size,index);
    }
    if(phen.dy > phen.dx) 
      ClockwiseRotation(phen);
    MinimalTilingRepresentation(phen.tiling);
    known_phenotypes[phen_size].emplace_back(phen);
    return std::make_pair(phen_size,known_phenotypes[phen_size].size()-1);
  }
};

namespace Stochastic
{
  //extern std::random_device rd;
  //extern std::mt19937 RNG_Generator;

  Phenotype_ID Analyse_Genotype_Outcome(Genotype genome, uint8_t N_Repeated_Checks, StochasticPhenotypeTable* pt,uint8_t seed);
  std::vector<int8_t> Stochastic_Polyomino_Builder(const Genotype& genome, uint8_t THRESHOLD_SIZE, uint8_t initial_Tile);
  void Interacting_Adjacency(std::vector<int8_t>& Interacting_Faces, uint8_t interacting_Face, uint8_t face_index, int8_t X, int8_t Y);
  Phenotype Generate_Spatial_Occupancy(std::vector<int8_t>& Placed_Tiles_Check, uint8_t generate_mode);  

}




