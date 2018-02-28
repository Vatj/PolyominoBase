#include <graph_analysis.hpp>
#include <core_phenotype.hpp>



struct StochasticPhenotypeTable;

namespace Stochastic
{
  //extern std::random_device rd;
  extern std::mt19937 RNG_Generator;
  std::vector<int> Stochastic_Polyomino_Builder(std::vector<int> genome, unsigned int THRESHOLD_SIZE,int initial_Tile,int initial_Rotation);
  void Stochastic_Interacting_Adjacency(std::vector<int>& Interacting_Faces, int interacting_Face, int face_index, int X, int Y);

  Phenotype_ID Analyse_Genotype_Outcome(std::vector<int> genome, int N_Repeated_Checks,StochasticPhenotypeTable* pt,int seed=0);
  bool Unbound_Deterministic_Check(std::vector<int>& Spatial_Occupation,int Delta_X,int Delta_Y);
  Phenotype Generate_Spatial_Occupancy(std::vector<int>& Placed_Tiles_Check,int generate_mode);
  

}


struct StochasticPhenotypeTable : PhenotypeTable {

  Phenotype_ID GetPhenotypeID(Phenotype phen) {
     uint8_t phen_size=std::count_if(phen.tiling.begin(),phen.tiling.end(),[](const uint8_t s) {return s>0;});
    for(size_t index=0;index<known_phenotypes[phen_size].size();++index) {
      Phenotype known=known_phenotypes[phen_size][index];
      if(ComparePolyominoes(phen,known))
        return std::make_pair(phen_size,index);
    }
    if(phen.dy > phen.dx) {
      ClockwiseRotation(phen);
    }
    MinimalTilingRepresentation(phen.tiling);
    known_phenotypes[phen_size].emplace_back(phen);
    return std::make_pair(phen_size,known_phenotypes[phen_size].size()-1);
  }

    
};

