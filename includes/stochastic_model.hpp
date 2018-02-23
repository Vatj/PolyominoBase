#include <graph_analysis.hpp>
#include <fstream>
struct Phenotype {
  int dx=1;
  int dy=1;
  std::vector<int> tiling{1};
};

struct PhenotypeTable;


namespace Stochastic
{
  //extern std::random_device rd;
  extern std::mt19937 RNG_Generator;
  std::vector<int> Stochastic_Polyomino_Builder(std::vector<int> genome, unsigned int THRESHOLD_SIZE,int initial_Tile,int initial_Rotation);
  void Stochastic_Interacting_Adjacency(std::vector<int>& Interacting_Faces, int interacting_Face, int face_index, int X, int Y);

  int Analyse_Genotype_Outcome(std::vector<int> genome, int N_Repeated_Checks,PhenotypeTable* pt,int seed=0);
  bool Unbound_Deterministic_Check(std::vector<int>& Spatial_Occupation,int Delta_X,int Delta_Y);
  std::vector<int> Generate_Spatial_Occupancy(std::vector<int>& Placed_Tiles_Check, int& DELTA_X_Check,int& DELTA_Y_Check,int generate_mode);
  

  bool Compare_Two_Polyominoes_Shapes(std::vector<int>& rs_lattice,int dx,int dy, std::vector<int>& rs_lattice_prime,int dx_prime,int dy_prime);
}

struct PhenotypeTable {
  std::vector<Phenotype> known_phenotypes;

  int GetPhenotypeID(Phenotype phen) {
    for(size_t index=0;index<known_phenotypes.size();++index) {
      Phenotype known=known_phenotypes[index];
      if(Stochastic::Compare_Two_Polyominoes_Shapes(phen.tiling,phen.dx,phen.dy,known.tiling,known.dx,known.dy))
        return index;
    }
    if(phen.dx > phen.dy) {
      Clockwise_Rotation(phen.tiling,phen.dx,phen.dy);
      std::swap(phen.dx,phen.dy);
    }
    known_phenotypes.push_back(phen);
    return known_phenotypes.size()-1;
  }
  void PrintTable(std::ofstream& fout) {
    for(Phenotype known : known_phenotypes) {
      fout<<known.dx<<" "<<known.dy<<" ";
      for(auto tile : known.tiling)
        fout<<tile<<" ";
      fout<<"\n";
    }
  }
    
};
