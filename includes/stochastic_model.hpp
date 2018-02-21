#include <graph_analysis.hpp>
namespace Brute_Force
{
  extern std::random_device rd;
  extern std::mt19937 RNG_Generator;
  std::vector<int> Brute_Force_Polyomino_Builder(std::vector<int> genome, unsigned int THRESHOLD_SIZE,int initial_Tile,int initial_Rotation);
  void Brute_Force_Interacting_Adjacency(std::vector<int>& Interacting_Faces, int interacting_Face, int face_index, int X, int Y);

  int Analyse_Genotype_Outcome(std::vector<int> genome, int N_Repeated_Checks);
  bool Unbound_Deterministic_Check(std::vector<int>& Spatial_Occupation,int Delta_X,int Delta_Y);
  std::vector<int> Generate_Spatial_Occupancy(std::vector<int>& Placed_Tiles_Check, int& DELTA_X_Check,int& DELTA_Y_Check,int generate_mode);
  

  bool Compare_Two_Polyominoes_Shapes(std::vector<int>& Spatial_Occupation_Check,int Delta_X_Check,int Delta_Y_Check, std::vector<int>& Spatial_Occupation_Compare,int Delta_X_Compare,int Delta_Y_Compare);

  int Compare_Two_Polyominoes_Tile_Details(std::vector<int>& Placed_Tiles_Check,std::vector<int>& Placed_Tiles_Compare);
}
