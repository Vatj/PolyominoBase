#include "xorshift.hpp"
#include "stdint.h"

#include <algorithm>
#include <vector>
#include <iostream>
#include <numeric>
#include <utility>
#include <tuple>

namespace interface_model
{

  template<typename T> inline T reverse_bits(T v);
  uint8_t SammingDistance(uint16_t face1,uint16_t face2);
  void MutateInterfaces(std::vector<uint16_t>& binary_genome,double mu_prob);


  int ProteinAssemblyOutcome(std::vector<uint16_t> binary_genome,uint8_t N_repeats);
  
  std::vector<int8_t> AssembleProtein(const std::vector<uint16_t>& binary_genome,double temperature);
  void PlaceNextUnit(const std::vector<uint16_t>& binary_genome,std::vector<uint8_t>& tile_types,std::vector<uint8_t>& faces,std::vector<int8_t>& growing_perimeter,std::vector<int8_t>& placed_tiles,double temperature);
  
  void PerimeterGrowth(int8_t x,int8_t y,uint8_t theta,uint8_t direction, uint8_t tile_type,std::vector<int8_t>& growing_perimeter,std::vector<int8_t>& placed_tiles);
  std::vector<uint8_t> SpatialGrid(std::vector<int8_t>& placed_tiles, uint8_t& dx,uint8_t& dy);

  bool ComparePolyominoes(std::vector<uint8_t>& Spatial_Occupation_Check,uint8_t& Delta_X_Check,uint8_t& Delta_Y_Check, std::vector<uint8_t>& Spatial_Occupation_Compare,uint8_t& Delta_X_Compare,uint8_t& Delta_Y_Compare);
  void ClockwiseRotation(std::vector<uint8_t>& Spatial_Occupation,uint8_t& DELTA_X,uint8_t& DELTA_Y);
  std::vector<uint8_t> ClockwisePiRotation(std::vector<uint8_t>& spatial_occupation);
}
