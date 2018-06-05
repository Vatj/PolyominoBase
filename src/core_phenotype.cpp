#include "core_phenotype.hpp"

bool ComparePolyominoes(Phenotype& phen1, const Phenotype& phen2,uint8_t GAUGE) {
  if(phen1.tiling.size()!=phen2.tiling.size() || std::count(phen1.tiling.begin(),phen1.tiling.end(),0)!=std::count(phen2.tiling.begin(),phen2.tiling.end(),0))
    return false;
  if(phen1.dx==phen2.dx && phen1.dy==phen2.dy && phen1.dx==phen2.dy) {
    if(phen1.tiling==phen2.tiling)
      return true;
    for(int rotation=0;rotation<3;++rotation) {
      ClockwiseRotation(phen1);
      MinimizePhenRep(phen1.tiling,GAUGE);
      if(phen1.tiling==phen2.tiling)
        return true;
    }
    return false;
  }
  if(phen1.dx==phen2.dx && phen1.dy==phen2.dy) {
    if(phen1.tiling==phen2.tiling)
      return true;
    ClockwisePiRotation(phen1);
    MinimizePhenRep(phen1.tiling,GAUGE);
    return phen1.tiling==phen2.tiling;
  }
  return false;
}

void ClockwiseRotation(Phenotype& phen) {
  std::vector<uint8_t> swapper;
  swapper.reserve(phen.tiling.size());
  for(uint8_t column=0;column<phen.dx;++column)
    for(uint8_t row=phen.dy;row!=0;--row)
      swapper.emplace_back(phen.tiling[(row-1)*phen.dx+column]);
  std::swap(phen.dx,phen.dy);
  phen.tiling=swapper;
}

void ClockwisePiRotation(Phenotype& phen) {
  std::reverse(phen.tiling.begin(),phen.tiling.end());
}

void MinimizePhenRep(std::vector<uint8_t>& tiling,uint8_t GAUGE) {
  if(tiling.size()==1)
    return;
  for(uint8_t& t:tiling)
    t+=128*(t!=0);
  uint8_t swap_count=1;
  for(std::vector<uint8_t>::iterator t_iter=std::find_if(tiling.begin(),tiling.end(),[](const int s) { return s>0; });t_iter!=tiling.end();) {
    const uint8_t static_swap=*t_iter;
    for(uint8_t cyclic = 0; cyclic<GAUGE; ++cyclic) {
      uint8_t pre_swap=(static_swap-(static_swap-1)%GAUGE)+((static_swap-1)%GAUGE+cyclic)%GAUGE;
      std::replace(tiling.begin(),tiling.end(),swap_count,uint8_t(255));
      std::replace(tiling.begin(),tiling.end(),pre_swap,swap_count);
      std::replace(tiling.begin(),tiling.end(),uint8_t(255),pre_swap);
      ++swap_count;
    }
    t_iter=std::find_if(tiling.begin(),tiling.end(),[swap_count,GAUGE](const int s) { return (s>0 && s-(s-1)%GAUGE >= swap_count); });
  }
}
