#include "core_genotype.hpp"


#include <fstream>
#include <unordered_map>

constexpr uint8_t FREE_POLYOMINO = true ? 2 : 1, GAUGE = true ? 4 : 1;

struct Phenotype {
  uint8_t dx=1;
  uint8_t dy=1;
  std::vector<uint8_t> tiling{1};
};

void ClockwiseRotation(Phenotype& phen);
void ClockwisePiRotation(Phenotype& phen);
void ChiralFlip(Phenotype& phen);

bool ComparePolyominoes(Phenotype& phen1, const Phenotype& phen2);
void MinimizePhenRep(std::vector<uint8_t>& tiling);
void GetMinPhenRepresentation(Phenotype& phen);

typedef std::pair<uint8_t,uint16_t> Phenotype_ID;

struct PhenotypeTable {
  std::unordered_map<uint8_t,std::vector<Phenotype> > known_phenotypes;
        
  void PrintTable(std::ofstream& fout) {
    for(auto known_phens : known_phenotypes) {
      uint16_t n_phen=0;
      for(Phenotype known : known_phens.second) {
	fout<<+known_phens.first<<" "<<+n_phen++<<" "<<+known.dx<<" "<<+known.dy<<" ";
	for(uint8_t tile : known.tiling)
	  fout<<+tile<<" ";
	fout<<"\n";
      }
    }
  }

};
