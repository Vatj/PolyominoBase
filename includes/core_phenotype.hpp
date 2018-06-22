#include "core_genotype.hpp"

#include <iostream>
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
typedef std::map<std::vector<Phenotype_ID>, std::vector<Genotype>> Set_to_Genome;

namespace simulation_params
{
  extern uint8_t n_genes, colours, metric_colours;
  extern std::mt19937 RNG_Engine;
  extern double UND_threshold;
  extern uint8_t phenotype_builds;
}

struct PhenotypeTable {
  std::unordered_map<uint8_t,std::vector<Phenotype> > known_phenotypes,undiscovered_phenotypes;
  std::unordered_map<uint8_t,std::vector<uint16_t> > new_phenotype_xfer;
  std::vector<uint16_t> undiscovered_phenotype_counts;

  Phenotype_ID GetPhenotypeID(Phenotype& phen) {
    uint8_t phenotype_size=std::count_if(phen.tiling.begin(),phen.tiling.end(),[](const int c){return c != 0;});
    for(uint16_t phenotype_index=0; phenotype_index != known_phenotypes[phenotype_size].size();++phenotype_index) {
      if(ComparePolyominoes(phen,known_phenotypes[phenotype_size][phenotype_index]))
	return std::make_pair(phenotype_size,phenotype_index);
    }
    uint8_t new_phenotype_index=0;
    for(Phenotype phen_p : undiscovered_phenotypes[phenotype_size]) {
      if(ComparePolyominoes(phen,phen_p)) {
	if(++undiscovered_phenotype_counts[new_phenotype_index]>=ceil(simulation_params::UND_threshold*simulation_params::phenotype_builds)) {
	  new_phenotype_xfer[phenotype_size].emplace_back(known_phenotypes[phenotype_size].size()+new_phenotype_index+simulation_params::phenotype_builds);
	  known_phenotypes[phenotype_size].emplace_back(phen);
	  new_phenotype_xfer[phenotype_size].emplace_back(known_phenotypes[phenotype_size].size()-1);

	  return std::make_pair(phenotype_size,known_phenotypes[phenotype_size].size()-1);
	}
	else
	  return std::make_pair(phenotype_size,known_phenotypes[phenotype_size].size()+new_phenotype_index+simulation_params::phenotype_builds);
      }
      ++new_phenotype_index;
    }
    GetMinPhenRepresentation(phen);
    undiscovered_phenotypes[phenotype_size].emplace_back(phen);
    undiscovered_phenotype_counts.emplace_back(1);
    return std::make_pair(phenotype_size,known_phenotypes[phenotype_size].size()+new_phenotype_index+simulation_params::phenotype_builds);
  }


  void RelabelPhenotypes(std::vector<Phenotype_ID >& pids) {
    for(std::unordered_map<uint8_t,std::vector<uint16_t> >::iterator x_iter=new_phenotype_xfer.begin();x_iter!=new_phenotype_xfer.end();++x_iter)
      for(std::vector<uint16_t>::iterator r_iter=x_iter->second.begin();r_iter!=x_iter->second.end();r_iter+=2)
	std::replace(pids.begin(),pids.end(),std::make_pair(x_iter->first,*(r_iter)),std::make_pair(x_iter->first,*(r_iter+1)));
    undiscovered_phenotypes.clear();
    undiscovered_phenotype_counts.clear();
    new_phenotype_xfer.clear();

  }

  /* Count each ID frequency */
  std::map<Phenotype_ID,uint8_t> PhenotypeFrequencies(std::vector<Phenotype_ID >& pids, bool& rare_phenotypes) {
    std::map<Phenotype_ID, uint8_t> ID_counter;
    for(std::vector<Phenotype_ID >::const_iterator ID_iter = pids.begin(); ID_iter!=pids.end(); ++ID_iter) {
      if(ID_iter->second < known_phenotypes[ID_iter->first].size())
	++ID_counter[std::make_pair(ID_iter->first,ID_iter->second)];
      else
	rare_phenotypes=true;
    }
    return ID_counter;
  }

  // void PrintTable(std::ofstream& fout)
  // {
  //   for(auto known_phens : known_phenotypes)
  //   {
  //     uint16_t n_phen=0;
  //     for(Phenotype known : known_phens.second)
  //     {
	//        fout<<+known_phens.first<<" "<<+n_phen++<<" "<<+known.dx<<" "<<+known.dy<<" ";
	//        for(uint8_t tile : known.tiling)
	//           fout<<+tile<<" ";
	//        fout<<"\n";
  //     }
  //   }
  // }

  void PrintTable(std::string phenotype_file)
  {
    std::cout << "Printing phenotypes to file : " << phenotype_file << "\n";
    std::ofstream fout(phenotype_file);
    for(auto known_phens : known_phenotypes)
    {
      uint16_t n_phen=0;
      for(Phenotype known : known_phens.second)
      {
	       fout<<+known_phens.first<<" "<<+n_phen++<<" "<<+known.dx<<" "<<+known.dy<<" ";
	       for(uint8_t tile : known.tiling)
	          fout<<+tile<<" ";
	       fout<<"\n";
      }
    }
  }

};
