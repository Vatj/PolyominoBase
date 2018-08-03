#include <graph_analysis.hpp>
#include <core_metrics.hpp>

namespace simulation_params
{
  extern std::mt19937 RNG_Engine;
  extern uint16_t phenotype_builds;
  extern double UND_threshold;
  extern bool STERIC_FORBIDDEN;
}

namespace Stochastic
{
  std::vector<Phenotype_ID> AssemblePlasticGenotype(Genotype genotype, PhenotypeTable* pt);
  std::map<Phenotype_ID, uint8_t> AssemblePlasticGenotypeFrequency(Genotype genotype, PhenotypeTable* pt);
  bool FastNoPIDs(Genotype genotype, PhenotypeTable* pt);

  Phenotype Generate_Spatial_Occupancy(std::vector<int8_t>& Placed_Tiles_Check);
  std::vector<int8_t> Stochastic_Polyomino_Builder(const Genotype& genome, uint8_t initial_Tile);
  void Interacting_Adjacency(std::vector<int8_t>& Interacting_Faces, uint8_t interacting_Face, uint8_t face_index, int8_t X, int8_t Y);
}
