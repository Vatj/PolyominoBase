#include <graph_analysis.hpp>
#include <core_metrics.hpp>

namespace Stochastic
{
  extern bool STERIC_FORBIDDEN;

  std::vector<Phenotype_ID> AssemblePlasticGenotype(Genotype genotype, PhenotypeTable* pt);
  
  Phenotype Generate_Spatial_Occupancy(std::vector<int8_t>& Placed_Tiles_Check, uint8_t generate_mode);
  Phenotype_ID Analyse_Genotype_Outcome(Genotype genome, uint8_t N_Repeated_Checks, PhenotypeTable* pt,uint8_t seed);
  std::vector<int8_t> Stochastic_Polyomino_Builder(const Genotype& genome, uint8_t THRESHOLD_SIZE, uint8_t initial_Tile);
  void Interacting_Adjacency(std::vector<int8_t>& Interacting_Faces, uint8_t interacting_Face, uint8_t face_index, int8_t X, int8_t Y);
}
