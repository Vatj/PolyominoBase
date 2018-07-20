#include <cstdint>
#include <vector>
#include <algorithm>
#include <utility>
#include <set>
#include <numeric>
#include <tuple>
#include <map>

typedef std::vector<uint8_t> Genotype;

uint8_t Interaction_Matrix(uint8_t input_face);
void Clean_Genome(Genotype& genome, bool Remove_Duplicates);
void Minimize_Tile_Set(Genotype& genome);
std::map<uint8_t,uint8_t> DuplicateGenes(Genotype& genome);
bool Disjointed_Check(Genotype& genome);
void Search_Next_Tile(Genotype& genome, std::vector<uint8_t>& Unvisited, std::vector<uint8_t>& Connected_Components,uint8_t tile);

/*Converting methods*/
uint64_t genotype_to_index(Genotype& genotype, uint16_t n_genes, uint16_t colours);
void index_to_genotype(uint64_t index, Genotype& genotype, uint16_t n_genes, uint16_t colours);
