#include <cstdint>
#include <vector>
#include <algorithm>
#include <utility>
#include <set>
#include <numeric>
#include <tuple>
#include <map>


uint8_t Interaction_Matrix(uint8_t input_face);
void Clean_Genome(std::vector<uint8_t>& genome,int secondNonInteracting,bool Remove_Duplicates);
void Minimize_Tile_Set(std::vector<uint8_t>& genome);
std::map<uint8_t,uint8_t> DuplicateGenes(std::vector<uint8_t>& genome);
bool Disjointed_Check(std::vector<uint8_t>& genome);
void Search_Next_Tile(std::vector<uint8_t>& genome, std::vector<uint8_t>& Unvisited, std::vector<uint8_t>& Connected_Components,uint8_t tile);
