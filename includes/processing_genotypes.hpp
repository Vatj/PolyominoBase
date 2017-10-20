#include "graph_analysis.hpp"
#include <sstream>
#include <iterator>

void DoShape(std::string sub_file,std::string folder_base,std::vector<int> runs);


void DeleteriousRobustness(int r);
bool CheckRobust(std::vector<int> target,int target_sets=0,int num_targets=0);
