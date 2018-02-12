#include <utility>
#include <algorithm>
#include <string>
#include <vector>
#include <tuple>
#include <map>
#include <iostream>
#include <set>
#include <numeric>
#include <tuple>

////////////
//Edge Set//
////////////
bool Produce_Edge_Sets(std::vector<int>& genome, std::vector<int>& E_Internal, std::vector<int>& E_External,bool Reject_Early=true);
int Interaction_Matrix(int input_face);
///////////////
//SIF RELATED//
///////////////
int SIF_Elimination(std::vector<int> &genome,bool erase_ends=true);
void Trim_Zero_Tiles(std::vector<int>& genome);
bool two_SIF_Check(std::vector<int>& genome);
bool Dangling_Tile(std::vector<int>& genome);

////////////////
//LOOP RELATED//
////////////////
void Trim_Topology(std::vector<int> genome,std::vector<int>& looping_genome);
int loop_Analysis(std::vector<int>& genome, int& Bound_Limit);
int FastLoopAnalysis(std::vector<int>& genome);
std::tuple<int,int,int> Take_Loop_Step(std::vector<int>& genome ,int tile,int face,int deltaFace,std::vector<int> loop_Path,std::vector<int>& loop_Histories,bool Passed_BP);

bool Cut_Rank_One_Loop(std::vector<int>& genome, std::vector<int>& loopR1_Path);
bool Check_If_Loop(int tile, std::vector<int>& loop_Path);
bool Loop_Rank_One_Check(std::vector<int>& loop_Histories);

bool Check_If_Disjointed_Loops(std::vector<int>& loop_Histories);

////////////////
//TRIM RELATED//
////////////////
void Clean_Genome(std::vector<int>& genome,int secondNonInteracting=-1,bool Remove_Duplicates=true);
void Minimize_Tile_Set(std::vector<int>& genome);
bool Disjointed_Check(std::vector<int>& genome);
void Order_2_Tile_Set(std::vector<int>& genome);
void Relabel_Genome(std::vector<int>& genome);

bool Disjointed_Check_Experimental(std::vector<int>& genome);
void Search_Next_Tile(std::vector<int>& genome, std::vector<int>& Unvisited, std::vector<int>& Connected_Components,int tile);

//////////////
//BP RELATED//
//////////////
bool BP_Check(std::vector<int>& genome);
bool BP_Check_2(std::vector<int>& genome);
bool Double_BP_Check(std::vector<int>& genome);
bool Symmetric_BP_Exception(std::vector<int>& genome);
bool BP_Disjointed_Check(std::vector<int> genome);



//Other
bool Take_Rank_One_Step(std::vector<int>& genome,std::vector<int>& genome_PURE,int tile,int face,int deltaFace,int x,int y,std::vector<int> loop_Path,std::vector<int>& loop_Histories);
