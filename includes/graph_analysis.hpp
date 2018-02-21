#include "graph_methods.hpp"
#include <fstream>


//Main//
int Graph_Analysis(std::vector<int> genotype);
int FastGraphAnalysis(std::vector<int> genotype);
bool Bounded_Loop_Growth_Check(std::vector<int>& genome, int Bound_Limit);

///////////
//STERICS//
///////////
double Steric_Check(std::vector<int>& genome,int Shape_Based_Fitness);
std::vector<int> Steric_Check_Occupation(std::vector<int>& genome);
int Steric_Check_Table(std::vector<int>& genome,std::vector<int>& Known_Shapes,int& Num_Shapes);

bool Traverse_Numbered_Sterics(std::vector<int>& Spatial_Grid,int DELTA_X,int DELTA_Y,int x,int y,int cameFrom);

//int SpiralCoordinate(int x, int y);

/////////////
//ROTATIONS//
/////////////

void Clockwise_Pi_Rotation(std::vector<int>& Spatial_Occupation,int DELTA_X,int DELTA_Y);

void Clockwise_Rotation(std::vector<int>& Spatial_Occupation,int DELTA_X,int DELTA_Y);



///////////
//FITNESS//
///////////
int Get_Phenotype_Fitness(std::vector<int> genome,int Shape_Matching_Fitness,bool zeroth_fitness=false);
bool GetMultiplePhenotypeFitness(std::vector<int> genome,std::vector<int> target_types,std::vector<double>& target_fitnesses,int active_targets);
double Shape_Matching_Fitness_Function(std::vector<int>& Spatial_Occupation,int& Delta_X,int& Delta_Y,int target_choice);



//defunct
int Phenotype_Fitness_Size_Check(std::vector<int> genome);
double Phenotype_Fitness_Shape_Check(std::vector<int> genome);


namespace Shape_Match {
  void Set_Shape_Data(std::vector<int> Targ_Occ,int Targ_X,int Targ_Y);
}
////////
//RUNS//
////////
void Run_3_Space();
void Sample_Space(int genomeLenght,int num_Sample);
void Run_2_Space(const int MAX_C);
void Read_Run_Space();
