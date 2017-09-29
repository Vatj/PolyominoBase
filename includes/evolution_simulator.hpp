#include "graph_analysis.hpp"


//////////////////////////////////////
//FITNESS FUNCTION RELATED FUNCTIONS//
//////////////////////////////////////
double Fitness_Function(double Phenotype_Size);
std::vector<int> Random_Selection(int Num_Genomes,int K_Samples);
std::vector<int> Stochastic_Acceptance_Selection(std::vector<double>& Fitness_Weights,int K_Samples);
std::vector<int> Roulette_Wheel_Selection(std::vector<double>& Fitness_Weights,int K_Samples);
int Find_Percentile_of_Vector(std::vector<int>& vec,double percentile);
///////////////////////////////////////
//GENOME EVOLTUTION RELATED FUNCTIONS//
///////////////////////////////////////
bool Evolve_Genomes(int& Discovery_Generation, int& Adaptation_Generation,std::bernoulli_distribution& Mutation_Chance,std::uniform_int_distribution<int>& Mutated_Colour);
void Evolution_Simulation(const double Mu, int& discovery_Median,int& adaptation_Median,int& Extinction_Events,int& Extinction_Average);
void Run_Evolution_Simulation_Over_Range(std::string inFileName);



void Evolve_Fitness3(std::string Run_Details,double Mu);
void Fitness_Evolution3();

void Evolve_Fitness4(std::string Run_Details,double Mu);
void Fitness_Evolution4();

void Evolve_Fitness5(std::string Run_Details,double Mu);
void Fitness_Evolution5();

/*
void Evolve_Fitness2(std::string Run_Details,double Mu);
void Fitness_Evolution2();
void Evolve_Fitness5(std::string Run_Details,double Mu);
void Fitness_Evolution5();

void Evolve_Fitness1(std::string Run_Details,double Mu);
void Fitness_Evolution1();

void Evolve_Fitness6(std::string Run_Details,double Mu);
void Fitness_Evolution6();
void Evolve_Fitness7(double Mu);
*/
/////////////////////
//UTILITY FUNCTIONS//
/////////////////////

void Set_Runtime_Configurations(int argc, char* argv[]);

void Run_Evolution_Simulation_Over_Range2(std::string inFileName);
void Evolution_Simulation2(const double Mu,std::vector<int>& adaptation);
