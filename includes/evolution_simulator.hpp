#include "evolution_methods.hpp"


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


void Run_Evolution_Simulation_Over_Range2(std::string inFileName);
void Evolution_Simulation2(const double Mu,std::vector<int>& adaptation);
