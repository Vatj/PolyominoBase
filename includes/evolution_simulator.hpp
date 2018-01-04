#include "evolution_methods.hpp"


///////////////////////////////////////
//GENOME EVOLTUTION RELATED FUNCTIONS//
///////////////////////////////////////
//bool Evolve_Genomes(int& Discovery_Generation, int& Adaptation_Generation,std::bernoulli_distribution& Mutation_Chance,std::uniform_int_distribution<int>& Mutated_Colour);

//void Evolution_Simulation(const double Mu, int& discovery_Median,int& adaptation_Median,int& Extinction_Events,int& Extinction_Average);
//void Run_Evolution_Simulation_Over_Range(std::string inFileName);


void SetEvolutionTargets(std::vector<int>& targets,int generation);
void EvolveFitnessDynamic(std::string Run_Details,double Mu);
void FitnessEvolutionDynamic();


/*
void EvolveFitnessStatic(std::string Run_Details,double Mu);
void FitnessEvolutionStatic();

void EvolveFitnessDynamicDoublet(std::string Run_Details,double Mu);
void FitnessEvolutionDynamicDoublet();

void EvolveFitnessDynamicSinglet(std::string Run_Details,double Mu);
void FitnessEvolutionDynamicSinglet();


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


void ManyEvolutionSimulations(std::string inFileName);
void EvolveRegulated(int& Discovery_Generation, int& Adaptation_Generation,double Mu);
void EvolveSimple(int& Discovery_Generation, int& Adaptation_Generation,double Mu);
  
void EvolutionSimulation(const double Mu,std::vector<int>& adaptation);
