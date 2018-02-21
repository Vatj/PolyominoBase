//#include "graph_analysis.hpp"
#include "stochastic_model.hpp"
#include "stdint.h"


namespace Brute_Force {

  
  

  //RUNNING COMMANDS
  void Test_2_Space(int MAX_C);

  void Test_Func(int C);

  void Brute_vs_Graph_Methods_Comparison(int MAX_C, int K_Repetions);

  void Brute_vs_Graph_Methods_Comparison_Random_Sample(unsigned int genome_Length,int MAX_C,int num_Samples, int K_Repetitions);
  void Brute_vs_Brute_Comparison_Random_Sample(unsigned int genome_Length,int MAX_C,int num_Samples, int K_Repetitions);

  void Enumerate_Topologies(int MAX_C);
  void Compare_Topology_Classes();
  void Topology_Multiplicity(int MAX_C);

  bool False_Fail_Detection(int Num_Trials);
    
}

/////////////////
//PAPER METHODS//
/////////////////

//Timing//
double get_wall_time();
double get_cpu_time();

//Generating
void Generate_Random_Genotypes(int genome_Length,std::vector<std::vector<int> >& GENOME_VECTOR, bool triangularly_random,int condition=-1, int C=-1);
bool Generating_Condition(std::vector<int>& genome,int condition=-1);

//Methods
void Run_Selection(int genome_length,int Test_Cases);
void Timing_Mode(std::vector<std::vector<int> >& GENOME_VECTOR, std::vector<int> k_repeats);
void Accuracy_Mode(std::vector<std::vector<int> >& GENOME_VECTOR, std::vector<int> k_repeats);

//other
void Run_Over_BDs(std::vector<int> genome_Lengths,std::vector<int> Cs, int N_Cases);
void BD_Fractional_Occupation(int genome_Length,int C,int N_Cases);


void Topology_Robustness(std::vector<int> genotype,int N_runs, std::vector<int> Colours,const int R_N);

void GenotypeSpaceSampler(int N_samples);
