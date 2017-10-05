#include "graph_analysis.hpp"

extern int Target_Fitness;
extern int Colour_Space; 
extern int GENERATION_LIMIT;
extern int Num_Runs;
extern int Num_Genomes;
extern int Num_Tiles;
extern bool Regulated;
extern int Shape_Matching_Fitness;
extern int Fitness_Mode;
extern int Fitness_Oscillation_Rate;
extern int burn_in_period;
extern int mutation_cofactor;
extern int RUN;

extern std::random_device rd;
extern xorshift RNG_Engine;
extern const double MINIMUM_FITNESS_THRESHOLD;

double Fitness_Function(double Phenotype_Size);

std::vector<int> Random_Selection(int Num_Genomes,int K_Samples);
std::vector<int> Stochastic_Acceptance_Selection(std::vector<double>& Fitness_Weights,int K_Samples);
std::vector<int> Roulette_Wheel_Selection(std::vector<double>& Fitness_Weights,int K_Samples);


int Find_Percentile_of_Vector(std::vector<int>& vec,double percentile);

void Set_Runtime_Configurations(int argc, char* argv[]);

/*
class GenotypeMutator {
   public:
  void mutate(std::vector<int> genotype);
  GenotypeMutator(int c, double m);
 
   private:
  std::uniform_int_distribution<int> Mutated_Colours;
  std::bernoulli_distribution Mutation_Chance;

};
*/
