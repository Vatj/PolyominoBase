#include "interface_simulator.hpp"

//./bin/ProteinEvolution -E -N 3 -P 50 -K 5 -B 50 -S 1 -D 1 -V 0 -F 1 -M 0.05 -T 0.09 -X 0.2
const uint16_t printing_resolution=100; /* NOTE RESOLUTION IS 1 */



void EvolvePopulation(std::string run_details) {
  /* Output files */
  std::string file_base_path="//scratch//asl47//Data_Runs//Bulk_Data//";    //"//rscratch//asl47//Bulk_Run//Interfaces//";
  std::string file_simulation_details="I"+std::to_string(model_params::interface_size)+"_T"+std::to_string(model_params::binding_threshold)+"_Mu"+std::to_string(model_params::mu_prob)+"_Gamma"+std::to_string(model_params::fitness_factor)+run_details+".txt";
    
  std::ofstream fout_strength(file_base_path+"Strengths_"+file_simulation_details, std::ios_base::out);
  std::ofstream fout_fitness(file_base_path+"Fitness_"+file_simulation_details, std::ios_base::out);
  std::ofstream fout_phenotype(file_base_path+"Phenotypes_"+file_simulation_details, std::ios_base::out);
  std::ofstream fout_genotype_history(file_base_path+"GenotypeHistory_"+file_simulation_details, std::ios_base::out);
  std::ofstream fout_phenotype_history(file_base_path+"PhenotypeHistory_"+file_simulation_details, std::ios_base::out);
  
  interface_model::InterfacePhenotypeTable pt = interface_model::InterfacePhenotypeTable();

  std::vector<double> population_fitnesses(simulation_params::population_size);
  std::vector<uint32_t> interface_counter(1.5*model_params::interface_size+2);


  bool record_strengths=false;

  std::vector<PopulationGenotype> evolving_population(simulation_params::population_size),reproduced_population(simulation_params::population_size);

  if(simulation_params::random_initilisation) {
    std::uniform_int_distribution<interface_type> dist;
    auto interface_filler = std::bind(dist, std::ref(interface_model::RNG_Engine));
    for(auto& species : evolving_population)
      std::generate(species.genotype.begin(),species.genotype.end(),interface_filler);

  }


  /* Median time to complete interface mutation, characteristic time for fitness re-assignment */
  std::poisson_distribution<uint16_t> landscape_changer(log(1-pow(2,-1./model_params::interface_size))/log(1-model_params::mu_prob/(4*simulation_params::n_tiles*model_params::interface_size)));
  uint16_t fitness_jiggle=landscape_changer(interface_model::RNG_Engine);

  /* Start main evolution loop */
  
  for(uint32_t generation=0;generation<simulation_params::generation_limit;++generation) {
    if(fitness_jiggle--==0) {
      //pt.ReassignFitness(); /* NOTE JIGGLER OFF */
      fitness_jiggle=landscape_changer(interface_model::RNG_Engine);
    }
    if(generation+100>=simulation_params::generation_limit)
      record_strengths=false;


    /* Start genotype loop */
    int nth_genotype=0;
    for(PopulationGenotype& evolving_genotype : evolving_population) {
      

      interface_model::MutateInterfaces(evolving_genotype.genotype);

      std::vector<interaction_pair> pid_interactions;
      population_fitnesses[nth_genotype]=interface_model::ProteinAssemblyOutcome(evolving_genotype.genotype,&pt,evolving_genotype.pid,pid_interactions);
      for(auto x : pid_interactions)
        fout_strength<<+x.first<<" "<<+x.second<<",";  
      fout_strength<<".";

      if(record_strengths)
        InterfaceStrengths(evolving_genotype.genotype,interface_counter);
      
      for(interface_type base_value : evolving_genotype.genotype)
        fout_genotype_history << +base_value << " ";
      fout_phenotype_history << +evolving_genotype.pid.first <<" "<<+evolving_genotype.pid.second<<" ";
      
      ++nth_genotype;
    } /* End genotype loop */
    fout_genotype_history<<"\n";
    fout_phenotype_history<<"\n";
    fout_strength<<"\n";
    

    /* Write data to file */
    if(record_strengths) {
      for(uint32_t count : interface_counter)
        fout_strength<<count<<" ";
      fout_strength<<"\n";
      std::fill(interface_counter.begin(),interface_counter.end(),0);
    }
    if(simulation_params::generation_limit>printing_resolution && generation%(simulation_params::generation_limit/printing_resolution)==0) {
      double mu=0,sigma=0;
      DistributionStatistics(population_fitnesses,mu,sigma);
      fout_fitness<<mu<<" "<<sigma<<"\n";
    }

    /* Start selection */
    std::vector<uint16_t> selection_indices=RouletteWheelSelection(population_fitnesses);
    for(uint16_t nth_reproduction = 0; nth_reproduction<simulation_params::population_size;++nth_reproduction)
      reproduced_population[nth_reproduction]=evolving_population[selection_indices[nth_reproduction]];
    evolving_population.assign(reproduced_population.begin(),reproduced_population.end());
    for(uint16_t selection_index : selection_indices)
      fout_phenotype_history << +selection_index << " ";
    fout_phenotype_history<<"\n";    
    /* End selection */
    
  } /* End main evolution loop */
  
  pt.PrintTable(fout_phenotype);  
  
  fout_fitness.close();
  fout_strength.close();
  fout_phenotype.close();
  fout_genotype_history.close();

  
}

  
std::vector<uint8_t> SequenceDifference(const std::vector<interface_type>& parent, const std::vector<interface_type>& child) {
  std::vector<uint8_t> divergence_indices;
  for(uint8_t base=0; base < parent.size(); ++base)
    if(parent[base]!=child[base])
      divergence_indices.emplace_back(base);
  return divergence_indices;
}

uint8_t ConjugateInterface(std::vector<interface_type>& genotype,uint8_t mutation_site) {
  for(uint8_t base =0; base < genotype.size(); ++base)
    if(interface_model::SammingDistance(genotype[base],genotype[mutation_site])==0) {
      return base;
    }
  return 255;
}

std::vector<uint16_t> RouletteWheelSelection(std::vector<double>& fitnesses) {
  std::partial_sum(fitnesses.begin(), fitnesses.end(), fitnesses.begin());
  std::vector<uint16_t> selected_indices(simulation_params::population_size);
  std::uniform_real_distribution<double> random_interval(0,fitnesses.back());
  for(uint16_t nth_selection=0; nth_selection<simulation_params::population_size; ++nth_selection) 
    selected_indices[nth_selection]=static_cast<uint16_t>(std::lower_bound(fitnesses.begin(),fitnesses.end(),random_interval(interface_model::RNG_Engine))-fitnesses.begin());
  std::sort(selected_indices.begin(),selected_indices.end());
  return selected_indices;
}


/********************/
/*******!MAIN!*******/
/********************/
void EvolutionRunner() {
#pragma omp parallel for schedule(dynamic)
  for(uint16_t r=0;r<simulation_params::independent_trials;++r) {
    std::string run_details="_Run"+std::to_string(r+simulation_params::run_offset);
    EvolvePopulation(run_details);
  }
}

int main(int argc, char* argv[]) {
  char run_option;
  if(argc<2) {
    std::cout<<"no Params"<<std::endl;
    run_option='H';
  }
  else {
    run_option=argv[1][1];
    SetRuntimeConfigurations(argc,argv);
  }
  
  //auto a = BindingProbabilities();
  BGenotype bg{15540537610786488263,3375271998381208445,0,2018349694661896724, 0,0,3375271998381208445,4697909892076606603};
  std::set<interaction_pair > aa;
  std::vector<int8_t> a;
  Phenotype phen;
  switch(run_option) {
  case 'E':
    EvolutionRunner();
    break;
  case 'X':
    a=interface_model::AssembleProteinNew(bg,aa);
    for(auto qqqq :a)
      std::cout<<+qqqq<<" ";
    std::cout<<std::endl;
    for(auto qq : aa)
      std::cout<<+qq.first<<"/"<<+qq.second<<std::endl;
    std::cout<<"Unused at this time"<<std::endl;
    phen=SpatialGrid(a);
    std::cout<<"phen: "<<+phen.dx<<", "<<+phen.dy<<", ";
    for(auto qqqqqq: phen.tiling)
      std::cout<<+qqqqqq<<", ";
    std::cout<<std::endl;
    break;
  case 'D':
    std::cout<<"ProteinEvolution -E -N 2 -P 100 -K 250 -B 25 -S 1 -R 0 -F 1 -X 1 -T 0.000001 -M 0.25 -D 1"<<std::endl;
    break;
  case '?':
    for(auto b : model_params::binding_probabilities)
      std::cout<<b<<std::endl;
    break;
  case 'H':
  default:
    std::cout<<"Protein interface model\n**Simulation Parameters**\nN: number of tiles\nP: population size\nK: generation limit\nB: number of phenotype builds\n";
    std::cout<<"\n**Model Parameters**\nU: mutation probability (per interface)\nT: temperature\nI: unbound size factor\nA: misbinding rate\nM: Fitness factor\n";
    std::cout<<"\n**Run options**\nR: evolution without fitness\nE: evolution with fitness\n";
    break;
  }
  return 0;
}

void SetRuntimeConfigurations(int argc, char* argv[]) {
  if(argc<3 && argv[1][1]!='H')
    std::cout<<"Invalid Parameters"<<std::endl;
  else {
    for(uint8_t arg=2;arg<argc;arg+=2) {
      switch(argv[arg][1]) {
      case 'N': simulation_params::n_tiles=std::stoi(argv[arg+1]);break;
      case 'P': simulation_params::population_size=std::stoi(argv[arg+1]);break;
      case 'K': simulation_params::generation_limit=std::stoi(argv[arg+1]);break;
      case 'B': simulation_params::phenotype_builds=std::stoi(argv[arg+1]);break;
	
      case 'D': simulation_params::independent_trials=std::stoi(argv[arg+1]);break;
      case 'V': simulation_params::run_offset=std::stoi(argv[arg+1]);break;
      case 'R': simulation_params::random_initilisation=std::stoi(argv[arg+1])>0;break; 
        
      case 'F': model_params::fitness_factor=std::stod(argv[arg+1]);break;

      case 'S': model_params::fixed_seed=std::stoi(argv[arg+1])>0;break;
      case 'M': model_params::mu_prob=std::stod(argv[arg+1]);break;
      case 'T': model_params::temperature=std::stod(argv[arg+1]);break;
      case 'Y': model_params::binding_threshold=double(static_cast<uint8_t>(ceil(model_params::interface_size*std::stod(argv[arg+1]))))/model_params::interface_size;break;
      case 'X': model_params::UND_threshold=std::stod(argv[arg+1]);break;
      case 'I': model_params::interface_threshold=std::stod(argv[arg+1]);break;
	
        
      default: std::cout<<"Unknown Parameter Flag: "<<argv[arg][1]<<std::endl;
      }
    }
    model_params::b_dist.param(std::binomial_distribution<uint8_t>::param_type(model_params::interface_size,model_params::mu_prob/(model_params::interface_size*4*simulation_params::n_tiles)));
    model_params::binding_probabilities=GenBindingProbsLUP();
  }
}





