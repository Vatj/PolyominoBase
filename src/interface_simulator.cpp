#include "interface_simulator.hpp"

//./bin/ProteinEvolution -E -N 3 -P 50 -K 5 -B 50 -S 1 -D 1 -V 0 -F 1 -M 0.05 -T 0.09 -X 0.2
const uint16_t printing_resolution=100;

std::vector<uint16_t> RouletteWheelSelection(std::vector<double>& fitnesses) {
  std::partial_sum(fitnesses.begin(), fitnesses.end(), fitnesses.begin());
  std::vector<uint16_t> selected_indices(simulation_params::population_size);
  std::uniform_real_distribution<double> random_interval(0,fitnesses.back());
  for(uint16_t nth_selection=0; nth_selection<simulation_params::population_size; ++nth_selection) 
    selected_indices[nth_selection]=static_cast<uint16_t>(std::lower_bound(fitnesses.begin(),fitnesses.end(),random_interval(interface_model::RNG_Engine))-fitnesses.begin());
  return selected_indices;
}

void EvolvePopulation(std::string run_details) {
  /* Output files */
  std::string file_base_path="//rscratch//asl47//Bulk_Run//Interfaces//";
  std::string file_simulation_details=std::string(simulation_params::fitness_selection? "S":"R")+"_T"+std::to_string(model_params::temperature)+"_Mu"+std::to_string(model_params::mu_prob)+"_Gamma"+std::to_string(model_params::fitness_factor)+run_details+".txt";
    
  std::ofstream fout_size(file_base_path+"Sizes_"+file_simulation_details, std::ios_base::out);
  std::ofstream fout_strength(file_base_path+"Strengths_"+file_simulation_details, std::ios_base::out);
  std::ofstream fout_fitness(file_base_path+"Fitness_"+file_simulation_details, std::ios_base::out);
  std::ofstream fout_phenotype(file_base_path+"Phenotypes_"+file_simulation_details, std::ios_base::out);
  std::ofstream fout_genotype_history(file_base_path+"GenotypeHistory_"+file_simulation_details, std::ios_base::out);

  interface_model::PhenotypeTable pt = interface_model::PhenotypeTable();
  std::vector<double> population_fitnesses(simulation_params::population_size);
  std::vector<uint32_t> interface_counter(1.5*model_params::interface_size+2);
  bool record_strengths=false;

  /* Genotype initilisation, either zero or random, + duplicate */
  std::vector< std::vector<interface_model::interface_type> > population_genotypes(simulation_params::population_size, std::vector<interface_model::interface_type>(simulation_params::n_tiles*4, 0));
  if(simulation_params::random_initilisation) {
    std::uniform_int_distribution<interface_model::interface_type> dist;
    auto interface_filler = std::bind(dist, interface_model::RNG_Engine);
    for(std::vector<interface_model::interface_type>& genotype : population_genotypes)
      std::generate(genotype.begin(),genotype.end(),interface_filler);   
  }
  std::vector< std::vector<interface_model::interface_type> > reproducing_genotypes(population_genotypes);

  /* Write initial population to file */
  for(std::vector<interface_model::interface_type>& genotype : population_genotypes)
    for(interface_model::interface_type base_value : genotype)
      fout_genotype_history << +base_value << " "; 
  fout_genotype_history<<"\n";
  
  
  

  /* Median time to complete interface mutation, characteristic time for fitness re-assignment */
  std::poisson_distribution<uint16_t> landscape_changer(log(1-pow(2,-1./model_params::interface_size))/log(1-model_params::mu_prob/(4*simulation_params::n_tiles*model_params::interface_size)));
  uint16_t fitness_jiggle=landscape_changer(interface_model::RNG_Engine);

  /* Start main evolution loop */
  for(uint32_t generation=0;generation<simulation_params::generation_limit;++generation) {
    if(fitness_jiggle--==0) {
      pt.ReassignFitness();
      fitness_jiggle=landscape_changer(interface_model::RNG_Engine);
    }
    if(generation+100>=simulation_params::generation_limit)
      record_strengths=true;


    /* Start genotype loop */
    int nth_genotype=0;
    for(std::vector<interface_model::interface_type>& evolving_genotype : population_genotypes) {
      population_fitnesses[nth_genotype++]=interface_model::ProteinAssemblyOutcome(evolving_genotype,&pt);
      if(record_strengths)
        InterfaceStrengths(evolving_genotype,interface_counter);
      interface_model::MutateInterfaces(evolving_genotype);
      for(interface_model::interface_type base_value : evolving_genotype)
        fout_genotype_history << +base_value << " "; 
    } 
    fout_genotype_history<<"\n";
    /* End genotype loop */
    
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
    if(simulation_params::fitness_selection) {
      std::vector<uint16_t> selection_indices=RouletteWheelSelection(population_fitnesses);
      for(uint16_t nth_reproduction = 0; nth_reproduction<simulation_params::population_size;++nth_reproduction)
        reproducing_genotypes[nth_reproduction]=population_genotypes[selection_indices[nth_reproduction]];
      population_genotypes.assign(reproducing_genotypes.begin(),reproducing_genotypes.end());
      for(uint16_t selection_index : selection_indices)
        fout_genotype_history << +selection_index << " ";
      fout_genotype_history<<"\n";
    }
    /* End selection */
    
  } /* End main evolution loop */
  

  
  for(std::unordered_map<uint8_t,std::vector<uint8_t> >::iterator phen_iter=pt.known_phenotypes.begin();phen_iter!=pt.known_phenotypes.end();++phen_iter) {
    uint32_t n_sized_phenotypes=0;
    for(std::vector<uint8_t>::iterator shape_iter=phen_iter->second.begin();shape_iter!=phen_iter->second.end();) {
      ++n_sized_phenotypes;
      fout_phenotype<<+*(shape_iter)<<" "<<+*(shape_iter+1)<<" ";
      for(std::vector<uint8_t>::iterator p_iter=shape_iter+2;p_iter!=shape_iter+*(shape_iter) * *(shape_iter+1)+2;++p_iter)
        fout_phenotype<<+*p_iter<<" ";
      fout_phenotype<<"\n";
      shape_iter+=*(shape_iter) * *(shape_iter+1)+2;
    }
    fout_size <<+phen_iter->first<<" "<<n_sized_phenotypes<<"\n";
  }
  
  
  
  fout_fitness.close();
  fout_size.close();
  fout_strength.close();
  fout_phenotype.close();
  fout_genotype_history.close();

  
}

  
 
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
    //std::vector<uint8_t> s{1,1,1, 0,1,0, 1,1,1};
    //uint8_t dx=3,dy=3;
    //std::cout<<+PhenotypeSymmetryFactor(s,dx,dy)<<std::endl;
    for(auto g : SequenceDifference({0,1,2,3,4,9,6,7},{19,1,2,3,4,5,6,7}))
      std::cout<<+g<<" ";
    std::cout<<std::endl;
  }
  else {
    run_option=argv[1][1];
    SetRuntimeConfigurations(argc,argv);
  }

  switch(run_option) {
  case 'E':
    EvolutionRunner();
    break;
  case 'X':
    std::cout<<"Unused at this time"<<std::endl;
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
      case 'S': simulation_params::fitness_selection=std::stoi(argv[arg+1])>0;break;
      case 'D': simulation_params::independent_trials=std::stoi(argv[arg+1]);break;
      case 'V': simulation_params::run_offset=std::stoi(argv[arg+1]);break;
      case 'R': simulation_params::random_initilisation=std::stoi(argv[arg+1])>0;break; 
        
      case 'F': model_params::fitness_factor=std::stod(argv[arg+1]);break;
      case 'A': model_params::misbinding_rate=std::stod(argv[arg+1]);break;
      case 'M': model_params::mu_prob=std::stod(argv[arg+1]);break;
      case 'T': model_params::temperature=std::stod(argv[arg+1]);break;
      case 'U': model_params::unbound_factor=std::stod(argv[arg+1]);break;
      case 'X': model_params::UND_threshold=std::stod(argv[arg+1]);break;
        
      default: std::cout<<"Unknown Parameter Flag: "<<argv[arg][1]<<std::endl;
      }
    }
    model_params::b_dist.param(std::binomial_distribution<uint8_t>::param_type(model_params::interface_size,model_params::mu_prob/(model_params::interface_size*4*simulation_params::n_tiles)));
  }
}

std::vector<uint8_t> SequenceDifference(std::vector<interface_model::interface_type> parent, std::vector<interface_model::interface_type> child) {
  std::vector<uint8_t> divergence_indices;
  for(uint8_t base=0; base < parent.size(); ++base) {
    if(parent[base]!=child[base])
      divergence_indices.emplace_back(base);
  }
  return divergence_indices;
}
 
