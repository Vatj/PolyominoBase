#include "interface_simulator.hpp"

namespace simulation_params
{
  population_size_type population_size=10;
  uint8_t phenotype_builds=10,n_tiles=2;
  uint32_t generation_limit=5,independent_trials=1,run_offset=0;
  bool fitness_selection=false;
}

std::vector<simulation_params::population_size_type> RouletteWheelSelection(std::vector<double>& fitnesses) {
  std::partial_sum(fitnesses.begin(), fitnesses.end(), fitnesses.begin());
  std::vector<simulation_params::population_size_type> selected_indices(simulation_params::population_size);
  std::uniform_real_distribution<double> random_interval(0,*(fitnesses.end()-1));
  for(simulation_params::population_size_type nth_selection=0; nth_selection<simulation_params::population_size; ++nth_selection) 
    selected_indices[nth_selection]=static_cast<uint16_t>(std::lower_bound(fitnesses.begin(),fitnesses.end(),random_interval(interface_model::RNG_Engine))-fitnesses.begin());
  return selected_indices;
}

void EvolvePopulation(std::string run_details) {
  std::string file_base_path="";//"//rscratch//asl47//Bulk_Run//Interfaces//";
  std::string out_name_strength=file_base_path+"Strengths_"+std::string(simulation_params::fitness_selection? "S":"R")+"_T"+std::to_string(model_params::temperature)+"_Mu"+std::to_string(model_params::mu_prob)+"_Gamma"+std::to_string(model_params::fitness_factor)+run_details+".txt";
  std::string out_name_size=file_base_path+"//rscratch//asl47//Bulk_Run//Interfaces//Sizes_"+std::string(simulation_params::fitness_selection? "S":"R")+"_T"+std::to_string(model_params::temperature)+"_Mu"+std::to_string(model_params::mu_prob)+"_Gamma"+std::to_string(model_params::fitness_factor)+run_details+".txt";
  std::string out_name_fitness=file_base_path+"//rscratch//asl47//Bulk_Run//Interfaces//Fitness_"+std::string(simulation_params::fitness_selection? "S":"R")+"_T"+std::to_string(model_params::temperature)+"_Mu"+std::to_string(model_params::mu_prob)+"_Gamma"+std::to_string(model_params::fitness_factor)+run_details+".txt";
  
  std::ofstream fout_size(out_name_size, std::ios_base::out);
  std::ofstream fout_strength(out_name_strength, std::ios_base::out);
  std::ofstream fout_fitness(out_name_fitness, std::ios_base::out);
  
  std::vector< std::vector<interface_model::interface_type> > population_genotypes(simulation_params::population_size, std::vector<interface_model::interface_type>(simulation_params::n_tiles*4, 0));
  std::vector< std::vector<interface_model::interface_type> > reproducing_genotypes(simulation_params::population_size, std::vector<interface_model::interface_type>(simulation_params::n_tiles*4, 0));
  std::vector<double> population_fitnesses(simulation_params::population_size);
  
  interface_model::PhenotypeTable pt = interface_model::PhenotypeTable();
  for(uint32_t generation=0;generation<simulation_params::generation_limit;++generation) {

    std::vector<uint32_t> interface_counter(1.5*model_params::interface_size+2);
    int nth_genotype=0;
    for(std::vector< std::vector<interface_model::interface_type> >::iterator evolving_genotype_iter=population_genotypes.begin(); evolving_genotype_iter!=population_genotypes.end();++evolving_genotype_iter) {

      std::transform(interface_counter.begin(), interface_counter.end(),InterfaceStrengths(*evolving_genotype_iter).begin() , interface_counter.begin(),std::plus<uint32_t>());
      population_fitnesses[nth_genotype++]=interface_model::ProteinAssemblyOutcome(*evolving_genotype_iter,simulation_params::phenotype_builds,&pt);
      interface_model::MutateInterfaces(*evolving_genotype_iter);

    }

    for(uint32_t count : interface_counter)
      fout_strength<<count<<" ";
    fout_strength<<"\n";

    
    double mu=0,sigma=0;
    DistributionStatistics(population_fitnesses,mu,sigma);
    fout_fitness<<mu<<" "<<sigma<<"\n";

    if(simulation_params::fitness_selection) {
      std::vector<simulation_params::population_size_type> selection_indices=RouletteWheelSelection(population_fitnesses);
      for(simulation_params::population_size_type nth_reproduction = 0; nth_reproduction<simulation_params::population_size;++nth_reproduction)
        reproducing_genotypes[nth_reproduction]=population_genotypes[selection_indices[nth_reproduction]];
      population_genotypes.assign(reproducing_genotypes.begin(),reproducing_genotypes.end());
    }
  }
  
  std::cout<<"printing table"<<std::endl;
  std::cout<<"N_P "<<pt.n_phenotypes<<std::endl;
  for(std::unordered_map<uint8_t,std::vector<uint8_t> >::iterator phen_iter=pt.known_phenotypes.begin();phen_iter!=pt.known_phenotypes.end();++phen_iter) {
    uint32_t n_sized_phenotypes=0;
    for(std::vector<uint8_t>::iterator shape_iter=phen_iter->second.begin();shape_iter!=phen_iter->second.end();) {
      ++n_sized_phenotypes;
      shape_iter+=*(shape_iter) * *(shape_iter+1)+2;
    }
    fout_size <<+phen_iter->first<<" "<<n_sized_phenotypes<<"\n";
  }

  /*
  for(std::vector<uint8_t>::iterator phen_iter = pt.known_phenotypes.begin();phen_iter!=pt.known_phenotypes.end();) {
  std::vector<uint8_t> temp_phenotype;
  for(std::vector<uint8_t>::iterator phen_iter = pt.known_phenotypes.begin();phen_iter!=pt.known_phenotypes.end();) {
    out_file_r<<static_cast<int>(*phen_iter)<<" "<<static_cast<int>(*(phen_iter+1))<<" ";
    temp_phenotype.assign(phen_iter+2,phen_iter+2+*phen_iter* *(phen_iter+1));
    for(uint8_t m :temp_phenotype)
      out_file_r<<static_cast<int>(m)<<" ";
    out_file_r<<"\n";
    phen_iter+=*phen_iter* *(phen_iter+1)+2;

  }
  */    
  
  
  fout_fitness.close();
  fout_size.close();
  fout_strength.close();
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
    std::vector<uint8_t> s{1,1,1, 0,1,0, 1,1,1};
    uint8_t dx=3,dy=3;
    std::cout<<+PhenotypeSymmetryFactor(s,dx,dy)<<std::endl;
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
    //RandomStrings();
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
        
      case 'F': model_params::fitness_factor=std::stod(argv[arg+1]);break;
      case 'A': model_params::misbinding_rate=std::stod(argv[arg+1]);break;
      case 'M': model_params::mu_prob=std::stod(argv[arg+1]);break;
      case 'T': model_params::temperature=std::stod(argv[arg+1]);break;
      case 'U': model_params::unbound_factor=std::stod(argv[arg+1]);break;
      case 'X': model_params::UND_threshold=std::stod(argv[arg+1]);break;
        
      default: std::cout<<"Unknown Parameter Flag: "<<argv[arg][1]<<std::endl;
      }
    }
    model_params::b_dist.param(std::binomial_distribution<uint8_t>::param_type(model_params::interface_size,model_params::mu_prob));
  }
}
