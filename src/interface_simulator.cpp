#include "interface_simulator.hpp"


#include <omp.h>

//typedef uint16_t interface_type;
namespace simulation_params
{
  
  population_size_type population_size=100;
  uint8_t phenotype_builds=10;
  

}

std::vector<simulation_params::population_size_type> RouletteWheelSelection(std::vector<double>& fitnesses) {
  std::partial_sum(fitnesses.begin(), fitnesses.end(), fitnesses.begin());
  std::vector<simulation_params::population_size_type> selected_indicies(simulation_params::population_size);
  std::uniform_real_distribution<double> random_interval(0,*(fitnesses.end()-1));
  for(simulation_params::population_size_type nth_selection=0; nth_selection<simulation_params::population_size; ++nth_selection)
    selected_indicies[nth_selection]=static_cast<uint16_t>(std::lower_bound(fitnesses.begin(),fitnesses.end(),random_interval(interface_model::RNG_Engine))-fitnesses.begin());
  return selected_indicies;
}

void EvolvePopulation() {

  //initialise pool of genomes (all zeroed)
  uint32_t GENERATION_LIMIT=5;
  
  std::vector< std::vector<interface_model::interface_type> > population_genotypes(simulation_params::population_size, std::vector<interface_model::interface_type>(8, 0));
  std::vector< std::vector<interface_model::interface_type> > reproducing_genotypes(simulation_params::population_size, std::vector<interface_model::interface_type>(8, 0));
  std::vector<double> population_fitnesses(simulation_params::population_size);

  
  /*
  for(std::vector<interface_model::interface_type>& g : v) {
    for(interface_model::interface_type b : g) {
      std::cout<<b<<" ";
    }
    std::cout<<std::endl;
    for(int x=0;x<2;x++) {
    interface_model::MutateInterfaces(g);
    for(interface_model::interface_type b : g) {
      std::cout<<b<<" ";
    }
    std::cout<<std::endl;
  }
  }
  */
  std::vector<uint16_t> g{0,0,0,0, 0,0,0,0};//{0,0,0,0, 65535,31,3,31, 31,31,100,16383, 31,31,31,55807};

  interface_model::PhenotypeTable pt = interface_model::PhenotypeTable();
  for(uint32_t generation=0;generation<GENERATION_LIMIT;++generation) {
    int nth_genotype=0;
    for(std::vector< std::vector<interface_model::interface_type> >::iterator evolving_genotype_iter=population_genotypes.begin(); evolving_genotype_iter!=population_genotypes.end();++evolving_genotype_iter) {
      population_fitnesses[nth_genotype++]=interface_model::ProteinAssemblyOutcome(*evolving_genotype_iter,simulation_params::phenotype_builds,&pt);
      interface_model::MutateInterfaces(*evolving_genotype_iter);

      

    }
    std::cout<<"fitness: ";
    for(double f : population_fitnesses)
      std::cout<<f<<" ";

    std::cout<<"\nselection: ";
    std::vector<simulation_params::population_size_type> selection_indices=RouletteWheelSelection(population_fitnesses);
    for(simulation_params::population_size_type nth_reproduction = 0; nth_reproduction<simulation_params::population_size;++nth_reproduction) {
      std::cout<<+selection_indices[nth_reproduction]<<" ";
      reproducing_genotypes[nth_reproduction]=population_genotypes[selection_indices[nth_reproduction]];
    }
    population_genotypes.assign(reproducing_genotypes.begin(),reproducing_genotypes.end());
    std::cout<<std::endl;
    //interface_model::ProteinAssemblyOutcome(g,10,&pt)
    //std::cout<<"fitness "<<<<std::endl;
    
    //interface_model::MutateInterfaces(g);
  
  }

  /* TO DO

     SELECTION
     
     SYMMETRY OF INTERFACE

   */
  
  

  std::cout<<"N_P "<<pt.n_phenotypes<<std::endl;

}
  
 



int main(int argc, char* argv[]) {
  model_params::temperature=std::stod(argv[1]);
  simulation_params::population_size=std::stoi(argv[2]);
  model_params::mu_prob=std::stod(argv[3]);
  model_params::b_dist.param(std::binomial_distribution<uint8_t>::param_type(model_params::interface_size,model_params::mu_prob));
  EvolvePopulation();
  uint16_t x=49164;
  //std::cout<<+interface_model::SymmetryFactor(x)<<std::endl;
  return 0;

  //std::vector<uint16_t> g{0,0,0,0, 65535,31,3,31, 31,31,100,16383, 31,31,31,55807};
  




  

  // 0 binds with 65535
  // 3 binds with 16383
  // 100 binds with 55807
  // 31 is fairly neutral
  /*
  if(argc>1)
    std::cout<<interface_model::ProteinAssemblyOutcome(g,10)<<std::endl;
  else
    std::cout<<interface_model::ProteinAssemblyOutcome(g,10)<<std::endl;
  

  int x=0,f1=0,f2=0,f3=0,f4=0;
#pragma omp parallel for schedule(static, 1000) default(none) firstprivate(g) num_threads(4) reduction(+:x,f1,f2,f3,f4)
  for(uint32_t j = 0; j < 40000; ++j) {

    //std::cout<<"j "<<j<<" t "<<omp_get_thread_num()<<" : ";
    //for(auto f: g)
    //   std::cout<<+f<<" ";
    // std::cout<<std::endl;
    
      if(interface_model::ProteinAssemblyOutcome(g,10)>.5) {
      ++x;
      //std::cout<<"j "<<j<<" f "<<omp_get_thread_num()<<" win"<<std::endl;
      }
    else {
      //std::cout<<"j "<<j<<" f "<<omp_get_thread_num()<<"\n";
      switch(omp_get_thread_num()) {
      case 0:
        f1++;
        break;
      case 1:
        f2++;
        break;
      case 2:
        f3++;
        break;
      case 3:
        f4++;
        break;
      }
      //for(auto f: model_params::faces)
      //  std::cout<<+f<<" ";
      //std::cout<<std::endl;
    }
  }
  std::cout<<x<<std::endl;
  std::cout<<"fails "<<"T0: "<<f1<<" T1: "<<f2<<" T2: "<<f3<<" T3: "<<f4<<std::endl;
  */
}
