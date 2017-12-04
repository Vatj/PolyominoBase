#include "interface_simulator.hpp"


#include <omp.h>

//typedef uint16_t interface_type;
namespace simulation_params
{
  
  population_size_type population_size=10;
  uint8_t phenotype_builds=10,n_tiles=2;
  uint32_t generation_limit=5;

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


  std::string out_name_f="//rscratch//asl47//Bulk_Run//Interfaces//Test3.txt";
  //std::string out_name_r="//rscratch//asl47//Bulk_Run//Modular//A"+std::to_string(active_targets)+Run_Details+"_Robust.txt";
  //std::ofstream out_file_g(out_name_g, std::ios_base::out);
  std::ofstream out_file_f(out_name_f, std::ios_base::out);
  
  //initialise pool of genomes (all zeroed)
  //std::vector<uint16_t> g{0,0,0,0, 0,0,0,0};
  
  std::vector< std::vector<interface_model::interface_type> > population_genotypes(simulation_params::population_size, std::vector<interface_model::interface_type>(simulation_params::n_tiles*4, 0));
  std::vector< std::vector<interface_model::interface_type> > reproducing_genotypes(simulation_params::population_size, std::vector<interface_model::interface_type>(simulation_params::n_tiles*4, 0));
  std::vector<double> population_fitnesses(simulation_params::population_size);
  std::vector<double> symmetries;//(simulation_params::generation_limit);
  
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
  //{0,0,0,0, 65535,31,3,31, 31,31,100,16383, 31,31,31,55807};

  interface_model::PhenotypeTable pt = interface_model::PhenotypeTable();
  for(uint32_t generation=0;generation<simulation_params::generation_limit;++generation) {
    //std::cout<<"generation: "<<generation<<std::endl;
    int nth_genotype=0;
    for(std::vector< std::vector<interface_model::interface_type> >::iterator evolving_genotype_iter=population_genotypes.begin(); evolving_genotype_iter!=population_genotypes.end();++evolving_genotype_iter) {
      population_fitnesses[nth_genotype++]=interface_model::ProteinAssemblyOutcome(*evolving_genotype_iter,simulation_params::phenotype_builds,&pt);
      interface_model::MutateInterfaces(*evolving_genotype_iter);

      

    }
    std::vector<double> syms;
    for(std::vector< std::vector<interface_model::interface_type> >::iterator evolving_genotype_iter=population_genotypes.begin(); evolving_genotype_iter!=population_genotypes.end();++evolving_genotype_iter) {
      for(std::vector<interface_model::interface_type>::iterator face = evolving_genotype_iter->begin();face!=evolving_genotype_iter->end();++face) {
        //std::cout<<*face<<", "<< interface_model::SymmetryFactor(*face)<<" ";//std::endl;
        syms.emplace_back(interface_model::SymmetryFactor(*face));
      }
      //std::cout<<std::endl;
      

      }
    symmetries.emplace_back(std::accumulate(syms.begin(),syms.end(),0.0)/(4*simulation_params::n_tiles*simulation_params::population_size));
    /*
    std::cout<<"fitness: ";
    for(double f : population_fitnesses)
      std::cout<<f<<" ";

    std::cout<<"\nselection: ";
    */
    std::vector<simulation_params::population_size_type> selection_indices=RouletteWheelSelection(population_fitnesses);
    for(simulation_params::population_size_type nth_reproduction = 0; nth_reproduction<simulation_params::population_size;++nth_reproduction) {
      //std::cout<<+selection_indices[nth_reproduction]<<" ";
      reproducing_genotypes[nth_reproduction]=population_genotypes[selection_indices[nth_reproduction]];
    }
    population_genotypes.assign(reproducing_genotypes.begin(),reproducing_genotypes.end());
    //std::cout<<std::endl;
    //interface_model::ProteinAssemblyOutcome(g,10,&pt)
    //std::cout<<"fitness "<<<<std::endl;
    
    //interface_model::MutateInterfaces(g);
  
  }

  /* TO DO

     SELECTION
     
     SYMMETRY OF INTERFACE

  
  for(double s:symmetries)
    out_file_f<<s<<" ";
  out_file_f.close();
  */
  
  std::cout<<"printng table"<<std::endl;
  std::cout<<"N_P "<<pt.n_phenotypes<<std::endl;
  std::vector<uint8_t> temp_phenotype;
  for(std::vector<uint8_t>::iterator phen_iter = pt.known_phenotypes.begin();phen_iter!=pt.known_phenotypes.end();) {
    out_file_f<<static_cast<int>(*phen_iter)<<" "<<static_cast<int>(*(phen_iter+1))<<" ";
    temp_phenotype.assign(phen_iter+2,phen_iter+2+*phen_iter* *(phen_iter+1));
    for(uint8_t m :temp_phenotype)
      out_file_f<<static_cast<int>(m)<<" ";
    out_file_f<<"\n";
    //std::cout<<std::endl;
    // if(ComparePolyominoes(phenotype,dx,dy,temp_phenotype,*phen_iter,*(phen_iter+1)))
    //     return phenotype_ID;
        phen_iter+=*phen_iter* *(phen_iter+1)+2;
        //++phenotype_ID;
      }    

  

  //std::cout<<"N_P "<<pt.n_phenotypes<<std::endl;
  out_file_f.close();
}

void EvolvePopulation_Off() {
  std::string out_name_f="Test3.txt";
  std::ofstream out_file_f(out_name_f, std::ios_base::out);
  
  std::vector< std::vector<interface_model::interface_type> > population_genotypes(simulation_params::population_size, std::vector<interface_model::interface_type>(simulation_params::n_tiles*4, 0));
  std::vector< std::vector<interface_model::interface_type> > reproducing_genotypes(simulation_params::population_size, std::vector<interface_model::interface_type>(simulation_params::n_tiles*4, 0));
  std::vector<double> population_fitnesses(simulation_params::population_size);
  std::vector<double> symmetries;
  

  interface_model::PhenotypeTable pt = interface_model::PhenotypeTable();
  for(uint32_t generation=0;generation<simulation_params::generation_limit;++generation) {
    //std::cout<<"generation: "<<generation<<std::endl;
    int nth_genotype=0;
    for(std::vector< std::vector<interface_model::interface_type> >::iterator evolving_genotype_iter=population_genotypes.begin(); evolving_genotype_iter!=population_genotypes.end();++evolving_genotype_iter) {
      population_fitnesses[nth_genotype++]=1;//
      interface_model::ProteinAssemblyOutcome(*evolving_genotype_iter,simulation_params::phenotype_builds,&pt);
      interface_model::MutateInterfaces(*evolving_genotype_iter);

      

    }
    std::vector<double> syms;
    for(std::vector< std::vector<interface_model::interface_type> >::iterator evolving_genotype_iter=population_genotypes.begin(); evolving_genotype_iter!=population_genotypes.end();++evolving_genotype_iter) {
      for(std::vector<interface_model::interface_type>::iterator face = evolving_genotype_iter->begin();face!=evolving_genotype_iter->end();++face) {
        //std::cout<<*face<<", "<< interface_model::SymmetryFactor(*face)<<" ";//std::endl;
        syms.emplace_back(interface_model::SymmetryFactor(*face));
      }
      //std::cout<<std::endl;
      

      }
    symmetries.emplace_back(std::accumulate(syms.begin(),syms.end(),0.0)/(4*simulation_params::n_tiles*simulation_params::population_size));

    std::vector<simulation_params::population_size_type> selection_indices=RouletteWheelSelection(population_fitnesses);
    for(simulation_params::population_size_type nth_reproduction = 0; nth_reproduction<simulation_params::population_size;++nth_reproduction) {
      reproducing_genotypes[nth_reproduction]=population_genotypes[selection_indices[nth_reproduction]];
    }
    population_genotypes.assign(reproducing_genotypes.begin(),reproducing_genotypes.end());
  
  }

  /* TO DO

     SELECTION
     
     SYMMETRY OF INTERFACE

  
  for(double s:symmetries)
    out_file_f<<s<<" ";
  out_file_f.close();
  */
  
  std::cout<<"printng table"<<std::endl;
  std::cout<<"N_P "<<pt.n_phenotypes<<std::endl;
  for(double i : pt.interface_strengths)
    std::cout<<i<<" ";
  std::cout<<std::endl;
  
  std::vector<uint8_t> temp_phenotype;
  for(std::vector<uint8_t>::iterator phen_iter = pt.known_phenotypes.begin();phen_iter!=pt.known_phenotypes.end();) {
    out_file_f<<static_cast<int>(*phen_iter)<<" "<<static_cast<int>(*(phen_iter+1))<<" ";
    temp_phenotype.assign(phen_iter+2,phen_iter+2+*phen_iter* *(phen_iter+1));
    for(uint8_t m :temp_phenotype)
      out_file_f<<static_cast<int>(m)<<" ";
    out_file_f<<"\n";
    //std::cout<<std::endl;
    // if(ComparePolyominoes(phenotype,dx,dy,temp_phenotype,*phen_iter,*(phen_iter+1)))
    //     return phenotype_ID;
        phen_iter+=*phen_iter* *(phen_iter+1)+2;
        //++phenotype_ID;
      }    

  

  //std::cout<<"N_P "<<pt.n_phenotypes<<std::endl;
  out_file_f.close();
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
  switch(run_option) {
  case 'R':
    EvolvePopulation_Off();
    break;
  case 'E':
    EvolvePopulation();
    break;
  case 'H':
  default:
    std::cout<<"Protein interface model\n**Simulation Parameters**\nN: number of tiles\nP: population size\nK: generation limit\nB: number of phenotype builds\n";
    std::cout<<"\n**Model Parameters**\nU: mutation probability (per interface)\nT: temperature\nI: unbound size factor\nA: misbinding rate\nM: Fitness factor\n";
    std::cout<<"\n**Run options**\nR: evolution without fitness\nE: evolution with fitness\n";
    break;

  }
  //model_params::temperature=std::stod(argv[1]);
  //simulation_params::population_size=std::stoi(argv[2]);
  //model_params::mu_prob=std::stod(argv[3]);
  
  //EvolvePopulation();
  //uint16_t x=49164;
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


void SetRuntimeConfigurations(int argc, char* argv[]) {
  if(argc<3 && argv[1][1]!='H') {
    std::cout<<"Invalid Parameters"<<std::endl;
  }
  else {
    for(uint8_t arg=2;arg<argc;arg+=2) {
      switch(argv[arg][1]) {
        //BASIC PARAMETERS//
      case 'N': simulation_params::n_tiles=std::stoi(argv[arg+1]);break;
      case 'P': simulation_params::population_size=std::stoi(argv[arg+1]);break;
      case 'K': simulation_params::generation_limit=std::stoi(argv[arg+1]);break;
      case 'B': simulation_params::phenotype_builds=std::stoi(argv[arg+1]);break;
        
      case 'M': model_params::fitness_factor=std::stod(argv[arg+1]);break;
      case 'A': model_params::misbinding_rate=std::stod(argv[arg+1]);break;
      case 'U': model_params::mu_prob=std::stod(argv[arg+1]);break;
      case 'T': model_params::temperature=std::stod(argv[arg+1]);break;
      case 'I': model_params::unbound_factor=std::stod(argv[arg+1]);break;


        //Default//
      default: std::cout<<"Unknown Parameter Flag: "<<argv[arg][1]<<std::endl;
      }
    }
    model_params::b_dist.param(std::binomial_distribution<uint8_t>::param_type(model_params::interface_size,model_params::mu_prob));
  }
}
