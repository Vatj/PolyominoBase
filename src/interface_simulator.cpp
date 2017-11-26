//#include "interface_simulator.hpp"
#include "interface_model.hpp"

#include <omp.h>

//typedef uint16_t interface_type;

void EvolvePopulation(const int pop_size) {

  //initialise pool of genomes (all zeroed)
  uint32_t GENERATION_LIMIT=10;
  
  std::vector< std::vector<interface_model::interface_type> > v(pop_size, std::vector<interface_model::interface_type>(16, 0));
  std::vector<double> population_fitnesses(pop_size);

  
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
    std::cout<<"fitness "<<interface_model::ProteinAssemblyOutcome(g,10,&pt)<<std::endl;
    interface_model::MutateInterfaces(g);
  
  }

  /* TO DO

     STORE FITNESS
     SELECTION
     
     SYMMETRY OF INTERFACE

   */
  
  

  std::cout<<"N_P "<<pt.n_phenotypes<<std::endl;

}
  
 



int main(int argc, char* argv[]) {
  params::temperature=std::stod(argv[1]);
  EvolvePopulation(std::stoi(argv[2]));
  return 0;

  std::vector<uint16_t> g{0,0,0,0, 65535,31,3,31, 31,31,100,16383, 31,31,31,55807};
  




  

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
      //for(auto f: params::faces)
      //  std::cout<<+f<<" ";
      //std::cout<<std::endl;
    }
  }
  std::cout<<x<<std::endl;
  std::cout<<"fails "<<"T0: "<<f1<<" T1: "<<f2<<" T2: "<<f3<<" T3: "<<f4<<std::endl;
  */
}
