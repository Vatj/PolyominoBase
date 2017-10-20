#include <evolution_simulator.hpp>

///////////////////////////////////////
//GENOME EVOLTUTION RELATED FUNCTIONS//
///////////////////////////////////////

void FitnessEvolutionDynamic() {
  double mu=1.0/(Num_Tiles*mutation_cofactor);
#pragma omp parallel for schedule(dynamic)
  for(int r=0;r<Num_Runs;++r) {
    std::string runN="_T"+std::to_string(Num_Tiles)+"_C"+std::to_string(Colour_Space+1)+"_N"+std::to_string(Num_Genomes)+"_Mu"+std::to_string(mu)+"_O"+std::to_string(Fitness_Oscillation_Rate)+"_K"+std::to_string(GENERATION_LIMIT)+"_I"+std::to_string(initial_condition)+"_Run"+std::to_string(r+RUN);
    EvolveFitnessDynamic(runN,mu);
  }
}

void EvolveFitnessDynamic(std::string Run_Details,double Mu) {
  std::string out_name_g="//rscratch//asl47//Bulk_Run//Modular//A"+std::to_string(active_targets)+Run_Details+"_Genotype.txt";
  std::string out_name_f="//rscratch//asl47//Bulk_Run//Modular//A"+std::to_string(active_targets)+Run_Details+"_Fitness.txt";
  //std::string out_name_r="//rscratch//asl47//Bulk_Run//Modular//A"+std::to_string(active_targets)+Run_Details+"_Robust.txt";
  std::ofstream out_file_g(out_name_g, std::ios_base::out);
  std::ofstream out_file_f(out_name_f, std::ios_base::out);
  //std::ofstream out_file_r(out_name_r, std::ios_base::out);

  const int FITNESS_CLIFF=2;
  const int TOTAL_TARGETS=3;
  
  std::uniform_int_distribution<int> Mutated_Colour(0,Colour_Space);
  std::bernoulli_distribution Mutation_Chance(Mu);
  std::vector<int> Initial_Genome(Num_Tiles*4,0),Evolving_Genome(Num_Tiles*4),target_types(TOTAL_TARGETS), Index_Selections(Num_Genomes);
  
  //InitialGenotypeConfiguration(Initial_Genome);
  std::vector<std::vector<int> > Genome_Pool(Num_Genomes,Initial_Genome), Temporary_Pool(Num_Genomes);
  std::vector<double> Phenotype_Fitness_Sizes(Num_Genomes), Phenotype_Fitness_Sizes_Full(Num_Genomes),target_fitnesses(TOTAL_TARGETS);
  std::set< std::vector<int> > unique_genotypes;

  int successive_generations=0;
  int writing_generation=GENERATION_LIMIT;
  
  for(int g=1;g<=GENERATION_LIMIT;++g) {
    switch(((g-1)/Fitness_Oscillation_Rate)%3) {
    case 0:
      target_types={4,5,6};
      break;
    case 1:
      target_types={5,6,4};
      break;
    case 2:
      target_types={6,4,5};
      break;
    }
    for(int n=0;n<Num_Genomes;++n) {
      Evolving_Genome=Genome_Pool[n];
      for(int t=0;t<Num_Tiles*4;++t) {
        if(Mutation_Chance(RNG_Engine)) {
          int previousFace=Evolving_Genome[t];
          do {
            Evolving_Genome[t]=Mutated_Colour(RNG_Engine);
          }while(Evolving_Genome[t]==previousFace);
        }
      }
      if(GetMultiplePhenotypeFitness(Evolving_Genome,target_types,target_fitnesses,active_targets)) {    
        Phenotype_Fitness_Sizes[n]=std::accumulate(target_fitnesses.begin(),target_fitnesses.begin()+active_targets,0.0,[](double cum,double nex){return cum+Fitness_Function(nex);})/(active_targets*pow(FITNESS_CLIFF,active_targets-std::count(target_fitnesses.begin(),target_fitnesses.begin()+active_targets,1.)));
        Phenotype_Fitness_Sizes_Full[n]=std::accumulate(target_fitnesses.begin(),target_fitnesses.end(),0.0);
        Phenotype_Fitness_Sizes_Full[n]=Phenotype_Fitness_Sizes_Full[n]>=3. ? 1. : 0;
      }
      else {
        Phenotype_Fitness_Sizes[n]=0;
        Phenotype_Fitness_Sizes_Full[n]=0;
      }
      Genome_Pool[n]=Evolving_Genome;
      
      if(Phenotype_Fitness_Sizes_Full[n]>=1.) {
        std::vector<int> duplicate_genome=(Evolving_Genome);
        Clean_Genome(duplicate_genome);
        if(Disjointed_Check(duplicate_genome))
          duplicate_genome.erase(std::find(duplicate_genome.begin(),duplicate_genome.end(),-1),duplicate_genome.end());
        unique_genotypes.insert(duplicate_genome);
      }
      
    } //END GENOME LOOP

    double max_fitness=*std::max_element(Phenotype_Fitness_Sizes_Full.begin(),Phenotype_Fitness_Sizes_Full.end());
    if(max_fitness>=1.) {
      ++successive_generations;
      out_file_f << max_fitness << " " << std::count(Phenotype_Fitness_Sizes_Full.begin(),Phenotype_Fitness_Sizes_Full.end(),max_fitness) << "\n";
    }
    else {
      out_file_f << "0 0\n";
      successive_generations=0;
    }
    /*
    if(g%5000==0) {
      out_file_g << "g: "<<g<<"\n";
      for(int nth_genome=0;nth_genome<Num_Genomes;++nth_genome) {
        if(Phenotype_Fitness_Sizes_Full[nth_genome]>1.0-MINIMUM_FITNESS_THRESHOLD) {
          std::vector<int> dup_g(Genome_Pool[nth_genome]);
          Clean_Genome(dup_g);
          if(Disjointed_Check(dup_g))
             dup_g.erase(std::find(dup_g.begin(),dup_g.end(),-1),dup_g.end());
          for(int base : dup_g) {
            out_file_g << base <<" ";
          }
          out_file_g << "\n";
        }
      }
    }
    */
    /*
    if(false && successive_generations>100 && writing_generation==GENERATION_LIMIT) {
      writing_generation=(g/Fitness_Oscillation_Rate+1)*Fitness_Oscillation_Rate;
      std::cout<<"fixed at "<<g<<" and writing at "<< writing_generation<<std::endl;
      successive_generations=0;
    }
    if(writing_generation!=GENERATION_LIMIT && g>=writing_generation && g<writing_generation+4*Fitness_Oscillation_Rate && (g-writing_generation)%(Fitness_Oscillation_Rate/5)==0) {
      for(int nth=0;nth<Num_Genomes;++nth) {
        if(Phenotype_Fitness_Sizes_Full[nth]>=1.0) {
          for(int base : Genome_Pool[nth]) {
          out_file_r <<  base << " ";
          }
          out_file_r << "\n";
        }
      }
      out_file_r << "x,"<<g<<"\n";
    }
    */
    //if(g==writing_generation+9*Fitness_Oscillation_Rate)
    //  return;
    
    
    Index_Selections=Roulette_Wheel_Selection(Phenotype_Fitness_Sizes,Num_Genomes);
    for(int nth_select=0;nth_select<Num_Genomes;++nth_select) {
      Temporary_Pool[nth_select]=Genome_Pool[Index_Selections[nth_select]];
    }   
    Genome_Pool.assign(Temporary_Pool.begin(),Temporary_Pool.end());
  }
  
  for(std::vector<int> genotype : unique_genotypes) {
    for(int base : genotype) {
      out_file_g << base <<" ";
    }
    out_file_g << "\n";
  }
  
  out_file_g.close();
  out_file_f.close();
}

void EvolveRegulated(int& Discovery_Generation, int& Adaptation_Generation,double Mu) {
  std::uniform_int_distribution<int> Mutated_Colour(0,Colour_Space);
  std::bernoulli_distribution Mutation_Chance(Mu);
  std::vector<int> Initial_Genome(Num_Tiles*5), Evolving_Genome(Num_Tiles*5), Index_Selections(Num_Genomes),regulated_Genome;
  

  InitialGenotypeConditions(Initial_Genome);

  std::vector<std::vector<int> > Genome_Pool(Num_Genomes,Initial_Genome), Temporary_Pool(Num_Genomes);
  std::vector<double> Phenotype_Fitness_Sizes(Num_Genomes);
  

  for(int g=1;g<=GENERATION_LIMIT;++g) {
    int Num_Maximally_Fit=0;
    for(int n=0;n<Num_Genomes;++n) {
      Evolving_Genome=Genome_Pool[n];
      for(int t=0;t<Num_Tiles*5;++t) {
          if(Mutation_Chance(RNG_Engine)) {
            if(t%5==0) //Transcription Factor
              Evolving_Genome[t]^=1;
            else {
              int previousFace=Evolving_Genome[t];
              do {
                Evolving_Genome[t]=Mutated_Colour(RNG_Engine);
              } while(Evolving_Genome[t]==previousFace);
            }
          }
        }
      regulated_Genome.clear();
      for(std::vector<int>::iterator evol_it=Evolving_Genome.begin();evol_it!=Evolving_Genome.end();evol_it+=5) {
        if(*evol_it==1) 
          regulated_Genome.insert(regulated_Genome.end(),evol_it+1,evol_it+5);
      }
      if(regulated_Genome.empty())
        Phenotype_Fitness_Sizes[n]=0;
      else {
        Phenotype_Fitness_Sizes[n]=Fitness_Function(Get_Phenotype_Fitness(regulated_Genome,0,true));          
      }
      if(Phenotype_Fitness_Sizes[n]>=1.-MINIMUM_FITNESS_THRESHOLD) {
        ++Num_Maximally_Fit;
        if(Discovery_Generation==-1) //DISCOVERY
          Discovery_Generation=g;
      }
      Genome_Pool[n]=Evolving_Genome;
    }
    if(Num_Maximally_Fit>=std::floor(Num_Genomes/2.)) { //ADAPTATION
      Adaptation_Generation=g;
      return;
    }
    Index_Selections=Roulette_Wheel_Selection(Phenotype_Fitness_Sizes,Num_Genomes);
    for(int nth_select=0;nth_select<Num_Genomes;++nth_select) {
      Temporary_Pool[nth_select]=Genome_Pool[Index_Selections[nth_select]];
    }   
    Genome_Pool.assign(Temporary_Pool.begin(),Temporary_Pool.end());
  }
}

void EvolveSimple(int& Discovery_Generation, int& Adaptation_Generation,double Mu) {
  std::uniform_int_distribution<int> Mutated_Colour(0,Colour_Space);
  std::bernoulli_distribution Mutation_Chance(Mu);
  std::vector<int> Initial_Genome(Num_Tiles*4), Evolving_Genome(Num_Tiles*4), Index_Selections(Num_Genomes);

  InitialGenotypeConditions(Initial_Genome);
  
  std::vector<std::vector<int> > Genome_Pool(Num_Genomes,Initial_Genome),Temporary_Pool(Num_Genomes);
  std::vector<double> Phenotype_Fitness_Sizes(Num_Genomes);
 
  for(int g=1;g<=GENERATION_LIMIT;++g) {
    int Num_Maximally_Fit=0;
    for(int n=0;n<Num_Genomes;++n) {
      Evolving_Genome=Genome_Pool[n];
      for(int t=0;t<Num_Tiles*4;++t) {
        if(Mutation_Chance(RNG_Engine)) {
          int previousFace=Evolving_Genome[t];
          do {
            Evolving_Genome[t]=Mutated_Colour(RNG_Engine);
          } while(Evolving_Genome[t]==previousFace);
        }
      }
      Phenotype_Fitness_Sizes[n]=Fitness_Function(Get_Phenotype_Fitness(Evolving_Genome,0,true));          
      if(Phenotype_Fitness_Sizes[n]>=1.-MINIMUM_FITNESS_THRESHOLD) {
        ++Num_Maximally_Fit;
        if(Discovery_Generation==-1) //DISCOVERY
          Discovery_Generation=g;
      }
      Genome_Pool[n]=Evolving_Genome;
    }
    if(Num_Maximally_Fit>=std::floor(Num_Genomes/2.)) { //ADAPTATION
      Adaptation_Generation=g;
      return;
    }
    Index_Selections=Roulette_Wheel_Selection(Phenotype_Fitness_Sizes,Num_Genomes);
    for(int nth_select=0;nth_select<Num_Genomes;++nth_select) {
      Temporary_Pool[nth_select]=Genome_Pool[Index_Selections[nth_select]];
    }   
    Genome_Pool.assign(Temporary_Pool.begin(),Temporary_Pool.end());
  }
}


void EvolutionSimulation(const double Mu,std::vector<int>& discovery,std::vector<int>& adaptation) {
#pragma omp parallel for schedule(dynamic) 
  for(int run =0;run<Num_Runs;++run) {
    int Discovery_Generation=-1, Adaptation_Generation=-1;
    if(Regulated)
      EvolveRegulated(Discovery_Generation,Adaptation_Generation,Mu);
    else
      EvolveSimple(Discovery_Generation,Adaptation_Generation,Mu);                               
    discovery[run]=Discovery_Generation>0 ? Discovery_Generation : GENERATION_LIMIT;
    adaptation[run]=Adaptation_Generation>0 ? Adaptation_Generation : GENERATION_LIMIT;;
  }
}

void ManyEvolutionSimulations() {
  std::string fileName="//rscratch//asl47//Bulk_Run//Configs//MuL_Values.txt";
  std::ifstream inFile(fileName);
  std::string outName;
  outName ="//rscratch//asl47//Bulk_Run//Regulation//Evolution_T"+std::to_string(Num_Tiles)+"_C"+std::to_string(Colour_Space+1)+"_N"+std::to_string(Num_Genomes)+"_K"+std::to_string(GENERATION_LIMIT)+"_M"+std::to_string(Fitness_Mode)+"_R"+std::to_string(Regulated)+"_I"+std::to_string(initial_condition)+".txt";
  std::ofstream outFile(outName,std::ios_base::out);
  double muL;
  while (inFile >>muL) {
    std::cout<<"On muL "<<muL<<std::endl;
    std::vector<int> discoveries(Num_Runs), adaptations(Num_Runs);
    double mu= muL/(Num_Tiles*4.);
    EvolutionSimulation(mu,discoveries,adaptations);
    outFile <<"muL: "<<muL<<"\nD: ";
    for(int dis: discoveries) 
      outFile << dis<<" ";
    outFile << "\nA: ";
    for(int ad: adaptations)
      outFile << ad<<" ";
    outFile<<"\n";
  }
  inFile.close();
  outFile.close();  
}



////////////////////////////
//        FLAGS           //
// -T = Tile Number       //
// -N = Genome Number     //
// -F = Fitness Target    //
// -C = Colour Space      //
// -R = Regulated         //
// -D = Number of Runs    //
// -K = Generation Limit  //
////////////////////////////

int main(int argc, char* argv[]) {
  if(argc>3) {
    Set_Runtime_Configurations(argc,argv);
  }

  
  if(argc>1) {
    switch(argv[1][1]) {
    case 'D':
      FitnessEvolutionDynamic();
      break;
    case 'R':
      ManyEvolutionSimulations();
      break;
      /*
    case '3':
      FitnessEvolutionStatic();
      break;
    case '4':
      FitnessEvolutionDynamicDoublet();
      break;
    case '5':
      FitnessEvolutionDynamicSinglet();
      break;   
      */   
    case 'H':
      std::cout<<"\n**Evolution running options**\n -Z for oscillating\n -X for summed\n -Q for sequential\n -A for 5\n -B for 6\n"<<std::endl;
      break;
    default:
      std::cout<<"Unknown Parameter Flag: "<<argv[1][1]<<std::endl;
      break;
    }
  }               
}
/*
bool Evolve_Genomes(int& Discovery_Generation, int& Adaptation_Generation,std::bernoulli_distribution& Mutation_Chance,std::uniform_int_distribution<int>& Mutated_Colour) {
  std::vector<int> Initial_Genome, Index_Selections(Num_Genomes), Evolving_Genome;
  std::vector<std::vector<int> > Genome_Pool(Num_Genomes,target_genotype),Temporary_Pool(Num_Genomes);
  std::vector<double> Phenotype_Fitness_Sizes(Num_Genomes);

  //ZERO SEED CONDITIONS
  if(Regulated) {
    Initial_Genome.reserve(Num_Tiles*5);
    for(int t=0;t<Num_Tiles;++t) {
      Initial_Genome.insert(Initial_Genome.end(),{1,0,0,0,0});
    }
  }  
  else {
    Initial_Genome.reserve(Num_Tiles*4);
    for(int t=0;t<Num_Tiles;++t) {
      Initial_Genome.insert(Initial_Genome.end(),{0,0,0,0});
    }
  }
  

  //START GENERATIONS LOOPS
  for(int g=1;g<=GENERATION_LIMIT;++g) {
    int Num_Maximally_Fit=0;
    for(int n=0;n<Num_Genomes;++n) {
      Evolving_Genome=Genome_Pool[n];
      
      /////////////
      //MUTATIONS//
      /////////////
      if(Regulated) {
        for(int t=0;t<Num_Tiles*5;++t) {
          if(Mutation_Chance(RNG_Engine)) {
            if(t%5==0) //Transcription Factor
              Evolving_Genome[t]^=1;
            else {
              int previousFace=Evolving_Genome[t];
              do {
                Evolving_Genome[t]=Mutated_Colour(RNG_Engine);
              } while(Evolving_Genome[t]==previousFace);
            }
          }
        }
      }
      else {
        for(int t=0;t<Num_Tiles*4;++t) {
          if(Mutation_Chance(RNG_Engine)) {
            int previousFace=Evolving_Genome[t];
            do {
              Evolving_Genome[t]=Mutated_Colour(RNG_Engine);
            } while(Evolving_Genome[t]==previousFace);
          }
        } 
      }
      /////////////////
      //END MUTATIONS//
      /////////////////
      double Phenotype_Size=0;
      if(Regulated) {
        std::vector<int> regulated_Genome;
        for(std::vector<int>::iterator it=Evolving_Genome.begin();it!=Evolving_Genome.end();it+=5) {
          if(*it==1) {
            regulated_Genome.insert(regulated_Genome.end(),it+1,it+5);
          } 
        }
        if(regulated_Genome.empty()) {
          Phenotype_Size=0;
        }
        else {
          Phenotype_Size=Get_Phenotype_Fitness(regulated_Genome,Shape_Matching_Fitness);          
          //G_L[g]+=regulated_Genome.size()/4;
        }
      }
      else {
        Phenotype_Size=Get_Phenotype_Fitness(Evolving_Genome,Shape_Matching_Fitness);
      }
      if(Phenotype_Size>=Target_Fitness-MINIMUM_FITNESS_THRESHOLD) {
        ++Num_Maximally_Fit;
        if(Discovery_Generation==-1) //DISCOVERY
          Discovery_Generation=g;
      }
      Phenotype_Fitness_Sizes[n]=Fitness_Function(Phenotype_Size);
      Genome_Pool[n]=Evolving_Genome;
    }
    //END GENOME LOOP

    
    if(Num_Maximally_Fit>=std::floor(Num_Genomes/2.)) { //ADAPTATION
      Adaptation_Generation=g;
      return true;
    }
    if(*std::max_element(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end())<MINIMUM_FITNESS_THRESHOLD) { //EXTINCTION
      Adaptation_Generation=-1;//g
      Discovery_Generation=-1;    
      return false;
    }

    //Bump_Fitness(Phenotype_Fitness_Sizes,false);
    Index_Selections=Roulette_Wheel_Selection(Phenotype_Fitness_Sizes,Num_Genomes);
    for(int nth_select=0;nth_select<Num_Genomes;++nth_select) {
      Temporary_Pool[nth_select]=Genome_Pool[Index_Selections[nth_select]];
    }   
    Genome_Pool.assign(Temporary_Pool.begin(),Temporary_Pool.end());
  }
  //END GENERATION LOOP
  if(Adaptation_Generation==-1) {
    Adaptation_Generation=GENERATION_LIMIT;
  }
  if(Discovery_Generation==-1) {
    Discovery_Generation=GENERATION_LIMIT;
  }
  return true;
}

void Evolution_Simulation(const double Mu, int& discovery_Median,int& adaptation_Median, int& Extinction_Events,int& Extinction_Average) {
  std::vector<int> Discoverys,Adaptations, Extinctions;
  std::bernoulli_distribution Mutation_Chance(Mu);
  std::uniform_int_distribution<int> Mutated_Colour(0,Colour_Space);
  
#pragma omp parallel for schedule(dynamic) 
  for(int run =0;run<Num_Runs;++run) {
    int Discovery_Generation=-1, Adaptation_Generation=-1;
    if(!Evolve_Genomes(Discovery_Generation,Adaptation_Generation,Mutation_Chance,Mutated_Colour)) {
#pragma omp critical(Extinctions)
        {
          Extinctions.emplace_back(Adaptation_Generation);
        }
    }
    else {
#pragma omp critical(Survivals)
        {
          Discoverys.emplace_back(Discovery_Generation);
          Adaptations.emplace_back(Adaptation_Generation);
        }
    }
  }
  Extinction_Events=Extinctions.size();
  if(Extinction_Events>0) {
    Extinction_Average=std::accumulate(Extinctions.begin(),Extinctions.end(),0.0)/Extinction_Events;
  }
  if(Extinction_Events<Num_Runs) {
    discovery_Median=Find_Percentile_of_Vector(Discoverys,50);
    adaptation_Median=Find_Percentile_of_Vector(Adaptations,50);
  }

}

void Run_Evolution_Simulation_Over_Range(std::string inFileName) {
  std::string fileName="//rscratch//asl47//Bulk_Run//"+inFileName+".txt";
  std::ifstream inFile(fileName);
  std::string outName;
  bool unique_File_Check=true;  
    
  std::string Fit_Type=Shape_Matching_Fitness? "Shape":"Size";
  do {
    outName ="//rscratch//asl47//Bulk_Run//Evolution_T"+std::to_string(Num_Tiles)+"_C"+std::to_string(Colour_Space+1)+"_N"+std::to_string(Num_Genomes)+"_R"+std::to_string(Num_Runs)+"_M"+std::to_string(Fitness_Mode)+"_Targ"+std::to_string(Target_Fitness)+"_Thresh"+std::to_string(GENERATION_LIMIT)+"_Reg_"+std::to_string(Regulated)+"_Run"+std::to_string(RUN)+"_"+Fit_Type+".txt";
    if (std::ifstream(outName)) {
      std::cout << "File already exists for this, increasing run number" << std::endl;
      ++RUN;
    }
    else {
      unique_File_Check=false;
    }
  }
  while(unique_File_Check);
  

  std::ofstream outFile(outName,std::ios_base::app | std::ios_base::out);
  double muL;
  while (inFile >>muL) {
    int discovery_Median=-1,adaptation_Median=-1,Extinction_Events=0,Extinction_Average=-1;
    double mu= Regulated? muL/(Num_Tiles*5.) : muL/(Num_Tiles*4.);
    //std::cout<<" mu is "<<mu<<" (muL is "<<mu*Num_Tiles*4<<")"<<std::endl;
    Evolution_Simulation(mu,discovery_Median,adaptation_Median,Extinction_Events,Extinction_Average);
    outFile <<"muL: "<<muL<<" D: "<<discovery_Median<<" A: "<<adaptation_Median<<" E: "<<Extinction_Events<<" Avg: "<<Extinction_Average<<"\n";
    if(adaptation_Median>=GENERATION_LIMIT)
      break;

  }
  inFile.close();
  outFile.close();

  
}
*/
//DEAD CODE

/*
void Evolve_Fitness(std::string Run_Details,double Mu) {
  int Burn_In_Period=0;//100000;
  std::vector<double> Average_Fitnesses;
  Average_Fitnesses.reserve(GENERATION_LIMIT);
  std::vector<double> Max_Fitnesses;
  Max_Fitnesses.reserve(GENERATION_LIMIT);
  std::vector<double> Min_Fitnesses;
  Min_Fitnesses.reserve(GENERATION_LIMIT);
  //std::vector<double> Avg_Size;
  //Avg_Size.reserve(GENERATION_LIMIT);

  std::uniform_int_distribution<int> Mutated_Colour(0,Colour_Space);
  std::bernoulli_distribution Mutation_Chance(Mu);
  std::vector<int> Initial_Genome;//{4,5,0,2, 0,1,0,3, 6,0,7,0, 14,14,14,14, 13,0,11,0, 8,9,0,9, 0,0,0,10};

  //Square shape
  
  std::vector<int> Target_1{0,0,0,0,0,0,0,0,0, 0,0,0,0,1,1,1,0,0, 0,1,0,0,0,1,0,0,0, 0,1,1,1,1,1,0,0,0, 0,1,0,1,0,1,0,1,0, 0,0,0,1,1,1,1,1,0, 0,0,0,1,0,0,0,1,0, 0,0,1,1,1,0,0,0,0, 0,0,0,0,0,0,0,0,0};
  int D_X_1=9,D_Y_1=9;
  //Cross shape
  std::vector<int> Target_2{0,0,0,0,0,0,0,0,0, 0,0,0,1,1,1,0,0,0, 0,0,0,0,1,0,0,0,0, 0,1,0,0,1,0,0,1,0, 0,1,1,1,1,1,1,1,0, 0,1,0,0,1,0,0,1,0, 0,0,0,0,1,0,0,0,0, 0,0,0,1,1,1,0,0,0, 0,0,0,0,0,0,0,0,0};
  int D_X_2=9,D_Y_2=9;
  bool Use_Target_1=false;

  //Shape_Match::Set_Shape_Data(Target_2,D_X_2,D_Y_2);
  
  
  Initial_Genome.reserve(Num_Tiles*4);
  for(int t=0;t<Num_Tiles;++t) {
    Initial_Genome.insert(Initial_Genome.end(),{0,0,0,0});
  }
  //  for(int t=0;t<Num_Tiles*4;++t) {
//	Initial_Genome.emplace_back(Mutated_Colour(RNG_Engine)
  //  }

  
  //}
  
  std::vector<std::vector<int> > Genome_Pool(Num_Genomes,Initial_Genome);
  std::vector<double> Phenotype_Fitness_Sizes(Num_Genomes);
  std::vector<std::vector<int> > Temporary_Pool(Num_Genomes);

  //std::string outName2="//rscratch//asl47//Bulk_Run//Fitness//"+Run_Details+"_Genomes.txt";
  //std::ofstream outFile2(outName2, std::ios_base::out);
  //std::set<std::vector<int> > genome_set;
  

  for(int g=0;g<GENERATION_LIMIT+Burn_In_Period;++g) {
    int Num_Maximally_Fit=0, Num_Minimally_Fit=0;
    if(g%Fitness_Oscillation_Rate==0) {
      Use_Target_1= !Use_Target_1;
      //if(Use_Target_1)
      //  Shape_Match::Set_Shape_Data(Target_1,D_X_1,D_Y_1);
      //else
      //  Shape_Match::Set_Shape_Data(Target_2,D_X_2,D_Y_2);
    }
    //std::vector<int> Sizess(Num_Genomes);
#pragma omp parallel for schedule(dynamic) reduction(+:Num_Maximally_Fit,Num_Minimally_Fit)
    //int temp_max=0;
    for(int n=0;n<Num_Genomes;++n) {
      
      std::vector<int> Evolving_Genome=Genome_Pool[n];
      for(int t=0;t<Num_Tiles*4;++t) {
        if(Mutation_Chance(RNG_Engine)) {
          int previousFace=Evolving_Genome[t];
          do {
            Evolving_Genome[t]=Mutated_Colour(RNG_Engine);
          }while(Evolving_Genome[t]==previousFace);
        }
      } 
      double Phenotype_Size=Fitness_Function(Get_Phenotype_Fitness(Evolving_Genome,Shape_Matching_Fitness));
      //}
      
      Phenotype_Fitness_Sizes[n]=Phenotype_Size;
      if(Phenotype_Size>1.0-MINIMUM_FITNESS_THRESHOLD) {
        ++Num_Maximally_Fit;
        //#pragma omp critical(setf)
        //{
        //  std::vector<int> cleaned_up(Evolving_Genome);
        //  Clean_Genome(cleaned_up);
        //  genome_set.insert(cleaned_up);          
        //}
        
      }
      if(Phenotype_Size<MINIMUM_FITNESS_THRESHOLD) {
        ++Num_Minimally_Fit;
      }
      Genome_Pool[n]=Evolving_Genome;
    }
    //END GENOME LOOP
    if(g>=Burn_In_Period) {
      Average_Fitnesses.push_back(std::accumulate(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end(),0.0)/Phenotype_Fitness_Sizes.size());
      Max_Fitnesses.push_back(Num_Maximally_Fit*1./Num_Genomes);
      Min_Fitnesses.push_back(Num_Minimally_Fit*1./Num_Genomes);
      //Avg_Size.push_back(std::accumulate(Sizess.begin(),Sizess.end(),0.0)/(Sizess.size()-std::count(Sizess.begin(),Sizess.end(),0)));
    }
    std::vector<int> Index_Selections=Roulette_Wheel_Selection(Phenotype_Fitness_Sizes,Num_Genomes);
#pragma omp parallel for
    for(int nth_select=0;nth_select<Num_Genomes;++nth_select) {
      Temporary_Pool[nth_select]=Genome_Pool[Index_Selections[nth_select]];
    }   
    Genome_Pool.assign(Temporary_Pool.begin(),Temporary_Pool.end());
  }
    
  
  //END GENERATION LOOP
  std::string outName="//rscratch//asl47//Bulk_Run//Fitness//"+Run_Details+".txt";
  std::ofstream outFile(outName, std::ios_base::out);
  for(int i=0;i<Average_Fitnesses.size();++i) {
    outFile << "g: "<<i+1<<" a: "<<Average_Fitnesses[i]<< " m: "<<Max_Fitnesses[i]<<" n: "<<Min_Fitnesses[i]<<"\n";//<<" s: "<<Avg_Size[i]<<"\n"; 
  }
  outFile.close();
}



void Fitness_Evolution() {
  std::cout<<"Evolving for "<<GENERATION_LIMIT<< " generations with oscillation rate of "<<Fitness_Oscillation_Rate<<std::endl;
  std::vector<double> Mu_List{1.0f/(Num_Tiles*4*2)};//,1./56,1./28}; //{0.0005,0.005,0.05,.5};
  for(int n=0;n<Mu_List.size();++n) {
    std::cout<<"Now running for Mu "<<Mu_List[n]<<std::endl;
    //#pragma omp parallel for schedule(dynamic)
    for(int r=0;r<Num_Runs;++r) { 
      std::string runN="T_"+std::to_string(Num_Tiles)+"C_"+std::to_string(Colour_Space)+"N_"+std::to_string(Num_Genomes)+"Mu_"+std::to_string(Mu_List[n])+"_Osc_"+std::to_string(Fitness_Oscillation_Rate)+"_Gen_"+std::to_string(GENERATION_LIMIT)+"_Run_"+std::to_string(r);
      Evolve_Fitness(runN,Mu_List[n]);
    }
  }
}

void Fitness_Evolution2() {
  double mu=1.0/(Num_Tiles*16);
#pragma omp parallel for schedule(dynamic)
  for(int r=0;r<Num_Runs;++r) {
    std::string runN="T"+std::to_string(Num_Tiles)+"_C"+std::to_string(Colour_Space)+"_N"+std::to_string(Num_Genomes)+"_Mu"+std::to_string(mu)+"_Osc"+std::to_string(Fitness_Oscillation_Rate)+"_B"+std::to_string(burn_in_period)+"_S"+std::to_string(Shape_Matching_Fitness)+"_Run"+std::to_string(r);
    Evolve_Fitness2(runN,mu);
  }
}

void Evolve_Fitness2(std::string Run_Details,double Mu) {
  std::string outName="//rscratch//asl47//Bulk_Run//Modular//Modular2_"+Run_Details+".txt";
  std::ofstream outFile(outName, std::ios_base::out);
  std::uniform_int_distribution<int> Mutated_Colour(0,Colour_Space);
  std::bernoulli_distribution Mutation_Chance(Mu);
  std::vector<int> Initial_Genome(Num_Tiles*4,0);
  std::vector<std::vector<int> > Genome_Pool(Num_Genomes,Initial_Genome);
  std::vector<double> Phenotype_Fitness_Sizes(Num_Genomes);
  std::vector<std::vector<int> > Temporary_Pool(Num_Genomes); 

  for(int g=0;g<GENERATION_LIMIT+burn_in_period;++g) {
    int Num_Maximally_Fit=0, Num_Minimally_Fit=0;
    if(g%Fitness_Oscillation_Rate==0) {
      ++Shape_Matching_Fitness;
      if(Shape_Matching_Fitness>3)
        Shape_Matching_Fitness=1;
    }
    for(int n=0;n<Num_Genomes;++n) {
      std::vector<int> Evolving_Genome=Genome_Pool[n];
      for(int t=0;t<Num_Tiles*4;++t) {
        if(Mutation_Chance(RNG_Engine)) {
          int previousFace=Evolving_Genome[t];
          do {
            Evolving_Genome[t]=Mutated_Colour(RNG_Engine);
          }while(Evolving_Genome[t]==previousFace);
        }
      }
      double Phenotype_Size=Fitness_Function(Get_Phenotype_Fitness(Evolving_Genome,Shape_Matching_Fitness));
      Phenotype_Fitness_Sizes[n]=Phenotype_Size;
      Genome_Pool[n]=Evolving_Genome;
    } //END GENOME LOOP
    if(g>=burn_in_period) {
      outFile << "g: "<<g-burn_in_period<<"\n";
      for(int nth_genome=0;nth_genome<Num_Genomes;++nth_genome) {
        if(Phenotype_Fitness_Sizes[nth_genome]>1.0-MINIMUM_FITNESS_THRESHOLD) {
          for(int base : Genome_Pool[nth_genome]) {
            outFile << base <<" ";
          }
          std::vector<int> dup_g(Genome_Pool[nth_genome]);
          Clean_Genome(dup_g);
          outFile << " x ";
          for(int base : dup_g) {
            outFile << base <<" ";
          }
          outFile << "\n";
        }
      }
    }
    std::vector<int> Index_Selections=Roulette_Wheel_Selection(Phenotype_Fitness_Sizes,Num_Genomes);
    for(int nth_select=0;nth_select<Num_Genomes;++nth_select) {
      Temporary_Pool[nth_select]=Genome_Pool[Index_Selections[nth_select]];
    }   
    Genome_Pool.assign(Temporary_Pool.begin(),Temporary_Pool.end());
  }

  outFile.close();
}


void Fitness_Evolution5() {
  double mu=1.0/(Num_Tiles*4*mutation_cofactor);
#pragma omp parallel for schedule(dynamic)
  for(int r=0;r<Num_Runs;++r) {
    std::string runN="T"+std::to_string(Num_Tiles)+"_C"+std::to_string(Colour_Space+1)+"_N"+std::to_string(Num_Genomes)+"_Mu"+std::to_string(mu)+"_O"+std::to_string(Fitness_Oscillation_Rate)+"_B"+std::to_string(burn_in_period)+"_Run"+std::to_string(r+RUN);
    Evolve_Fitness5(runN,mu);
  }
}

void Evolve_Fitness5(std::string Run_Details,double Mu) {
  std::string out_name_g="//rscratch//asl47//Bulk_Run//Modular//Modular5_"+Run_Details+"_Genotype.txt";
  std::string out_name_f="//rscratch//asl47//Bulk_Run//Modular//Modular5_"+Run_Details+"_Fitness.txt";
  std::ofstream out_file_g(out_name_g, std::ios_base::out);
  std::ofstream out_file_f(out_name_f, std::ios_base::out);

  const double COMPLETED_FITNESS=5;
  const int NUM_TARGETS=1;
  
  std::uniform_int_distribution<int> Mutated_Colour(0,Colour_Space);
  std::bernoulli_distribution Mutation_Chance(Mu);
  std::vector<int> Initial_Genome(Num_Tiles*4,0),target_types(NUM_TARGETS+1);
  std::vector<std::vector<int> > Genome_Pool(Num_Genomes,Initial_Genome), Temporary_Pool(Num_Genomes);
  std::vector<double> Phenotype_Fitness_Sizes(Num_Genomes), Phenotype_Fitness_Sizes_Full(Num_Genomes),target_fitnesses(NUM_TARGETS+1);
  bool discovered=false;

  for(int g=0;g<GENERATION_LIMIT+burn_in_period;++g) {
    int Num_Maximally_Fit=0, Num_Minimally_Fit=0;
    switch((g/Fitness_Oscillation_Rate)%3) {
    case 0:
      target_types={4,5,6};
      break;
    case 1:
      target_types={5,6,4};
      break;
    case 2:
      target_types={6,4,5};
      break;
    }
    for(int n=0;n<Num_Genomes;++n) {
      std::vector<int> Evolving_Genome=Genome_Pool[n];
      for(int t=0;t<Num_Tiles*4;++t) {
        if(Mutation_Chance(RNG_Engine)) {
          int previousFace=Evolving_Genome[t];
          do {
            Evolving_Genome[t]=Mutated_Colour(RNG_Engine);
          }while(Evolving_Genome[t]==previousFace);
        }
      }      
      
      if(GetMultiplePhenotypeFitness(Evolving_Genome,target_types,target_fitnesses)) {    
        Phenotype_Fitness_Sizes[n]=std::accumulate(target_fitnesses.begin(),target_fitnesses.begin()+NUM_TARGETS,0.0,[](double cum,double nex){return cum+Fitness_Function(nex);})/(NUM_TARGETS*2);
        Phenotype_Fitness_Sizes_Full[n]=std::accumulate(target_fitnesses.begin(),target_fitnesses.end(),0.0,[](double cum,double nex){return cum+Fitness_Function(nex);})/(NUM_TARGETS+2);
        for(int nth=0;nth<NUM_TARGETS;++nth) {
         if(target_fitnesses[nth]>1.0-MINIMUM_FITNESS_THRESHOLD)
           Phenotype_Fitness_Sizes[n]*=2;
        }
      }
      else {
        Phenotype_Fitness_Sizes[n]=0;
        Phenotype_Fitness_Sizes_Full[n]=0;
      }
      //if(!discovered && Phenotype_Fitness_Sizes[n]>1.0-MINIMUM_FITNESS_THRESHOLD) {
      //  std::cout<<"Discovered at generation "<<g<<" run("<<Run_Details.back()<<")"<<std::endl;
      //  discovered=true;
      //}      
      Genome_Pool[n]=Evolving_Genome;
    } //END GENOME LOOP
    
    out_file_f << std::accumulate(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end(),0.0)/Num_Genomes<<" "<<*std::max_element(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end())<<" "<<std::accumulate(Phenotype_Fitness_Sizes_Full.begin(),Phenotype_Fitness_Sizes_Full.end(),0.0)/Num_Genomes<<" "<<*std::max_element(Phenotype_Fitness_Sizes_Full.begin(),Phenotype_Fitness_Sizes_Full.end())   <<"\n";
    if(g>=burn_in_period) {
      //std::cout<<"Avg fitness: "<<std::accumulate(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end(),0.0)/Num_Genomes<<" Max: "<<*std::max_element(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end())<<std::endl;
      //for(int b: Genome_Pool[std::max_element(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end())-Phenotype_Fitness_Sizes.begin()])
      //  std::cout<<b<<" ";
      //std::cout<<std::endl;
      out_file_g << "g: "<<g-burn_in_period<<"\n";
      for(int nth_genome=0;nth_genome<Num_Genomes;++nth_genome) {
        if(Phenotype_Fitness_Sizes[nth_genome]>1.0-MINIMUM_FITNESS_THRESHOLD) {
          std::vector<int> dup_g(Genome_Pool[nth_genome]);
          Clean_Genome(dup_g);
          for(int base : dup_g) {
            out_file_g << base <<" ";
          }
          out_file_g << "\n";
        }
      }
    }
    std::vector<int> Index_Selections=Roulette_Wheel_Selection(Phenotype_Fitness_Sizes,Num_Genomes);
    for(int nth_select=0;nth_select<Num_Genomes;++nth_select) {
      Temporary_Pool[nth_select]=Genome_Pool[Index_Selections[nth_select]];
    }   
    Genome_Pool.assign(Temporary_Pool.begin(),Temporary_Pool.end());
  }
  out_file_g.close();
  out_file_f.close();
}


void Fitness_Evolution1() {
  double mu=1.0/(Num_Tiles*4*mutation_cofactor);
#pragma omp parallel for schedule(dynamic)
  for(int r=0;r<Num_Runs;++r) {
    std::string runN="T"+std::to_string(Num_Tiles)+"_C"+std::to_string(Colour_Space+1)+"_N"+std::to_string(Num_Genomes)+"_Mu"+std::to_string(mu)+"_O"+std::to_string(Fitness_Oscillation_Rate)+"_B"+std::to_string(burn_in_period)+"_Run"+std::to_string(r+RUN);
    Evolve_Fitness1(runN,mu);
  }
}

void Evolve_Fitness1(std::string Run_Details,double Mu) {
  std::string out_name_g="//rscratch//asl47//Bulk_Run//Modular//Modular5_"+Run_Details+"_Genotype.txt";
  std::string out_name_f="//rscratch//asl47//Bulk_Run//Modular//Modular5_"+Run_Details+"_Fitness.txt";
  std::ofstream out_file_g(out_name_g, std::ios_base::out);
  std::ofstream out_file_f(out_name_f, std::ios_base::out);

  const double COMPLETED_FITNESS=5;
  const int NUM_TARGETS=1;
  
  std::uniform_int_distribution<int> Mutated_Colour(0,Colour_Space);
  std::bernoulli_distribution Mutation_Chance(Mu);
  std::vector<int> Initial_Genome(Num_Tiles*4,0),target_types(NUM_TARGETS+1);
  std::vector<std::vector<int> > Genome_Pool(Num_Genomes,Initial_Genome), Temporary_Pool(Num_Genomes);
  std::vector<double> Phenotype_Fitness_Sizes(Num_Genomes), Phenotype_Fitness_Sizes_Full(Num_Genomes),target_fitnesses(NUM_TARGETS+1);
  bool discovered=false;

  for(int g=0;g<GENERATION_LIMIT+burn_in_period;++g) {
    int Num_Maximally_Fit=0, Num_Minimally_Fit=0;
    switch((g/Fitness_Oscillation_Rate)%3) {
    case 0:
      target_types={4,5,6};
      break;
    case 1:
      target_types={5,6,4};
      break;
    case 2:
      target_types={6,4,5};
      break;
    }
    for(int n=0;n<Num_Genomes;++n) {
      std::vector<int> Evolving_Genome=Genome_Pool[n];
      for(int t=0;t<Num_Tiles*4;++t) {
        if(Mutation_Chance(RNG_Engine)) {
          int previousFace=Evolving_Genome[t];
          do {
            Evolving_Genome[t]=Mutated_Colour(RNG_Engine);
          }while(Evolving_Genome[t]==previousFace);
        }
      }      
      
      if(GetMultiplePhenotypeFitness(Evolving_Genome,target_types,target_fitnesses)) {    
        Phenotype_Fitness_Sizes[n]=std::accumulate(target_fitnesses.begin(),target_fitnesses.begin()+NUM_TARGETS,0.0,[](double cum,double nex){return cum+Fitness_Function(nex);})/(NUM_TARGETS*2);
        Phenotype_Fitness_Sizes_Full[n]=std::accumulate(target_fitnesses.begin(),target_fitnesses.end(),0.0,[](double cum,double nex){return cum+Fitness_Function(nex);})/(NUM_TARGETS+2);
        for(int nth=0;nth<NUM_TARGETS;++nth) {
         if(target_fitnesses[nth]>1.0-MINIMUM_FITNESS_THRESHOLD)
           Phenotype_Fitness_Sizes[n]*=2;
        }
      }
      else {
        Phenotype_Fitness_Sizes[n]=0;
        Phenotype_Fitness_Sizes_Full[n]=0;
      }

      if(Phenotype_Fitness_Sizes[n]>0) {
        //clean and get first CC
        //shape table check
        //record which shape that was 

      }

      Genome_Pool[n]=Evolving_Genome;
    } //END GENOME LOOP
    
    out_file_f << std::accumulate(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end(),0.0)/Num_Genomes<<" "<<*std::max_element(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end())<<" "<<std::accumulate(Phenotype_Fitness_Sizes_Full.begin(),Phenotype_Fitness_Sizes_Full.end(),0.0)/Num_Genomes<<" "<<*std::max_element(Phenotype_Fitness_Sizes_Full.begin(),Phenotype_Fitness_Sizes_Full.end())   <<"\n";
    if(g>=burn_in_period) {

      out_file_g << "g: "<<g-burn_in_period<<"\n";
      for(int nth_genome=0;nth_genome<Num_Genomes;++nth_genome) {
        if(Phenotype_Fitness_Sizes[nth_genome]>1.0-MINIMUM_FITNESS_THRESHOLD) {
          std::vector<int> dup_g(Genome_Pool[nth_genome]);
          Clean_Genome(dup_g);
          for(int base : dup_g) {
            out_file_g << base <<" ";
          }
          out_file_g << "\n";
        }
      }
    }
    std::vector<int> Index_Selections=Roulette_Wheel_Selection(Phenotype_Fitness_Sizes,Num_Genomes);
    for(int nth_select=0;nth_select<Num_Genomes;++nth_select) {
      Temporary_Pool[nth_select]=Genome_Pool[Index_Selections[nth_select]];
    }   
    Genome_Pool.assign(Temporary_Pool.begin(),Temporary_Pool.end());
  }
  out_file_g.close();
  out_file_f.close();
}


void Fitness_Evolution51() {
  double mu=1.0/(Num_Tiles*40);
#pragma omp parallel for schedule(dynamic)
  for(int r=0;r<Num_Runs;++r) {
    std::string runN="T"+std::to_string(Num_Tiles)+"_C"+std::to_string(Colour_Space)+"_N"+std::to_string(Num_Genomes)+"_Mu"+std::to_string(mu)+"_Osc"+std::to_string(Fitness_Oscillation_Rate)+"_B"+std::to_string(burn_in_period)+"_S"+std::to_string(Shape_Matching_Fitness)+"_Run"+std::to_string(r);
    Evolve_Fitness51(runN,mu);
  }
}

void Evolve_Fitness51(std::string Run_Details,double Mu) {
  std::string outName="//rscratch//asl47//Bulk_Run//Modular//Modular5_"+Run_Details+".txt";
  std::ofstream outFile(outName, std::ios_base::out);
  std::uniform_int_distribution<int> Mutated_Colour(0,Colour_Space);
  std::bernoulli_distribution Mutation_Chance(Mu);
  std::vector<int> Initial_Genome(Num_Tiles*4,0);
  std::vector<std::vector<int> > Genome_Pool(Num_Genomes,Initial_Genome);
  std::vector<double> Phenotype_Fitness_Sizes(Num_Genomes);
  std::vector<double> Phenotype_Fitness_Sizes1(Num_Genomes);
  std::vector<double> Phenotype_Fitness_Sizes2(Num_Genomes);
  std::vector<double> Phenotype_Fitness_Sizes3(Num_Genomes);
  std::vector<std::vector<int> > Temporary_Pool(Num_Genomes); 

  for(int g=0;g<GENERATION_LIMIT+burn_in_period;++g) {
    int Num_Maximally_Fit=0, Num_Minimally_Fit=0;
    for(int n=0;n<Num_Genomes;++n) {
      std::vector<int> Evolving_Genome=Genome_Pool[n];
      for(int t=0;t<Num_Tiles*4;++t) {
        if(Mutation_Chance(RNG_Engine)) {
          int previousFace=Evolving_Genome[t];
          do {
            Evolving_Genome[t]=Mutated_Colour(RNG_Engine);
          }while(Evolving_Genome[t]==previousFace);
        }
      }
      Phenotype_Fitness_Sizes1[n]=Fitness_Function(Get_Phenotype_Fitness(Evolving_Genome,1,true));
      Phenotype_Fitness_Sizes2[n]=Fitness_Function(Get_Phenotype_Fitness(Evolving_Genome,4,true));
      Phenotype_Fitness_Sizes[n]=(Phenotype_Fitness_Sizes1[n]+Phenotype_Fitness_Sizes2[n])/2.;
      
      Genome_Pool[n]=Evolving_Genome;
    } //END GENOME LOOP
    if(g>=burn_in_period) {
      std::cout<<"Avg fitness: "<<std::accumulate(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end(),0.0)/Num_Genomes<<" Max: "<<*std::max_element(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end())<<std::endl;
      for(int b: Genome_Pool[std::max_element(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end())-Phenotype_Fitness_Sizes.begin()])
        std::cout<<b<<" ";
      std::cout<<std::endl;
      outFile << "g: "<<g-burn_in_period<<"\n";
      for(int nth_genome=0;nth_genome<Num_Genomes;++nth_genome) {
        if(Phenotype_Fitness_Sizes[nth_genome]>1.0-MINIMUM_FITNESS_THRESHOLD) {
          for(int base : Genome_Pool[nth_genome]) {
            outFile << base <<" ";
          }
          std::vector<int> dup_g(Genome_Pool[nth_genome]);
          Clean_Genome(dup_g);
          outFile << " x ";
          for(int base : dup_g) {
            outFile << base <<" ";
          }
          outFile << "\n";
        }
      }
    }
    std::vector<int> Index_Selections=Roulette_Wheel_Selection(Phenotype_Fitness_Sizes,Num_Genomes);
    for(int nth_select=0;nth_select<Num_Genomes;++nth_select) {
      Temporary_Pool[nth_select]=Genome_Pool[Index_Selections[nth_select]];
    }   
    Genome_Pool.assign(Temporary_Pool.begin(),Temporary_Pool.end());
  }
  outFile.close();
}

void Fitness_Evolution6() {
  double mu=1.0/(Num_Tiles*40);
#pragma omp parallel for schedule(dynamic)
  for(int r=0;r<Num_Runs;++r) {
    std::string runN="T"+std::to_string(Num_Tiles)+"_C"+std::to_string(Colour_Space)+"_N"+std::to_string(Num_Genomes)+"_Mu"+std::to_string(mu)+"_Osc"+std::to_string(Fitness_Oscillation_Rate)+"_B"+std::to_string(burn_in_period)+"_S"+std::to_string(Shape_Matching_Fitness)+"_Run"+std::to_string(r);
    Evolve_Fitness6(runN,mu);
  }
}

void Evolve_Fitness6(std::string Run_Details,double Mu) {
  std::string outName="//rscratch//asl47//Bulk_Run//Modular//Modular6_"+Run_Details+".txt";
  std::ofstream outFile(outName, std::ios_base::out);
  std::uniform_int_distribution<int> Mutated_Colour(0,Colour_Space);
  std::bernoulli_distribution Mutation_Chance(Mu);
  std::vector<int> Initial_Genome(Num_Tiles*4,0);
  std::vector<std::vector<int> > Genome_Pool(Num_Genomes,Initial_Genome);
  std::vector<double> Phenotype_Fitness_Sizes(Num_Genomes);
  std::vector<double> Phenotype_Fitness_Sizes1(Num_Genomes);
  std::vector<double> Phenotype_Fitness_Sizes2(Num_Genomes);
  std::vector<std::vector<int> > Temporary_Pool(Num_Genomes); 

  for(int g=0;g<GENERATION_LIMIT+burn_in_period;++g) {
    int Num_Maximally_Fit=0, Num_Minimally_Fit=0;
    for(int n=0;n<Num_Genomes;++n) {
      std::vector<int> Evolving_Genome=Genome_Pool[n];
      for(int t=0;t<Num_Tiles*4;++t) {
        if(Mutation_Chance(RNG_Engine)) {
          int previousFace=Evolving_Genome[t];
          do {
            Evolving_Genome[t]=Mutated_Colour(RNG_Engine);
          }while(Evolving_Genome[t]==previousFace);
        }
      }
      Phenotype_Fitness_Sizes1[n]=Fitness_Function(Get_Phenotype_Fitness(Evolving_Genome,0,true));
      Phenotype_Fitness_Sizes2[n]=Fitness_Function(Get_Phenotype_Fitness(Evolving_Genome,4,true));
      Phenotype_Fitness_Sizes[n]=(Phenotype_Fitness_Sizes1[n]+Phenotype_Fitness_Sizes2[n])/2.;
      
      Genome_Pool[n]=Evolving_Genome;
    } //END GENOME LOOP
    if(g>=burn_in_period) {
      std::cout<<"Avg fitness: "<<std::accumulate(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end(),0.0)/Num_Genomes<<" Max: "<<*std::max_element(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end())<<std::endl;
      for(int b: Genome_Pool[std::max_element(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end())-Phenotype_Fitness_Sizes.begin()])
        std::cout<<b<<" ";
      std::cout<<std::endl;
      outFile << "g: "<<g-burn_in_period<<"\n";
      for(int nth_genome=0;nth_genome<Num_Genomes;++nth_genome) {
        if(Phenotype_Fitness_Sizes[nth_genome]>1.0-MINIMUM_FITNESS_THRESHOLD) {
          for(int base : Genome_Pool[nth_genome]) {
            outFile << base <<" ";
          }
          std::vector<int> dup_g(Genome_Pool[nth_genome]);
          Clean_Genome(dup_g);
          outFile << " x ";
          for(int base : dup_g) {
            outFile << base <<" ";
          }
          outFile << "\n";
        }
      }
    }
    std::vector<int> Index_Selections=Roulette_Wheel_Selection(Phenotype_Fitness_Sizes1,Num_Genomes);
    for(int nth_select=0;nth_select<Num_Genomes;++nth_select) {
      Temporary_Pool[nth_select]=Genome_Pool[Index_Selections[nth_select]];
      Phenotype_Fitness_Sizes[nth_select]=Phenotype_Fitness_Sizes2[Index_Selections[nth_select]];
    }   
    Genome_Pool.assign(Temporary_Pool.begin(),Temporary_Pool.end());
    
    Index_Selections=Roulette_Wheel_Selection(Phenotype_Fitness_Sizes,Num_Genomes);
    for(int nth_select=0;nth_select<Num_Genomes;++nth_select) {
      Temporary_Pool[nth_select]=Genome_Pool[Index_Selections[nth_select]];
    }   
    Genome_Pool.assign(Temporary_Pool.begin(),Temporary_Pool.end());


  }
  outFile.close();
}

void Evolve_Fitness7(double Mu) {
  std::uniform_int_distribution<int> Mutated_Colour(0,Colour_Space);
  std::bernoulli_distribution Mutation_Chance(Mu);
  std::vector<int> Initial_Genome(Num_Tiles*4,0);
  std::vector<std::vector<int> > Genome_Pool(Num_Genomes,Initial_Genome);
  std::vector<double> Phenotype_Fitness_Sizes(Num_Genomes);
  std::vector<std::vector<int> > Temporary_Pool(Num_Genomes); 

  for(int g=0;g<GENERATION_LIMIT+burn_in_period;++g) {
    int Num_Maximally_Fit=0, Num_Minimally_Fit=0;
#pragma omp parallel for schedule(dynamic)
    for(int n=0;n<Num_Genomes;++n) {
      std::vector<int> Evolving_Genome=Genome_Pool[n];
      for(int t=0;t<Num_Tiles*4;++t) {
        if(Mutation_Chance(RNG_Engine)) {
          int previousFace=Evolving_Genome[t];
          do {
            Evolving_Genome[t]=Mutated_Colour(RNG_Engine);
          }while(Evolving_Genome[t]==previousFace);
        }
      }
      Phenotype_Fitness_Sizes[n]=Fitness_Function(Get_Phenotype_Fitness(Evolving_Genome,Shape_Matching_Fitness,true));

      
      Genome_Pool[n]=Evolving_Genome;
    } //END GENOME LOOP
    if(g>=burn_in_period) {
      std::cout<<"Avg fitness: "<<std::accumulate(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end(),0.0)/Num_Genomes<<" Max: "<<*std::max_element(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end())<<std::endl;
      for(int b: Genome_Pool[std::max_element(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end())-Phenotype_Fitness_Sizes.begin()])
        std::cout<<b<<" ";
      std::cout<<std::endl;
    }
    std::vector<int> Index_Selections=Roulette_Wheel_Selection(Phenotype_Fitness_Sizes,Num_Genomes);
    for(int nth_select=0;nth_select<Num_Genomes;++nth_select) {
      Temporary_Pool[nth_select]=Genome_Pool[Index_Selections[nth_select]];
    }   
    Genome_Pool.assign(Temporary_Pool.begin(),Temporary_Pool.end());

  }

}

void FitnessEvolutionStatic() {
  double mu=1.0/(Num_Tiles*mutation_cofactor);
#pragma omp parallel for schedule(dynamic)
  for(int r=0;r<Num_Runs;++r) {
    std::string runN="T"+std::to_string(Num_Tiles)+"_C"+std::to_string(Colour_Space+1)+"_N"+std::to_string(Num_Genomes)+"_Mu"+std::to_string(mu)+"_K"+std::to_string(GENERATION_LIMIT)+"_Run"+std::to_string(r+RUN);
    EvolveFitnessStatic(runN,mu);
  }
}

void EvolveFitnessStatic(std::string Run_Details,double Mu) {
  std::string out_name_g="//rscratch//asl47//Bulk_Run//Modular//Static_"+Run_Details+"_Genotype.txt";
  std::string out_name_f="//rscratch//asl47//Bulk_Run//Modular//Static_"+Run_Details+"_Fitness.txt";
  std::ofstream out_file_g(out_name_g, std::ios_base::out);
  std::ofstream out_file_f(out_name_f, std::ios_base::out);

  //const double COMPLETED_FITNESS=5;

  std::uniform_int_distribution<int> Mutated_Colour(0,Colour_Space);
  std::bernoulli_distribution Mutation_Chance(Mu);
  std::vector<int> Initial_Genome(Num_Tiles*4,0),target_types{4,5,6}, Index_Selections(Num_Genomes);
  std::vector<std::vector<int> > Genome_Pool(Num_Genomes,Initial_Genome), Temporary_Pool(Num_Genomes);
  std::vector<double> Phenotype_Fitness_Sizes(Num_Genomes),target_fitnesses{0,0,0};

  for(int g=1;g<=GENERATION_LIMIT;++g) {
    for(int n=0;n<Num_Genomes;++n) {
      std::vector<int> Evolving_Genome=Genome_Pool[n];
      for(int t=0;t<Num_Tiles*4;++t) {
        if(Mutation_Chance(RNG_Engine)) {
          int previousFace=Evolving_Genome[t];
          do {
            Evolving_Genome[t]=Mutated_Colour(RNG_Engine);
          }while(Evolving_Genome[t]==previousFace);
        }
      }      
      
      if(GetMultiplePhenotypeFitness(Evolving_Genome,target_types,target_fitnesses)) {    
        Phenotype_Fitness_Sizes[n]=std::accumulate(target_fitnesses.begin(),target_fitnesses.end(),0.0,[](double cum,double nex){return cum+Fitness_Function(nex);})/(target_types.size()*8);
        for(unsigned int nth=0;nth<target_types.size();++nth) {
          if(target_fitnesses[nth]>1.0-MINIMUM_FITNESS_THRESHOLD)
            Phenotype_Fitness_Sizes[n]*=2.;
        }
      }
      else 
        Phenotype_Fitness_Sizes[n]=0;
         
      Genome_Pool[n]=Evolving_Genome;
    } //END GENOME LOOP

    double max_fitness=*std::max_element(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end());
    out_file_f << max_fitness << " " << std::count(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end(),max_fitness) << "\n";
    if(g%5000==0) {
      out_file_g << "g: "<<g<<"\n";
      for(int nth_genome=0;nth_genome<Num_Genomes;++nth_genome) {
        if(Phenotype_Fitness_Sizes[nth_genome]>1.0-MINIMUM_FITNESS_THRESHOLD) {
          std::vector<int> dup_g(Genome_Pool[nth_genome]);
          Clean_Genome(dup_g);
          if(Disjointed_Check(dup_g))
             dup_g.erase(std::find(dup_g.begin(),dup_g.end(),-1),dup_g.end());
          for(int base : dup_g) {
            out_file_g << base <<" ";
          }
          out_file_g << "\n";
        }
      }
    }
    Index_Selections=Roulette_Wheel_Selection(Phenotype_Fitness_Sizes,Num_Genomes);
    for(int nth_select=0;nth_select<Num_Genomes;++nth_select) {
      Temporary_Pool[nth_select]=Genome_Pool[Index_Selections[nth_select]];
    }   
    Genome_Pool.assign(Temporary_Pool.begin(),Temporary_Pool.end());
  }
  out_file_g.close();
  out_file_f.close();
}

void FitnessEvolutionDynamicDoublet() {
  double mu=1.0/(Num_Tiles*mutation_cofactor);
#pragma omp parallel for schedule(dynamic)
  for(int r=0;r<Num_Runs;++r) {
    std::string runN="T"+std::to_string(Num_Tiles)+"_C"+std::to_string(Colour_Space+1)+"_N"+std::to_string(Num_Genomes)+"_Mu"+std::to_string(mu)+"_O"+std::to_string(Fitness_Oscillation_Rate)+"_K"+std::to_string(GENERATION_LIMIT)+"_Run"+std::to_string(r+RUN);
    EvolveFitnessDynamicDoublet(runN,mu);
  }
}

void EvolveFitnessDynamicDoublet(std::string Run_Details,double Mu) {
  std::string out_name_g="//rscratch//asl47//Bulk_Run//Modular//Modular4_"+Run_Details+"_Genotype.txt";
  std::string out_name_f="//rscratch//asl47//Bulk_Run//Modular//Modular4_"+Run_Details+"_Fitness.txt";
  std::ofstream out_file_g(out_name_g, std::ios_base::out);
  std::ofstream out_file_f(out_name_f, std::ios_base::out);

  //const double COMPLETED_FITNESS=5;
  //const int NUM_TARGETS=2;
  
  std::uniform_int_distribution<int> Mutated_Colour(0,Colour_Space);
  std::bernoulli_distribution Mutation_Chance(Mu);
  std::vector<int> Initial_Genome(Num_Tiles*4,0),target_types(3), Index_Selections(Num_Genomes);
  std::vector<std::vector<int> > Genome_Pool(Num_Genomes,Initial_Genome), Temporary_Pool(Num_Genomes);
  std::vector<double> Phenotype_Fitness_Sizes(Num_Genomes), Phenotype_Fitness_Sizes_Full(Num_Genomes),target_fitnesses(3);

  for(int g=1;g<=GENERATION_LIMIT;++g) {
    switch((g/Fitness_Oscillation_Rate)%3) {
    case 0:
      target_types={4,5,6};
      break;
    case 1:
      target_types={5,6,4};
      break;
    case 2:
      target_types={6,4,5};
      break;
    }
    for(int n=0;n<Num_Genomes;++n) {
      std::vector<int> Evolving_Genome=Genome_Pool[n];
      for(int t=0;t<Num_Tiles*4;++t) {
        if(Mutation_Chance(RNG_Engine)) {
          int previousFace=Evolving_Genome[t];
          do {
            Evolving_Genome[t]=Mutated_Colour(RNG_Engine);
          }while(Evolving_Genome[t]==previousFace);
        }
      }      
      
      if(GetMultiplePhenotypeFitness(Evolving_Genome,target_types,target_fitnesses)) {    
        Phenotype_Fitness_Sizes[n]=std::accumulate(target_fitnesses.begin(),target_fitnesses.begin()+2,0.0,[](double cum,double nex){return cum+Fitness_Function(nex);})/(2*4.);
        Phenotype_Fitness_Sizes_Full[n]=std::accumulate(target_fitnesses.begin(),target_fitnesses.end(),0.0,[](double cum,double nex){return cum+Fitness_Function(nex);})/(3.);
        for(int nth=0;nth<2;++nth) {
         if(target_fitnesses[nth]>1.0-MINIMUM_FITNESS_THRESHOLD)
           Phenotype_Fitness_Sizes[n]*=2;
        }
      }
      else {
        Phenotype_Fitness_Sizes[n]=0;
        Phenotype_Fitness_Sizes_Full[n]=0;
      }
 
      Genome_Pool[n]=Evolving_Genome;
    } //END GENOME LOOP


    double max_fitness=*std::max_element(Phenotype_Fitness_Sizes_Full.begin(),Phenotype_Fitness_Sizes_Full.end());
    out_file_f << max_fitness << " " << std::count(Phenotype_Fitness_Sizes_Full.begin(),Phenotype_Fitness_Sizes_Full.end(),max_fitness) << "\n";
    if(g%5000==0) {
      out_file_g << "g: "<<g<<"\n";
      for(int nth_genome=0;nth_genome<Num_Genomes;++nth_genome) {
        if(Phenotype_Fitness_Sizes_Full[nth_genome]>1.0-MINIMUM_FITNESS_THRESHOLD) {
          std::vector<int> dup_g(Genome_Pool[nth_genome]);
          Clean_Genome(dup_g);
          if(Disjointed_Check(dup_g))
             dup_g.erase(std::find(dup_g.begin(),dup_g.end(),-1),dup_g.end());
          for(int base : dup_g) {
            out_file_g << base <<" ";
          }
          out_file_g << "\n";
        }
      }
    }
    Index_Selections=Roulette_Wheel_Selection(Phenotype_Fitness_Sizes,Num_Genomes);
    for(int nth_select=0;nth_select<Num_Genomes;++nth_select) {
      Temporary_Pool[nth_select]=Genome_Pool[Index_Selections[nth_select]];
    }   
    Genome_Pool.assign(Temporary_Pool.begin(),Temporary_Pool.end());
  }
  out_file_g.close();
  out_file_f.close();
}

void FitnessEvolutionDynamicSinglet() {
  double mu=1.0/(Num_Tiles*mutation_cofactor);
#pragma omp parallel for schedule(dynamic)
  for(int r=0;r<Num_Runs;++r) {
    std::string runN="T"+std::to_string(Num_Tiles)+"_C"+std::to_string(Colour_Space+1)+"_N"+std::to_string(Num_Genomes)+"_Mu"+std::to_string(mu)+"_O"+std::to_string(Fitness_Oscillation_Rate)+"_K"+std::to_string(GENERATION_LIMIT)+"_Run"+std::to_string(r+RUN);
    EvolveFitnessDynamicSinglet(runN,mu);
  }
}

void EvolveFitnessDynamicSinglet(std::string Run_Details,double Mu) {
  std::string out_name_g="//rscratch//asl47//Bulk_Run//Modular//Modular5_"+Run_Details+"_Genotype.txt";
  std::string out_name_f="//rscratch//asl47//Bulk_Run//Modular//Modular5_"+Run_Details+"_Fitness.txt";
  std::ofstream out_file_g(out_name_g, std::ios_base::out);
  std::ofstream out_file_f(out_name_f, std::ios_base::out);

  //const double COMPLETED_FITNESS=5;
  //const int NUM_TARGETS=1;
  
  std::uniform_int_distribution<int> Mutated_Colour(0,Colour_Space);
  std::bernoulli_distribution Mutation_Chance(Mu);
  std::vector<int> Initial_Genome(Num_Tiles*4,0),target_types(3), Index_Selections(Num_Genomes);
  std::vector<std::vector<int> > Genome_Pool(Num_Genomes,Initial_Genome), Temporary_Pool(Num_Genomes);
  std::vector<double> Phenotype_Fitness_Sizes(Num_Genomes), Phenotype_Fitness_Sizes_Full(Num_Genomes),target_fitnesses(3);

  for(int g=1;g<=GENERATION_LIMIT;++g) {
    switch((g/Fitness_Oscillation_Rate)%3) {
    case 0:
      target_types={4,5,6};
      break;
    case 1:
      target_types={5,6,4};
      break;
    case 2:
      target_types={6,4,5};
      break;
    }
    for(int n=0;n<Num_Genomes;++n) {
      std::vector<int> Evolving_Genome=Genome_Pool[n];
      for(int t=0;t<Num_Tiles*4;++t) {
        if(Mutation_Chance(RNG_Engine)) {
          int previousFace=Evolving_Genome[t];
          do {
            Evolving_Genome[t]=Mutated_Colour(RNG_Engine);
          }while(Evolving_Genome[t]==previousFace);
        }
      }      
      
      if(GetMultiplePhenotypeFitness(Evolving_Genome,target_types,target_fitnesses)) {    
        Phenotype_Fitness_Sizes[n]=std::accumulate(target_fitnesses.begin(),target_fitnesses.begin()+1,0.0,[](double cum,double nex){return cum+Fitness_Function(nex);});
        Phenotype_Fitness_Sizes_Full[n]=std::accumulate(target_fitnesses.begin(),target_fitnesses.end(),0.0,[](double cum,double nex){return cum+Fitness_Function(nex);})/(3.);
      }
      else {
        Phenotype_Fitness_Sizes[n]=0;
        Phenotype_Fitness_Sizes_Full[n]=0;
      }
       Genome_Pool[n]=Evolving_Genome;
    } //END GENOME LOOP


    double max_fitness=*std::max_element(Phenotype_Fitness_Sizes_Full.begin(),Phenotype_Fitness_Sizes_Full.end());
    out_file_f << max_fitness << " " << std::count(Phenotype_Fitness_Sizes_Full.begin(),Phenotype_Fitness_Sizes_Full.end(),max_fitness) << "\n";
    if(g%5000==0) {
      out_file_g << "g: "<<g<<"\n";
      for(int nth_genome=0;nth_genome<Num_Genomes;++nth_genome) {
        if(Phenotype_Fitness_Sizes_Full[nth_genome]>1.0-MINIMUM_FITNESS_THRESHOLD) {
          std::vector<int> dup_g(Genome_Pool[nth_genome]);
          Clean_Genome(dup_g);
          if(Disjointed_Check(dup_g))
             dup_g.erase(std::find(dup_g.begin(),dup_g.end(),-1),dup_g.end());
          for(int base : dup_g) {
            out_file_g << base <<" ";
          }
          out_file_g << "\n";
        }
      }
    }
    Index_Selections=Roulette_Wheel_Selection(Phenotype_Fitness_Sizes,Num_Genomes);
    for(int nth_select=0;nth_select<Num_Genomes;++nth_select) {
      Temporary_Pool[nth_select]=Genome_Pool[Index_Selections[nth_select]];
    }   
    Genome_Pool.assign(Temporary_Pool.begin(),Temporary_Pool.end());
  }
  out_file_g.close();
  out_file_f.close();
}
*/
