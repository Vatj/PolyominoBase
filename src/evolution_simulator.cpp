#include <evolution_simulator.hpp>
//#include <stochastic_model.hpp>

///////////////////////////////////////
//GENOME EVOLTUTION RELATED FUNCTIONS//
///////////////////////////////////////

void FitnessEvolutionDynamic() {
  double mu=1.0/(Num_Tiles*mutation_cofactor);
#pragma omp parallel for schedule(dynamic)
  for(int r=0;r<Num_Runs;++r) {
    std::string runN="_T"+std::to_string(Num_Tiles)+"_C"+std::to_string(Colour_Space+1)+"_N"+std::to_string(Num_Genomes)+"_Mu"+std::to_string(mu)+"_O"+std::to_string(Fitness_Oscillation_Rate)+"_K"+std::to_string(GENERATION_LIMIT)+"_P"+std::to_string(periodic_changes)+"_Run"+std::to_string(r+RUN);
    EvolveFitnessDynamic(runN,mu);
  }
}

void SetEvolutionTargets(std::vector<int>& targets,int generation) {
  if(periodic_changes) {
    switch(((generation-1)/Fitness_Oscillation_Rate)%3) {
    case 0:
      targets={4,5,6};
      break;
    case 1:
      targets={5,6,4};
      break;
    case 2:
      targets={6,4,5};
      break;
    }
  }
  else {
    std::shuffle(targets.begin(),targets.end(),RNG_Engine);
  }
    

  
}
void EvolveFitnessDynamic(std::string Run_Details,double Mu) {
  std::string out_name_g="//rscratch//asl47//Bulk_Run//Modular//A"+std::to_string(active_targets)+Run_Details+"_Genotype.txt";
  std::string out_name_f="//rscratch//asl47//Bulk_Run//Modular//A"+std::to_string(active_targets)+Run_Details+"_Fitness.txt";
  //std::string out_name_r="//rscratch//asl47//Bulk_Run//Modular//A"+std::to_string(active_targets)+Run_Details+"_Robust.txt";
  //std::ofstream out_file_g(out_name_g, std::ios_base::out);
  std::ofstream out_file_f(out_name_f, std::ios_base::out);
  //std::ofstream out_file_r(out_name_r, std::ios_base::out);

  const int FITNESS_CLIFF=5;
  const int TOTAL_TARGETS=3;
  
  std::uniform_int_distribution<int> Mutated_Colour(0,Colour_Space);
  std::bernoulli_distribution Mutation_Chance(Mu);
  std::vector<int> Initial_Genome(Num_Tiles*4,0),Evolving_Genome(Num_Tiles*4),target_types, Index_Selections(Num_Genomes);
  target_types={4,5,6};
  
  //InitialGenotypeConfiguration(Initial_Genome);
  std::vector<std::vector<int> > Genome_Pool(Num_Genomes,Initial_Genome), Temporary_Pool(Num_Genomes);
  std::vector<double> Phenotype_Fitness_Sizes(Num_Genomes), Phenotype_Fitness_Sizes_Full(Num_Genomes),target_fitnesses(TOTAL_TARGETS);
  std::set< std::vector<int> > unique_genotypes;

  
  for(int g=1;g<=GENERATION_LIMIT;++g) {
    if((g-1)%Fitness_Oscillation_Rate==0)
      SetEvolutionTargets(target_types,g);
    
    int a=0,b=0,c=0,ab=0,bc=0,ac=0,abc=0;
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
        //Phenotype_Fitness_Sizes_Full[n]=std::accumulate(target_fitnesses.begin(),target_fitnesses.end(),0.0)/.;
        Phenotype_Fitness_Sizes_Full[n]=std::accumulate(target_fitnesses.begin(),target_fitnesses.end(),0.0)>=static_cast<double>(TOTAL_TARGETS) ? 1. : 0;

	if(target_fitnesses[std::find(target_types.begin(),target_types.end(),4)-target_types.begin()]==1) {
	  if(target_fitnesses[std::find(target_types.begin(),target_types.end(),5)-target_types.begin()]==1) {
	    if(target_fitnesses[std::find(target_types.begin(),target_types.end(),6)-target_types.begin()]==1) {
              ++abc;
            }
            else {
              ++ab;
            }
          }
          else {
	    if(target_fitnesses[std::find(target_types.begin(),target_types.end(),6)-target_types.begin()]==1) {
	      ++ac;
	    }
	    else {
	      ++a;
	    }
	  }
	}
	else {
	  if(target_fitnesses[std::find(target_types.begin(),target_types.end(),5)-target_types.begin()]==1) {
	    if(target_fitnesses[std::find(target_types.begin(),target_types.end(),6)-target_types.begin()]==1) {
	      ++bc;
	    }
	    else {
	      ++b;
	    }
	  }
	  else {
	    if(target_fitnesses[std::find(target_types.begin(),target_types.end(),6)-target_types.begin()]==1) {
	      ++c;
	    }
	  }
        
	}
      }
      
      else {
	Phenotype_Fitness_Sizes[n]=0;
	Phenotype_Fitness_Sizes_Full[n]=0;
      }
      Genome_Pool[n]=Evolving_Genome;
	/*
	  if(Phenotype_Fitness_Sizes_Full[n]>=1.) {
	  std::vector<int> duplicate_genome=(Evolving_Genome);
	  Clean_Genome(duplicate_genome);
	  if(Disjointed_Check(duplicate_genome))
          duplicate_genome.erase(std::find(duplicate_genome.begin(),duplicate_genome.end(),-1),duplicate_genome.end());
	  unique_genotypes.insert(duplicate_genome);
	  }
	*/
      
    } //END GENOME LOOP

    
    double max_fitness=*std::max_element(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end());
    if(true || max_fitness>=1.) {
      //out_file_f << max_fitness << " " << std::count(Phenotype_Fitness_Sizes_Full.begin(),Phenotype_Fitness_Sizes_Full.end(),max_fitness) << "\n";
      //out_file_f << std::count(Phenotype_Fitness_Sizes.begin(),Phenotype_Fitness_Sizes.end(),1) << " " << std::count(Phenotype_Fitness_Sizes_Full.begin(),Phenotype_Fitness_Sizes_Full.end(),1) << "\n";
      out_file_f <<abc<<" "<<ab<<" "<<bc<<" "<<ac<<" "<<a<<" "<<b<<" "<<c<<"\n"; 
    }
    else {
      out_file_f << "0 0 0 0 0 0 0\n";
    }
    
    

    if(abc>=Num_Genomes/2)
      return;
    Index_Selections=Roulette_Wheel_Selection(Phenotype_Fitness_Sizes);
    for(int nth_select=0;nth_select<Num_Genomes;++nth_select) {
      Temporary_Pool[nth_select]=Genome_Pool[Index_Selections[nth_select]];
    }    
    Genome_Pool.assign(Temporary_Pool.begin(),Temporary_Pool.end());
  }
  /*
  for(std::vector<int> genotype : unique_genotypes) {
    for(int base : genotype) {
      out_file_g << base <<" ";
    }
    out_file_g << "\n";
  }
  */
  
  
  //out_file_g.close();
  out_file_f.close();
}
/*
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
    Index_Selections=Roulette_Wheel_Selection(Phenotype_Fitness_Sizes);
    for(int nth_select=0;nth_select<Num_Genomes;++nth_select) {
      Temporary_Pool[nth_select]=Genome_Pool[Index_Selections[nth_select]];
    }   
    Genome_Pool.assign(Temporary_Pool.begin(),Temporary_Pool.end());
  }
}
*/
void EvolveSimple(int& Discovery_Generation, int& Adaptation_Generation,int& misclass, double Mu) {
  std::uniform_int_distribution<int> Mutated_Colour(0,Colour_Space);
  std::bernoulli_distribution Mutation_Chance(Mu);
  std::vector<int> Initial_Genome(Num_Tiles*4), Evolving_Genome(Num_Tiles*4), Index_Selections(Num_Genomes);

  //InitialGenotypeConditions(Initial_Genome);
  
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
      int bfit=0;
      if(initial_condition>0) {
        int cf_fit=Fitness_Function(Get_Phenotype_Fitness(Evolving_Genome,-1,true));//Fitness_Function(Get_Phenotype_Fitness(Evolving_Genome,50,false));
      
        switch(initial_condition) {
        case 1:
          bfit=cf_fit;
          break;
        case 2:
          bfit=Fitness_Function(Get_Phenotype_Fitness(Evolving_Genome,10,false));
          break;
        default:
          break;
        }

        if(bfit!=cf_fit) {
          misclass++;
        }
      }

      
      Phenotype_Fitness_Sizes[n]=bfit;

      if(Phenotype_Fitness_Sizes[n]>=Target_Fitness) {
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
    Index_Selections=Roulette_Wheel_Selection(Phenotype_Fitness_Sizes);
    for(int nth_select=0;nth_select<Num_Genomes;++nth_select) {
      Temporary_Pool[nth_select]=Genome_Pool[Index_Selections[nth_select]];
    }   
    Genome_Pool.assign(Temporary_Pool.begin(),Temporary_Pool.end());
  }
}


void EvolutionSimulation(const double Mu,std::vector<int>& discovery,std::vector<int>& adaptation,int& misclass) {
#pragma omp parallel for schedule(dynamic) reduction(+:misclass )
  for(int run =0;run<Num_Runs;++run) {
    int Discovery_Generation=-1, Adaptation_Generation=-1;
    if(false && Regulated)
      continue;
      //EvolveRegulated(Discovery_Generation,Adaptation_Generation,Mu);
    else
      EvolveSimple(Discovery_Generation,Adaptation_Generation,misclass,Mu);                               
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
    int misclass=0;
    EvolutionSimulation(mu,discoveries,adaptations,misclass);
    outFile <<"muL: "<<muL<<" D: ";
    for(int dis: discoveries) 
      outFile << dis<<" ";
    outFile << " A: ";
    for(int ad: adaptations)
      outFile << ad<<" ";
    outFile<<"\nMc: "<<misclass<<"\n";
    
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

  std::vector<int> g{0, 60, 19, 70, 24, 0, 52, 35, 26, 37, 0, 85, 0, 0, 43, 0, 57, 86, 38, 47, 0, 0, 89, 32, 76, 50, 0, 0, 23, 84, 76, 14, 0, 39, 14, 5, 0, 0, 0, 49, 17, 68, 81, 99, 0, 54, 61, 4, 0, 0, 42, 0, 99, 62, 51, 24, 0, 0, 44, 0, 0, 13, 0, 78, 4, 39, 13, 95, 63, 69, 0, 0, 0, 18, 66, 36, 48, 60, 82, 52};
 
  std::vector<int> g2;
  if(argc>1) {
    switch(argv[1][1]) {
    case 'D':
      FitnessEvolutionDynamic();
      break;
    case 'R':
      ManyEvolutionSimulations();
      break;
    case 'X':
      Clean_Genome(g);
      Disjointed_Check(g);
      for(auto x: g)
        std::cout<<x<<" ";
      std::cout<<std::endl;
      g2.assign(g.begin(),std::find(g.begin(),g.end(),-1));
      //std::cout<<Graph_Analysis(g2)<<std::endl;
      std::cout<<Brute_Force::Analyse_Genotype_Outcome(g2,2)<<std::endl;;
      break;
    case 'H':
      std::cout<<"\n**Evolution running options**\n -Z for oscillating\n -X for summed\n -Q for sequential\n -A for 5\n -B for 6\n"<<std::endl;
      break;
    default:
      std::cout<<"Unknown Parameter Flag: "<<argv[1][1]<<std::endl;
      break;
    }
  }               
}

