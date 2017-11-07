#include "brute_simulator.hpp"

#include <time.h>
#include <sys/time.h>

namespace Brute_Force
{
  std::random_device rd;
  xorshift RNG_Generator(rd());
  
  std::vector<int> Brute_Force_Polyomino_Builder(std::vector<int> genome, unsigned int THRESHOLD_SIZE, int initial_Tile,int initial_Rotation) {
    Clean_Genome(genome,-1); 
    std::vector<int> Placed_Tiles{0,0,initial_Tile,initial_Rotation}; //DEFINED AS (X,Y,Tile Type Number, Tile Rotation[in CW rotation])
    std::vector<int> Interacting_Faces; //DEFINED AS (X,Y,Target Face colour, Target Face Index) 
    for(int face=0; face<4; ++face) {
      if(genome[initial_Tile*4+ (face-initial_Rotation+4)%4]!=0)
        Brute_Force_Interacting_Adjacency(Interacting_Faces,genome[initial_Tile*4+ (face-initial_Rotation+4)%4],face,0,0);
    }
    
    while(!Interacting_Faces.empty()) {  
      for(int number_possibilities=(Interacting_Faces.size()/4)-1; number_possibilities >=0; --number_possibilities) {
        bool collision_Detected=false;
        for(std::vector<int>::iterator occupied_iter = Placed_Tiles.begin(); occupied_iter!=Placed_Tiles.end(); occupied_iter+=4) {
          if(Interacting_Faces[number_possibilities*4]==*occupied_iter && Interacting_Faces[number_possibilities*4+1]==*(occupied_iter+1)) {
            if(Interacting_Faces[number_possibilities*4+2]!=genome[*(occupied_iter+2)*4+(Interacting_Faces[number_possibilities*4+3]-*(occupied_iter+3)+4)%4])
              return {}; //Steric mismatch, reject empty vector as sign
            Interacting_Faces.erase(Interacting_Faces.end()-4,Interacting_Faces.end()); 
            collision_Detected=true;
            break;
          }
        }
        if(collision_Detected)
          continue;
        else {
          int conjugate_count=std::count(genome.begin(),genome.end(),Interacting_Faces[number_possibilities*4+2]);
          std::uniform_int_distribution<int> Random_Count(0,conjugate_count-1);
          int nth_conjugate=Random_Count(RNG_Generator);
          std::vector<int>::iterator current_conjugate=std::find(genome.begin(),genome.end(),Interacting_Faces[number_possibilities*4+2]);
          for(int conj_cnt=1;conj_cnt<=nth_conjugate;++conj_cnt) {
            current_conjugate=std::find(current_conjugate+1,genome.end(),Interacting_Faces[number_possibilities*4+2]);
          }
          int new_Tile=(current_conjugate-genome.begin())/4;
          int new_Face=(current_conjugate-genome.begin())%4;
          int rotation= (Interacting_Faces[number_possibilities*4+3]-new_Face+4)%4;
          Placed_Tiles.insert(Placed_Tiles.end(),{Interacting_Faces[number_possibilities*4],Interacting_Faces[number_possibilities*4+1],new_Tile,rotation});
          if(Placed_Tiles.size()/4 >THRESHOLD_SIZE) 
            return Placed_Tiles;
          int placed_X=Interacting_Faces[number_possibilities*4];
          int placed_Y=Interacting_Faces[number_possibilities*4+1];
          Interacting_Faces.erase(Interacting_Faces.begin()+number_possibilities*4,Interacting_Faces.begin()+number_possibilities*4+4);
          for(int face=1;face<4;++face) {
            int temp_Face=genome[new_Tile*4+(new_Face+face)%4];
            if(temp_Face!=0)
              Brute_Force_Interacting_Adjacency(Interacting_Faces,temp_Face,(new_Face+face+rotation)%4,placed_X,placed_Y);
          }
          break;
        }
      }
    }
    return Placed_Tiles;
  }

  void Brute_Force_Interacting_Adjacency(std::vector<int>& Interacting_Faces, int interacting_Face, int face_index, int X, int Y) {
    int X_OFFSET=0,Y_OFFSET=0;
    switch(face_index) {
    case 0:
      X_OFFSET=0;
      Y_OFFSET=1;
      break;
    case 1:
      X_OFFSET=1;
      Y_OFFSET=0;
      break;
    case 2:
      X_OFFSET=0;
      Y_OFFSET=-1;
      break;
    case 3:
      X_OFFSET=-1;
      Y_OFFSET=0;
      break;
    }
    int conjugate_Face=Interaction_Matrix(interacting_Face);//(1-interacting_Face%2)*(interacting_Face-1)+(interacting_Face%2)*(interacting_Face+1); 
    std::uniform_int_distribution<int> Random_Insertion(0,Interacting_Faces.size()/4);
    Interacting_Faces.insert(Interacting_Faces.begin()+Random_Insertion(RNG_Generator)*4,{X+X_OFFSET,Y+Y_OFFSET,conjugate_Face,(face_index+2)%4});
  }

  int Analyse_Genotype_Outcome(std::vector<int> genome, int N_Repeated_Checks) {
    //////////////////
    // RETURN CODES //
    // -5 = UNBD    //
    // -4 = BND     //
    // -2 = disjoint//
    // -1 = steric  //
    // 4 = BD       //
    //////////////////
   
    const unsigned int THRESHOLD_SIZE=(genome.size()*genome.size())/2;
    bool UND_Check=false;
    std::vector<int> Placed_Tiles_Check=Brute_Force_Polyomino_Builder(genome,THRESHOLD_SIZE,0,0);
    
    if(Placed_Tiles_Check.empty()) //STERIC
      return -1;
    if(Placed_Tiles_Check.size()/4 > THRESHOLD_SIZE)
      return -1;
    
    int Delta_X_Check=-1,Delta_Y_Check=-1;
    std::vector<int> Spatial_Occupation_Check=Generate_Spatial_Occupancy(Placed_Tiles_Check,Delta_X_Check,Delta_Y_Check,0);    
    std::vector<int> Spatial_Tile_Check=Generate_Spatial_Occupancy(Placed_Tiles_Check,Delta_X_Check,Delta_Y_Check,1);
    //Declare but not initialise
    std::vector<int> Placed_Tiles_Compare, Spatial_Occupation_Compare, Spatial_Tile_Compare;
    
    bool broken_Comparison=false;
    for(int nth_repeat=1;nth_repeat<N_Repeated_Checks;++nth_repeat) {
      Placed_Tiles_Compare=Brute_Force_Polyomino_Builder(genome,THRESHOLD_SIZE,nth_repeat%(genome.size()/4),0);
      if(Placed_Tiles_Compare.empty()) //STERIC
        return -1;
      if(Placed_Tiles_Compare.size()/4 > THRESHOLD_SIZE)
        return -1;
      if(Placed_Tiles_Check.size()!=Placed_Tiles_Compare.size())
        return -5;
      
      int Delta_X_Compare=-1,Delta_Y_Compare=-1;
      Spatial_Occupation_Compare=Generate_Spatial_Occupancy(Placed_Tiles_Compare,Delta_X_Compare,Delta_Y_Compare,0);
      if(!Compare_Two_Polyominoes_Shapes(Spatial_Occupation_Check,Delta_X_Check,Delta_Y_Check,Spatial_Occupation_Compare,Delta_X_Compare,Delta_Y_Compare)) {
        broken_Comparison=true;
        break;
      }
      else {
        Spatial_Tile_Compare=Generate_Spatial_Occupancy(Placed_Tiles_Compare,Delta_X_Compare,Delta_Y_Compare,1);
        if(!Compare_Two_Polyominoes_Shapes(Spatial_Tile_Check,Delta_X_Check,Delta_Y_Check,Spatial_Tile_Compare,Delta_X_Compare,Delta_Y_Compare)) {                      
          broken_Comparison=true;
          break;
        }
      }
      if(N_Repeated_Checks-nth_repeat>1) {
        Placed_Tiles_Check=Placed_Tiles_Compare;
        Spatial_Occupation_Check=Spatial_Occupation_Compare;
        Spatial_Tile_Check=Spatial_Tile_Compare;
        Delta_X_Check=Delta_X_Compare;
        Delta_Y_Check=Delta_Y_Compare;
      }
    }    
    if(UND_Check) {
      return -5; 
    }
    else { //Bound
      if(broken_Comparison) { //BND
        return -5;
      }
      else {
        return Placed_Tiles_Check.size()/4;//Steric_Check(genome);
      }
    }
  }
  
  std::vector<int> Generate_Spatial_Occupancy(std::vector<int>& Placed_Tiles_Check, int& DELTA_X_Check,int& DELTA_Y_Check, int generate_mode) {
    std::vector<int> X_Locs_Check, Y_Locs_Check, Tile_Type_Check,Tile_Orientation_Check;
    for(std::vector<int>::iterator check_iter = Placed_Tiles_Check.begin();check_iter!=Placed_Tiles_Check.end();check_iter+=4) {
      X_Locs_Check.emplace_back(*check_iter);
      Y_Locs_Check.emplace_back(*(check_iter+1));
      Tile_Type_Check.emplace_back(*(check_iter+2));
      Tile_Orientation_Check.emplace_back(*(check_iter+3));
    }
    std::vector<int>::iterator LEFT_X_Check,RIGHT_X_Check,TOP_Y_Check,BOTTOM_Y_Check;
    std::tie(LEFT_X_Check,RIGHT_X_Check)=std::minmax_element(X_Locs_Check.begin(),X_Locs_Check.end());
    std::tie(BOTTOM_Y_Check,TOP_Y_Check)=std::minmax_element(Y_Locs_Check.begin(),Y_Locs_Check.end());
    DELTA_X_Check=*RIGHT_X_Check-*LEFT_X_Check+1;
    DELTA_Y_Check=*TOP_Y_Check-*BOTTOM_Y_Check+1;
    std::vector<int> Spatial_Occupancy_Check(DELTA_X_Check*DELTA_Y_Check);

    for(unsigned int tileIndex=0;tileIndex<X_Locs_Check.size();++tileIndex) {
      if(generate_mode==0) {
        Spatial_Occupancy_Check[(*TOP_Y_Check-Y_Locs_Check[tileIndex])*DELTA_X_Check + (X_Locs_Check[tileIndex]-*LEFT_X_Check)]=1;
      }
      if(generate_mode==1) {
      Spatial_Occupancy_Check[(*TOP_Y_Check-Y_Locs_Check[tileIndex])*DELTA_X_Check + (X_Locs_Check[tileIndex]-*LEFT_X_Check)]=Tile_Type_Check[tileIndex];
      }
      if(generate_mode==2) {
      Spatial_Occupancy_Check[(*TOP_Y_Check-Y_Locs_Check[tileIndex])*DELTA_X_Check + (X_Locs_Check[tileIndex]-*LEFT_X_Check)]=Tile_Orientation_Check[tileIndex];
      } 
    }
    return Spatial_Occupancy_Check;
  }
  int Compare_Two_Polyominoes_Tile_Details(std::vector<int>& Placed_Tiles_Check, std::vector<int>& Placed_Tiles_Compare) {
    //return Codes
    //1 - Not tile deterministic
    //2 - Tile, but not orientation determinsitic
    //3 - Tile and orientation determinstic 
    int return_code=3;
    for(std::vector<int>::iterator check_iter = Placed_Tiles_Check.begin(); check_iter!=Placed_Tiles_Check.end(); check_iter+=4) {
      for(std::vector<int>::iterator compare_iter = Placed_Tiles_Compare.begin(); compare_iter!=Placed_Tiles_Compare.end(); compare_iter+=4) {
        if(*check_iter==*compare_iter && *(check_iter+1) == *(compare_iter+1)) {
          if(*(check_iter+2)!=*(compare_iter+2)) {
            return 1;
          }
          if(*(check_iter+3)!=*(compare_iter+3)) {
            return_code=2;
          }
        }
      }
    }
    return return_code;
  }
  bool Compare_Two_Polyominoes_Shapes(std::vector<int>& Spatial_Occupation_Check,int Delta_X_Check,int Delta_Y_Check, std::vector<int>& Spatial_Occupation_Compare,int Delta_X_Compare,int Delta_Y_Compare) {
    if(std::accumulate(Spatial_Occupation_Check.begin(),Spatial_Occupation_Check.end(),0)!=std::accumulate(Spatial_Occupation_Compare.begin(),Spatial_Occupation_Compare.end(),0)) {
      return false;
    }
    if(Delta_X_Check==Delta_X_Compare && Delta_Y_Check==Delta_Y_Compare && Delta_X_Check==Delta_Y_Check) { //bounding boxes match, symmetric
      if(Spatial_Occupation_Check==Spatial_Occupation_Compare) {
        return true;
      }
      else {
        for(int rotation=0;rotation<3;++rotation) {
          Clockwise_Rotation(Spatial_Occupation_Check,Delta_X_Check,Delta_Y_Check);
          if(Spatial_Occupation_Check==Spatial_Occupation_Compare) {
            return true;
          }
        }
        return false;
      }
    }
    if(Delta_X_Check==Delta_X_Compare && Delta_Y_Check==Delta_Y_Compare) { //bounding boxes match, asymmetric
      if(Spatial_Occupation_Check==Spatial_Occupation_Compare) {
        return true;
      }
      else {
        Clockwise_Pi_Rotation(Spatial_Occupation_Check,Delta_X_Check,Delta_Y_Check);
        if(Spatial_Occupation_Check==Spatial_Occupation_Compare) {
          return true;
        }
        return false;
      }
    }
    if(Delta_X_Check==Delta_Y_Compare && Delta_Y_Check==Delta_X_Compare) { //bounding boxes pi/2 off, asymmetric
      Clockwise_Rotation(Spatial_Occupation_Check,Delta_X_Check,Delta_Y_Check);
      std::swap(Delta_X_Check,Delta_Y_Check);
      if(Spatial_Occupation_Check==Spatial_Occupation_Compare) {
        return true;
      }
      else {
        Clockwise_Pi_Rotation(Spatial_Occupation_Check,Delta_X_Check,Delta_Y_Check);
        if(Spatial_Occupation_Check==Spatial_Occupation_Compare) {
          return true;
        }
        return false;
      }
    }
    //nomiminally else, completely mismatch
    return false;
  }



  
  void Brute_vs_Graph_Methods_Comparison_Random_Sample(unsigned int genome_Length,int MAX_C,int num_Samples, int K_Repetitions) {
    std::cout<<"Running for tile size of "+std::to_string(genome_Length/4)+" for "+std::to_string(num_Samples)+" iterations"<<std::endl;
    std::string outFile_Title="/rscratch/asl47/Method_Comparisons_Sample_T"+std::to_string(genome_Length/4)+"_Trial.txt";
    std::ofstream Mismatch_outFile(outFile_Title,std::ios_base::app | std::ios_base::out);
    Mismatch_outFile<<"N_t:"<<genome_Length/4<<", C="<<MAX_C<<" ||  K="<<K_Repetitions<<", N:"<<num_Samples<<"\n";
    bool print_detailed=true;
    int numRun=0, numDelB=0, numDelG=0;
    std::uniform_int_distribution<int> Color(0,MAX_C-1);
#pragma omp parallel for  schedule(dynamic)
    for(int n=0;n<num_Samples;++n) { 
      std::vector<int> genome(genome_Length);
      do {
        genome.assign(genome_Length,0);
        for(unsigned int i=0;i<genome_Length;++i)
          genome[i]=Color(RNG_Generator);
  
        Clean_Genome(genome);
      }while(genome.size()!=genome_Length);
      
      if(false &&n>14600000) {
        for(auto& g:genome) {
          std::cout<<g<<" ";
        }
        std::cout<<std::endl;
      }
      
      std::vector<int> Bgenome(genome);
      int B_result=Brute_Force::Analyse_Genotype_Outcome(genome,K_Repetitions);
      int G_result;

      if(Disjointed_Check(genome)) {
        G_result=-5;
      }
      else {
        G_result=Graph_Analysis(genome);
        //G_result=FastGraphAnalysis(genome);
        if(G_result>0) {
          G_result=Steric_Check(genome,-1);
        }
      }
#pragma omp critical(mapping)
      {
        ++numRun;
        if(numRun%(num_Samples/100)==0) {
          std::cout<<"on "<<numRun<<std::endl;
        }
        if(B_result!=G_result) {
          if(B_result>0 && G_result<=0) {
            ++numDelB;
            if(print_detailed) {
              Mismatch_outFile<<"Mismatch on: ";
              int tinc=0;
              for(auto& g:Bgenome) {
                Mismatch_outFile<<g<<", ";
                ++tinc;
                if(tinc==4) {
                  Mismatch_outFile<<"  ";
                  tinc=0;
                }
              }
              Mismatch_outFile<<"|| with B:"<<B_result<<" and G:"<<G_result<<"\n";
            }
          }
          if(B_result<0 && G_result>0) {
            ++numDelG;
            if(print_detailed) {
              Mismatch_outFile<<"Mismatch on: ";
              int tinc=0;
              for(auto& g:Bgenome) {
                Mismatch_outFile<<g<<", ";
                ++tinc;
                if(tinc==4) {
                  Mismatch_outFile<<"  ";
                  tinc=0;
                }
              }
              Mismatch_outFile<<"|| with B:"<<B_result<<" and G:"<<G_result<<"\n";
            }
          }
        }
      }
    }
    Mismatch_outFile<<"Fail_Rate_B>G:"<<numDelB<<"     Fail_Rate_B<G:"<<numDelG<<"\n";
    Mismatch_outFile.close();
  }

  void Brute_vs_Brute_Comparison_Random_Sample(unsigned int genome_Length,int MAX_C,int num_Samples, int K_Repetitions) {

    std::cout<<"Running BvB for tile size of "+std::to_string(genome_Length/4)+" for "+std::to_string(num_Samples)+" iterations"<<std::endl;
    
    std::string outFile_Title="/rscratch/asl47/Method_Comparisons_BvB_T"+std::to_string(genome_Length/4)+"_Run_1.txt";
    std::ofstream Mismatch_outFile(outFile_Title,std::ios_base::app | std::ios_base::out);
    
    Mismatch_outFile<<"N_t:"<<genome_Length/4<<", C="<<MAX_C<<" ||  K="<<K_Repetitions<<", N:"<<num_Samples<<"\n";
   
    int numRun=0, numDelB=0, numDelG=0;
#pragma omp parallel for  schedule(dynamic)
    for(int n=0;n<num_Samples;++n) { 
      std::vector<int> genome(genome_Length);
      do {
        genome.assign(genome_Length,0);
        int upper=1;
        for(unsigned int i=0;i<genome_Length;++i) { 
          if(*std::max_element(genome.begin(),genome.begin()+i)==upper) 
            upper=std::min(upper+2,MAX_C-1);
          std::uniform_int_distribution<int> Color(0,upper);
          genome[i]=Color(RNG_Generator);
        }
        Clean_Genome(genome);
      }while(genome.size()!=genome_Length);

      std::vector<int> Bgenome(genome);
      int B_Perfect_result=Brute_Force::Analyse_Genotype_Outcome(genome,250);
      int B_Test_result=Brute_Force::Analyse_Genotype_Outcome(Bgenome,K_Repetitions);

#pragma omp critical(mapping)
      {
        ++numRun;
        if(numRun%(num_Samples/100)==0) {
          std::cout<<"on "<<numRun<<std::endl;
        }
        if((B_Test_result>0 && B_Perfect_result<0) || (B_Test_result<0 && B_Perfect_result>0 )) {
          if(B_Test_result>0 && B_Perfect_result<0) {
            ++numDelG;
          }
          else {
            ++numDelB;
          }
          int tinc=0;
          for(auto& g:Bgenome) {
            Mismatch_outFile<<g<<", ";
            ++tinc;
            if(tinc==4) {
              Mismatch_outFile<<"  ";
              tinc=0;
            }
          }
          Mismatch_outFile<<"|| B_P:"<<B_Perfect_result<<" B_T:"<<B_Test_result<<"\n";          
        }
      }
    }
    
    Mismatch_outFile<<"FR P>T:"<<numDelB<<" FR P<T:"<<numDelG<<"\n";
    Mismatch_outFile.close();
  }


}

void GenotypeSpaceSampler(int N_samples) {
  std::uniform_int_distribution<int> Color(0,199);
  std::random_device rd;
  xorshift RNG_source(rd());

  int a=0,b=0,c=0,ab=0,bc=0,ac=0,abc=0;
#pragma omp parallel for schedule(dynamic),reduction(+:a,b,c,ab,bc,ac,abc)
  for(int n=0;n<N_samples;++n) {
    std::vector<int> targets={4,5,6};
    std::vector<double> fitnesses(3);
    std::vector<int> genotype(20);
    for(int b=0;b<20;++b)
      genotype[b]=Color(RNG_source);
    
    if(GetMultiplePhenotypeFitness(genotype,targets,fitnesses,3)) {    
      if(fitnesses[0]==1) {
        if(fitnesses[1]==1) {
          if(fitnesses[2]==1) {
            ++abc;
          }
          else {
            ++ab;
          }
        }
        else {
          if(fitnesses[2]==1) {
            ++ac;
          }
          else {
            ++a;
          }
        }
      }
      else {
        if(fitnesses[1]==1) {
          if(fitnesses[2]==1) {
            ++bc;
          }
          else {
            ++b;
            
          }
        }
        else {
          if(fitnesses[2]==1) {
            ++c;
          }
        }
      }
    } 

  }
  std::cout<<"G Sample\nABC: "<<abc<<"\nAB: "<<ab<<"\nBC: "<<bc<<"\nAC: "<<ac<<"\nA: "<<a<<"\nB: "<<b<<"\nC: "<<c<<std::endl; 
}

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

                  


int main(int argc, char* argv[]) {

  std::vector<std::vector<int> > GENOME_VECTOR(1);
  std::vector<int> test_genome;
  int graph_result;
  if(argc>1) {
    switch(argv[1][1]) {
    case 'G':
      GenotypeSpaceSampler(std::stoi(argv[2]));
      break;
    case 'R':
      //Topology_Robustness({0,0,0,1,2,2,3,4},std::stoi(argv[2]),{6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 38, 42, 50, 60, 70, 84, 100, 150, 250, 500},std::stoi(argv[3]));
      Topology_Robustness({0,0,0,1,2,2,2,2},std::stoi(argv[2]),{6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 38, 42, 50, 60, 70, 84, 100, 150, 250, 500},std::stoi(argv[3]));
      break;
    case 'C':
      Brute_Force::Brute_vs_Graph_Methods_Comparison_Random_Sample(std::stoi(argv[2])*4,std::stoi(argv[2])*4+2,std::stoi(argv[3]),25);
      break;
    case 'S':
      for(int i=0;i<10;++i) {
        //std::cout<<"STARTING RUN "<<std::to_string(i)<<std::endl;
        Run_Selection(std::stoi(argv[2])*4,std::stoi(argv[3]));
        //std::cout<<"Safe exit"<<std::endl;
      }
      break;
    case 'X':
      Generate_Random_Genotypes(std::stoi(argv[2])*4,GENOME_VECTOR,false,2);
      Timing_Mode(GENOME_VECTOR,{-1,-2});
      break;
    case 'T':
      test_genome={0, 0, 2, 19,   0, 14, 20, 18,   0, 0, 17, 5,   0, 0, 6, 1,   0, 20, 18, 13};
        //{0, 4, 8, 7,   0, 9, 3, 2,   0, 10, 0, 13,   0, 0, 1, 14};
        //{0, 11, 6, 4, 0, 0, 3, 7,      0, 0, 8, 14,   0, 13, 5, 12};
        //{0, 1, 3, 16,   0, 0, 6, 9,   0, 0, 15, 5,   0, 2, 10, 4};
        //{27, 61, 93, 90, 0, 0, 28, 176, 0, 175, 109, 184, 0, 110, 89, 197, 0, 198, 3, 183, 0, 0, 4, 94, 0, 38, 0, 62, 0, 37, 124, 123};
      graph_result=Graph_Analysis(test_genome);
      std::cout<<"GR "<<graph_result<<std::endl;
      break;
    case 'H':
      std::cout<<"Brute Tile running options\n -R [# runs] [run #] for topology robustness\n -C [# tiles] [# runs] for brute v graph comparison\n -S [# tiles] [# runs] for run selection\n"<<std::endl;
    }
  }

  

  
  
}


void Run_Selection(int genome_length,int Test_Cases) {
  std::vector<int> K_s{0};//,10};//,5,10,20};
  //std::cout<<"Running for T "<<std::to_string(genome_length/4)<<" for "<<std::to_string(Test_Cases)<<" cases, with Ks: ";
  //for(int k:K_s) {
  //  std::cout<<k<<" ";
  //}
  //std::cout<<std::endl;

  double cpu_start,cpu_finish;
  for(int G_C=0;G_C<1; ++G_C) {
    if(genome_length>28 && G_C==0) //skip B/D for T>7
      continue;
    int modified_cases= G_C==0 ? 1000000 : 5000000;//Test_Cases; //* pow(5,G_C);
    std::vector<std::vector<int> > GENOME_VECTOR(modified_cases);
    //std::cout<<"Generating "+std::to_string(modified_cases)+" with condition "+std::to_string(G_C)+" ...";
    cpu_start= get_cpu_time();
    Generate_Random_Genotypes(genome_length,GENOME_VECTOR,false,G_C);
    cpu_finish  = get_cpu_time();
    double LOAD_TIME=cpu_finish-cpu_start;
    std::cout<<"Generated after "<<LOAD_TIME<<" seconds"<<std::endl;
    Timing_Mode(GENOME_VECTOR,K_s);
    continue;
    if(G_C==0) {
      for(auto& k:K_s) {
        std::cout<<"K: "<<k<<" Over: "<<0<<" Under: "<<0<<std::endl;
      }
    }
    else {
    Accuracy_Mode(GENOME_VECTOR,K_s);
    }
  }
}

void Generate_Random_Genotypes(int genome_Length,std::vector<std::vector<int> >& GENOME_VECTOR, bool triangularly_random,int condition, int C) {
  //std::cout<<"Note, running this for a fixed seed"<<std::endl;
  const int COLOUR_SPACE= C>0 ? C : genome_Length+1;
#pragma omp parallel for schedule(dynamic) 
  for(unsigned int t=0;t<GENOME_VECTOR.size();++t) {
    std::random_device rd;
    xorshift RNG_Genotype(rd());
    std::vector<int> genome(genome_Length);
    do{
      do {
        genome.assign(genome_Length,0);
        if(triangularly_random) {
          int upper=1;
          for(int i=0;i<genome_Length;++i) { 
            if(*std::max_element(genome.begin(),genome.begin()+i)==upper) 
              upper=std::min(upper+2,genome_Length+1);
            std::uniform_int_distribution<int> Color(0,upper);
            genome[i]=Color(RNG_Genotype);
          }
        }
        else {
          for(int i=0;i<genome_Length;++i) { 
            std::uniform_int_distribution<int> Color(0,COLOUR_SPACE);
            genome[i]=Color(RNG_Genotype);
          }
        }
        Clean_Genome(genome);
      }while(genome.size()!=genome_Length);
    }while(!Generating_Condition(genome,condition));
    GENOME_VECTOR[t]=genome;    
  }
}

bool Generating_Condition(std::vector<int>& genome,int condition) {
  switch(condition) {
  case 0: //Want B/D topologies only
    return (!Disjointed_Check(genome) && Graph_Analysis(genome)>0 && Steric_Check(genome,false)>0); //Brute_Force::Analyse_Genotype_Outcome(genome,50)>0;
  case 1: //Want connected topologies only
    return !Disjointed_Check(genome);
  default: //Want any topology
    return true;
  }
}


void Timing_Mode(std::vector<std::vector<int> >& GENOME_VECTOR, std::vector<int> k_repeats) {  
  double cpu_start,cpu_finish;
  //std::cout<<"Running \"Timing\""<<std::endl;
  cpu_start= get_cpu_time();
#pragma omp parallel for  schedule(dynamic)
  for(unsigned int t=0;t<GENOME_VECTOR.size();++t) {
    std::vector<int> genome=GENOME_VECTOR[t];
  }
  cpu_finish  = get_cpu_time();
  double LOAD_TIME=cpu_finish-cpu_start;
  //std::cout<<"Load time was "<<LOAD_TIME<<std::endl;
  for(int k:k_repeats) {
    int BD=0;
    //std::cout<<"Running for K: "<<k<<" ------> ";
    cpu_start= get_cpu_time();
#pragma omp parallel for  schedule(dynamic,100) reduction(+:BD)
    for(unsigned int t=0;t<GENOME_VECTOR.size();++t) {
      std::vector<int> genome=GENOME_VECTOR[t];
      Steric_Check(genome,-1);
      
      continue; /*TEMPORARY*/
      
      if(k>0) {
        if(Brute_Force::Analyse_Genotype_Outcome(genome,k)>0)
          ++BD;
      }
      else {
        if(k==-1) {
          if(!Disjointed_Check(genome) && Graph_Analysis(genome)>0 && Steric_Check(genome,false)>0) 
            ++BD;
        }
        else {
          if(!Disjointed_Check(genome) && FastGraphAnalysis(genome)>0 && Steric_Check(genome,false)>0) 
            ++BD;
        }
          
      }
    }
    cpu_finish  = get_cpu_time();
    std::cout<<"Runtime taken: "<<cpu_finish-cpu_start-LOAD_TIME<<" seconds"<<std::endl;
    //std::cout<<"B/D: "<<BD<<" and UND: "<<GENOME_VECTOR.size()-BD<<std::endl;
  }
}

void Accuracy_Mode(std::vector<std::vector<int> >& GENOME_VECTOR, std::vector<int> k_repeats) {
  std::cout<<"Running \"Accuracy\""<<std::endl;
  std::map<int,int> Over_Results;
  std::map<int,int> Under_Results;
  const bool PRINTER=false;
#pragma omp parallel for  schedule(dynamic)
  for(unsigned int t=0;t<GENOME_VECTOR.size();++t) {
    std::vector<int> genome=GENOME_VECTOR[t];
    int TRUTH_RESULT=0;//Brute_Force::Analyse_Genotype_Outcome(genome,150);
    std::vector<int> G_genome2(genome);
    if(Disjointed_Check(G_genome2))
      TRUTH_RESULT=-5;
    else {
      TRUTH_RESULT=Graph_Analysis(G_genome2);
      if(TRUTH_RESULT>0)
        TRUTH_RESULT=Steric_Check(G_genome2,false);
    }
    for(int k:k_repeats) {
      int Test_Result;
      if(k>0)
        Test_Result=Brute_Force::Analyse_Genotype_Outcome(genome,k);
      else {
        std::vector<int> G_genome(genome);
        if(Disjointed_Check(G_genome))
          Test_Result=-5;
        else {
          Test_Result=Graph_Analysis(G_genome);
          if(Test_Result>0)
            Test_Result=Steric_Check(G_genome,false);
        }
      }
#pragma omp critical(Mapping)
      {
        if(Test_Result>0 && TRUTH_RESULT<0) {
          if(PRINTER) {
          for(auto& g:genome) {
                std::cout<<g<<" ";
              }
          std::cout<<"Test "<<Test_Result<<" Truth "<<TRUTH_RESULT<<std::endl;
          }
          ++Over_Results[k];
        }
        if(Test_Result<0 && TRUTH_RESULT>0) {
          if(PRINTER) {
          for(auto& g:genome) {
                std::cout<<g<<" ";
              }
          }
          std::cout<<"Test "<<Test_Result<<" Truth "<<TRUTH_RESULT<<std::endl;
          ++Under_Results[k];
        }
        if(Test_Result>0 && TRUTH_RESULT>0) {
          if(Test_Result>TRUTH_RESULT) {
            if(PRINTER) {
            for(auto& g:genome) {
                std::cout<<g<<" ";
              }
            }
          std::cout<<"Test "<<Test_Result<<" Truth "<<TRUTH_RESULT<<std::endl;
            ++Over_Results[k];
          }
          if(Test_Result<TRUTH_RESULT) {
            if(PRINTER) {
            for(auto& g:genome) {
                std::cout<<g<<" ";
              }
          std::cout<<"Test "<<Test_Result<<" Truth "<<TRUTH_RESULT<<std::endl;
            }
            ++Under_Results[k];
          }
        }
      }
    }       
  }
  for(int k: k_repeats) {
    std::cout<<"K: "<<k<<" Over: "<<Over_Results[k]<<" Under: "<<Under_Results[k]<<std::endl;
  }
}

void Run_Over_BDs(std::vector<int> genome_Lengths,std::vector<int> Cs, int N_Cases) {
  for(auto& T:genome_Lengths) {    
    for(auto& C:Cs) {
      BD_Fractional_Occupation(T*4,C,N_Cases);
    }
  }
}

void BD_Fractional_Occupation(int genome_Length,int C,int N_Cases) {
  std::ofstream outFile("/rscratch/asl47/Frac_Oc/Fractional_Occupation_T"+std::to_string(genome_Length/4)+"_Run_3.txt",std::ios_base::app | std::ios_base::out);
  std::vector<std::vector<int> > GENOME_VECTOR(N_Cases);
  Generate_Random_Genotypes(genome_Length,GENOME_VECTOR,false,-1,C);
  int BD=0, BD2=0;
#pragma omp parallel for  schedule(dynamic,100) reduction(+:BD,BD2)
  for(unsigned int t=0;t<GENOME_VECTOR.size();++t) {
    std::vector<int> genome=GENOME_VECTOR[t];
    if(Get_Phenotype_Fitness(genome,false)>0)
      ++BD2;
    if(!Disjointed_Check(genome) && Graph_Analysis(genome)>0 && Steric_Check(genome,false)>0) 
      ++BD;
  }
  outFile << "C: " << std::to_string(C)<< " "<<BD <<" "<<BD2<<"\n";
  outFile.close(); 
}

void Topology_Robustness(std::vector<int> genotype,int N_runs, std::vector<int> Colours,const int R_N) {
  std::string outName="/rscratch/asl47/Topology_Robustness_SC_R"+std::to_string(R_N)+".txt";
  std::ofstream outFile(outName,std::ios_base::app | std::ios_base::out);
  outFile<<"Running topology robustness with "<<N_runs<<" runs\n";
  const int FITNESS_TARGET=Steric_Check(genotype,false);
  static int faces[8]={0,1,2,3,4,5,6,7};
#pragma omp threadprivate(faces)
  std::random_device rd;
  xorshift RNG_Shuffler(rd());

  for(int N_mut=0;N_mut<=8;++N_mut) {
    if(N_mut>0)
      outFile<<"\n";
    outFile<< "N_mutations: "<<N_mut<<"\n";

    for(int C : Colours) {
      outFile<<"C: "<<C<<" ";
      if(N_mut==0) {
        outFile <<"R: "<<N_runs<<" ";
        continue;
      }
      std::uniform_int_distribution<int> Random_Colour(0,C-1);
      int robust=0;
#pragma omp parallel for schedule(dynamic) reduction(+:robust)
      for(int n=0;n<N_runs;++n) {
	std::shuffle(&faces[0],&faces[8],RNG_Shuffler);
	std::vector<int> mutated_genotype(genotype);
	for(int mut=0;mut<N_mut;++mut) {
	  int prev_c=mutated_genotype[faces[mut]];
	  do {
	    mutated_genotype[faces[mut]]=Random_Colour(RNG_Shuffler);
	  }while(mutated_genotype[faces[mut]]==prev_c); 
	}
        Clean_Genome(mutated_genotype);
        if(!Disjointed_Check(mutated_genotype) && Graph_Analysis(mutated_genotype)>0 && Steric_Check(mutated_genotype,false)==FITNESS_TARGET)
          ++robust;
      }
      outFile <<"R: "<<robust<<" ";
    }
  }
  outFile<<"\n";
  outFile.close();
}


void Purify_Topologies() {
  std::ifstream infile("//scratch//asl47//Search_Spaces//SearchSpace_2_8.txt");//Topology_3_14_RAW.txt");
  std::ofstream outFile("//scratch//asl47//Search_Spaces//SearchSpace_2_8_P.txt");//Topology_3_14_Purified3.txt");
  int a,b,c,d,e,f,g,h;//,i,j,k,l;
  std::vector<std::vector<int> > BUFFER_READ;
  std::vector<std::vector<int> > Uniques;
  std::cout<<"reading"<<std::endl;
  while (infile >>a>>b>>c>>d>>e>>f>>g>>h) {//>>i>>j>>k>>l) {
    BUFFER_READ.push_back({a,b,c,d,e,f,g,h});//,i,j,k,l});
  }
  infile.close();
  std::cout<<"Read "<<BUFFER_READ.size()<<" genotypes"<<std::endl;
  std::cout<<"running"<<std::endl;
  for(unsigned int n=0;n<BUFFER_READ.size();++n) {
    if(n%100000==0) {
      std::cout<<"on "<<n<<std::endl;
    }
    std::vector<int> genome=BUFFER_READ[n];
    std::vector<int> genome2(genome);
    Clean_Genome(genome,13,false);
    if(genome==genome2) {
      Uniques.push_back(genome2);
      for(auto& g:genome2) {
        outFile<<g<<" ";
      }
      outFile<<"\n";
    }
  }
  outFile.close();
}

