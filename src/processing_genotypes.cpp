#include "processing_genotypes.hpp"

void PrintShapes(std::string fileName) {
  std::cout<<"hi"<<std::endl;
  std::ifstream input( "filename.ext" );

}

void DoShape(std::string sub_file,std::string folder_base,std::vector<int> runs) {
  std::vector<int> shape_table;
  int num_shapes=0;

  //std::cout<<sub_file<<std::endl;
  //std::cout<<folder_base<<std::endl;

  std::vector<int> component_lengths;

  std::string modularity_file_name="//rscratch//asl47//Processed/Dynamic//"+folder_base+"_Modular_Data.txt";
  std::ofstream mod_file_out(modularity_file_name);

  //std::cout<<"//scratch//asl47//Data_Runs//Dynamic_2//"+folder_base+"//"+sub_file+"_Run"+std::to_string(0)+"_Genotype.txt"<<std::endl;
  for(int run : runs) {
    //std::cout<<"Merging for run "<<run<<std::endl;
    std::string in_file_name="//scratch//asl47//Data_Runs//Dynamic_2//"+folder_base+"//"+sub_file+"_Run"+std::to_string(run)+"_Genotype.txt";
    std::ifstream in_file(in_file_name);
    std::string line_input;
    bool collect_input=false;
    while (std::getline(in_file, line_input)) {
      if(line_input[0]=='g') {
        if(line_input.find("10000")!=std::string::npos)
          collect_input=true;
        else
          collect_input=false;
        continue;
      }
      if(!collect_input)  
        continue;
      
      std::istringstream string_stream(line_input);
      std::vector<int> genome{std::istream_iterator<int>(string_stream), std::istream_iterator<int>()};
      if(Disjointed_Check(genome)) {
        std::cout<<"disjoined"<<std::endl;
        std::vector<int> shape_keys;
        std::vector<int>::iterator found_at=std::find(genome.begin(),genome.end(),-1);
        std::vector<int> partial_genome(genome.begin(),found_at);
        int genome_length=partial_genome.size();
        if(Graph_Analysis(partial_genome)<=0)
          continue;
        int steric_result=Steric_Check(partial_genome,-1);
        if(steric_result<=0)
          continue;

        mod_file_out << genome_length << " " <<steric_result << "\n";
        Steric_Check_Table(partial_genome,shape_table,num_shapes);
        component_lengths.emplace_back(partial_genome.size());
      }
      else {
        mod_file_out << genome.size()/4 << " " <<Steric_Check(genome,-1) << "\n";
        Steric_Check_Table(genome,shape_table,num_shapes);
        component_lengths.emplace_back(genome.size());
      }
    }
  }
  mod_file_out.close();


  std::string out_file_name;
  std::string out_file_name_lengths;
  out_file_name ="//rscratch//asl47//Processed/Dynamic//"+folder_base+"Shape_Data";
  out_file_name_lengths=out_file_name+"_L.txt";
  std::ofstream out_file_lengths(out_file_name_lengths);
  for(int length :component_lengths)
    out_file_lengths<<length<<" ";
  out_file_lengths.close();
  out_file_name+=".txt";
 
  
   
  std::ofstream out_file(out_file_name);
  //std::cout<<"** RESULTS **"<<std::endl;
  int shape_N=0;
  for(std::vector<int>::iterator shape_it = shape_table.begin(); shape_it != shape_table.end();) {
    out_file <<"*Code* "<<shape_N<<" *Deltas* ";//<<*it<<" "<<*(it+1)<<" Shape:";
    for(int p=0;p<*shape_it* *(shape_it+1)+2;++p) {
      out_file <<*(shape_it+p)<<" ";
      if(p==1)
        out_file<<" *Shape* ";
    }

    int dx=*shape_it,dy=*(shape_it+1);
    std::vector<int> check_shape_data(shape_it+2,shape_it+2+*shape_it* *(shape_it+1));
    std::vector<int> phenotype_Shape;
    phenotype_Shape.assign(((dx+2)*(dy+2)),0);
    for(int x_shift=0;x_shift<dx;++x_shift) {
      for(int y_shift=0;y_shift<dy;++y_shift) {
        phenotype_Shape[(dx+2)+1 + x_shift+y_shift*(dx+2)]=check_shape_data[x_shift+y_shift*dx];
      }
    }
    dx+=2;
    dy+=2;
    //std::cout<<dx<<","<<dy<<": ";
    //for(int s:phenotype_Shape)
    //  std::cout<<s<<" ";
    //std::cout<<std::endl;
    
    bool total_match=true;
    for(int sh=4;sh<7;++sh) {
      if(Shape_Matching_Fitness_Function(phenotype_Shape,dx,dy,sh)<1) {
        total_match=false;
        break;
      }
    }
    int symmetry_factor=-1;
    std::vector<int> phenotype_Shape_spare(phenotype_Shape);
    int dx2=dx,dy2=dy;

    Clockwise_Pi_Rotation(phenotype_Shape_spare,dx2,dy2);
    std::swap(dx2,dy2);
    if(phenotype_Shape==phenotype_Shape_spare)
        symmetry_factor=2;
    else
        symmetry_factor=0;

    if(symmetry_factor==2 && dx2==dy2) {
      Clockwise_Rotation(phenotype_Shape_spare,dx2,dy2);
      if(phenotype_Shape==phenotype_Shape_spare)
        symmetry_factor=4;
    }
 
    
    if(total_match)
      out_file<<" *total* 1";
    else
      out_file<<" *total* 0";

    switch(symmetry_factor) {
    case 0:
      out_file<<" *symmetry* 0";
      break;
    case 2:
      out_file<<" *symmetry* 2";
      break;
    case 4:
      out_file<<" *symmetry* 4";
      break;
    default:
      std::cout<<"symmetry factor is impossible"<<std::endl;
      break;      
    }

    
    out_file<<"\n";
    shape_it+=*shape_it* *(shape_it+1)+2;
    ++shape_N;
  }
  out_file.close();

}

std::random_device rd;
xorshift RNG_ENGINE(rd());

bool CheckRobust(std::vector<int> target,int target_sets,int num_targets) {
  std::uniform_int_distribution<int> mutated_base(0, 79);
  std::uniform_int_distribution<int> mutated_colours(0, 199);

  std::vector<int> targets;
  switch(target_sets) {
  case 0:
    targets={4,5,6};
    break;
  case 1:
    targets={5,6,4};
    break;
  case 2:
    targets={4,6,5};
    break;
  }
  std::vector<double> fitnesses(3);

  int chosen_base=mutated_base(RNG_ENGINE);
  int old_base_colour=target[chosen_base];
  do {
    target[chosen_base]=mutated_colours(RNG_ENGINE);
  }while(target[chosen_base]==old_base_colour);

  if(GetMultiplePhenotypeFitness(target,targets,fitnesses,num_targets)) {
    return true;
    if(std::accumulate(fitnesses.begin(),fitnesses.begin()+num_targets,0.0)==num_targets)
      return true;
    else
      return false;
  }
  else
    return false;
}



void DeleteriousRobustness(int r) {
  int N_REPS=1000;
  std::vector< std::vector<double> > complete_robustnesses(45);
  std::vector<int> robs(N_REPS);
  int slice_count=0;
  
  std::string in_file_name="//rscratch//asl47//Bulk_Run//Modular//A2_T20_C200_N5000_Mu0.003125_O125_K5000_Run"+std::to_string(r)+"_Robust.txt";
  std::ifstream in_file(in_file_name);
  std::string line_input;
  while (std::getline(in_file, line_input)) {
    if(line_input[0]=='x') {
      ++slice_count;
      continue;
    }
    else {
      std::istringstream string_stream(line_input);
      std::vector<int> genotype{std::istream_iterator<int>(string_stream), std::istream_iterator<int>()};
#pragma omp parallel for schedule(dynamic)
      for(int n=0;n<N_REPS;++n) {
        robs[n]=CheckRobust(genotype);        
      }
      complete_robustnesses[slice_count].push_back(std::accumulate(robs.begin(),robs.end(),0.0)/N_REPS);

    }
  }
  std::ofstream out_file("//rscratch//asl47//Processed//Dynamic//TN_rob_O125_r"+std::to_string(r)+".txt");

  for(int n=0;n<slice_count;++n) {
    out_file <<"s: " << n<<"\n";
    for(double rob_val : complete_robustnesses[n]) {
      out_file << rob_val<< " ";
    }
    out_file<< "\n";

  }
  out_file.close();
    
  
  

}





int main(int argc, char* argv[]) {
  int N=1000;
  if(argc>2) {
    N=std::stoi(argv[2]);
  }
  
  std::vector<int> runs(N);
  std::iota(runs.begin(),runs.end(),0);
  std::string file_base, folder_base,mu_t, mu,run_t, temp_input;
  std::map<int,std::string> Mu_map ={{1, "0.050000"}, {4, "0.012500"}, {8, "0.006250"}, {16, "0.003125"}, {32,"0.001563"}};

  std::vector<int> genotype;
  int num_rob=0,num_rob2=0,num_rob3=0,num_robx=0,num_robx2=0,num_robx3=0;;
  if(argc>1) {
    switch(argv[1][1]) {
    case '3':
      mu=Mu_map[std::stoi(argv[3])];
      mu_t=argv[3];
      file_base="Modular3_T20_C200_N500_Mu"+mu+"_K15000";
      folder_base="T3Mu"+mu_t;
      DoShape(file_base,folder_base,runs);
      break;
    case '4':
    case '5':
      run_t=argv[1][1];
      mu=Mu_map[std::stoi(argv[3])];
      mu_t=argv[3];
      temp_input=argv[4];
      file_base="Modular"+run_t+"_T20_C200_N500_Mu"+mu+"_O"+temp_input+"_K15000";
      folder_base="T"+run_t+"Mu"+mu_t+"O"+temp_input;
      //std::cout<<file_base<<std::endl;
      DoShape(file_base,folder_base,runs);
      break;

      break;
    case 'R':
      DeleteriousRobustness(N);
      break;
    case 'C':
      for(int i=0;i<10;++i) {
        num_rob=0;num_rob2=0;num_rob3=0;
        num_robx=0;num_robx2=0;num_robx3=0;
        switch(i) {
        case 0:
        default: //SIF zeroed
          genotype={0,0,0,1, 0,0,0,3, 0,0,0,5, 2,0,7,0, 8,6,9,4, 10,0,11,0, 12,0,13,0, 14,15,0,0, 0,19,18,0, 17,0,0,16, 0,21,0,20, 0,0,23,22, 24,25,0,0, 0,0,0,26, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
          break;   
        case 1: //SIF non-interacting
          genotype={27, 29, 31, 1, 33, 35, 37, 3, 39, 41, 43, 5, 2, 45, 7, 47, 8, 6, 9, 4, 10, 49, 11, 51, 12, 53, 13, 55, 14, 15, 57, 59, 61, 19, 18, 63, 17, 65, 67, 16, 69, 21, 71, 20, 73, 75, 23, 22, 24, 25, 77, 79, 81, 83, 85, 26, 87, 89, 91, 93, 95, 97, 99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 125, 127, 129, 131, 133};
          break;
        case 2: //compact zeroed
          genotype={0,0,1,0, 2,0,3,0, 4,2,5,2, 6,0,7,0, 8,9,2,0, 0,11,2,10, 0,13,0,12, 0,0,15,14, 16,2,0,0, 0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0};
          break;
        case 3: //compact non-interacting
          genotype={17, 19, 1, 21, 2, 23, 3, 25, 4, 2, 5, 2, 6, 27, 7, 29, 8, 9, 2, 31, 33, 11, 2, 10, 35, 13, 37, 12, 39, 41, 15, 14, 16, 2, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 125, 127, 129, 131, 133};
          break;
        case 4: //4 fold complex Z
          genotype={1,2,0,3, 0,4,0,5, 0,6,7,15, 8,0,9,0, 10,11,13,11, 0,12,0,0, 14,0,11,0, 0,16,0,17, 19,18,0,0, 0,0,20,11,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0};
          break;
        case 5: //4fold non-int
          genotype={1, 2, 21, 3, 23, 4, 25, 5, 27, 6, 7, 15, 8, 29, 9, 31, 10, 11, 13, 11, 33, 12, 35, 37, 14, 39, 11, 41, 43, 16, 45, 17, 19, 18, 47, 49, 51, 53, 20, 11, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 125, 127, 129, 131, 133};
          break;
        case 6: //single seeded
          genotype={0,1,0,0, 0,0,3,2, 4,5,0,0, 0,27,0,6, 0,19,7,28, 8,0,9,0, 10,13,15,11, 0,12,0,0, 0,0,0,14, 16,0,17,0, 18,0,0,0, 21,0,0,20, 0,23,22,0, 0,0,25,24, 26,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
          break;
        case 7: //single seeded
          genotype={29,1,31,33, 35,37,3,2, 4,5,39,41, 43,27,45,6, 47,19,7,28, 8,49,9,51, 10,13,15,11, 53,12,55,57, 59,61,63,14, 16,65,17,67, 18,69,71,73, 21,75,77,20, 79,23,22,81, 83,85,25,24, 26,87,89,91, 93,95,97,99, 101,103,105,107, 109,111,113,115, 117,119,121,123, 125,127,129,131};
          break;
        case 8: // sif with rank 1
          genotype={0,0,0,1, 0,0,0,3, 0,0,0,5, 2,0,7,0, 8,6,9,4, 10,0,11,0, 12,27,13,0, 14,15,0,0, 0,19,18,28, 17,0,0,16, 0,21,0,20, 0,0,23,22, 24,25,0,0, 0,0,0,26, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
          break;
        case 9: //SIF non-interacting
          genotype={27, 29, 31, 1, 33, 35, 37, 3, 39, 41, 43, 5, 2, 45, 7, 47, 8, 6, 9, 4, 10, 49, 11, 51, 12, 135, 13, 55, 14, 15, 57, 59, 61, 19, 18, 136, 17, 65, 67, 16, 69, 21, 71, 20, 73, 75, 23, 22, 24, 25, 77, 79, 81, 83, 85, 26, 87, 89, 91, 93, 95, 97, 99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 125, 127, 129, 131, 133};
          break;
        }
     
      
      
#pragma omp parallel for schedule(dynamic),reduction(+:num_rob,num_rob2,num_rob3)
        for(int n =0; n<N;++n) {
          if(CheckRobust(genotype,0,2))
            ++num_rob;
          //if(CheckRobust(genotype,0,3))
          //  ++num_robx;
          //if(CheckRobust(genotype,1,2))
          //  ++num_rob2;
          //if(CheckRobust(genotype,1,3))
          //  ++num_robx2;
          //if(CheckRobust(genotype,2,2))
          //  ++num_rob3;
          //if(CheckRobust(genotype,2,3))
          //  ++num_robx3;
        }
        
        std::cout<<"**Robustnesses for genotype "<<i<<"**\n Targets 4,5 -> "<<static_cast<float>(num_rob)/N << " || "<<static_cast<float>(num_robx)/N <<"\n Targets 5,6 -> " <<static_cast<float>(num_rob2)/N <<  " || "<<static_cast<float>(num_robx2)/N <<"\n Targets 6,4 -> "<<static_cast<float>(num_rob3)/N <<" || "<<static_cast<float>(num_robx3)/N << "\n";
        }
        break;
       
    default:
      std::cout<<"unknown"<<std::endl;
    }
  }
  else 
    std::cout<<"Wrong input"<<std::endl;
  
  



  
}
