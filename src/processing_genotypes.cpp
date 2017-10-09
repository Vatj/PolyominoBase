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

  //std::cout<<"//scratch//asl47//Data_Runs//Dynamic_2//"+folder_base+"//"+sub_file+"_Run"+std::to_string(0)+"_Genotype.txt"<<std::endl;
  for(int run : runs) {
    //std::cout<<"Merging for run "<<run<<std::endl;
    std::string in_file_name="//scratch//asl47//Data_Runs//Dynamic_2//"+folder_base+"//"+sub_file+"_Run"+std::to_string(run)+"_Genotype.txt";
    std::ifstream in_file(in_file_name);
    std::string line_input;
    bool collect_input=false;
    while (std::getline(in_file, line_input)) {
      if(line_input[0]=='g') {
        if(line_input.find("10000")!=std::string::npos) {
            collect_input=true;
          }
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
        if(Graph_Analysis(partial_genome)<=0)
          continue;
        if(Steric_Check(partial_genome,-1)<=0)
          continue;
        Steric_Check_Table(partial_genome,shape_table,num_shapes);
        component_lengths.emplace_back(partial_genome.size());
      }
      else {
        Steric_Check_Table(genome,shape_table,num_shapes);
        component_lengths.emplace_back(genome.size());
      }
    }
  }


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
      std::cout<<"symmetry factor is weird?"<<std::endl;
      break;      
    }

    
    out_file<<"\n";
    shape_it+=*shape_it* *(shape_it+1)+2;
    ++shape_N;
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
    default:
      std::cout<<"unknown"<<std::endl;
    }
  }
  else 
    std::cout<<"Wrong input"<<std::endl;
  
  



  
}
