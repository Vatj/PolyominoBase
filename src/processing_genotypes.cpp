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

  
  for(int run : runs) {
    std::cout<<"Merging for run "<<run<<std::endl;
    std::string in_file_name="//scratch//asl47//Data_Runs//Dynamic//"+folder_base+"//"+sub_file+"_Run"+std::to_string(run)+"_Genotype.txt";
    std::ifstream in_file(in_file_name);
    std::string line_input;
    bool first_line=true;
    while (std::getline(in_file, line_input)) {
      if(first_line) {
        first_line=false;
        continue;
      }
      std::string line_buffer(line_input,line_input.find('x')+1,std::string::npos);
      std::istringstream string_stream(line_buffer);
      std::vector<int> genome{std::istream_iterator<int>(string_stream), std::istream_iterator<int>()};
      if(Disjointed_Check(genome)) {
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
        int steric_result=Steric_Check_Table(genome,shape_table,num_shapes);
      }
    }
  }


  std::string out_file_name;
  std::string out_file_name_lengths;
  if(folder_base=="T3")
    out_file_name ="//rscratch//asl47//Bulk_Run//Modular//Shape_Data_T3";
  else {
    out_file_name ="//rscratch//asl47//Bulk_Run//Modular//Shape_Data_T4_O";
    for(int s=3;s<folder_base.size();++s)
      out_file_name+=folder_base[s];
  }
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
  std::string file_base;
  std::string folder_base;
  std::string temp_input;
  if(argc>1) {
    switch(argv[1][1]) {
    case '3':
      file_base="Modular3_T20_C200_N1000_Mu0.003125_B5000";
      folder_base="T3";
      DoShape(file_base,folder_base,runs);
      break;
    case '4':
      temp_input=argv[3];
      file_base="Modular4_T20_C200_N1000_Mu0.003125_O"+temp_input+"_B5000";
      folder_base="T4O"+temp_input;
      //std::cout<<file_base<<std::endl;
      DoShape(file_base,folder_base,runs);
      break;
    case '5':
      //DoShape("Modular4_T20_C200_N1000_Mu0.003125_O1_B5000",runs);
      break;
    case '6':
      //DoShape("Modular4_T20_C200_N1000_Mu0.003125_O10_B5000",runs);
      break;
    default:
      std::cout<<"unknown"<<std::endl;
    }
  }
  else 
    std::cout<<"Wrong input"<<std::endl;
  
  



  
}
