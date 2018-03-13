#include "genotype_phenotype.hpp"
#include <sstream>
#include <iterator>

ulong *f;  // data in f[1..N],  f[0] = 0
ulong N;  // word length
ulong K;  // K-ary

struct NecklaceFactory {
  uint8_t colours=1;
  std::vector<std::vector<uint8_t> > necklaces;
  std::vector<uint8_t> necklace_grower;

  NecklaceFactory()  {necklace_grower.assign(5,0);}

  void GenNecklaces(uint8_t c) {
    colours=c;
    crsms_gen(1,1);
    std::cout<<necklaces.size()<<" necklaces generated"<<std::endl;
  }
  void is_necklace(uint8_t j) {
    if(4%j==0)
      necklaces.emplace_back(std::vector<uint8_t>{necklace_grower.begin()+1,necklace_grower.end()});
  }
  void crsms_gen(int n, int j) {
    if(n>4)
      is_necklace(j);
    else {
      necklace_grower[n]=necklace_grower[n-j];
      crsms_gen(n+1,j);
      for(uint8_t i=necklace_grower[n-j]+1;i<colours;++i) {
        necklace_grower[n]=i;
        crsms_gen(n+1,n);
      }
    }
  }

};

struct GenotypeGenerator {
  bool is_done=false;
  uint8_t n_genes,colours;
  std::vector<uint8_t> necklace_states;
  NecklaceFactory necks;

  void init() {
    necks=NecklaceFactory();
    necks.GenNecklaces(colours);
  }
  
  GenotypeGenerator(uint8_t a,uint8_t b) {n_genes=a;colours=b; necklace_states.assign(a,0);}

  std::vector<uint8_t> operator() () {
    return !is_done ? next_genotype() : std::vector<uint8_t>{};
  }

  bool valid_growing_faces(std::vector<uint8_t>& genotype) {
    uint8_t max_face=1;
    for(auto face : genotype) {
      if(face>max_face)
        return false;
      if(face==max_face)
        max_face+=2;
    }
    return true;
  }

  bool valid_bindings(std::vector<uint8_t>& genotype) {
    for(uint8_t interface=1;interface<=*std::max_element(genotype.begin(),genotype.end());interface+=2) {
      if(std::find(genotype.begin(),genotype.end(),interface)!=genotype.end()) { //is present
        if(std::find(genotype.begin(),genotype.end(),interface+1)==genotype.end()) //is not present
          return false;
      }
      else {
        if(std::find(genotype.begin(),genotype.end(),interface+1)!=genotype.end()) 
          return false;
      }
    }
    return true;
  }
  bool valid_genotype(std::vector<uint8_t>& genotype) {
    //std::cout<<"genotpe is ";
    //for(auto x : genotype)
    //  std::cout<<+x<<" ";
    //std::cout<<std::endl;
    if(!valid_growing_faces(genotype))
      return false;
    //std::cout<<"valid growing"<<std::endl;
    if(!valid_bindings(genotype))
      return false;
    //std::cout<<"valid binding"<<std::endl;
    return true;
  }
  void increment_states(std::vector<uint8_t>& states) {
    ++states.back();
    for(int8_t rind=states.size()-1;rind>=0;--rind) {
      if(states[rind]>=necks.necklaces.size()) {
        if(rind==0) {
          is_done=true;
          return;
        }
        else { 
          states[rind]=0;
          ++states[rind-1];
        }
      }

    }
   
    auto max_iter=std::max_element(states.begin(),states.end());
    std::replace(max_iter,states.end(),uint8_t(0),*max_iter);
  }

  std::vector<uint8_t> next_genotype() {
    std::vector<uint8_t> genotype;

    while(!is_done) {
      genotype.clear();
      genotype.reserve(n_genes*4);
      for(auto index : necklace_states)
        genotype.insert(genotype.end(),necks.necklaces[index].begin(),necks.necklaces[index].end());
      //std::cout<<"pre ";
      //for(auto x: necklace_states)
      //  std::cout<<+x<<" ";
      //std::cout<<std::endl;
      increment_states(necklace_states);
      //std::cout<<"post ";
      //for(auto x: necklace_states)
      //  std::cout<<+x<<" ";
      //std::cout<<std::endl;
      if(valid_genotype(genotype)) 
        return genotype;
        
    }
    genotype.clear();
    return genotype;
    
  }
  
  
  
};



/*
std::vector<int> GetPhenotypeID(std::vector<int>& genome,std::vector<int>& Known_Shapes, int Num_Shapes) {
  Clean_Genome(genome,-1);
  std::vector<int> IDs;
  if(Disjointed_Check(genome)) {
    while(std::find(genome.begin(),genome.end(),-1)!=genome.end()) {
      std::vector<int>::iterator foundAt=std::find(genome.begin(),genome.end(),-1);
      std::vector<int> partialGenome(genome.begin(),foundAt);
      genome.erase(genome.begin(),foundAt+1);
      if(Graph_Analysis(partialGenome)>0) {
        int steric_result=Steric_Check_Table(partialGenome,Known_Shapes,Num_Shapes);
        IDs.emplace_back(steric_result);
      }
      else 
        IDs.emplace_back(0);
    }
    
  }
  else  { //Single connected component
    if(Graph_Analysis(genome)>0) {
      int steric_result=Steric_Check_Table(genome,Known_Shapes,Num_Shapes);
      IDs.emplace_back(steric_result);
    }
    IDs.emplace_back(0);
  }
  return IDs;

}
*/
/*
extern "C" void WrappedGetPhenotypeID(int g_size, int* genotype,int k_size,int* K_Shapes, int Num_Shapes, int* IDs_p) {
  std::vector<int> genome(genotype,genotype+g_size);
  std::vector<int> Known_Shapes(K_Shapes,K_Shapes+k_size-(g_size*g_size)/4-2);
  Clean_Genome(genome,-1);
  std::vector<int> IDs;
  if(Disjointed_Check(genome)) {
    while(std::find(genome.begin(),genome.end(),-1)!=genome.end()) {
      std::vector<int>::iterator foundAt=std::find(genome.begin(),genome.end(),-1);
      std::vector<int> partialGenome(genome.begin(),foundAt);
      genome.erase(genome.begin(),foundAt+1);
      if(Graph_Analysis(partialGenome)>0) {
        int steric_result=Steric_Check_Table(partialGenome,Known_Shapes,Num_Shapes);
        IDs.emplace_back(steric_result);
      }
      else 
        IDs.emplace_back(-1);
    }
    
  }
  else  { //Single connected component
    if(Graph_Analysis(genome)>0) {
      int steric_result=Steric_Check_Table(genome,Known_Shapes,Num_Shapes);
      IDs.emplace_back(steric_result);
    }
    else
      IDs.emplace_back(-1);
  }



  while(Known_Shapes.size()<k_size) 
    Known_Shapes.emplace_back(-1);
  std::copy(Known_Shapes.begin(), Known_Shapes.end(), K_Shapes);
  
  while(IDs.size()<g_size/4)
    IDs.emplace_back(-2);
  std::copy(IDs.begin(), IDs.end(), IDs_p);
  
  //return IDs[0];

}
*/

std::map<uint8_t,uint8_t> DuplicateGenes(std::vector<int>& genome) {
  std::map<uint8_t,uint8_t> dups;
  for(int check_index=genome.size()/4-1;check_index>0;--check_index) {
    for(int compare_index=0;compare_index<check_index;++compare_index) {
      if(std::equal(genome.begin()+check_index*4,genome.begin()+check_index*4+4,genome.begin()+compare_index*4)) {
        genome.erase(genome.begin()+check_index*4,genome.begin()+check_index*4+4);
        dups[check_index]=compare_index;
        break;
      }
    }
  }
  return dups;
}

extern "C" uint64_t genotype_to_index(uint8_t* genotype, uint8_t colours, uint8_t n_genes) {
  uint64_t count=0;
  for(uint8_t index=0;index<n_genes*4;++index)
    count+= genotype[index] * pow(colours,n_genes*4-index-1);
  return count;
}

extern "C" void index_to_genotype(int index, int* genotype, int colours, int n_genes) {
  for(uint8_t count=0;count<n_genes*4;++count) {
    int value=index/pow(colours,n_genes*4-count-1);
    genotype[count]=value;
    index-= value * pow(colours,n_genes*4-count-1);
  }
}

extern "C" void GGenerator(const char* a,bool file_of_genotypes,uint8_t colours=0,uint8_t n_genes=0) {
  
  std::string filename(a);
  std::ofstream fout(filename);
   
  std::vector<uint8_t> geno,nullg;
  GenotypeGenerator ggenerator(n_genes,colours);
  ggenerator.init();
  while((geno=ggenerator())!=nullg) {
    if(file_of_genotypes) {
      for(auto x : geno)
        fout<<+x<<" ";
      fout<<"\n";
    }
    else {
      fout<<genotype_to_index(geno.data(),colours,n_genes)<<"\n";
    }
  }
}

extern "C" void WrappedGetPhenotypesID(const char* a,const char* b,bool file_of_genotypes,int colours=0,int n_genes=0) {

  std::string file_path(b);
  std::ofstream fout(file_path+"Genotype_Codes_N"+std::to_string(n_genes)+"_C"+std::to_string(colours)+".txt", std::ios_base::out);
  std::ofstream fout2(file_path+"Phenotype_Table_N"+std::to_string(n_genes)+"_C"+std::to_string(colours)+".txt", std::ios_base::out);

  StochasticPhenotypeTable pt;
  int k_builds=10;
  std::vector<int> genotype(n_genes*4);

 
  std::string filename(a);
  std::ifstream file(filename);
  std::string str; 
  while (std::getline(file, str)) {
    if(file_of_genotypes) {
      std::istringstream is( str );
      genotype.assign( std::istream_iterator<int>( is ), std::istream_iterator<int>() );
    }
    else {
      index_to_genotype(std::stoi(str),genotype.data(),colours,n_genes);
    }
    Clean_Genome(genotype,-1,false);
    std::map<uint8_t,uint8_t> dups=DuplicateGenes(genotype);
    std::map<uint8_t,Phenotype_ID> seed_m;
    for(uint8_t seed=0;seed<genotype.size()/4;++seed)
      seed_m[seed]=Stochastic::Analyse_Genotype_Outcome(genotype,k_builds,&pt,seed);
    uint8_t index=0;
    for(uint8_t gene=0;gene<n_genes;++gene) {
      if(dups.count(gene)) {
        Phenotype_ID phen_id=seed_m[dups[gene]];
        fout<<+phen_id.first<<" "<<+phen_id.second<<" ";
      }
      else {
        Phenotype_ID phen_id=seed_m[index];
        fout<<+phen_id.first<<" "<<+phen_id.second<<" ";
        ++index;
      }
    }
    fout<<"\n";
  }
  pt.PrintTable(fout2);

}

void ExhaustivePhen() {
  StochasticPhenotypeTable pt;
  int k_builds=3;

  std::ofstream fout("TestG.txt", std::ios_base::out);
  std::ofstream fout2("TestP.txt", std::ios_base::out);
  //int max_col1=5;
  int max_col2=7;
  
  for(int i=0;i<2;++i) {
    for(int j=0;j<4;++j) {
      for(int k=0;k<6;++k) {
        for(int l=0;l<max_col2;++l) {
          
          for(int q=0;q<max_col2;++q) {
            for(int w=0;w<max_col2;++w) {
              for(int e=0;e<max_col2;++e) {
                for(int r=0;r<max_col2;++r) {
                  
                  std::vector<int> genotype{i,j,k,l, q,w,e,r};
                  for(uint8_t seed=0;seed<genotype.size()/4;++seed) {
		    Phenotype_ID phen_id=Stochastic::Analyse_Genotype_Outcome(genotype,k_builds,&pt,0);
		    fout<<+phen_id.first<<" "<<+phen_id.second<<" ";
		  }
                  fout<<"\n";
                  
                }
              }
            }
          }
        }
      }
    }
  }
  pt.PrintTable(fout2);
  



}

/*

extern "C" void PhenTuples(int g_size, int* genotype_p,int k_builds, int* IDs_p) {
  std::vector<int> genotype(genotype_p,genotype_p+g_size);
 
  std::vector<int> outcomes;
  for(uint8_t seed=0;seed<genotype.size()/4;++seed) {
    *(IDs_p+seed)=Stochastic::Analyse_Genotype_Outcome(genotype,k_builds,seed);
  }
  
  for(auto x: outcomes)
    std::cout<<x<<" ";
  std::cout<<std::endl;
  



}
*/

int main(int argc, char* argv[]) {
  if(argc>2) {
    std::cout<<argv[1][0]<<std::endl;
  }
  auto nf =NecklaceFactory();
  nf.GenNecklaces(std::stoi(argv[2]));
  std::cout<<"necklaces"<<std::endl;
  for(auto x: nf.necklaces) {
    for(auto y: x)
      std::cout<<+y<<" ";
    std::cout<<std::endl;
  }
  std::vector<uint8_t> geno,nullg;
  GenotypeGenerator ggenerator(std::stoi(argv[1]),std::stoi(argv[2]));
  ggenerator.init();
  while((geno=ggenerator())!=nullg) {
    for(auto x : geno)
      std::cout<<+x<<" ";
    std::cout<<std::endl;
  }

  
  return 0;


}




