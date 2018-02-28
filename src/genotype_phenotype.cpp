#include "genotype_phenotype.hpp"
#include <sstream>
#include <iterator>

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

extern "C" void WrappedGetPhenotypesID(char* a) {

  std::ofstream fout("Genotype_Codes.txt", std::ios_base::out);
  std::ofstream fout2("Phenotype_Table.txt", std::ios_base::out);

  StochasticPhenotypeTable pt;
  int k_builds=10;
  std::vector<int> genotype;
  
  std::string filename(a);
  std::ifstream file(filename);
  std::string str; 
  while (std::getline(file, str)) {
    std::istringstream is( str );
    genotype.assign( std::istream_iterator<int>( is ), std::istream_iterator<int>() );
    for(uint8_t seed=0;seed<genotype.size()/4;++seed) {
      Phenotype_ID phen_id=Stochastic::Analyse_Genotype_Outcome(genotype,k_builds,&pt,0);
      fout<<+phen_id.first<<" "<<+phen_id.second<<" ";
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
  int max_col1=5;
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
                  for(int seed=0;seed<genotype.size()/4;++seed) {
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
  std::vector<int> g{0,0,1,0, 2,0,0,0, 2,0,3,4};
  //PhenTuples(g);
  ExhaustivePhen();
  return 0;


}




