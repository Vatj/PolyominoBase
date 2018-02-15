#include "genotype_phenotype.hpp"
#include <fstream>

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


int main(int argc, char* argv[]) {
  std::vector<int> Known_Shapes;
  int Num_Shapes=0;
  
  std::ofstream fout("TestG.txt", std::ios_base::out);
  std::ofstream fout2("TestP.txt", std::ios_base::out);

  //TEST ALL GENOTYPES
  for(int i=0;i<7;++i) {
    for(int j=0;j<7;++j) {
      for(int k=0;k<7;++k) {
        for(int l=0;l<7;++l) {
          
          for(int q=0;q<7;++q) {
            for(int w=0;w<7;++w) {
              for(int e=0;e<7;++e) {
                for(int r=0;r<7;++r) {
                  std::vector<int> genotype={i,j,k,l, q,w,e,r};
                  for(auto ID : GetPhenotypeID(genotype,Known_Shapes,Num_Shapes))
                    fout<<ID<<" ";
                  fout<<"\n";
                }
              }
            }
          }
        }
      }
    }
  }

  //WRITE SHAPES TO FILE
  for(auto it=Known_Shapes.begin();it!=Known_Shapes.end();) {
    fout2<<*it<<" "<<*(it+1)<<" ";
    for(auto it2=it+2;it2<it+2+*(it)* (*(it+1));++it2)
      fout2<<*it2<<" ";
    fout2<<"\n";
    it=it+2+*(it)* (*(it+1));
  }
  
}
