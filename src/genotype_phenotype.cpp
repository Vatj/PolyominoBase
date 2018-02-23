#include "genotype_phenotype.hpp"

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
void ExhaustivePhen() {
  PhenotypeTable pt;
  int k_builds=3;

  std::ofstream fout("TestG.txt", std::ios_base::out);
  std::ofstream fout2("TestP.txt", std::ios_base::out);
  /*
  std::vector<int> genotype2={0,0,1,0, 2,0,0,0, 2,3,4,0};
  for(uint8_t seed=0;seed<genotype2.size()/4;++seed)
    fout<<Stochastic::Analyse_Genotype_Outcome(genotype2,k_builds,&pt,seed)<<" ";
  fout<<"\n";
  pt.PrintTable(fout2);
  */
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
                  
                  std::vector<int> genotype={i,j,k,l, q,w,e,r};
                  for(int seed=0;seed<genotype.size()/4;++seed)
                    fout<<Stochastic::Analyse_Genotype_Outcome(genotype,k_builds,&pt,0)<<" ";
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

  /*
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
  */
}





/*
int gen_neck_next(std::vector<uint8_t> neck, const uint8_t m, const uint8_t n)
{
  int j; //index
  int i; //help index
  int r; //temporary remainder, for hand made modulo computation only

 SKIP: //previous prenecklace skipped

  //find rightmost element to increase
  j = n - 1;
  while(j >= 0 && neck[j] == m - 1)
    j--;

  //terminate if all elements are m - 1
  if(j < 0)
    return 1;

  //increase
  neck[j]++;

  //set right-hand elements
  for(i = j + 1; i < n; i++)
    neck[i] = neck[i - j - 1];

  //necklaces only
  r = n;
  j++;

  while(r >= j)
    r -= j;

  if(r != 0)
    goto SKIP; //skip this prenecklace

  return 0;
}
*/
