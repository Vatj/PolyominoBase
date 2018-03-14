#include "stochastic_model.hpp"
#include <iostream>

extern "C" void WrappedGetPhenotypesID(const char* a,const char* b,bool file_of_genotypes,uint8_t colours,uint8_t n_genes);
extern "C" void GGenerator(const char* a,bool file_of_genotypes,uint8_t colours,uint8_t n_genes);
extern "C" void SampleMinimalGenotypes(uint8_t n_genes, uint8_t colours,const uint32_t N_SAMPLES,bool allow_duplicates);

uint64_t genotype_to_index(uint8_t* genotype, uint8_t colours, uint8_t n_genes);
void index_to_genotype(uint64_t index, uint8_t* genotype, uint8_t colours, uint8_t n_genes);


struct NecklaceFactory {
  uint8_t colours=1;
  std::vector<std::vector<uint8_t> > necklaces;
  std::vector<uint8_t> necklace_grower;

  NecklaceFactory()  {necklace_grower.assign(5,0);}

  void GenNecklaces(uint8_t c) {
    colours=c;
    crsms_gen(1,1);
  }
  
  bool is_finite_necklace(std::vector<uint8_t>& neck) {
    if((neck[1] && Interaction_Matrix(neck[1])==neck[3]) || (neck[2] && Interaction_Matrix(neck[2])==neck[4]))
      return false;
    return true;
  }
  
  void is_necklace(uint64_t j) {
    if(4%j==0)
      if(is_finite_necklace(necklace_grower))
        necklaces.emplace_back(std::vector<uint8_t>{necklace_grower.begin()+1,necklace_grower.end()});
  }
  
  void crsms_gen(uint64_t n, uint64_t j) {
    if(n>4)
      is_necklace(j);
    else {
      necklace_grower[n]=necklace_grower[n-j];
      crsms_gen(n+1,j);
      for(uint64_t i=necklace_grower[n-j]+1;i<colours;++i) {
        necklace_grower[n]=i;
        crsms_gen(n+1,n);
      }
    }
  }
};

struct GenotypeGenerator {
  bool is_done=false;
  uint8_t n_genes,colours;
  std::vector<uint32_t> necklace_states;
  std::vector<std::vector<uint8_t> > necklaces;
  uint32_t n_necklaces;

  void init() {
    NecklaceFactory necks=NecklaceFactory();
    necks.GenNecklaces(colours);
    necklaces=necks.necklaces;
    n_necklaces=necklaces.size();
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
    //std::vector<uint8_t> int_pairs((colours-1)/2);
    //for(uint8_t g : genotype)
    //++int_pairs[(g-1)/2]
    
    return true;
  }
  bool valid_genotype(std::vector<uint8_t>& genotype) {
    if(!valid_growing_faces(genotype))
      return false;
    if(!valid_bindings(genotype))
      return false;
    return true;
  }
  
  void increment_states(std::vector<uint32_t>& states) {
    ++states.back();
    for(uint32_t rind=states.size();rind>0;--rind) {
      if(states[rind-1]>=n_necklaces) {
        if(rind==1) {
          is_done=true;
          return;
        }
        else { 
          states[rind-1]=0;
          ++states[rind-2];
        }
      }
    }
    auto max_iter=std::max_element(states.begin(),states.end());
    std::replace(max_iter,states.end(),static_cast<uint32_t>(0),*max_iter);
  }

  std::vector<uint8_t> next_genotype() {
    std::vector<uint8_t> genotype;
    while(!is_done) {
      genotype.clear();
      genotype.reserve(n_genes*4);
      for(auto index : necklace_states)
        genotype.insert(genotype.end(),necklaces[index].begin(),necklaces[index].end());
      increment_states(necklace_states);
      if(valid_genotype(genotype)) 
        return genotype;
    }   
    genotype.clear();
    return genotype;
  }
};
