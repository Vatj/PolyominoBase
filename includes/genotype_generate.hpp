// #include "stochastic_model.hpp"
#include "genotype_phenotype.hpp"
#include <iostream>

/*External wrappers for python integration */
// extern "C"
// {
//   void ExhaustiveMinimalGenotypes(const char* file_path_c, uint8_t n_genes, uint8_t colours);
//   void SampleMinimalGenotypes(const char* file_path_c, uint8_t n_genes, uint8_t colours, const uint32_t N_SAMPLES, bool allow_duplicates);
// }

std::vector<Genotype> ExhaustiveMinimalGenotypes(uint8_t n_genes, uint8_t colours, StochasticPhenotypeTable* pt);
std::vector<Genotype> SampleMinimalGenotypes(uint8_t n_genes, uint8_t colours,
  const uint32_t N_SAMPLES, bool allow_duplicates, StochasticPhenotypeTable* pt);
std::vector<Genotype> ExhaustiveFullGenotypes2(uint8_t colours, StochasticPhenotypeTable* pt);

/*Minimal genotype methods*/
struct NecklaceFactory
{
  uint8_t colours=1;
  std::vector<std::vector<uint8_t> > necklaces;
  std::vector<uint8_t> necklace_grower;

  NecklaceFactory()  {necklace_grower.assign(5,0);}

  void GenNecklaces(uint8_t c);

  bool is_finite_necklace(std::vector<uint8_t>& neck);

  void is_necklace(uint64_t j);

  void crsms_gen(uint64_t n, uint64_t j);
};

struct GenotypeGenerator
{
  bool is_done=false;
  uint8_t n_genes,colours;
  std::vector<uint32_t> necklace_states;
  std::vector<std::vector<uint8_t> > necklaces;
  uint32_t n_necklaces;
  Genotype nullg;

  void init()
  {
    NecklaceFactory necks=NecklaceFactory();
    necks.GenNecklaces(colours);
    necklaces=necks.necklaces;
    n_necklaces=necklaces.size();
    necklace_states[0]=1;
  }

  GenotypeGenerator(uint8_t a,uint8_t b) {n_genes= a; colours=b; necklace_states.assign(a,0);}

  Genotype operator() ()
  {
    return !is_done ? next_genotype() : Genotype{};
  }

  bool valid_growing_faces(Genotype& genotype,uint8_t max_face)
  {
    for(auto face : genotype)
    {
      if(face>max_face)
        return false;
      if(face==max_face)
        max_face+=2;
    }
    return true;
  }

  bool valid_bindings(Genotype& genotype)
  {
    for(uint8_t interface=1; interface<=*std::max_element(genotype.begin(), genotype.end()); interface+=2)
    {
      if(std::find(genotype.begin(),genotype.end(),interface)!=genotype.end())
      { //is present
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

  bool valid_genotype(Genotype& genotype)
  {
    if(!valid_growing_faces(genotype,1))
      return false;
    if(!valid_bindings(genotype))
      return false;
    return true;
  }

  void increment_states(std::vector<uint32_t>& states)
  {
    ++states.back();
    uint32_t zero_state_init=states[0];
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
    if(zero_state_init==1 && states[0]!=1)
      states[1]=colours+2;
    if(zero_state_init==static_cast<uint32_t>(colours+2) && states[0]!=static_cast<uint32_t>(colours+2)) {
      is_done=true;
      return;
    }
    auto max_iter=std::max_element(states.begin(),states.end());
    std::replace(max_iter,states.end(),static_cast<uint32_t>(0),*max_iter);
  }

  Genotype next_genotype() {
    Genotype genotype;
    while(!is_done) {
      inc_lab:
      genotype.clear();
      genotype.reserve(n_genes*4);
      for(uint8_t index=0; index!= n_genes;++index) {
	uint8_t input_face=*std::max_element(genotype.begin(),genotype.end());
	uint8_t next_max_face= input_face>0 ? ((input_face-1)/2)*2+3 : 1;

	while(!valid_growing_faces(necklaces[necklace_states[index]],next_max_face)) {
	  const uint32_t post_inc = necklace_states[index]+1;
	  if(post_inc==necklaces.size()) {
	    increment_states(necklace_states);
	    goto inc_lab;
	  }
	  std::fill(necklace_states.begin()+index,necklace_states.end(),post_inc);
	}

        genotype.insert(genotype.end(),necklaces[necklace_states[index]].begin(),necklaces[necklace_states[index]].end());
      }
      increment_states(necklace_states);

      if(valid_genotype(genotype))
        return genotype;
    }
    genotype.clear();
    return genotype;
  }
};
