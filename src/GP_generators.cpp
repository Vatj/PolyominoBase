#include "GP_generators.hpp"

std::mt19937 RNG_Engine(std::random_device{}());

uint64_t genotype_to_index(Genotype& genotype, uint8_t n_genes, uint8_t colours) {
  uint64_t count=0;
  for(uint8_t index=0;index<n_genes*4;++index)
    count+= genotype[index] * pow(colours,n_genes*4-index-1);
  return count;
}

void index_to_genotype(uint64_t index, Genotype& genotype, uint8_t n_genes, uint8_t colours) {
  for(uint8_t count=0;count<n_genes*4;++count) {
    uint64_t value=index/pow(colours,n_genes*4-count-1);
    genotype[count]=value;
    index-= value * pow(colours,n_genes*4-count-1);
  }
}

void SampleMinimalGenotypes(const char* file_path_c, uint8_t n_genes, uint8_t colours,const uint32_t N_SAMPLES,bool allow_duplicates, bool file_of_genotypes) {
  uint32_t good_genotypes=0;
  std::string file_path(file_path_c),file_details="_N"+std::to_string(n_genes)+"_C"+std::to_string(colours)+".txt";
  std::ofstream fout(file_path+"SampledGenotypes"+file_details);
  
  GenotypeGenerator ggenerator = GenotypeGenerator(n_genes,colours);
  ggenerator.init();
  StochasticPhenotypeTable pt;
  uint8_t k_builds=5;

#pragma omp parallel
  while(true) {
    uint32_t local_counter;
     Genotype genotype;
     std::vector<uint32_t> states;
#pragma omp atomic read
    local_counter = good_genotypes;
    if(local_counter >= N_SAMPLES)
      break;
    
    uint32_t lbound_neck=0;
    for(uint8_t neck_ind=0;neck_ind<n_genes;++neck_ind) {
      lbound_neck=std::uniform_int_distribution<uint32_t>{lbound_neck,ggenerator.n_necklaces-1}(RNG_Engine);
      states.emplace_back(lbound_neck);
      if(!allow_duplicates) {
        if((lbound_neck+n_genes-neck_ind)>ggenerator.n_necklaces)
          continue;
        ++lbound_neck;
      }
    }
    for(auto state : states)
      genotype.insert(genotype.end(),ggenerator.necklaces[state].begin(),ggenerator.necklaces[state].end());
    
    if(!ggenerator.valid_genotype(genotype))
      continue;

    bool all_zero_phen=true;
    for(uint8_t seed=0;seed<n_genes;++seed) {
      uint8_t phen_size = Stochastic::Analyse_Genotype_Outcome(genotype,k_builds,&pt,seed).first;
      if(phen_size==255) 
        continue;
      if(phen_size>0)
        all_zero_phen=false;
    }
    if(all_zero_phen)
      continue;
#pragma omp atomic update
    ++good_genotypes;
#pragma omp critical
    {
      if(file_of_genotypes) {
	for(auto g : genotype)
	  fout<<+g<<" ";
      }
      else
	fout<<genotype_to_index(genotype,n_genes,colours)<<" ";
      fout<<"\n";
    }
  }
}

void ExhaustiveMinimalGenotypes(const char* file_path_c, uint8_t n_genes, uint8_t colours, bool file_of_genotypes) {
  std::string file_path(file_path_c),file_details="_N"+std::to_string(n_genes)+"_C"+std::to_string(colours)+".txt";
  std::ofstream fout(file_path+"ExhaustiveGenotypes"+file_details);
  
  GenotypeGenerator ggenerator = GenotypeGenerator(n_genes,colours);
  ggenerator.init();
  StochasticPhenotypeTable pt;
  uint8_t k_builds=5;
  Genotype genotype,nullg;

  while((genotype=ggenerator())!=nullg) {
    bool all_zero_phen=true;
    for(uint8_t seed=0;seed<n_genes;++seed) {
      uint8_t phen_size = Stochastic::Analyse_Genotype_Outcome(genotype,k_builds,&pt,seed).first;
      if(phen_size==255) 
        continue;
      if(phen_size>0)
        all_zero_phen=false;
    }
    if(all_zero_phen)
      continue;
    if(file_of_genotypes)
      for(auto g : genotype)
        fout<<+g<<" ";
    else
      fout<<genotype_to_index(genotype,n_genes,colours)<<" ";
    fout<<"\n";
  }
}
