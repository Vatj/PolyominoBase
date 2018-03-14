#include "genotype_phenotype.hpp"
#include <sstream>
#include <iterator>
#include <functional>


extern "C" void SampleMinimalGenotypes(uint8_t n_genes, uint8_t colours,const uint32_t N_SAMPLES,bool allow_duplicates) {
  uint32_t good_genotypes=0;
 

  std::ofstream fout("temp_wri.txt");
  
  GenotypeGenerator ggenerator = GenotypeGenerator(n_genes,colours);
  ggenerator.init();
  StochasticPhenotypeTable pt;
  uint8_t k_builds=3;

  std::mt19937 RNG_Engine(std::random_device{}());
  //#pragma omp threadprivate(RNG_Engine)

#pragma omp parallel //firstprivate(RNG_Engine)
  while(true) {
    uint32_t local_counter;
     std::vector<uint8_t> genotype;
     std::vector<uint32_t> states;
#pragma omp atomic read
    local_counter = good_genotypes;
    if(local_counter >= N_SAMPLES)
      break;
    
    //std::cout<<"samplign"<<std::endl;
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
    for(auto g : genotype)
      fout<<+g<<" ";
    fout<<"\n";
    }
  }
}

extern "C" void SampleMinimalGenotypes2(uint8_t n_genes, uint8_t colours,const uint32_t N_SAMPLES,bool allow_duplicates) {
  uint32_t good_genotypes=0;
 

  std::ofstream fout("temp_wri2.txt");

  StochasticPhenotypeTable pt;
  uint8_t k_builds=3;
  GenotypeGenerator ggenerator = GenotypeGenerator(n_genes,colours);

  std::mt19937 RNG_Engine(std::random_device{}());

  std::uniform_int_distribution<int> dist(0, colours-1);
  auto gen = std::bind(dist, std::ref(RNG_Engine));

#pragma omp parallel //firstprivate(RNG_Engine)
  while(true) {
    uint32_t local_counter;
    std::vector<uint8_t> genotype(n_genes*4);
     std::generate(genotype.begin(), genotype.end(), gen);
     std::vector<uint32_t> states;
#pragma omp atomic read
    local_counter = good_genotypes;
    if(local_counter >= N_SAMPLES)
      break;
    
       
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
    for(auto g : genotype)
      fout<<+g<<" ";
    fout<<"\n";
    }
  }
}

extern "C" void GGenerator(const char* a,bool file_of_genotypes,uint8_t colours,uint8_t n_genes) {
  
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


extern "C" void WrappedGetPhenotypesID(const char* a,const char* b,bool file_of_genotypes,uint8_t colours,uint8_t n_genes) {
  std::string file_path(b),filename(a),str,details="_N"+std::to_string(n_genes)+"_C"+std::to_string(colours)+".txt";
  std::ifstream file(filename);
  std::ofstream fout(file_path+"Genotype_Codes"+details, std::ios_base::out);
  std::ofstream fout2(file_path+"Phenotype_Table_N"+details, std::ios_base::out);

  StochasticPhenotypeTable pt;
  uint8_t k_builds=10;
  std::vector<uint8_t> genotype(n_genes*4);
  
  while (std::getline(file, str)) {
    if(file_of_genotypes) {
      std::istringstream is( str );
      genotype.assign( std::istream_iterator<uint8_t>( is ), std::istream_iterator<uint8_t>() );
    }
    else {
      index_to_genotype(std::stoull(str,nullptr),genotype.data(),colours,n_genes);
    }
    Clean_Genome(genotype,0,false);
    std::map<uint8_t,uint8_t> dups=DuplicateGenes(genotype);
    std::map<uint8_t,Phenotype_ID> seed_m;
    for(uint8_t seed=0;seed<n_genes;++seed)
      seed_m[seed]=Stochastic::Analyse_Genotype_Outcome(genotype,k_builds,&pt,seed);
    uint8_t index=0;
    Phenotype_ID phen_id;
  
    for(uint8_t gene=0;gene<n_genes;++gene) {
      phen_id = dups.count(gene) ? seed_m[dups[gene]] : phen_id=seed_m[index++];
      fout<<+phen_id.first<<" "<<+phen_id.second<<" ";
    }
    fout<<"\n";
  }
  pt.PrintTable(fout2);

}

int main(int argc, char* argv[]) {
  if(argc>2) {
    std::cout<<argv[1][0]<<std::endl;
    if(argv[1][0]=='1') {
      std::cout<<"1"<<std::endl;
      SampleMinimalGenotypes(3,7,std::stoi(argv[2]),true);
    }
    else {
      std::cout<<"2"<<std::endl;
      SampleMinimalGenotypes2(3,7,std::stoi(argv[2]),true);
    }
  }
  /*
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
  */
  

  return 0;


}

uint64_t genotype_to_index(uint8_t* genotype, uint8_t colours, uint8_t n_genes) {
  uint64_t count=0;
  for(uint8_t index=0;index<n_genes*4;++index)
    count+= genotype[index] * pow(colours,n_genes*4-index-1);
  return count;
}

void index_to_genotype(uint64_t index, uint8_t* genotype, uint8_t colours, uint8_t n_genes) {
  for(uint8_t count=0;count<n_genes*4;++count) {
    uint64_t value=index/pow(colours,n_genes*4-count-1);
    genotype[count]=value;
    index-= value * pow(colours,n_genes*4-count-1);
  }
}
