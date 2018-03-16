#include "genotype_phenotype.hpp"
#include <sstream>
#include <iterator>
#include <functional>

std::mt19937 RNG_Engine(std::random_device{}());

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


void GP_MapSampler(const char* file_path_c,uint8_t n_genes, uint8_t rcolours,uint8_t colours, bool file_of_genotypes) {
  
  const uint8_t k_builds=10;
  const uint32_t N_JIGGLE=100;
  StochasticPhenotypeTable pt;
  Genotype genotype(n_genes*4);

  std::string file_path(file_path_c),str,file_details="_N"+std::to_string(n_genes)+"_C"+std::to_string(rcolours)+".txt";
  std::ifstream file_in(file_path+"SampledGenotypes"+file_details);
  std::ofstream fout(file_path+"GP_Map_Cx"+std::to_string(colours)+file_details, std::ios_base::out);
 
  while (std::getline(file_in, str)) {
    if(file_of_genotypes) {
      std::istringstream is( str );
      genotype.assign( std::istream_iterator<int>( is ), std::istream_iterator<int>() );
    }
    else 
      index_to_genotype(std::stoull(str,nullptr),genotype,n_genes,colours);
    std::vector<Phenotype_ID> genotype_pIDs=GetPhenotypeIDs(genotype,k_builds,&pt);
#pragma omp parallel for schedule(dynamic) firstprivate(genotype)
    for(uint32_t nth_jiggle=0; nth_jiggle<N_JIGGLE;++nth_jiggle) {
      Clean_Genome(genotype,0,false);
      JiggleGenotype(genotype,colours);
      uint32_t robust=0,evolve=0;
      for(Genotype neighbour : genotype_neighbourhood(genotype,n_genes,colours)) {
	std::vector<Phenotype_ID> neighbour_pIDs=GetPhenotypeIDs(neighbour,k_builds,&pt);
	if(neighbour_pIDs!=genotype_pIDs)
	  ++evolve;
	else
	  ++robust;
      }
#pragma omp critical(file_writing)
      {
	if(file_of_genotypes) 
	  for(auto g : genotype)
	    fout<<+g<<" ";
        else
	  fout<<genotype_to_index(genotype,n_genes,colours)<<" ";
	fout<<": R "<<robust<<" E "<<evolve<<"\n";
      }
    }
  }  
    
}

std::vector<Phenotype_ID> GetPhenotypeIDs(Genotype& genotype, uint8_t k_builds, StochasticPhenotypeTable* pt_it) {
  std::vector<Phenotype_ID> pIDs;
  Clean_Genome(genotype,0,false);
  const uint8_t n_genes=genotype.size()/4;
  std::map<uint8_t,uint8_t> dups=DuplicateGenes(genotype);
  std::map<uint8_t,Phenotype_ID> seed_m;
  for(uint8_t seed=0;seed<genotype.size()/4;++seed)
    seed_m[seed]=Stochastic::Analyse_Genotype_Outcome(genotype,k_builds,pt_it,seed);
  uint8_t index=0;
  Phenotype_ID phen_id;
  for(uint8_t gene=0;gene<n_genes; ++gene) {
    pIDs.emplace_back(dups.count(gene) ? seed_m[dups[gene]] : phen_id=seed_m[index++]);
  }
  std::sort(pIDs.begin(),pIDs.end());
  return pIDs;
}

/****************/
/***** MAIN *****/
/****************/
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
  /*
  
  Genotype geno,nullg;
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

void JiggleGenotype(Genotype& genotype, uint8_t max_colour) {
  uint8_t min_colour=*std::max_element(genotype.begin(),genotype.end());
  if(min_colour+1==max_colour)
    return;
  std::vector<uint8_t> neutral_colours(1+(max_colour-min_colour)/2);
  std::generate(neutral_colours.begin()+1, neutral_colours.end(), [n = min_colour-1] () mutable { return n+=2; });
  std::uniform_int_distribution<size_t> jiggle_index(0,neutral_colours.size()-1);

  for(auto& base : genotype)
    base= (base==0) ? neutral_colours[jiggle_index(RNG_Engine)] : base;
}
std::vector<Genotype> genotype_neighbourhood(const Genotype& genome, uint8_t ngenes, uint8_t colours) {
  std::vector<Genotype> neighbours(1,genome);
  Genotype neighbour;
  std::vector<uint8_t> mutants(colours);
  std::iota(mutants.begin(), mutants.end(), 0);
  for(uint8_t index=0; index<ngenes*4; ++index) {
    std::swap(*std::find(mutants.begin(), mutants.end(),genome[index]), mutants.back());
    neighbour = genome;
    for(int j=0; j<colours-1; ++j) {
      neighbour[index] = mutants[j];
      neighbours.emplace_back(neighbour);
    }
  }
  return neighbours;
}

void GetPhenotypeIDs(const char* file_path_c,const char* file_name_c, uint8_t n_genes, uint8_t colours, bool file_of_genotypes) {

  std::string file_path(file_path_c),file_name(file_name_c),str,details="_N"+std::to_string(n_genes)+"_C"+std::to_string(colours)+".txt";
  std::ifstream file_in(file_path+file_name);
  std::ofstream gfout(file_path+"Genotype_Codes"+details, std::ios_base::out);
  std::ofstream pfout(file_path+"Phenotype_Table"+details, std::ios_base::out);

  StochasticPhenotypeTable pt;
  uint8_t k_builds=10;
  Genotype genotype(n_genes*4);
  
  while (std::getline(file_in, str)) {
    if(file_of_genotypes) {
      std::istringstream is( str );
      genotype.assign( std::istream_iterator<int>( is ), std::istream_iterator<int>() );
    }
    else 
      index_to_genotype(std::stoull(str,nullptr),genotype,n_genes,colours);
    for(auto x: genotype)
      std::cout<<+x<<" ";
    std::cout<<std::endl;
    for(auto phen_id : GetPhenotypeIDs(genotype,k_builds,&pt)) 
      gfout<<+phen_id.first<<" "<<+phen_id.second<<" ";
    gfout<<"\n";
  }
  pt.PrintTable(pfout);

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

 
