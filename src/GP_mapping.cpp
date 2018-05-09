#include "GP_mapping.hpp"



void GP_MapSampler(const char* file_path_c,uint8_t n_genes, uint8_t rcolours,uint8_t colours, bool file_of_genotypes) {
  
  const uint8_t k_builds=10;
  const uint32_t N_JIGGLE=100;
  StochasticPhenotypeTable pt;
  Genotype genotype(n_genes*4);

  std::string file_path(file_path_c),str,file_details="_N"+std::to_string(n_genes)+"_C"+std::to_string(rcolours);
  std::ifstream file_in(file_path+"SampledGenotypes"+file_details+"_Processed.txt");
  std::ofstream fout(file_path+"GP_Map_Cx"+std::to_string(colours)+file_details+".txt", std::ios_base::out);
  std::set<Phenotype_ID> evolvable_union;
  std::vector<double> robustnesses,deaths;
  std::vector<uint64_t> neutral_weightings;
  Phenotype_ID pid_section;
  
  while (std::getline(file_in, str)) {
    if(str=="x") {
      //next phenotype, dump everything

      double total_neutral_size=std::accumulate(neutral_weightings.begin(),neutral_weightings.end(),uint64_t(0));
      double avg_R=0,avg_D=0;
      for(size_t ind=0;ind<robustnesses.size();++ind) {
	avg_R+=robustnesses[ind]*neutral_weightings[ind]/total_neutral_size;
	avg_D+=deaths[ind]*neutral_weightings[ind]/total_neutral_size;
      }
      
      
      fout<<+pid_section.first<<" "<<pid_section.second<<": ";
      fout<<avg_R<<" "<<avg_D<<" "<<evolvable_union.size()<<"\n";


      //clear everything
      evolvable_union.clear();
      robustnesses.clear();
      deaths.clear();
      neutral_weightings.clear();
      continue;
    }
    if(file_of_genotypes) {
      std::istringstream is( str );
      genotype.assign( std::istream_iterator<int>( is ), std::istream_iterator<int>() );
    }
    else 
      index_to_genotype(std::stoull(str,nullptr),genotype,n_genes,colours);
    
    std::vector<Phenotype_ID> genotype_pIDs=GetPhenotypeIDs(genotype,k_builds,&pt);
    pid_section=genotype_pIDs.front();
#pragma omp parallel for schedule(dynamic) firstprivate(genotype)
    for(uint32_t nth_jiggle=0; nth_jiggle<N_JIGGLE;++nth_jiggle) {
      Clean_Genome(genotype,0,false);
      JiggleGenotype(genotype,colours);
      uint32_t robust=0,death=0;//,evolve=0;
      for(Genotype neighbour : genotype_neighbourhood(genotype,n_genes,colours)) {
	std::vector<Phenotype_ID> neighbour_pIDs=GetPhenotypeIDs(neighbour,k_builds,&pt);
	if(neighbour_pIDs.front()!=genotype_pIDs.front()) {
	  if(neighbour_pIDs.front().first < 255 && neighbour_pIDs.front().first > 0) {
#pragma omp critical(set_insertion)
	    {
	      evolvable_union.insert(neighbour_pIDs.front());
	    }
	  }
	  else
	    ++death;  
	}
	else
	  ++robust;
      }
#pragma omp critical(vector_emplacing)
      {
      robustnesses.emplace_back(static_cast<double>(robust)/N_JIGGLE);
      deaths.emplace_back(static_cast<double>(death)/N_JIGGLE);
      neutral_weightings.emplace_back(NeutralSize(genotype,1,colours-1)); //DEFAULT SET TO 1 NEUTRAL
      }
    }
  }  
    
}

uint64_t nChoosek(uint8_t n, uint8_t k ) {
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;
    uint64_t result = n;
    for(uint8_t i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}
uint64_t combination_with_repetiton(uint8_t space_size , uint8_t sample_size) {
  if(sample_size==0)
    return 1;
  std::vector<uint8_t> v(sample_size,1); 
  uint64_t comb_sum=0;
  while (true){ 
    for (uint8_t i = 0; i < sample_size; ++i){               
      if (v[i] > space_size){ 
        v[i + 1] += 1; 
        for (int16_t k = i; k >= 0; --k)
          v[k] = v[i + 1]; 
        v[i] = v[i + 1]; 
      } 
    } 
    if (v[sample_size] > 0) 
      break;
    uint64_t comb_prod=1;
    for(auto x: v)
      comb_prod*=x;
    comb_sum+=comb_prod;
    v[0] += 1; 
  } 
  return comb_sum; 
}


uint64_t NeutralSize(Genotype genotype,uint32_t N_neutral_colours,uint32_t N_possible_interacting_colours) {
  uint8_t neutral_faces=std::count(genotype.begin(),genotype.end(),0);
  Clean_Genome(genotype,0,false);
  std::set<uint8_t> unique_cols(genotype.begin(),genotype.end());
  uint32_t N_interacting_colours=N_possible_interacting_colours,N_interacting_pairs=(unique_cols.size()-1)/2;
  uint64_t neutral_interacting=1;
  for(uint8_t n=0;n<N_interacting_pairs;++n)
    neutral_interacting*=(N_interacting_colours-2*n);
  
  uint32_t N_noninteracting_colours=N_possible_interacting_colours-(unique_cols.size()-1);
  uint64_t neutral_noninteracting=0;
  for(uint8_t f=0;f<=neutral_faces;++f) {
    uint64_t pre_sum=nChoosek(neutral_faces,f)*pow(N_neutral_colours,neutral_faces-f);

    uint64_t sum_term=0;
    for(uint8_t U=0;U<=f;++U) {
      uint64_t pre_prod=1;
      for(uint8_t Un=0;Un<U;++Un)
        pre_prod*=(N_noninteracting_colours-2*Un);
      sum_term+=pre_prod*combination_with_repetiton(U,f-U);
    }
    neutral_noninteracting+=pre_sum*sum_term;   
  }
  return neutral_noninteracting*neutral_interacting;
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
  std::vector<Genotype> neighbours;
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
