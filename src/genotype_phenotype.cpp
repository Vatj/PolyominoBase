#include "genotype_phenotype.hpp"
#include <sstream>
#include <iterator>
#include <functional>
#include <set>
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


void GP_MapSampler(const char* file_path_c,uint8_t n_genes, uint8_t rcolours,uint8_t colours, bool file_of_genotypes)
{

  const uint8_t k_builds=10;
  const uint32_t N_JIGGLE=100;
  StochasticPhenotypeTable pt;
  Genotype genotype(n_genes*4);

  std::string file_path(file_path_c), str, file_details="_N"+std::to_string(n_genes)+"_C"+std::to_string(rcolours);
  std::ifstream file_in(file_path+"SampledGenotypes"+file_details+"_Processed.txt");
  std::ofstream fout(file_path+"GP_Map_Cx"+std::to_string(colours)+file_details+".txt", std::ios_base::out);

  std::string pmetric_str = "pmetrics", ending = ".txt";
  std::ofstream pmetric_out(file_path + pmetric_str + file_details + ending);

  std::string gmetric_str = "gmetrics";
  std::ofstream gmetric_out(file_path + gmetric_str + file_details + ending);

  Phenotype_Metrics pmetrics(n_genes, colours);
  // std::set<Phenotype_ID> evolvable_union;
  // std::vector<double> robustnesses, deaths;
  // std::vector<uint64_t> neutral_weightings;
  // Phenotype_ID pid_section;

  while (std::getline(file_in, str))
  {
    if(str=="x")
    {
      //next phenotype, dump everything

      // double total_neutral_size=std::accumulate(pmetrics.neutral_weightings.begin(), pmetrics.neutral_weightings.end(), uint64_t(0));
      // double avg_R=0, avg_D=0;
      //
      // for(size_t ind=0; ind<robustnesses.size(); ++ind)
      // {
	    //    avg_R += robustness[ind] * neutral_weightings[ind] / total_neutral_size;
	    //    avg_D += deaths[ind] * neutral_weightings[ind] / total_neutral_size;
      // }

      // //clear everything
      // evolvable_union.clear();
      // robustnesses.clear();
      // deaths.clear();
      // neutral_weightings.clear();
      pmetrics.save_to_file(pmetric_out);
      pmetrics.clear();

      continue;
    }
    if(file_of_genotypes)
    {
      std::istringstream is( str );
      genotype.assign( std::istream_iterator<int>( is ), std::istream_iterator<int>() );
    }
    else
      index_to_genotype(std::stoull(str,nullptr),genotype,n_genes,colours);

    std::vector<Phenotype_ID> genotype_pIDs = GetPhenotypeIDs(genotype,k_builds,&pt);
    pmetrics.pIDs = genotype_pIDs;

    // pid_section=genotype_pIDs.front();

    Genotype_Metrics gmetrics(n_genes);

#pragma omp parallel for schedule(dynamic) firstprivate(genotype, gmetrics, genotype_pIDs)
    for(uint32_t nth_jiggle=0; nth_jiggle<N_JIGGLE;++nth_jiggle)
    {
      Clean_Genome(genotype, 0, false);
      JiggleGenotype(genotype, colours);

      gmetrics.set_reference(genotype, genotype_pIDs);

      for(Genotype neighbour : genotype_neighbourhood(genotype, n_genes, colours))
      {
	       std::vector<Phenotype_ID> neighbour_pIDs = GetPhenotypeIDs(neighbour, k_builds, &pt);
         gmetrics.analyse_pIDs(neighbour_pIDs);
       }

      #pragma omp critical
      {
        gmetrics.save_to_file(gmetric_out);
        pmetrics.add_genotype_metrics(gmetrics);
      }

      gmetrics.clear();
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



  StochasticPhenotypeTable pt;
  uint8_t k_builds=5;
  double UND_frac=.1
  for(auto g: Stochastic::AssemblePlasticGenotype({0,0,0,1, 2,1,0,0, 0,0,2,2, 2,2,2,2}, k_builds, &pt,UND_frac)) {
      std::cout<<+g.first<<" "<<g.second<<std::endl;
    }
 
  return 0;

  if(argc>2) {
    std::cout<<argv[1][0]<<std::endl;

  }
  Genotype genotype{0,0,0,1,0,0,0,2};
  int x=0;
  for(Genotype neighbour : genotype_neighbourhood(genotype,2,std::stoi(argv[1]))) {
    ++x;
    for(auto q:neighbour)
      std::cout<<+q<<" ";
    std::cout<<std::endl;
  }
  std::cout<<x<<std::endl;
  return 0;


  Genotype geno,nullg;
  GenotypeGenerator ggenerator(std::stoi(argv[1]),std::stoi(argv[2]));
  ggenerator.init();
  while((geno=ggenerator())!=nullg) {
    for(auto x : geno)
      std::cout<<+x<<" ";
    std::cout<<std::endl;
  }



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

void GetPhenotypesIDs(const char* file_path_c,const char* file_name_c, uint8_t n_genes, uint8_t colours, bool file_of_genotypes) {
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

void PreProcessGenotypes(const char* file_path_c, uint8_t n_genes, uint8_t colours,bool file_of_genotypes) {
  std::string str,file_path(file_path_c),file_details="_N"+std::to_string(n_genes)+"_C"+std::to_string(colours);
  std::ifstream fin(file_path+"SampledGenotypes"+file_details+"_Iso.txt");
  std::ofstream fout(file_path+"SampledGenotypes"+file_details+"_Processed.txt");

  StochasticPhenotypeTable pt;
  uint8_t k_builds=5;
  std::map<Phenotype_ID,std::vector<Genotype> > phen_sets;
  Genotype genotype;

  while (std::getline(fin, str)) {
    if(file_of_genotypes) {
      std::istringstream is( str );
      genotype.assign( std::istream_iterator<int>( is ), std::istream_iterator<int>() );
    }
    const Genotype g_write(genotype);
    phen_sets[GetPhenotypeIDs(genotype, k_builds, &pt).front()].emplace_back(g_write);

  }

  for(auto kv : phen_sets) {
    for(Genotype& gen : kv.second) {
      for(uint8_t base : gen) {
	fout<<+base<<" ";
      }
      fout<<"\n";
    }
    fout<<"x\n";
  }
}
