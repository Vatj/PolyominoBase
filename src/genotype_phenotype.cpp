#include "genotype_phenotype.hpp"
#include <sstream>
#include <iterator>
#include <functional>
#include <set>

std::mt19937 RNG_Engine(std::random_device{}());

void SampleMinimalGenotypes(const char* file_path_c, uint8_t n_genes, uint8_t colours,const uint32_t N_SAMPLES,bool allow_duplicates, bool file_of_genotypes)
{
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

std::vector<Phenotype_ID> GetPhenotypeIDs(Genotype& genotype, uint8_t k_builds, StochasticPhenotypeTable* pt_it)
{
  std::vector<Phenotype_ID> pIDs;
  Clean_Genome(genotype, 0,false);
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
int main()
{

  StochasticPhenotypeTable pt;
  uint8_t k_builds=5;
  double UND_frac=.1;

  for(auto g: Stochastic::AssemblePlasticGenotype({0,0,0,1, 2,1,0,0, 0,0,2,2, 2,2,2,2}, k_builds, &pt,UND_frac))
  {
      std::cout<<+g.first<<" "<<g.second<<std::endl;
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

void GetPhenotypesIDs(const char* file_path_c,const char* file_name_c, uint8_t n_genes, uint8_t colours, bool file_of_genotypes)
{
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

// Member functions of the NecklaceFactory structure
void NecklaceFactory::GenNecklaces(uint8_t c)
{
  colours=c;
  crsms_gen(1,1);
}


bool NecklaceFactory::is_finite_necklace(std::vector<uint8_t>& neck)
{
  //internal infinite loop
  if((neck[1] && Interaction_Matrix(neck[1])==neck[3]) || (neck[2] && Interaction_Matrix(neck[2])==neck[4]))
    return false;
  bool a_pair=false;
  for(uint8_t base=1; base<4;++base) {
    if(neck[base]==0)
continue;
    uint8_t b1=std::count(neck.begin()+1,neck.end(),neck[base]), b2=std::count(neck.begin()+1,neck.end(),Interaction_Matrix(neck[base]));
    if(b1 && b2) {
//internal branching point/degenerate double loop
if(b1+b2>2)
  return false;
else {
  //internal unique double loop
  if(a_pair)
    return false;
  else
    a_pair=true;
}
    }
  }
  return true;
}

void NecklaceFactory::is_necklace(uint64_t j)
{
  if(4%j==0)
    if(is_finite_necklace(necklace_grower))
      necklaces.emplace_back(std::vector<uint8_t>{necklace_grower.begin()+1,necklace_grower.end()});
}

void NecklaceFactory::crsms_gen(uint64_t n, uint64_t j)
{
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
