#include "genotype_generate.hpp"
#include <sstream>
#include <iterator>
#include <functional>
#include <set>

std::mt19937 RNG_Engine(std::random_device{}());

// void SampleMinimalGenotypes_old(const char* file_path_c, uint8_t n_genes, uint8_t colours,const uint32_t N_SAMPLES,bool allow_duplicates)
// {
//   uint32_t good_genotypes=0;
//   std::string file_path(file_path_c),file_details="_N"+std::to_string(n_genes)+"_C"+std::to_string(colours)+".txt";
//   std::ofstream fout(file_path+"SampledGenotypes"+file_details);
//
//   GenotypeGenerator ggenerator = GenotypeGenerator(n_genes, colours);
//   ggenerator.init();
//   StochasticPhenotypeTable pt;
//   uint8_t k_builds=5;
//
//
// #pragma omp parallel
//   while(true) {
//     uint32_t local_counter;
//      Genotype genotype;
//      std::vector<uint32_t> states;
// #pragma omp atomic read
//     local_counter = good_genotypes;
//     if(local_counter >= N_SAMPLES)
//       break;
//
//     uint32_t lbound_neck=0;
//     for(uint8_t neck_ind=0;neck_ind<n_genes;++neck_ind) {
//       lbound_neck=std::uniform_int_distribution<uint32_t>{lbound_neck,ggenerator.n_necklaces-1}(RNG_Engine);
//       states.emplace_back(lbound_neck);
//       if(!allow_duplicates) {
//         if((lbound_neck+n_genes-neck_ind)>ggenerator.n_necklaces)
//           continue;
//         ++lbound_neck;
//       }
//     }
//     for(auto state : states)
//       genotype.insert(genotype.end(),ggenerator.necklaces[state].begin(),ggenerator.necklaces[state].end());
//
//     if(!ggenerator.valid_genotype(genotype))
//       continue;
//
//     bool all_zero_phen=true;
//     for(uint8_t seed=0;seed<n_genes;++seed) {
//       uint8_t phen_size = Stochastic::Analyse_Genotype_Outcome(genotype,k_builds,&pt,seed).first;
//       if(phen_size==255)
//         continue;
//       if(phen_size>0)
//         all_zero_phen=false;
//     }
//     if(all_zero_phen)
//       continue;
// #pragma omp atomic update
//     ++good_genotypes;
// #pragma omp critical
//     {
// 	     for(auto g : genotype)
// 	      fout<<+g<<" ";
//        fout<<"\n";
//     }
//   }
// }

void SampleMinimalGenotypes(const char* file_path_c, uint8_t n_genes, uint8_t colours,const uint32_t N_SAMPLES,bool allow_duplicates)
{
  uint64_t good_genotypes=0, generated_genotypes=0;
  std::string file_path(file_path_c),file_details="_N"+std::to_string(n_genes)+"_C"+std::to_string(colours)+".txt";
  std::ofstream fout(file_path+"SampledGenotypes"+file_details);
  std::ofstream logut(file_path+"log"+file_details);

  GenotypeGenerator ggenerator = GenotypeGenerator(n_genes, colours);
  ggenerator.init();
  StochasticPhenotypeTable pt;
  uint8_t k_builds = 5;
  Phenotype_ID death_pID = {0, 0}, loop_pID = {255, 0};
  std::vector<Phenotype_ID> pIDs;
  Genotype nullg;

#pragma omp parallel firstprivate(pIDs)
  while(true)
  {
    uint64_t local_counter;
    Genotype genotype;
    std::vector<uint32_t> states;

#pragma omp atomic read
    local_counter = good_genotypes;

    if(local_counter >= N_SAMPLES)
      break;

#pragma omp atomic update
  ++generated_genotypes;

    uint32_t lbound_neck=0;
    for(uint8_t neck_ind=0; neck_ind < n_genes; ++neck_ind)
    {
      lbound_neck = std::uniform_int_distribution<uint32_t>{lbound_neck, ggenerator.n_necklaces-1}(RNG_Engine);
      states.emplace_back(lbound_neck);
      if(!allow_duplicates)
      {
        if((lbound_neck + n_genes - neck_ind) > ggenerator.n_necklaces)
          continue;
        ++lbound_neck;
      }
    }

    for(auto state : states)
      genotype.insert(genotype.end(),ggenerator.necklaces[state].begin(), ggenerator.necklaces[state].end());

    // if(!ggenerator.valid_genotype(genotype))
    //   continue;

    pIDs = GetSetPIDs(genotype, k_builds, &pt);

    if(pIDs.size() > 0)
      if ((pIDs.front() == death_pID) || (pIDs.back() == loop_pID))
        continue;

    #pragma omp atomic update
      ++good_genotypes;

    #pragma omp critical
    {
	     for(auto g : genotype)
	      fout<<+g<<" ";
       fout << "\n";

       logut <<+ good_genotypes << " " << generated_genotypes << "\n";
    }
  }
}



/****************/
/***** MAIN *****/
/****************/
int main()
{

  StochasticPhenotypeTable pt;
  uint8_t k_builds=5;
  double UND_frac=.1;

  for(auto g: Stochastic::AssemblePlasticGenotype({0,0,0,1, 2,1,0,0, 0,0,2,2, 2,2,2,2}, k_builds, &pt, UND_frac))
  {
      std::cout<<+g.first<<" "<<g.second<<std::endl;
  }

  return 0;
}



void ExhaustiveMinimalGenotypes(const char* file_path_c, uint8_t n_genes, uint8_t colours)
{
  std::string file_path(file_path_c), file_details="_N" + std::to_string(n_genes) + "_C" + std::to_string(colours)+".txt";
  std::ofstream fout(file_path + "ExhaustiveGenotypes" + file_details);

  GenotypeGenerator ggenerator = GenotypeGenerator(n_genes, colours);
  ggenerator.init();
  StochasticPhenotypeTable pt;
  uint8_t k_builds=5;
  Genotype genotype, nullg;


  while((genotype=ggenerator())!=nullg)
  {
    bool all_zero_phen=true;
    for(uint8_t seed=0;seed<n_genes;++seed)
    {
      uint8_t phen_size = Stochastic::Analyse_Genotype_Outcome(genotype,k_builds,&pt,seed).first;
      if(phen_size==255)
        continue;
      if(phen_size>0)
        all_zero_phen=false;
    }
    if(all_zero_phen)
      continue;
    for(auto g : genotype)
      fout<<+g<<" ";
    fout<<"\n";
  }
}

void PreProcessGenotypesTopology(const char* file_path_c, uint8_t n_genes, uint8_t colours)
{
  std::string str,file_path(file_path_c),file_details="_N"+std::to_string(n_genes)+"_C"+std::to_string(colours);
  std::ifstream fin(file_path+"SampledGenotypes"+file_details+"_Iso.txt");
  std::ofstream fout(file_path+"SampledGenotypes"+file_details+"_Processed.txt");

  StochasticPhenotypeTable pt;
  uint8_t k_builds=5;
  std::map<std::vector<Phenotype_ID>, std::vector<Genotype>> phen_sets;
  Genotype genotype;

  while (std::getline(fin, str))
  {
    std::istringstream is( str );
    genotype.assign( std::istream_iterator<int>( is ), std::istream_iterator<int>() );

    const Genotype g_write(genotype);
    std::vector<Phenotype_ID> pIDs = GetSetPIDs(genotype, k_builds, &pt);
    phen_sets[pIDs].emplace_back(g_write);
  }

  for(auto kv : phen_sets)
  {
    for(Genotype& gen : kv.second)
    {
      for(uint8_t base : gen)
      {
	       fout<<+base<<" ";
      }
      fout<<"\n";
    }

    fout << "x ";
    for(auto pID: kv.first)
      fout <<+ pID.first << " " <<+ pID.second << " ";
    fout << "\n";

  }
}

void PreProcessGenotypes(const char* file_path_c, uint8_t n_genes, uint8_t colours)
{
  std::string str,file_path(file_path_c),file_details="_N"+std::to_string(n_genes)+"_C"+std::to_string(colours);
  std::ifstream fin(file_path+"SampledGenotypes"+file_details+".txt");
  std::ofstream fout(file_path+"SampledGenotypes"+file_details+"_Processed.txt");

  StochasticPhenotypeTable pt;
  uint8_t k_builds=5;
  std::map<std::vector<Phenotype_ID>, std::vector<Genotype>> phen_sets;
  Genotype genotype;


  while (std::getline(fin, str))
  {
    std::istringstream is( str );
    genotype.assign( std::istream_iterator<int>( is ), std::istream_iterator<int>() );

    const Genotype g_write(genotype);
    std::vector<Phenotype_ID> pIDs = GetSetPIDs(genotype, k_builds, &pt);
    phen_sets[pIDs].emplace_back(g_write);
  }

  for(auto kv : phen_sets)
  {
    for(Genotype& gen : kv.second)
    {
      for(uint8_t base : gen)
      {
	       fout<<+base<<" ";
      }
      fout<<"\n";
    }

    fout << "x ";
    for(auto pID: kv.first)
      fout <<+ pID.first << " " <<+ pID.second << " ";
    fout << "\n";

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
