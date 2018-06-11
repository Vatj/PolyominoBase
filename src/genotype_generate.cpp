#include "genotype_generate.hpp"
#include <sstream>
#include <iterator>
#include <functional>
#include <set>
#include <algorithm>

std::mt19937 RNG_Engine(std::random_device{}());

void SampleMinimalGenotypes(const char* file_path_c, uint8_t n_genes, uint8_t colours,const uint32_t N_SAMPLES,bool allow_duplicates)
{
  uint64_t good_genotypes=0, generated_genotypes=0;
  std::string file_path(file_path_c),file_details="_N"+std::to_string(n_genes)+"_C"+std::to_string(colours)+".txt";
  std::ofstream fout(file_path + "SampledGenotypes" + file_details);
  std::ofstream logut(file_path+"log"+file_details);

  GenotypeGenerator ggenerator = GenotypeGenerator(n_genes, colours);
  ggenerator.init();
  StochasticPhenotypeTable pt;
  uint8_t k_builds = 20;
  std::vector<Phenotype_ID> loop_pID = {{255, 0}};
  std::vector<Phenotype_ID> pIDs;
  Genotype nullg;

  std::ifstream pheno_in(file_path + "PhenotypeTable" + file_details);
  // if(pheno_in)
  //   LoadExistingTable(pheno_in, &pt);


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

    if(!ggenerator.valid_genotype(genotype))
      continue;

    // pIDs = GetSetPIDs(genotype, k_builds, &pt);
    if(GetSetPIDs(genotype, k_builds, &pt) == loop_pID)
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

std::vector<Genotype> SampleMinimalGenotypes(uint8_t n_genes, uint8_t colours,
  const uint32_t N_SAMPLES, bool allow_duplicates, StochasticPhenotypeTable* pt)
{
  uint64_t good_genotypes=0, generated_genotypes=0;
  std::vector<Genotype> genomes(N_SAMPLES + 4);

  std::cout << "Generating " <<+ N_SAMPLES << " samples \n";

  GenotypeGenerator ggenerator = GenotypeGenerator(n_genes, colours);
  ggenerator.init();

  uint8_t k_builds = 20;
  Phenotype_ID rare_pID = {0, 0}, loop_pID = {255, 0};
  std::vector<Phenotype_ID> pIDs;

#pragma omp parallel firstprivate(pIDs)
  while(true)
  {
    uint64_t local_counter, local_generated;
    Genotype genotype;
    std::vector<uint32_t> states;

#pragma omp atomic read
    local_counter = good_genotypes;
#pragma omp atomic read
    local_generated = generated_genotypes;

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
      genotype.insert(genotype.end(), ggenerator.necklaces[state].begin(), ggenerator.necklaces[state].end());

    if(!ggenerator.valid_genotype(genotype))
      continue;

    pIDs = GetSetPIDs(genotype, k_builds, pt);
    if(pIDs.front() == rare_pID || pIDs.back() == loop_pID)
      continue;

      #pragma omp critical
        genomes[good_genotypes] = genotype;

      #pragma omp atomic update
        ++good_genotypes;

      local_counter++;
      if(local_counter % 100 == 0)
        std::cout << "Found " <<+ local_counter << " out of " <<+ local_generated << " generated \n";
    }

  // std::sort(std::begin(genomes), std::end(genomes));
  // auto last = std::unique(std::begin(genomes), std::end(genomes));
  // genomes.erase(last, std::end(genomes));

  while(genomes.size() > N_SAMPLES)
    genomes.pop_back();

  std::cout << "Final Values : Found " <<+ genomes.size() << " out of " <<+ generated_genotypes << " generated \n";

  return genomes;
}


/****************/
/***** MAIN *****/
/****************/
// int main()
// {
//
//   StochasticPhenotypeTable pt;
//   uint8_t k_builds=5;
//   double UND_frac=.1;
//
//   for(auto g: Stochastic::AssemblePlasticGenotype({0,0,0,1, 2,1,0,0, 0,0,2,2, 2,2,2,2}, k_builds, &pt, UND_frac))
//   {
//       std::cout<<+g.first<<" "<<g.second<<std::endl;
//   }
//
//   return 0;
// }



std::vector<Genotype> ExhaustiveMinimalGenotypes(uint8_t n_genes, uint8_t colours, StochasticPhenotypeTable* pt)
{
  std::vector<Genotype> genomes;

  std::cout << "Generating all minimal samples\n";

  GenotypeGenerator ggenerator = GenotypeGenerator(n_genes, colours);
  ggenerator.init();
  uint8_t k_builds=100;
  Genotype genotype, nullg;
  Phenotype_ID rare_pID = {0, 0}, loop_pID = {255, 0};
  std::vector<Phenotype_ID> pIDs;
  uint64_t good_genotypes = 0, generated_genotypes = 0;

  while((genotype=ggenerator())!=nullg)
  {
    generated_genotypes++;
    if(!ggenerator.valid_genotype(genotype))
      continue;

    pIDs = GetSetPIDs(genotype, k_builds, pt);
    if(pIDs.front() == rare_pID || pIDs.back() == loop_pID)
      continue;

    good_genotypes++;
    genomes.emplace_back(genotype);

    if(good_genotypes % 100 == 0)
      std::cout << "Found " <<+ good_genotypes << " out of " <<+ generated_genotypes << " generated \n";
  }

  std::cout << "Final Values : Found " <<+ genomes.size() << " out of " <<+ generated_genotypes << " generated \n";
  return genomes;
}

std::vector<Genotype> ExhaustiveFullGenotypes2(uint8_t colours, StochasticPhenotypeTable* pt)
{
  std::vector<Genotype> genomes;
  //
  // // GenotypeGenerator ggenerator = GenotypeGenerator(n_genes, colours);
  // // ggenerator.init();
  // uint8_t k_builds=100;
  // // Genotype genotype(2 * 4);
  // Phenotype_ID rare_pID = {0, 0}, loop_pID = {255, 0};
  // std::vector<Phenotype_ID> pIDs;
  // uint64_t good_genotypes = 0, generated_genotypes = 0;

  // for(int i=0;i<colours;++i) {
  //   for(int j=0;j<colours;++j) {
  //     for(int k=0;k<colours;++k) {
  //       for(int l=0;l<colours;++l) {
  //         for(int q=0;q<colours;++q) {
  //           for(int w=0;w<colours;++w) {
  //             for(int e=0;e<colours;++e) {
  //               for(int r=0;r<colours;++r) {
  //
  //                 Genotype genotype = {i,j,k,l,q,w,e,r};
  //                 Genotype copy = {i,j,k,l,q,w,e,r};
  //                 generated_genotypes++;
  //
  //                 pIDs = GetSetPIDs(genotype, k_builds, pt);
  //
  //                 if (pIDs.front() != rare_pID && pIDs.back() != loop_pID)
  //                 {
  //                   good_genotypes++;
  //                   genomes.emplace_back(copy);
  //
  //                   if(good_genotypes % 10000 == 0)
  //                     std::cout << "Found " <<+ good_genotypes << " out of " <<+ generated_genotypes << " generated \n";
  //                 }
  //               }
  //             }
  //           }
  //         }
  //       }
  //     }
  //   }
  // }

  // std::cout << "Final Values : Found " <<+ genomes.size() << " out of " <<+ generated_genotypes << " generated \n";
  return genomes;
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
