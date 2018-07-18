#include "genotype_generate.hpp"
#include <sstream>
#include <iterator>
#include <functional>
#include <set>
#include <algorithm>

std::vector<Genotype> SampleMinimalGenotypes(PhenotypeTable* pt)
{
  uint64_t good_genotypes=0, generated_genotypes=0;
  std::vector<Genotype> genomes(simulation_params::n_samples + 4);
  Phenotype_ID rare_pID = {0, 0}, unbound_pID = {255, 0};
  std::vector<Phenotype_ID> pIDs;

  std::cout << "Generating " <<+ simulation_params::n_samples << " samples \n";

  GenotypeGenerator ggenerator = GenotypeGenerator(simulation_params::n_genes, simulation_params::colours);
  ggenerator.init();

  std::cout << "Threshold is : " << (ceil(simulation_params::phenotype_builds * simulation_params::UND_threshold));
  std::cout << " out of " <<+ simulation_params::phenotype_builds << " builds \n";

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

    if(local_counter >= simulation_params::n_samples)
      break;

#pragma omp atomic update
  ++generated_genotypes;

    uint32_t lbound_neck=0;
    for(uint8_t neck_ind=0; neck_ind < simulation_params::n_genes; ++neck_ind)
    {
      lbound_neck = std::uniform_int_distribution<uint32_t>{lbound_neck, ggenerator.n_necklaces-1}(simulation_params::RNG_Engine);
      states.emplace_back(lbound_neck);
      if(!simulation_params::allow_duplicates)
      {
        if((lbound_neck + simulation_params::n_genes - neck_ind) > ggenerator.n_necklaces)
          continue;
        ++lbound_neck;
      }
    }

    for(auto state : states)
      genotype.insert(genotype.end(), ggenerator.necklaces[state].begin(), ggenerator.necklaces[state].end());

    if(!ggenerator.valid_genotype(genotype))
      continue;

    pIDs = GetSetPIDs(genotype, pt);
    if(pIDs.front() == rare_pID || pIDs.back() == unbound_pID)
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

  while(genomes.size() > simulation_params::n_samples)
    genomes.pop_back();

  std::cout << "Final Values : Found " <<+ genomes.size() << " out of " <<+ generated_genotypes << " generated \n";

  return genomes;
}

std::vector<Genotype> ExhaustiveMinimalGenotypesIL(PhenotypeTable* pt)
{
  std::vector<Genotype> genomes;

  std::cout << "Generating all minimal samples \n";

  GenotypeGenerator ggenerator = GenotypeGenerator(simulation_params::n_genes, simulation_params::colours);
  ggenerator.init();
  Genotype genotype, nullg;
  Phenotype_ID unbound_pID = {255, 0};
  std::vector<Phenotype_ID> pIDs;
  uint64_t good_genotypes = 0, generated_genotypes = 0;

  std::cout << "Threshold is : " << (ceil(simulation_params::phenotype_builds * simulation_params::UND_threshold));
  std::cout << " out of " <<+ simulation_params::phenotype_builds << " builds \n";

  while((genotype=ggenerator())!=nullg)
  {
    generated_genotypes++;
    if(!ggenerator.valid_genotype(genotype))
      continue;

    pIDs = GetSetPIDs(genotype, pt);
    if(pIDs.back() == unbound_pID)
      continue;

    good_genotypes++;
    genomes.emplace_back(genotype);

    if(good_genotypes % 100 == 0)
      std::cout << "Found " <<+ good_genotypes << " out of " <<+ generated_genotypes << " generated \n";
  }

  std::cout << "Final Values : Found " <<+ genomes.size() << " out of " <<+ generated_genotypes << " generated \n";
  return genomes;
}

std::vector<Genotype> ExhaustiveMinimalGenotypesFiltered(PhenotypeTable* pt)
{
  std::vector<Genotype> genomes;

  std::cout << "Generating all minimal samples\n";

  GenotypeGenerator ggenerator = GenotypeGenerator(simulation_params::n_genes, simulation_params::colours);
  ggenerator.init();
  Genotype genotype, nullg;
  Phenotype_ID rare_pID = {0, 0}, unbound_pID = {255, 0};
  std::vector<Phenotype_ID> pIDs;
  uint64_t good_genotypes = 0, generated_genotypes = 0;

  std::cout << "Threshold is : " << (ceil(simulation_params::phenotype_builds * simulation_params::UND_threshold));
  std::cout << " out of " <<+ simulation_params::phenotype_builds << " builds \n";

  while((genotype=ggenerator())!=nullg)
  {
    generated_genotypes++;
    if(!ggenerator.valid_genotype(genotype))
      continue;

    pIDs = GetSetPIDs(genotype, pt);
    if(pIDs.front() == rare_pID || pIDs.back() == unbound_pID)
      continue;

    good_genotypes++;
    genomes.emplace_back(genotype);

    if(good_genotypes % 100 == 0)
      std::cout << "Found " <<+ good_genotypes << " out of " <<+ generated_genotypes << " generated \n";
  }

  std::cout << "Final Values : Found " <<+ genomes.size() << " out of " <<+ generated_genotypes << " generated \n";
  return genomes;
}

// Temporary Fix
std::vector<Genotype> ExhaustiveMinimalGenotypesFilteredDuplicate(std::vector<Genotype>& genomes, PhenotypeTable* pt)
{
  std::vector<Genotype> duplicates, dups;

  std::cout << "Generating all minimal samples\n";

  GenotypeGenerator ggenerator = GenotypeGenerator(simulation_params::n_genes - 1, simulation_params::colours);
  ggenerator.init();
  Genotype genotype, nullg;
  Phenotype_ID rare_pID = {0, 0}, unbound_pID = {255, 0};
  std::vector<Phenotype_ID> pIDs;
  uint64_t good_genotypes = 0, generated_genotypes = 0;
  bool keep_original;

  std::cout << "Threshold is : " << (ceil(simulation_params::phenotype_builds * simulation_params::UND_threshold));
  std::cout << " out of " <<+ simulation_params::phenotype_builds << " builds \n";

  while((genotype=ggenerator())!=nullg)
  {
    generated_genotypes++;
    if(!ggenerator.valid_genotype(genotype))
      continue;

    dups = GeneDuplication(genotype);
    keep_original = false;
    for(auto genome: dups)
    {
      pIDs = GetSetPIDs(genome, pt);
      if(pIDs.front() == rare_pID || pIDs.back() == unbound_pID)
        continue;
      good_genotypes++;
      duplicates.emplace_back(genome);
      keep_original = true;
      if(good_genotypes % 100 == 0)
        std::cout << "Found " <<+ good_genotypes << " out of " <<+ generated_genotypes << " generated \n";
    }
    if(keep_original)
      genomes.emplace_back(genotype);
  }

  std::cout << "Final Values : Found " <<+ genomes.size() << " out of " <<+ generated_genotypes << " generated \n";
  return duplicates;
}


// std::vector<Genotype> ExhaustiveFullGenotypes2(uint8_t colours, InterfacePhenotypeTable* pt)
// {
//   std::vector<Genotype> genomes;
//   //
//   // // GenotypeGenerator ggenerator = GenotypeGenerator(n_genes, colours);
//   // // ggenerator.init();
//   // // Genotype genotype(2 * 4);
//   // Phenotype_ID rare_pID = {0, 0}, unbound_pID = {255, 0};
//   // std::vector<Phenotype_ID> pIDs;
//   // uint64_t good_genotypes = 0, generated_genotypes = 0;
//
//   // for(int i=0;i<colours;++i) {
//   //   for(int j=0;j<colours;++j) {
//   //     for(int k=0;k<colours;++k) {
//   //       for(int l=0;l<colours;++l) {
//   //         for(int q=0;q<colours;++q) {
//   //           for(int w=0;w<colours;++w) {
//   //             for(int e=0;e<colours;++e) {
//   //               for(int r=0;r<colours;++r) {
//   //
//   //                 Genotype genotype = {i,j,k,l,q,w,e,r};
//   //                 Genotype copy = {i,j,k,l,q,w,e,r};
//   //                 generated_genotypes++;
//   //
//   //                 pIDs = GetSetPIDs(genotype, pt);
//   //
//   //                 if (pIDs.front() != rare_pID && pIDs.back() != unbound_pID)
//   //                 {
//   //                   good_genotypes++;
//   //                   genomes.emplace_back(copy);
//   //
//   //                   if(good_genotypes % 10000 == 0)
//   //                     std::cout << "Found " <<+ good_genotypes << " out of " <<+ generated_genotypes << " generated \n";
//   //                 }
//   //               }
//   //             }
//   //           }
//   //         }
//   //       }
//   //     }
//   //   }
//   // }
//
//   // std::cout << "Final Values : Found " <<+ genomes.size() << " out of " <<+ generated_genotypes << " generated \n";
//   return genomes;
// }




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
