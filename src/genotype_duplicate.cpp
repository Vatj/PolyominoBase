#include "genotype_duplicate.hpp"
#include <sstream>
#include <iostream>
#include <iterator>


std::vector<Genotype> GenomesDuplication(std::vector<Genotype> genomes)
{
  std::vector<Genotype> genomes_dup;
  std::cout << "Adding duplicate genes to " <<+ genomes.size() << " genomes \n";

  for(auto genome: genomes)
    for(auto duplicate: GeneDuplication(genome))
      genomes_dup.emplace_back(duplicate);

  return genomes_dup;
}

std::vector<Genotype> GeneDuplication(Genotype& genotype)
{
  std::vector <Genotype> duplicates;

  for(uint8_t index=0; index < simulation_params::n_genes; ++index)
  {
    Genotype duplicate(4 * simulation_params::n_genes);
    std::copy(std::begin(genotype), std::end(genotype), std::begin(duplicate));

    for(uint8_t tail=0; tail < 4; tail++)
      duplicate.emplace_back(genotype[(4 * index) + tail]);

    duplicates.emplace_back(duplicate);
  }

  return duplicates;
}

// void GenomesJiggleDuplication(std::vector<Genotype>& genomes, std::vector<Genotype>& jiggle_genomes, std::vector<Genotype>& duplicates)
// {
//   std::cout << "Adding duplicate genes to " <<+ genomes.size() << " genomes \n";
//
//   for(auto genome: genomes)
//   {
//     for(uint32_t nth_jiggle=0; nth_jiggle<simulation_params::n_jiggle; ++nth_jiggle)
//     {
//       // Clean_Genome(genome, false);
//       JiggleGenotype(genome);
//       jiggle_genomes.emplace_back(genome);
//
//       for(auto duplicate: GeneDuplication(genome))
//         duplicates.emplace_back(duplicate);
//     }
//   }
// }
