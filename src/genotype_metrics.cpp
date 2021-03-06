#include "genotype_metrics.hpp"
#include <sstream>
#include <iterator>
#include <functional>
#include <set>
#include <algorithm>
#include <cstring>

// Main Sampling function

void GP_MapSampler(Set_to_Genome& set_to_genome, PhenotypeTable* pt)
{
  Phenotype_ID unbound_pID = {255, 0}, rare_pID = {0, 0};
  double neutral_weight = 0;

  uint32_t number_of_genomes = 0;
  for(Set_to_Genome::iterator iter = std::begin(set_to_genome); iter != std::end(set_to_genome); iter++)
    number_of_genomes += (iter->second).size();
  std::cout << "There are " <<+ number_of_genomes << " genomes to jiggle a ";
  std::cout <<+ simulation_params::n_jiggle << "times! \n";

  // Create new files and open them for writing (erase previous data)
  std::ofstream set_metric_out(io_params::set_metric_file);
  std::ofstream genome_metric_out(io_params::genome_metric_file);
  header_metric_files(set_metric_out, genome_metric_out);

  for(Set_to_Genome::iterator iter = std::begin(set_to_genome); iter != std::end(set_to_genome); iter++)
  {
    std::cout << "Currently processing " <<+ (iter->second).size() << " genomes for {";
    for(auto pID: iter->first)
      std::cout << "(" <<+ pID.first << ", " <<+ pID.second << "), ";
    number_of_genomes -= (iter->second).size();
    std::cout << "}. Only " <<+ number_of_genomes << " left! \n";

    if(((iter->first).front() == rare_pID) || ((iter->first).back() == unbound_pID))
      continue;

    Set_Metrics set_metrics(simulation_params::n_genes, simulation_params::metric_colours);
    set_metrics.ref_pIDs = iter->first;

    for(auto genotype: iter->second)
    {
      neutral_weight = ((double) NeutralSize(genotype, 1, simulation_params::metric_colours - 1)) / simulation_params::n_jiggle; // This weight will be counted n_jiggle time when added to form the total neutral weight
      set_metrics.originals.emplace_back(genotype);

      #pragma omp parallel for schedule(dynamic) firstprivate(genotype, neutral_weight)
      for(uint32_t nth_jiggle=0; nth_jiggle<simulation_params::n_jiggle; ++nth_jiggle)
      {
        Clean_Genome(genotype, false);
        JiggleGenotype(genotype);

        Genotype_Metrics genome_metric(simulation_params::n_genes, simulation_params::metric_colours);
        genome_metric.set_reference(genotype, iter->first, neutral_weight);

        std::map<Phenotype_ID, uint8_t> pID_counter = GetPIDCounter(genotype, pt);
        genome_metric.pID_counter.insert(std::begin(pID_counter), std::end(pID_counter));

        std::vector<Phenotype_ID> pIDs;
        for(auto pID: pID_counter)
          pIDs.emplace_back(pID.first);

        if(pIDs != iter->first)
          set_metrics.misclassified[genotype] = pIDs;

        for(Genotype neighbour : genotype_neighbourhood(genotype))
        {
           std::vector<Phenotype_ID> neighbour_pIDs = GetSetPIDs(neighbour, pt);
           genome_metric.analyse_pIDs(neighbour_pIDs);
         }
         #pragma omp critical
         {
           set_metrics.add_genotype_metrics(genome_metric);
         }
      }
    }
    // metrics.emplace_back(set_metrics);
    set_metrics.save_to_file(set_metric_out, genome_metric_out);
  }
}

// Subroutine of the GP_MapSampler

void JiggleGenotype(Genotype& genotype)
{
  uint8_t max_colour = simulation_params::metric_colours;
  uint8_t min_colour =* std::max_element(genotype.begin(), genotype.end());

  if(min_colour + 1 == max_colour)
    return;

  std::vector<uint8_t> neutral_colours( 1+ (max_colour - min_colour) / 2);

  std::generate(neutral_colours.begin() + 1, neutral_colours.end(), [n = min_colour-1] () mutable { return n+=2; });

  std::uniform_int_distribution<size_t> jiggle_index(0, neutral_colours.size() - 1);

  std::map <uint8_t, uint8_t> dups;

  if(simulation_params::dup_aware)
    dups = DuplicateGenes(genotype);

  for(auto& base : genotype)
    base= (base==0) ? neutral_colours[jiggle_index(simulation_params::RNG_Engine)] : base;

  if(simulation_params::dup_aware)
  {
    for(auto dup: dups)
    {
      Genotype::const_iterator first = genotype.begin() + (4 * dup.second);
      Genotype::const_iterator last = genotype.begin() + (4 * (dup.second + 1));
      std::copy(first, last, std::back_inserter(genotype));
    }
  }
}

std::vector<Genotype> genotype_neighbourhood(const Genotype& genome)
{
  std::vector<Genotype> neighbours;
  Genotype neighbour;

  std::vector<uint8_t> mutants(simulation_params::metric_colours);
  std::iota(mutants.begin(), mutants.end(), 0);

  for(uint8_t index=0; index<simulation_params::n_genes*4; ++index)
  {
    std::swap(*std::find(mutants.begin(), mutants.end(),genome[index]), mutants.back());

    neighbour = genome;
    for(int j=0; j<simulation_params::metric_colours-1; ++j)
    {
      neighbour[index] = mutants[j];
      neighbours.emplace_back(neighbour);
    }
  }
  return neighbours;
}

uint64_t nChoosek(uint8_t n, uint8_t k)
{
    if (k > n) return 0;
    if (k * 2 > n) k = n - k;
    if (k == 0) return 1;

    uint64_t result = n;

    for(uint8_t i = 2; i <= k; ++i)
    {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

uint64_t combination_with_repetiton(uint8_t space_size , uint8_t sample_size)
{
  if(sample_size==0)
    return 1;
  std::vector<uint8_t> v(sample_size,1);
  uint64_t comb_sum=0;
  while (true)
  {
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

uint64_t NeutralSize(Genotype genotype, uint32_t N_neutral_colours, uint32_t N_possible_interacting_colours)
{
  uint8_t neutral_faces = std::count(genotype.begin(),genotype.end(),0);
  Clean_Genome(genotype, false);
  std::set<uint8_t> unique_cols(genotype.begin(),genotype.end());

  uint32_t N_interacting_colours= N_possible_interacting_colours, N_interacting_pairs = (unique_cols.size() - 1) / 2;
  uint64_t neutral_interacting=1;

  for(uint8_t n=0; n<N_interacting_pairs; ++n)
    neutral_interacting *= (N_interacting_colours - (2 * n));

  uint32_t N_noninteracting_colours = N_possible_interacting_colours - (unique_cols.size() - 1);
  uint64_t neutral_noninteracting = 0;

  for(uint8_t f=0; f <= neutral_faces; ++f)
  {
    uint64_t pre_sum = nChoosek(neutral_faces, f) * pow(N_neutral_colours, neutral_faces-f);
    uint64_t sum_term=0;

    for(uint8_t U = 0; U <= f; ++U)
    {
      uint64_t pre_prod=1;

      for(uint8_t Un=0;Un<U;++Un)
        pre_prod *= (N_noninteracting_colours - (2 * Un));

      sum_term += pre_prod*combination_with_repetiton(U, f - U);
    }
    neutral_noninteracting += pre_sum * sum_term;
  }
  return neutral_noninteracting*neutral_interacting;
}
