#include "genotype_metrics.hpp"
#include <sstream>
#include <iterator>
#include <functional>
#include <set>
#include <algorithm>
#include <cstring>

std::mt19937 RNG_Engine1(std::random_device{}());

// Constructor of the Genotype_Metrics structure
Genotype_Metrics::Genotype_Metrics(uint8_t n_genes, uint8_t colours):
n_genes(n_genes), colours(colours)
{
  number_of_neighbours = (colours - 1) * n_genes * 4.;
}

void Genotype_Metrics::set_reference(Genotype& genotype, std::vector <Phenotype_ID> pIDs, double neutral)
{
  ref_genotype = genotype;
  ref_pIDs = pIDs;

  neutral_weight = neutral;

  // for (auto pID: pIDs)
  //   shapes.emplace_back(Shape_Metrics(pID));
}

void Genotype_Metrics::analyse_pIDs(std::vector <Phenotype_ID>& pIDs)
{
  if (pIDs == ref_pIDs)
    strict_robustness+= 1.;

  std::vector <Phenotype_ID> intersection, union_set;

  std::set_intersection(std::begin(pIDs), std::end(pIDs), std::begin(ref_pIDs), std::end(ref_pIDs), std::back_inserter(intersection));
  std::set_union(std::begin(pIDs), std::end(pIDs), std::begin(ref_pIDs), std::end(ref_pIDs), std::back_inserter(union_set));

  if(ref_pIDs.size() > 0)
    intersection_robustness += (double) intersection.size() / (double) ref_pIDs.size();
  else
    intersection_robustness = 0;

  union_evolvability += (double) (union_set.size() - ref_pIDs.size());

  if (std::find(std::begin(pIDs), std::end(pIDs), rare_pID) != std::end(pIDs))
    rare += 1.;

  if (std::find(std::begin(pIDs), std::end(pIDs), loop_pID) != std::end(pIDs))
    loop += 1.;

  // for (auto shape: shapes)
  //   shape.robust_pID(pIDs);

  for (auto pID: pIDs)
    diversity.insert(pID);
}

void Genotype_Metrics::save_to_file(std::ofstream& fout)
{
  fout << "(";
  for (auto face: ref_genotype)
    fout <<+ face << ",";
  fout.seekp((long) fout.tellp() - 1);
  fout << ") ";

  fout <<+ strict_robustness / number_of_neighbours << " ";
  fout <<+ intersection_robustness / number_of_neighbours << " ";
  fout <<+ union_evolvability / number_of_neighbours << " ";
  fout <<+ rare / number_of_neighbours << " ";
  fout <<+ loop / number_of_neighbours << " ";
  fout <<+ diversity.size() << " ";
  fout <<+ neutral_weight << " ";

  if(ref_pIDs.size() > 0)
  {
    fout << "{";
    for (auto pID: ref_pIDs)
      fout <<+ "(" <<+ pID.first << "," <<+ pID.second << "),";

    fout.seekp((long) fout.tellp() - 1);
    fout << "}\n";
  }
  else
    fout << "{}\n";
}

void Genotype_Metrics::clear()
{
  strict_robustness = 0, union_evolvability = 0, rare = 0;
  intersection_robustness = 0, loop = 0;

  // shapes.clear();
  diversity.clear();
}

Shape_Metrics::Shape_Metrics(Phenotype_ID pID): pID(pID), robustness(0)
{}

void Shape_Metrics::robust_pID(std::vector <Phenotype_ID> pIDs)
{
  auto presence = std::find(std::begin(pIDs), std::end(pIDs), pID);

  if (presence != std::end(pIDs))
    robustness++;
}


// Constructor of the Set_Metrics structure
Set_Metrics::Set_Metrics(uint8_t n_genes, uint8_t colours):
n_genes(n_genes), colours(colours), analysed(0)
{}

// Member functions of the Set_Metrics structure

void Set_Metrics::add_genotype_metrics(Genotype_Metrics& genome_metric)
{
  double number_of_neighbours = (colours - 1) * n_genes * 4.;

  analysed++;

  strict_robustnesses.emplace_back(genome_metric.strict_robustness / number_of_neighbours);
  intersection_robustnesses.emplace_back(genome_metric.intersection_robustness / number_of_neighbours);
  union_evolvabilities.emplace_back(genome_metric.union_evolvability / number_of_neighbours);
  rares.emplace_back(genome_metric.rare / number_of_neighbours);
  loops.emplace_back(genome_metric.loop / number_of_neighbours);

  neutral_weightings.emplace_back(genome_metric.neutral_weight);

  for (auto pID: genome_metric.diversity)
    diversity.insert(pID);

  genome_metrics.emplace_back(genome_metric);
}

void Set_Metrics::save_to_file(std::ofstream& set_out, std::ofstream& genome_out)
{
  for(auto metric: genome_metrics)
    metric.save_to_file(genome_out);

  double total_neutral_size = std::accumulate(neutral_weightings.begin(), neutral_weightings.end(), uint64_t(0));

  double average_strict_robustness = std::inner_product(std::begin(strict_robustnesses), std::end(strict_robustnesses), std::begin(neutral_weightings), 0) / total_neutral_size;
  double average_intersection_robustness = std::inner_product(std::begin(intersection_robustnesses), std::end(intersection_robustnesses), std::begin(neutral_weightings), 0) / total_neutral_size;
  double average_union_evolvability = std::inner_product(std::begin(union_evolvabilities), std::end(union_evolvabilities), std::begin(neutral_weightings), 0) / total_neutral_size;
  double average_rare = std::inner_product(std::begin(rares), std::end(rares), std::begin(neutral_weightings), 0) / total_neutral_size;
  double average_loop = std::inner_product(std::begin(loops), std::end(loops), std::begin(neutral_weightings), 0) / total_neutral_size;

  set_out <<+ average_strict_robustness << " " <<+ average_intersection_robustness << " ";
  set_out <<+ average_union_evolvability << " " <<+ average_rare << " ";
  set_out <<+ average_loop << " " <<+ analysed << " ";
  set_out <<+ total_neutral_size << " " << diversity.size() << " ";

  set_out << "{";
  for (auto pID: ref_pIDs)
    set_out <<+ "(" <<+ pID.first << "," <<+ pID.second << "),";

  set_out.seekp((long) set_out.tellp() - 1);
  set_out << "}" << std::endl;
}

void Set_Metrics::clear()
{
  strict_robustnesses.clear(), union_evolvabilities.clear(), rares.clear();
  intersection_robustnesses.clear(), neutral_weightings.clear(), loops.clear();

  genome_metrics.clear();
  diversity.clear();
}

// Sampling functions

void GP_MapSampler(uint8_t n_genes, uint8_t metric_colours, std::vector<Set_Metrics>& metrics,
   Set_to_Genome& set_to_genome, StochasticPhenotypeTable* pt)
{
  const uint8_t k_builds = 20;
  const uint32_t N_JIGGLE = 10;
  std::vector<Phenotype_ID> loop_pID = {{255, 0}};

  // metrics.resize(set_to_genome.size(), Set_Metrics(n_genes, metric_colours));

  uint16_t number_of_genomes = 0;
  for(Set_to_Genome::iterator iter = std::begin(set_to_genome); iter != std::end(set_to_genome); iter++)
    number_of_genomes += (iter->second).size();

  std::cout << "There are " <<+ number_of_genomes << " genomes to jiggle! \n";

  for(Set_to_Genome::iterator iter = std::begin(set_to_genome); iter != std::end(set_to_genome); iter++)
  {
    std::cout << "Currently processing " <<+ (iter->second).size() << " genomes for {";
    for(auto pID: iter->first)
      std::cout << "(" <<+ pID.first << ", " <<+ pID.second << "), ";
    number_of_genomes -= (iter->second).size();
    std::cout << "}. Only " <<+ number_of_genomes << " left! \n";

    if(iter->first == loop_pID)
      continue;

    Set_Metrics set_metrics(n_genes, metric_colours);
    set_metrics.ref_pIDs = iter->first;

    for(auto genotype: iter->second)
    {
      double neutral_weight = (double) NeutralSize(genotype, 1, metric_colours) / N_JIGGLE; // This weight will be counted N_JIGGLE time when added to form the total neutral weight

      #pragma omp parallel for schedule(dynamic) firstprivate(genotype, neutral_weight)
      for(uint32_t nth_jiggle=0; nth_jiggle<N_JIGGLE; ++nth_jiggle)
      {
        Clean_Genome(genotype, 0, false);
        JiggleGenotype(genotype, metric_colours);

        Genotype_Metrics genome_metric(n_genes, metric_colours);
        genome_metric.set_reference(genotype, iter->first, neutral_weight);

        for(Genotype neighbour : genotype_neighbourhood(genotype, n_genes, metric_colours))
        {
           std::vector<Phenotype_ID> neighbour_pIDs = GetSetPIDs(neighbour, k_builds, pt);
           genome_metric.analyse_pIDs(neighbour_pIDs);
         }
         #pragma omp critical
         {
           set_metrics.add_genotype_metrics(genome_metric);
         }
      }
    }
    metrics.emplace_back(set_metrics);
  }


  // #pragma omp parallel for schedule(dynamic)
  // for(uint32_t index=0; index < set_to_genome.size(); ++index)
  // {
  //   if ((index % 10) == 0)
  //     std::cout << "iteration : " <<+ index << " out of " <<+ set_to_genome.size() << "\n";
  //
  //   Set_to_Genome::iterator iter = std::begin(set_to_genome);
  //   std::advance(iter, index);
  //
  //   if(iter->first == loop_pID)
  //     continue;
  //
  //   Set_Metrics set_metrics(n_genes, metric_colours);
  //   set_metrics.ref_pIDs = iter->first;
  //
  //   for(auto genotype: iter->second)
  //   {
  //     for(uint32_t nth_jiggle=0; nth_jiggle<N_JIGGLE; ++nth_jiggle)
  //     {
  //       Clean_Genome(genotype, 0, false);
  //       JiggleGenotype(genotype, metric_colours);
  //
  //       Genotype_Metrics genome_metric(n_genes, metric_colours);
  //       genome_metric.set_reference(genotype, iter->first);
  //
  //       for(Genotype neighbour : genotype_neighbourhood(genotype, n_genes, metric_colours))
  //       {
  //          std::vector<Phenotype_ID> neighbour_pIDs = GetSetPIDs(neighbour, k_builds, pt);
  //          genome_metric.analyse_pIDs(neighbour_pIDs);
  //        }
  //        set_metrics.add_genotype_metrics(genome_metric);
  //      }
  //    }
  //   metrics[index] = set_metrics;
  // }
}

// Subroutine of the GP_MapSampler

void JiggleGenotype(Genotype& genotype, uint8_t max_colour)
{
  uint8_t min_colour =* std::max_element(genotype.begin(), genotype.end());

  if(min_colour + 1 == max_colour)
    return;

  std::vector<uint8_t> neutral_colours( 1+ (max_colour - min_colour) / 2);

  std::generate(neutral_colours.begin() + 1, neutral_colours.end(), [n = min_colour-1] () mutable { return n+=2; });

  std::uniform_int_distribution<size_t> jiggle_index(0, neutral_colours.size() - 1);

  for(auto& base : genotype)
    base= (base==0) ? neutral_colours[jiggle_index(RNG_Engine1)] : base;
}

std::vector<Genotype> genotype_neighbourhood(const Genotype& genome, uint8_t ngenes, uint8_t colours)
{
  std::vector<Genotype> neighbours;
  Genotype neighbour;

  std::vector<uint8_t> mutants(colours);
  std::iota(mutants.begin(), mutants.end(), 0);

  for(uint8_t index=0; index<ngenes*4; ++index)
  {
    std::swap(*std::find(mutants.begin(), mutants.end(),genome[index]), mutants.back());

    neighbour = genome;
    for(int j=0; j<colours-1; ++j)
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
  Clean_Genome(genotype, 0, false);
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
