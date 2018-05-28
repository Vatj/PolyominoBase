#include "genotype_metrics.hpp"
#include <sstream>
#include <iterator>
#include <functional>
#include <set>
#include <cstring>

std::mt19937 RNG_Engine1(std::random_device{}());

// Constructor of the Genotype_Metrics structure
Genotype_Metrics::Genotype_Metrics(uint8_t n_genes, uint8_t colours):
n_genes(n_genes), colours(colours)
{
  number_of_neighbours = (colours - 1) * n_genes * 4.;
}

void Genotype_Metrics::set_reference(Genotype& genotype, std::vector <Phenotype_ID>& pIDs)
{
  ref_genotype = genotype;
  ref_pIDs = pIDs;

  neutral_weight = NeutralSize(genotype, 1, colours - 1);

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

  if (std::find(std::begin(pIDs), std::end(pIDs), death_pID) != std::end(pIDs))
    death += 1.;

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
  fout <<+ death / number_of_neighbours << " ";
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
  strict_robustness = 0, union_evolvability = 0, death = 0;
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
n_genes(n_genes), colours(colours)
{}

// Member functions of the Set_Metrics structure

void Set_Metrics::add_genotype_metrics(Genotype_Metrics& genome_metrics)
{
  double number_of_neighbours = (colours - 1) * n_genes * 4.;

  strict_robustnesses.emplace_back(genome_metrics.strict_robustness / number_of_neighbours);
  intersection_robustnesses.emplace_back(genome_metrics.intersection_robustness / number_of_neighbours);
  union_evolvabilities.emplace_back(genome_metrics.union_evolvability / number_of_neighbours);
  deaths.emplace_back(genome_metrics.death / number_of_neighbours);
  loops.emplace_back(genome_metrics.loop / number_of_neighbours);

  neutral_weightings.emplace_back(genome_metrics.neutral_weight);

  for (auto pID: genome_metrics.diversity)
  {
    diversity.insert(pID);
  }
}

void Set_Metrics::save_to_file(std::ofstream& fout)
{
  double total_neutral_size = std::accumulate(neutral_weightings.begin(), neutral_weightings.end(), uint64_t(0));

  double average_strict_robustness = 0, average_intersection_robustness = 0;
  double average_union_evolvability = 0, average_death = 0, average_loop = 0;

  for(size_t index = 0; index < strict_robustnesses.size(); ++index)
  {
     average_strict_robustness += strict_robustnesses[index] * (double) neutral_weightings[index] / total_neutral_size;
     average_intersection_robustness += intersection_robustnesses[index] * (double) neutral_weightings[index] / total_neutral_size;
     average_union_evolvability += union_evolvabilities[index] * (double) neutral_weightings[index] / total_neutral_size;
     average_death += deaths[index] * (double) neutral_weightings[index] / total_neutral_size;
     average_loop += loops[index] * (double) neutral_weightings[index] / total_neutral_size;
  }

  fout <<+ average_strict_robustness << " ";
  fout <<+ average_intersection_robustness << " ";
  fout <<+ average_union_evolvability << " ";
  fout <<+ average_death << " ";
  fout <<+ average_loop << " ";
  fout <<+ diversity.size() << " ";

  if(ref_pIDs.size() > 0)
  {
    fout << "{";
    for (auto pID: ref_pIDs)
      fout <<+ "(" <<+ pID.first << "," <<+ pID.second << "),";

    fout.seekp((long) fout.tellp() - 1);
    fout << "}" << std::endl;
  }
  else
    fout << "{}" << std::endl;
}

void Set_Metrics::clear()
{
  strict_robustnesses.clear(), union_evolvabilities.clear(), deaths.clear();
  intersection_robustnesses.clear(), neutral_weightings.clear(), loops.clear();

  diversity.clear();
}

// Sampling functions

void GP_MapSampler(const char* file_path_c, uint8_t n_genes, uint8_t rcolours, uint8_t colours)
{

  const uint8_t k_builds = 10;
  const uint32_t N_JIGGLE = 10;
  StochasticPhenotypeTable pt;
  Genotype genotype(n_genes * 4);

  std::string file_path(file_path_c), file_details="_N"+std::to_string(n_genes)+"_C"+std::to_string(rcolours);
  std::ifstream file_in(file_path + "SampledGenotypes" + file_details + "_Processed.txt");
  file_details += "_Cx" + std::to_string(colours);

  // std::ofstream fout(file_path+"GP_Map" + std::to_string(colours)+file_details+".txt", std::ios_base::out);

  std::string shape_metrics_str = "shape_metrics", ending = ".txt";
  std::ofstream shape_metrics_out(file_path + shape_metrics_str + file_details + ending);

  std::string genome_metrics_str = "genome_metrics";
  std::ofstream genome_metrics_out(file_path + genome_metrics_str + file_details + ending);

  std::string set_metrics_str = "set_metrics";
  std::ofstream set_metrics_out(file_path + set_metrics_str + file_details + ending);

  Set_Metrics set_metrics(n_genes, colours);

  std::string str, separator="x";

  while (std::getline(file_in, str))
  {
    if(str.compare(0, 1, separator) == 0)
    {
      set_metrics.save_to_file(set_metrics_out);
      set_metrics.clear();
      continue;
    }

    std::istringstream is( str );
    genotype.assign( std::istream_iterator<int>( is ), std::istream_iterator<int>() );

    // std::vector<Phenotype_ID> genotype_pIDs = GetPhenotypeIDs(genotype, k_builds, &pt);
    std::vector<Phenotype_ID> genotype_pIDs = GetSetPIDs(genotype, k_builds, &pt);
    set_metrics.ref_pIDs = genotype_pIDs;

    Genotype_Metrics genome_metrics(n_genes, colours);

    #pragma omp parallel for schedule(dynamic) firstprivate(genotype, genome_metrics)
    for(uint32_t nth_jiggle=0; nth_jiggle<N_JIGGLE; ++nth_jiggle)
    {
      Clean_Genome(genotype, 0, false);
      JiggleGenotype(genotype, colours);

      genome_metrics.set_reference(genotype, genotype_pIDs);

      for(Genotype neighbour : genotype_neighbourhood(genotype, n_genes, colours))
      {
	       // std::vector<Phenotype_ID> neighbour_pIDs = GetPhenotypeIDs(neighbour, k_builds, &pt);
         std::vector<Phenotype_ID> neighbour_pIDs = GetSetPIDs(neighbour, k_builds, &pt);
         genome_metrics.analyse_pIDs(neighbour_pIDs);

       }

      #pragma omp critical
      {
        genome_metrics.save_to_file(genome_metrics_out);
        set_metrics.add_genotype_metrics(genome_metrics);
      }

      genome_metrics.clear();
    }
  }
  shape_metrics_out << "Done C++" << std::endl;
}

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
