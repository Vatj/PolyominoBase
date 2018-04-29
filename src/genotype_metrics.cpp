#include "genotype_metrics.hpp"
#include <sstream>
#include <iterator>
#include <functional>
#include <set>

std::mt19937 RNG_Engine1(std::random_device{}());

// Constructor of the Genotype_Metrics structure
Genotype_Metrics::Genotype_Metrics(uint8_t n_genes):
n_genes(n_genes), robustness(n_genes, 0.0f), new_evolvability(n_genes, 0.0f),
death_InfiniteLoop(n_genes, 0.0f), death_NonDeterministic(n_genes, 0.0f), death(n_genes, 0.0f),
diversity(n_genes)
{}

void Genotype_Metrics::set_reference(Genotype& genotype, std::vector <Phenotype_ID>& pIDs)
{
  ref_genotype = genotype;
  ref_pIDs = pIDs;
}

void Genotype_Metrics::analyse_pIDs(std::vector <Phenotype_ID>& pIDs)
{
  for (unsigned index=0; index < n_genes; ++index)
  {
    if (pIDs[index] == ref_pIDs[index])
      robustness[index] += 1.0f;

    else if (pIDs[index] == pID_InfiniteLoop)
      death_InfiniteLoop[index] += 1.0f;

    else if (pIDs[index] == pID_NonDeterministic)
      death_NonDeterministic[index]  += 1.0f;

    else
    {
      new_evolvability[index] += 1.0f;
      diversity[index].insert(pIDs[index]);
    }
  }
  std::transform(death_InfiniteLoop.begin( ), death_InfiniteLoop.end( ), death_NonDeterministic.begin( ), death.begin( ), std::plus<double>( ));
}

void Genotype_Metrics::save_to_file(std::ofstream& fout)
{
  fout << "genotype : ";

  for (auto face: ref_genotype)
    fout <<+ face << " ";

  fout << "; pIDs : ";

  for (auto pID: ref_pIDs)
    fout <<+ pID.first << " " <<+ pID.second << " ";

  fout << "; robustness : ";

  for (auto value: robustness)
    fout <<+ value << " ";

  fout << "; evolvability : ";

  for (auto value: new_evolvability)
    fout <<+ value << " ";

  fout << "; death : ";

  // Possible to look at Infinite Loop and Non-Deterministic later
  for (auto value: death)
    fout <<+ value << " ";

  fout << "; diversity : ";

  // Just the length of the set of different accessible Phenotype_ID
  for (auto set: diversity)
    fout << set.size() << " ";

  fout << std::endl;
}

void Genotype_Metrics::clear()
{
  std::fill(robustness.begin(), robustness.end(), 0.0f);
  std::fill(new_evolvability.begin(), new_evolvability.end(), 0.0f);
  std::fill(death_InfiniteLoop.begin(), death_InfiniteLoop.end(), 0.0f);
  std::fill(death_NonDeterministic.begin(), death_NonDeterministic.end(), 0.0f);
  std::fill(death.begin(), death.end(), 0.0f);

  for (auto set: diversity)
    set.clear();
}

// Constructor of the Phenotype_Metrics structure
Phenotype_Metrics::Phenotype_Metrics(uint8_t n_genes, uint8_t colours):
n_genes(n_genes), colours(colours), diversity(n_genes)
{}

// Member functions of the Phenotype_Metrics structure

void Phenotype_Metrics::add_genotype_metrics(Genotype_Metrics& gmetrics)
{

  robustnesses.emplace_back(gmetrics.robustness);
  new_evolvabilities.emplace_back(gmetrics.new_evolvability);
  deaths_InfiniteLoop.emplace_back(gmetrics.death_InfiniteLoop);
  deaths_NonDeterministic.emplace_back(gmetrics.death_NonDeterministic);
  deaths.emplace_back(gmetrics.death);

  for (unsigned index=0; index < n_genes; ++index)
  {
    for (auto pID: gmetrics.diversity[index])
    {
      diversity[index].insert(pID);
    }
  }

  neutral_weightings.emplace_back(NeutralSize(gmetrics.ref_genotype, 1, colours - 1));
}

void Phenotype_Metrics::save_to_file(std::ofstream& fout)
{

  double total_neutral_size=std::accumulate(neutral_weightings.begin(), neutral_weightings.end(), uint64_t(0));
  double avg_R=0, avg_D=0;

  for(size_t ind=0; ind<robustnesses.size(); ++ind)
  {
     avg_R += robustnesses[ind][0] * (double) neutral_weightings[ind] / total_neutral_size;
     avg_D += deaths[ind][0] * (double) neutral_weightings[ind] / total_neutral_size;
  }

  fout << "pIDs : ";

  for (auto pID: pIDs)
    fout <<+ pID.first << " " <<+ pID.second << " ";

  fout << "; robustness : " <<+ avg_R << " ; death : " <<+ avg_D << " ; diversity : " << diversity[0].size() << "\n";

  // for (auto value: robustness)
  //   fout << value << " ";
  //
  // for (auto value: new_evolvability)
  //   fout << value << " ";
  //
  // // Possible to look at Infinite Loop and Non-Deterministic later
  // for (auto value: death)
  //   fout << value << " ";
  //
  // // Just the length of the set of different accessible Phenotype_ID
  // fout << diversity.size() << "\n";
}

void Phenotype_Metrics::clear()
{
  robustnesses.clear();
  new_evolvabilities.clear();
  deaths_InfiniteLoop.clear();
  deaths_NonDeterministic.clear();
  deaths.clear();

  for (auto set: diversity)
    set.clear();
}

void GP_MapSampler(const char* file_path_c,uint8_t n_genes, uint8_t rcolours, uint8_t colours, bool file_of_genotypes)
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

  while (std::getline(file_in, str))
  {
    if(str=="x")
    {
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

uint64_t NeutralSize(Genotype genotype,uint32_t N_neutral_colours,uint32_t N_possible_interacting_colours)
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
