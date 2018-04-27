#include "stochastic_model.hpp"
#include <iostream>

/*External wrappers for python integration */
extern "C"
{
  void GetPhenotypesIDs(const char* file_path_c,const char* file_name_c, uint8_t n_genes, uint8_t colours, bool file_of_genotypes);
  void ExhaustiveMinimalGenotypes(const char* file_path_c, uint8_t n_genes, uint8_t colours, bool file_of_genotypes);
  void SampleMinimalGenotypes(const char* file_path_c, uint8_t n_genes, uint8_t colours,const uint32_t N_SAMPLES,bool allow_duplicates, bool file_of_genotypes);
  void GP_MapSampler(const char* file_path_c,uint8_t n_genes, uint8_t rcolours,uint8_t colours,bool file_of_genotypes);
  void PreProcessGenotypes(const char* file_path_c, uint8_t n_genes, uint8_t colours,bool file_of_genotypes);
}

/*Converting methods*/
uint64_t genotype_to_index(Genotype& genotype, uint8_t n_genes, uint8_t colours);
void index_to_genotype(uint64_t index, Genotype& genotype, uint8_t n_genes, uint8_t colours);

/*GP map calculations*/
std::vector<Phenotype_ID> GetPhenotypeIDs(Genotype& genotype, uint8_t k_builds, StochasticPhenotypeTable* pt_it);
std::vector<Genotype> genotype_neighbourhood(const Genotype& genome, uint8_t ngenes, uint8_t colours);
void JiggleGenotype(Genotype& genotype, uint8_t max_colour);


/*Neutral size calculations*/
uint64_t NeutralSize(Genotype genotype,uint32_t N_neutral_colours,uint32_t N_possible_interacting_colours);
uint64_t combination_with_repetiton(uint8_t space_size , uint8_t sample_size);
uint64_t nChoosek(uint8_t n, uint8_t k);

/*Minimal genotype methods*/
struct NecklaceFactory {
  uint8_t colours=1;
  std::vector<std::vector<uint8_t> > necklaces;
  std::vector<uint8_t> necklace_grower;

  NecklaceFactory()  {necklace_grower.assign(5,0);}

  void GenNecklaces(uint8_t c) {
    colours=c;
    crsms_gen(1,1);
  }

  bool is_finite_necklace(std::vector<uint8_t>& neck) {
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

  void is_necklace(uint64_t j) {
    if(4%j==0)
      if(is_finite_necklace(necklace_grower))
        necklaces.emplace_back(std::vector<uint8_t>{necklace_grower.begin()+1,necklace_grower.end()});
  }

  void crsms_gen(uint64_t n, uint64_t j) {
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
};

struct GenotypeGenerator {
  bool is_done=false;
  uint8_t n_genes,colours;
  std::vector<uint32_t> necklace_states;
  std::vector<std::vector<uint8_t> > necklaces;
  uint32_t n_necklaces;

  void init() {
    NecklaceFactory necks=NecklaceFactory();
    necks.GenNecklaces(colours);
    necklaces=necks.necklaces;
    n_necklaces=necklaces.size();
    necklace_states[0]=1;
  }

  GenotypeGenerator(uint8_t a,uint8_t b) {n_genes=a;colours=b; necklace_states.assign(a,0);}

  Genotype operator() () {
    return !is_done ? next_genotype() : Genotype{};
  }

  bool valid_growing_faces(Genotype& genotype,uint8_t max_face) {
    for(auto face : genotype) {
      if(face>max_face)
        return false;
      if(face==max_face)
        max_face+=2;
    }
    return true;
  }

  bool valid_bindings(Genotype& genotype) {
    for(uint8_t interface=1;interface<=*std::max_element(genotype.begin(),genotype.end());interface+=2) {
      if(std::find(genotype.begin(),genotype.end(),interface)!=genotype.end()) { //is present
        if(std::find(genotype.begin(),genotype.end(),interface+1)==genotype.end()) //is not present
          return false;
      }
      else {
        if(std::find(genotype.begin(),genotype.end(),interface+1)!=genotype.end())
          return false;
      }
    }

    return true;
  }
  bool valid_genotype(Genotype& genotype) {
    if(!valid_growing_faces(genotype,1))
      return false;
    if(!valid_bindings(genotype))
      return false;
    return true;
  }

  void increment_states(std::vector<uint32_t>& states) {
    ++states.back();
    uint32_t zero_state_init=states[0];
    for(uint32_t rind=states.size();rind>0;--rind) {
      if(states[rind-1]>=n_necklaces) {
        if(rind==1) {
          is_done=true;
          return;
        }
        else {
          states[rind-1]=0;
          ++states[rind-2];
        }
      }
    }
    if(zero_state_init==1 && states[0]!=1)
      states[1]=colours+2;
    if(zero_state_init==static_cast<uint32_t>(colours+2) && states[0]!=static_cast<uint32_t>(colours+2)) {
      is_done=true;
      return;
    }
    auto max_iter=std::max_element(states.begin(),states.end());
    std::replace(max_iter,states.end(),static_cast<uint32_t>(0),*max_iter);
  }

  Genotype next_genotype() {
    Genotype genotype;
    while(!is_done) {
      inc_lab:
      genotype.clear();
      genotype.reserve(n_genes*4);
      for(uint8_t index=0; index!= n_genes;++index) {
	uint8_t input_face=*std::max_element(genotype.begin(),genotype.end());
	uint8_t next_max_face= input_face>0 ? ((input_face-1)/2)*2+3 : 1;

	while(!valid_growing_faces(necklaces[necklace_states[index]],next_max_face)) {
	  const uint32_t post_inc = necklace_states[index]+1;
	  if(post_inc==necklaces.size()) {
	    increment_states(necklace_states);
	    goto inc_lab;
	  }
	  std::fill(necklace_states.begin()+index,necklace_states.end(),post_inc);
	}

        genotype.insert(genotype.end(),necklaces[necklace_states[index]].begin(),necklaces[necklace_states[index]].end());
      }
      increment_states(necklace_states);

      if(valid_genotype(genotype))
        return genotype;
    }
    genotype.clear();
    return genotype;
  }
};

struct Genotype_Metrics
{
  uint8_t n_genes;
  double neutral_size;

  Genotype ref_genotype;
  std::vector <Phenotype_ID> ref_pIDs;

  Phenotype_ID pID_InfiniteLoop{0, 255}, pID_NonDeterministic{0, 0};

  std::vector <double> robustness;
  std::vector <double> new_evolvability;
  std::vector <double> death_InfiniteLoop;
  std::vector <double> death_NonDeterministic;
  std::vector <double> death;

  std::vector < std::set <Phenotype_ID> > diversity;

  Genotype_Metrics(uint8_t ngenes);

  void set_reference(Genotype& ref_genotype, std::vector <Phenotype_ID>& ref_pIDs);

  void clear();

  void analyse_pIDs(std::vector <Phenotype_ID>& pIDs);

  void save_to_file(std::ofstream& fout);
};

typedef struct Genotype_Metrics Genotype_Metrics;

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
  fout << diversity[0].size() << "\n";
}

struct Phenotype_Metrics
{
  uint8_t n_genes, colours;
  uint32_t N_JIGGLE=100;

  std::vector <Phenotype_ID> pIDs;

  std::vector <std::vector <double>> robustnesses;
  std::vector <std::vector <double>> new_evolvabilities;
  std::vector <std::vector <double>> deaths_InfiniteLoop;
  std::vector <std::vector <double>> deaths_NonDeterministic;
  std::vector <std::vector <double>> deaths;

  std::vector<uint64_t> neutral_weightings;

  std::vector <std::set <Phenotype_ID>> diversity;

  Phenotype_Metrics(uint8_t n_genes, uint8_t colours);

  void add_genotype_metrics(Genotype_Metrics& gmetrics);

  void save_to_file(std::ofstream& fout);

  void clear();

};

typedef struct Phenotype_Metrics Phenotype_Metrics;

Phenotype_Metrics::Phenotype_Metrics(uint8_t n_genes, uint8_t colours):
n_genes(n_genes), colours(colours), diversity(n_genes)
{}

void Phenotype_Metrics::add_genotype_metrics(Genotype_Metrics& gmetrics)
{
  // std::transform(robustness.begin( ), robustness.end( ), gmetrics.robustness.begin( ), robustness.begin( ), std::plus<double>( ));

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

  fout << "; robustness : ";

  fout <<+ avg_R << " ; death : " <<+ avg_D << " ; diversity : " << diversity[0].size() << "\n";

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
