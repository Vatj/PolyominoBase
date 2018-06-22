#include "interface_model.hpp"

namespace simulation_params
{
  uint16_t population_size=10;
  uint8_t phenotype_builds=10,n_tiles=2;
  uint32_t generation_limit=5,independent_trials=1,run_offset=0;
  bool random_initilisation=false;
}

namespace model_params
{
  std::mt19937 RNG_Engine(std::random_device{}());
  bool fixed_seed=true;
  double temperature=1,binding_threshold=1,mu_prob=0.2,fitness_factor=1,UND_threshold=0.2,interface_threshold=0.2;
  std::binomial_distribution<uint8_t> b_dist(interface_size,mu_prob);
  std::uniform_real_distribution<double> real_dist(0, 1);
  std::array<double,model_params::interface_size+1> binding_probabilities;
}


namespace interface_model
{

  inline interface_type reverse_bits(interface_type v) {
    interface_type s(model_params::interface_size);
    interface_type mask = ~0;
    while ((s >>= 1) > 0) {
      mask ^= (mask << s);
      v = ((v >> s) & mask) | ((v << s) & ~mask);
    }
    return v;
  }
  inline uint8_t ArbitraryPopcount(interface_type face) {
    uint8_t c;
    for(c = 0; face; c++)
      face &= face - 1;
    return c;
  }
  uint8_t SammingDistance(interface_type face1,interface_type face2) {
    return static_cast<uint8_t>(ArbitraryPopcount(face1 ^ reverse_bits(~face2)));
  }

  void MutateInterfaces(BGenotype& binary_genome) {
    std::vector<uint8_t> interface_indices(model_params::interface_size);
    std::iota(interface_indices.begin(),interface_indices.end(),0);
    for(interface_type& base : binary_genome) {
      std::shuffle(interface_indices.begin(), interface_indices.end(), model_params::RNG_Engine);
      for(uint8_t nth=0;nth<model_params::b_dist(model_params::RNG_Engine);++nth)
        base ^= (interface_type(1) << interface_indices[nth]);
    }
  }

  double ProteinAssemblyOutcome(BGenotype binary_genome,InterfacePhenotypeTable* pt,Phenotype_ID& pid,std::vector<interaction_pair>& pid_interactions) {
    std::vector<int8_t> assembly_information;
    Phenotype phen;
    std::vector<Phenotype_ID> Phenotype_IDs;Phenotype_IDs.reserve(simulation_params::phenotype_builds);
    std::set<interaction_pair > interacting_indices;
    std::map<Phenotype_ID, std::map<interaction_pair, uint8_t> > phenotype_interactions;
    for(uint8_t nth=0;nth<simulation_params::phenotype_builds;++nth) {
      assembly_information=AssembleProteinNew(binary_genome,interacting_indices);

      if(assembly_information.size()>0) {
        for(auto interacting_pair : interacting_indices)
	  ++phenotype_interactions[Phenotype_IDs.back()][interacting_pair];
        phen=SpatialGrid(assembly_information);
        Phenotype_IDs.emplace_back(pt->GetPhenotypeID(phen));
      }
      else
        Phenotype_IDs.emplace_back(0,0);
      interacting_indices.clear();
    }
    pt->RelabelPhenotypes(Phenotype_IDs,phenotype_interactions);
    std::map<Phenotype_ID,uint8_t> ID_counter=pt->PhenotypeFrequencies(Phenotype_IDs);
    pid=std::max_element(ID_counter.begin(),ID_counter.end(),[] (const auto & p1, const auto & p2) {return p1.second < p2.second;})->first;
    for(auto& imap : phenotype_interactions[pid])
      if(imap.second>=static_cast<uint8_t>(simulation_params::phenotype_builds*model_params::interface_threshold))
	pid_interactions.emplace_back(imap.first);
    pt->AssignInitialFitnesses();
    return pt->GenotypeFitness(ID_counter);
  }

  std::vector<int8_t> AssembleProtein(const BGenotype& binary_genome,std::set<interaction_pair>& interacting_indices) {

    std::vector<int8_t> placed_tiles{0,0,1},growing_perimeter;
    PerimeterGrowth(0,0,0,-1,0,growing_perimeter,placed_tiles);
    int8_t orientation, tile, direction, x, y;
    std::vector<uint8_t> genome_bases(binary_genome.size());
    std::iota(genome_bases.begin(), genome_bases.end(), 0);

    while(!growing_perimeter.empty()) {
      orientation=growing_perimeter.back();growing_perimeter.pop_back();
      tile=growing_perimeter.back();growing_perimeter.pop_back();
      direction=growing_perimeter.back();growing_perimeter.pop_back();
      y=growing_perimeter.back();growing_perimeter.pop_back();
      x=growing_perimeter.back();growing_perimeter.pop_back();

      std::shuffle(genome_bases.begin(), genome_bases.end(), model_params::RNG_Engine);

      for(uint8_t base : genome_bases) {
        if(model_params::real_dist(model_params::RNG_Engine)<model_params::binding_probabilities[SammingDistance(binary_genome[tile*4+orientation],binary_genome[base])]) {
	  interacting_indices.insert(std::minmax(static_cast<uint8_t>(tile*4+orientation),base));
          placed_tiles.insert(placed_tiles.end(),{x,y,static_cast<int8_t>(base%4+((GAUGE==4)*(-2*(base%4)+(direction-base%4+4)%4+1)))});
	  PerimeterGrowth(x,y,(4+direction-(base%4))%4,direction,base/4,growing_perimeter,placed_tiles);
	  break;
	}
      }
      if(placed_tiles.size()>(static_cast<uint8_t>(0.75*binary_genome.size()*binary_genome.size())))
        return {};
    }

    return placed_tiles;
  }

  void PerimeterGrowth(int8_t x,int8_t y,int8_t theta,int8_t direction, int8_t tile_type,std::vector<int8_t>& growing_perimeter,std::vector<int8_t>& placed_tiles) {
    for(std::vector<int8_t>::iterator tile_info=growing_perimeter.begin();tile_info!=growing_perimeter.end();) {
      if(x==*(tile_info) && y==*(tile_info+1)) {
	//if(tile_type!=*(tile_info+3))
	tile_info=growing_perimeter.erase(tile_info,tile_info+5);
      }
      else
	tile_info+=5;
    }
    int8_t dx=0,dy=0;
    for(uint8_t f=0;f<4;++f) {
      if(direction==f)
        continue;
      switch(f) {
      case 0:dx=0;dy=1;break;
      case 1:dx=1;dy=0;break;
      case 2:dx=0;dy=-1;break;
      case 3:dx=-1;dy=0;break;
      }
      bool occupied_site=false;
      std::uniform_int_distribution<int> index_randomizer(0, growing_perimeter.size()/5);
      for(std::vector<int8_t>::reverse_iterator tile_info=placed_tiles.rbegin();tile_info!=placed_tiles.rend();tile_info+=3) {
        if(x+dx==*(tile_info+2) && y+dy==*(tile_info+1)) {
          occupied_site=true;
          break;
        }
      }
      if(!occupied_site)
        growing_perimeter.insert(growing_perimeter.begin()+5*index_randomizer(model_params::RNG_Engine),{static_cast<int8_t>(x+dx),static_cast<int8_t>(y+dy),static_cast<int8_t>((f+2)%4),tile_type,static_cast<int8_t>((4+f-theta)%4)});
    }
  }

  std::vector<int8_t> AssembleProteinNew(const BGenotype& binary_genome,std::set<interaction_pair>& interacting_indices) {
    int8_t seed=1;
    if(!model_params::fixed_seed)
      seed=1+4*(std::uniform_int_distribution<uint8_t>(1,binary_genome.size()/4)(model_params::RNG_Engine)-1);

    std::vector<int8_t> placed_tiles{0,0,seed},growing_perimeter;
    std::vector<double> strengths,strengths_cdf;
    std::vector<interaction_pair> interaction_pairs;
    ExtendPerimeter(binary_genome,seed,0,0,placed_tiles,growing_perimeter,strengths,interaction_pairs);

    while(!strengths.empty()) {
      //select new site proportional to binding strength
      strengths_cdf.resize(strengths.size());
      std::partial_sum(strengths.begin(), strengths.end(), strengths_cdf.begin());
      std::uniform_real_distribution<double> random_interval(0,strengths_cdf.back());
      size_t selected_choice=static_cast<size_t>(std::lower_bound(strengths_cdf.begin(),strengths_cdf.end(),random_interval(model_params::RNG_Engine))-strengths_cdf.begin());
      //place new tile
      const int8_t f_x=*(growing_perimeter.begin()+selected_choice*3),f_y=*(growing_perimeter.begin()+selected_choice*3+1),f_t=*(growing_perimeter.begin()+selected_choice*3+2);
      placed_tiles.insert(placed_tiles.end(),growing_perimeter.begin()+selected_choice*3,growing_perimeter.begin()+selected_choice*3+3);
      interacting_indices.insert(interaction_pairs[selected_choice]);
      //remove all further options in same tile location
      for(size_t cut_index=0;cut_index<strengths.size();) {
	if(f_x==*(growing_perimeter.begin()+cut_index*3) && f_y==*(growing_perimeter.begin()+cut_index*3+1)) {
	  strengths.erase(strengths.begin()+cut_index);

          interaction_pairs.erase(interaction_pairs.begin()+cut_index);
	  growing_perimeter.erase(growing_perimeter.begin()+cut_index*3,growing_perimeter.begin()+(cut_index+1)*3);
	}
	else
	  ++cut_index;
      }

      //find new potential sites
      ExtendPerimeter(binary_genome,f_t,f_x,f_y,placed_tiles,growing_perimeter,strengths,interaction_pairs);
      //unbound limit reached
      if(placed_tiles.size()>static_cast<size_t>(0.75*binary_genome.size()*binary_genome.size()))
        return {};
    }
    return placed_tiles;
  }

  void ExtendPerimeter(const BGenotype& binary_genome,uint8_t tile_detail, int8_t x,int8_t y, std::vector<int8_t>& placed_tiles,std::vector<int8_t>& new_sites,std::vector<double>& binding_strengths,std::vector<interaction_pair>& interaction_pairs) {
    int8_t dx=0,dy=0,tile=(tile_detail-1)/4,theta=(tile_detail-1)%4;
    for(uint8_t f=0;f<4;++f) {
      switch(f) {
      case 0:dx=0;dy=1;break;
      case 1:dx=1;dy=0;break;
      case 2:dx=0;dy=-1;break;
      case 3:dx=-1;dy=0;break;
      }
      //site already occupied, move on
      for(std::vector<int8_t>::reverse_iterator tile_info=placed_tiles.rbegin();tile_info!=placed_tiles.rend();tile_info+=3)
        if((x+dx)==*(tile_info+2) && (y+dy)==*(tile_info+1))
	  goto nextloop;

      //find all above threshold bindings and add to new sites
      for(uint8_t base=0;base<binary_genome.size();++base) {
	uint8_t SD=SammingDistance(binary_genome[base],binary_genome[tile*4+(f-theta+4)%4]);
	if(SD<=static_cast<uint8_t>(model_params::interface_size*(1-model_params::binding_threshold))){
	  binding_strengths.emplace_back(model_params::binding_probabilities[SD]);
	  new_sites.insert(new_sites.end(),{static_cast<int8_t>(x+dx),static_cast<int8_t>(y+dy),static_cast<int8_t>(base-base%4+((f+2)%4-base%4+4)%4+1)});
          interaction_pairs.emplace_back(std::minmax(static_cast<uint8_t>(tile*4+(f-theta+4)%4),base));
	}
      }
    nextloop: ; //continue if site already occupied
    }

  }

}//end interface_model namespace

Phenotype SpatialGrid(std::vector<int8_t>& placed_tiles) {
  std::vector<int8_t> x_locs, y_locs,tile_vals;
  x_locs.reserve(placed_tiles.size()/3);y_locs.reserve(placed_tiles.size()/3);tile_vals.reserve(placed_tiles.size()/3);

  for(std::vector<int8_t>::iterator check_iter = placed_tiles.begin();check_iter!=placed_tiles.end();check_iter+=3) {
    x_locs.emplace_back(*check_iter);
    y_locs.emplace_back(*(check_iter+1));
    tile_vals.emplace_back(*(check_iter+2));
  }
  std::vector<int8_t>::iterator x_left,x_right,y_top,y_bottom;
  std::tie(x_left,x_right)=std::minmax_element(x_locs.begin(),x_locs.end());
  std::tie(y_bottom,y_top)=std::minmax_element(y_locs.begin(),y_locs.end());
  uint8_t dx=*x_right-*x_left+1,dy=*y_top-*y_bottom+1;
  std::vector<uint8_t> spatial_grid(dx*dy);
  for(uint16_t tileIndex=0;tileIndex<x_locs.size();++tileIndex)
    spatial_grid[(*y_top-y_locs[tileIndex])*dx + (x_locs[tileIndex]-*x_left)]=tile_vals[tileIndex];

  Phenotype phen{dx,dy,spatial_grid};
  if(dy>dx) {
    ClockwiseRotation(phen);
  }
  MinimizePhenRep(phen.tiling);
  return phen;
}


void PrintShape(Phenotype phen) {
  for(uint8_t y=0;y<phen.dy;++y) {
    for(uint8_t x=0;x<phen.dx;++x)
      std::cout<<+phen.tiling[y*phen.dx+x]<<" ";
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
}

uint8_t PhenotypeSymmetryFactor(std::vector<uint8_t>& original_shape, uint8_t dx, uint8_t dy) {
  std::vector<uint8_t> rotated_shape(original_shape);
  std::reverse(original_shape.begin(),original_shape.end());
  if(original_shape!=rotated_shape)
    return 1;
  if(dx==dy) {
    //ClockwiseRotation(Phenotyrotated_shape,dx,dy);
    if(original_shape==rotated_shape)
      return 4;
  }
  return 2;
}

void DistributionStatistics(std::vector<double>& intf, double& mean, double& variance) {
  uint16_t N = 0;
  double Mprev = 0;
  for(double x : intf) {
    ++N;
    Mprev = mean;
    mean += (x - Mprev) / N;
    variance += (x - Mprev) * (x - mean);
  }
  variance /=(N-1);
}

void InterfaceStrengths(BGenotype& interfaces, std::vector<uint32_t>& strengths) {
  for(BGenotype::const_iterator outer_face=interfaces.begin(); outer_face!=interfaces.end(); ++outer_face) {
    ++strengths[static_cast<uint8_t>(1.5*model_params::interface_size)+1-interface_model::SammingDistance(*outer_face,*outer_face)/2];
    for(BGenotype::const_iterator inner_face=outer_face+1; inner_face!=interfaces.end(); ++inner_face)
      ++strengths[model_params::interface_size-interface_model::SammingDistance(*outer_face,*inner_face)];
  }
}

std::array<double,model_params::interface_size+1> GenBindingProbsLUP() {
  std::array<double,model_params::interface_size+1> probs;
  size_t N_active_strengths=static_cast<size_t>((1-model_params::binding_threshold)*model_params::interface_size)+1;
  std::fill(probs.begin(),probs.begin()+N_active_strengths,1);
  std::fill(probs.begin()+N_active_strengths,probs.end(),0);
  if(model_params::temperature<0 || static_cast<uint32_t>(model_params::binding_threshold)==1) {
    model_params::temperature=0;
    return probs;
  }
  for(size_t i=0;i<probs.size();++i) {

    probs[i]=(i<=static_cast<size_t>((1.-model_params::binding_threshold)*model_params::interface_size)?1:0)*exp((double(i)/model_params::interface_size)/((model_params::binding_threshold-1.)*model_params::temperature));

  }
  return probs;
}
