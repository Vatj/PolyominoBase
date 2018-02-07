#include "interface_model.hpp"

namespace simulation_params
{
  uint16_t population_size=10;
  uint8_t phenotype_builds=10,n_tiles=2;
  uint32_t generation_limit=5,independent_trials=1,run_offset=0;
  bool fitness_selection=false,random_initilisation=false;
}

namespace model_params
{
  const uint8_t interface_size=CHAR_BIT*sizeof(interface_model::interface_type);
  
  double temperature=1,mu_prob=0.2,unbound_factor=1,misbinding_rate=0,fitness_factor=1,UND_threshold=0.2;
 
  std::binomial_distribution<uint8_t> b_dist(interface_size,mu_prob);
  std::uniform_real_distribution<double> real_dist(0, 1);
  

}

namespace interface_model
{
  std::random_device rd;
  std::mt19937 RNG_Engine(rd());  
  //std::mt19937 RNG_Engine(276358710);
  
  inline interface_type reverse_bits(interface_type v) {
    interface_type s = sizeof(v) * CHAR_BIT; // bit size; must be power of 2
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
  inline uint8_t SammingDistance(interface_type face1,interface_type face2) {
    return static_cast<uint8_t>(ArbitraryPopcount(face1 ^ reverse_bits(~face2)));
  }

  void MutateInterfaces(std::vector<interface_type>& binary_genome) {
    std::vector<uint8_t> interface_indices(model_params::interface_size);
    std::iota(interface_indices.begin(),interface_indices.end(),0);
    for(interface_type& base : binary_genome) {
      std::shuffle(interface_indices.begin(), interface_indices.end(), RNG_Engine);
      for(uint8_t nth=0;nth<model_params::b_dist(RNG_Engine);++nth) 
        base ^= (static_cast<interface_type>(1) << interface_indices[nth]);
    }
  }

  double ProteinAssemblyOutcome(std::vector<interface_type> binary_genome,PhenotypeTable* pt) {
    uint8_t dx=0,dy=0;
    std::vector<int8_t> assembly_information;
    std::vector<uint8_t> spatial_information;
    std::vector<std::pair<uint8_t,uint32_t> > phenotype_IDs;phenotype_IDs.reserve(simulation_params::phenotype_builds);
    for(uint8_t nth=0;nth<simulation_params::phenotype_builds;++nth) {
      assembly_information=AssembleProtein(binary_genome);
      spatial_information=SpatialGrid(assembly_information,dx,dy);
      phenotype_IDs.emplace_back(std::accumulate(spatial_information.begin(),spatial_information.end(),0),pt->PhenotypeCheck(spatial_information,dx,dy));
    }
    return pt->GenotypeFitness(phenotype_IDs);
  }  

  std::vector<int8_t> AssembleProtein(const std::vector<interface_type>& binary_genome) {
    std::vector<int8_t> placed_tiles{0,0};
    std::vector<int8_t> growing_perimeter; 
    PerimeterGrowth(0,0,0,-1,0,growing_perimeter,placed_tiles);
    int8_t current_orientation, current_tile, current_direction, current_x, current_y;
    std::vector<uint8_t> genome_bases(binary_genome.size());
    std::iota(genome_bases.begin(), genome_bases.end(), 0);
    
    while(!growing_perimeter.empty()) {
      current_orientation=growing_perimeter.back();growing_perimeter.pop_back();
      current_tile=growing_perimeter.back();growing_perimeter.pop_back();
      current_direction=growing_perimeter.back();growing_perimeter.pop_back();
      current_y=growing_perimeter.back();growing_perimeter.pop_back();
      current_x=growing_perimeter.back();growing_perimeter.pop_back();
      std::shuffle(genome_bases.begin(), genome_bases.end(), RNG_Engine);
      for(uint8_t base : genome_bases) {
        if(model_params::real_dist(RNG_Engine)<std::exp(-1*static_cast<double>(SammingDistance(binary_genome[current_tile*4+current_orientation],binary_genome[base]))/(model_params::interface_size*model_params::temperature))) {
            placed_tiles.insert(placed_tiles.end(),{current_x,current_y});
            PerimeterGrowth(current_x,current_y,(4+current_direction-(base%4))%4,current_direction,base/4,growing_perimeter,placed_tiles);
            break;
          }
      }    
      if(placed_tiles.size()>(model_params::unbound_factor*binary_genome.size()*binary_genome.size()/2.))
        return {};
    }
    return placed_tiles;   
  }
  
  void PerimeterGrowth(int8_t x,int8_t y,int8_t theta,int8_t direction, int8_t tile_type,std::vector<int8_t>& growing_perimeter,std::vector<int8_t>& placed_tiles) {
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
      for(std::vector<int8_t>::reverse_iterator tile_info=placed_tiles.rbegin();tile_info!=placed_tiles.rend();tile_info+=2) {
        if(x+dx==*(tile_info+1) && y+dy==*(tile_info)) {
          occupied_site=true;
          break;
        }
      }
      if(!occupied_site)
        growing_perimeter.insert(growing_perimeter.begin()+5*index_randomizer(RNG_Engine),{static_cast<int8_t>(x+dx),static_cast<int8_t>(y+dy),static_cast<int8_t>((f+2)%4),tile_type,static_cast<int8_t>((4+f-theta)%4)});      
    }
  }
  
  
  
}//end interface_model namespace

std::vector<uint8_t> SpatialGrid(std::vector<int8_t>& placed_tiles, uint8_t& dx,uint8_t& dy) {
  if(placed_tiles.empty())
    return {};
  std::vector<int8_t> x_locs, y_locs;
  x_locs.reserve(placed_tiles.size()/2);
  y_locs.reserve(placed_tiles.size()/2);
  for(std::vector<int8_t>::iterator check_iter = placed_tiles.begin();check_iter!=placed_tiles.end();check_iter+=2) {
    x_locs.emplace_back(*check_iter);
    y_locs.emplace_back(*(check_iter+1));
  }
  std::vector<int8_t>::iterator x_left,x_right,y_top,y_bottom;
  std::tie(x_left,x_right)=std::minmax_element(x_locs.begin(),x_locs.end());
  std::tie(y_bottom,y_top)=std::minmax_element(y_locs.begin(),y_locs.end());
  dx=*x_right-*x_left+1;
  dy=*y_top-*y_bottom+1;
  std::vector<uint8_t> spatial_grid(dx*dy); 
  for(uint16_t tileIndex=0;tileIndex<x_locs.size();++tileIndex)
    spatial_grid[(*y_top-y_locs[tileIndex])*dx + (x_locs[tileIndex]-*x_left)]=1;
  return spatial_grid;
}

bool ComparePolyominoes(std::vector<uint8_t>& phenotype,uint8_t& dx,uint8_t& dy, std::vector<uint8_t>& phenotype_prime,uint8_t dx_prime,uint8_t dy_prime) {
  if(std::accumulate(phenotype.begin(),phenotype.end(),0)!=std::accumulate(phenotype_prime.begin(),phenotype_prime.end(),0))
    return false;
  if(dx==dx_prime && dy==dy_prime && dx==dy) { //bounding boxes match, symmetric
    if(phenotype==phenotype_prime) 
      return true;
    else { //SYMMETRY ARGUMENT TO OPTIMIZE
      for(int rotation=0;rotation<3;++rotation) {
	ClockwiseRotation(phenotype,dx,dy);
	if(phenotype==phenotype_prime)
	  return true;
      }
      return false;
    }
  }
  if(dx==dx_prime && dy==dy_prime)  //bounding boxes match, asymmetric
    return phenotype==phenotype_prime ? true : phenotype==ClockwisePiRotation(phenotype_prime);
  if(dx==dy_prime && dy==dx_prime) { //bounding boxes pi/2 off, asymmetric
    ClockwiseRotation(phenotype,dx,dy);
    return phenotype==phenotype_prime ? true : phenotype==ClockwisePiRotation(phenotype_prime);
  }
  return false;
}
void ClockwiseRotation(std::vector<uint8_t>& phenotype,uint8_t& dx,uint8_t& dy) {
  std::vector<uint8_t> swapper;
  swapper.reserve(phenotype.size());
  for(uint8_t column=0;column<dx;++column) 
    for(uint8_t row=dy;row!=0;--row) 
      swapper.emplace_back(phenotype[(row-1)*dx+column]);
  std::swap(dx,dy);
  phenotype=swapper;
}
std::vector<uint8_t> ClockwisePiRotation(std::vector<uint8_t>& phenotype) {
  std::reverse(phenotype.begin(),phenotype.end());
  return phenotype;
}
void PrintShape(std::vector<uint8_t>& spatial_information,uint8_t dx,uint8_t dy) {
  for(uint8_t y=0;y<dy;++y) {
    for(uint8_t x=0;x<dx;++x)
      std::cout<<+spatial_information[y*dx+x]<<" ";
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
    ClockwiseRotation(rotated_shape,dx,dy);
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

std::vector<uint16_t> InterfaceStrengths(std::vector<interface_model::interface_type>& interfaces) {
  std::vector<uint16_t> strengths(static_cast<uint8_t>(model_params::interface_size*1.5)+2);
  for(std::vector<interface_model::interface_type>::const_iterator outer_face=interfaces.begin(); outer_face!=interfaces.end(); ++outer_face) {
    ++strengths[static_cast<uint8_t>(1.5*model_params::interface_size)+1-interface_model::SammingDistance(*outer_face,*outer_face)/2];
    for(std::vector<interface_model::interface_type>::const_iterator inner_face=outer_face+1; inner_face!=interfaces.end(); ++inner_face)
      ++strengths[model_params::interface_size-interface_model::SammingDistance(*outer_face,*inner_face)];
  }
  return strengths;
}

