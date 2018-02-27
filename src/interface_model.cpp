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
  const uint8_t interface_size=CHAR_BIT*sizeof(interface_type);
  double temperature=1,mu_prob=0.2,unbound_factor=1,misbinding_rate=0,fitness_factor=1,UND_threshold=0.2;
  std::binomial_distribution<uint8_t> b_dist(interface_size,mu_prob);
  std::uniform_real_distribution<double> real_dist(0, 1);
  

}

namespace interface_model
{
  std::random_device rd;
  //std::mt19937 RNG_Engine(rd());  
  std::mt19937 RNG_Engine(276358710);
  
  inline interface_type reverse_bits(interface_type v) {
    interface_type s = sizeof(v) * CHAR_BIT;
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

  void MutateInterfaces(std::vector<interface_type>& binary_genome) {
    std::vector<uint8_t> interface_indices(model_params::interface_size);
    std::iota(interface_indices.begin(),interface_indices.end(),0);
    for(interface_type& base : binary_genome) {
      std::shuffle(interface_indices.begin(), interface_indices.end(), RNG_Engine);
      for(uint8_t nth=0;nth<model_params::b_dist(RNG_Engine);++nth) 
        base ^= (static_cast<interface_type>(1) << interface_indices[nth]);
    }
  }

  double ProteinAssemblyOutcome(std::vector<interface_type> binary_genome,PhenotypeTable* pt,phenotype_ID& pid) {
    uint8_t dx=0,dy=0;
    std::vector<int8_t> assembly_information;
    Phenotype phen;
    std::vector<phenotype_ID> phenotype_IDs;phenotype_IDs.reserve(simulation_params::phenotype_builds);
    for(uint8_t nth=0;nth<simulation_params::phenotype_builds;++nth) {
      assembly_information=AssembleProtein(binary_genome);
      //std::cout<<"ass size "<<assembly_information.size()<<std::endl;
      if(assembly_information.size()>0) {
        phen=SpatialGrid(assembly_information,dx,dy);
	std::vector<uint8_t> bb8{1,5,7,3};
	if(false && phen.dx==phen.dy && phen.tiling==bb8) {
	  
	  std::cout<<"bb8d: "<<assembly_information.size()<<"@ ";
	  for(auto x : binary_genome)
	    std::cout<<+x<<" ";
	  std::cout<<std::endl;
	  return 1;
	}
        phenotype_IDs.emplace_back(std::count_if(phen.tiling.begin(),phen.tiling.end(),[](const int c){return c != 0;}),pt->PhenotypeCheck(phen));
      }
      else
        phenotype_IDs.emplace_back(0,0);
    }
    auto ID_counter=pt->PhenotypeFrequencies(phenotype_IDs);
    pid=std::max_element(ID_counter.begin(),ID_counter.end(),[] (const auto & p1, const auto & p2) {return p1.second < p2.second;})->first;
    return pt->GenotypeFitness(ID_counter);
  }


  std::vector<int8_t> AssembleProtein(const std::vector<interface_type>& binary_genome) {
    std::vector<int8_t> placed_tiles{0,0,1};
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

	  
          placed_tiles.insert(placed_tiles.end(),{current_x,current_y,static_cast<int8_t>(base-base%4+(current_direction-base%4+4)%4+1)});
	  //std::cout<<base/4+(current_direction-base%4+4)%4+1<<std::endl;
	  PerimeterGrowth(current_x,current_y,(4+current_direction-(base%4))%4,current_direction,base/4,growing_perimeter,placed_tiles);
	  
	  break;
	}
      }    
      if(placed_tiles.size()>(static_cast<uint8_t>(model_params::unbound_factor*0.75*binary_genome.size()*binary_genome.size())))
	
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
        growing_perimeter.insert(growing_perimeter.begin()+5*index_randomizer(RNG_Engine),{static_cast<int8_t>(x+dx),static_cast<int8_t>(y+dy),static_cast<int8_t>((f+2)%4),tile_type,static_cast<int8_t>((4+f-theta)%4)});      
    }
  }

}//end interface_model namespace

Phenotype SpatialGrid(std::vector<int8_t>& placed_tiles, uint8_t& dx,uint8_t& dy) {
  //if(placed_tiles.empty())
  //  return {};
  std::vector<int8_t> x_locs, y_locs,tile_vals;
  x_locs.reserve(placed_tiles.size()/3);
  y_locs.reserve(placed_tiles.size()/3);
  tile_vals.reserve(placed_tiles.size()/3);
  for(std::vector<int8_t>::iterator check_iter = placed_tiles.begin();check_iter!=placed_tiles.end();check_iter+=3) {
    x_locs.emplace_back(*check_iter);
    y_locs.emplace_back(*(check_iter+1));
    tile_vals.emplace_back(*(check_iter+2));//*4+*(check_iter+3));
  }
  std::vector<int8_t>::iterator x_left,x_right,y_top,y_bottom;
  std::tie(x_left,x_right)=std::minmax_element(x_locs.begin(),x_locs.end());
  std::tie(y_bottom,y_top)=std::minmax_element(y_locs.begin(),y_locs.end());
  dx=*x_right-*x_left+1;
  dy=*y_top-*y_bottom+1;
  std::vector<uint8_t> spatial_grid(dx*dy); 
  for(uint16_t tileIndex=0;tileIndex<x_locs.size();++tileIndex)
    spatial_grid[(*y_top-y_locs[tileIndex])*dx + (x_locs[tileIndex]-*x_left)]=tile_vals[tileIndex];

  Phenotype phen{dx,dy,spatial_grid};
  if(dx<dy) {
    ClockwiseRotation(phen);
  }
    
  MinimalTilingRepresentation(phen.tiling);
  return phen;
}

bool ComparePolyominoes(Phenotype& phen1, const Phenotype& phen2) {
  if(std::count(phen1.tiling.begin(),phen1.tiling.end(),0)!=std::count(phen2.tiling.begin(),phen2.tiling.end(),0))
    return false;
  if(phen1.dx==phen2.dx && phen1.dy==phen2.dy && phen1.dx==phen2.dy) { //bounding boxes match, symmetric
    if(phen1.tiling==phen2.tiling) 
      return true;
    else { //SYMMETRY ARGUMENT TO OPTIMIZE
      for(int rotation=0;rotation<3;++rotation) {
	ClockwiseRotation(phen1);
        if(phen1.tiling==phen2.tiling) 
	  return true;
      }
      return false;
    }
  }
  if(phen1.dx==phen2.dx && phen1.dy==phen2.dy) {  //bounding boxes match, asymmetric
    if(phen1.tiling==phen2.tiling)
      return true;
    else {
      ClockwisePiRotation(phen1);
      return phen1.tiling==phen2.tiling;
    }
  }
  /*
  if(phen1.dx==phen2.dy && phen1.dy==phen2.dx) { //bounding boxes pi/2 off, asymmetric
    ClockwiseRotation(phen1);
    if(phen1.tiling==phen2.tiling)
      return true;
    else {
      ClockwisePiRotation(phen1);
      return phen1.tiling==phen2.tiling;
    }
  }
  */
  return false;
}

void ClockwiseRotation(Phenotype& phen) {
  std::vector<uint8_t> swapper;
  swapper.reserve(phen.tiling.size());

  for(uint8_t column=0;column<phen.dx;++column) 
    for(uint8_t row=phen.dy;row!=0;--row)  {
      if(phen.tiling[(row-1)*phen.dx+column])
        swapper.emplace_back(phen.tiling[(row-1)*phen.dx+column]);//+(phen.tiling[(row-1)*phen.dx+column])%4+((phen.tiling[(row-1)*phen.dx+column]-1)/4) *4);
      else
        swapper.emplace_back(0);                            
    }
  std::swap(phen.dx,phen.dy);
  phen.tiling=swapper;
}
void ClockwisePiRotation(Phenotype& phen) {
  std::reverse(phen.tiling.begin(),phen.tiling.end());
  //for(auto& base : phen.tiling)
  //  base=((base-1)+2)%4+1;
}

void MinimalTilingRepresentation(std::vector<uint8_t>& tiling) {
  if(tiling.size()==1)
    return;
  for(uint8_t& t:tiling)
    t+=128*(t!=0);
  uint8_t swap_count=1;
  std::set<std::pair<uint8_t,uint8_t> > swapped_already;
  //std::unordered_map<uint8_t,uint8_t> xfer;
  //std::
  bool prin=false;
 
  if(tiling.size()==0 && std::count(tiling.begin(),tiling.end(),0)==0) {// &&std::accumulate(tiling.begin(),tiling.end(),0)>5) {
    std::cout<<"pre ";
    for(auto x:tiling)
      std::cout<<+x<<" ";
    std::cout<<std::endl;
    prin=true;
  }
  
  for(std::vector<uint8_t>::iterator t_iter=std::find_if(tiling.begin(),tiling.end(),[](const int s) { return s>0; });t_iter!=tiling.end();) {
    /*
      if(*t_iter-(*t_iter-1)%4==swap_count) {
	if(*t_iter==swap_count) { //exact replica, just skip
	  if(prin)
	    std::cout<<"not swapping "<<+swap_count<<std::endl;
	  swap_count+=4;
	
	}
	else {
	  const uint8_t static_swap=*t_iter;
	  const uint8_t delta_theta=(*t_iter-swap_count+4)%4;
	  for(uint8_t cyclic = 0; cyclic<2; ++cyclic) {
	    uint8_t pre_swap=(static_swap-(static_swap-1)%4)+((static_swap-1)%4+cyclic)%4;
	    if(prin)
	      std::cout<<"swapping same"<<+pre_swap<<" and "<<+swap_count<<std::endl;
	    std::replace(tiling.begin(),tiling.end(),swap_count,uint8_t(255));
	    std::replace(tiling.begin(),tiling.end(),pre_swap,swap_count);
	    std::replace(tiling.begin(),tiling.end(),uint8_t(255),pre_swap);
	    ++swap_count;
	  }
	  swap_count=+2;
	}
      }
      */
    //else {
	const uint8_t static_swap=*t_iter;
	for(uint8_t cyclic = 0; cyclic<4; ++cyclic) {
	  uint8_t pre_swap=(static_swap-(static_swap-1)%4)+((static_swap-1)%4+cyclic)%4;
	  if(prin)
	    std::cout<<"swapping other"<<+pre_swap<<" and "<<+swap_count<<std::endl;
	  std::replace(tiling.begin(),tiling.end(),swap_count,uint8_t(255));
	  std::replace(tiling.begin(),tiling.end(),pre_swap,swap_count);
	  std::replace(tiling.begin(),tiling.end(),uint8_t(255),pre_swap);
	  ++swap_count;
	}
	//}
    
      t_iter=std::find_if(tiling.begin(),tiling.end(),[swap_count](const int s) { return (s>0 && s-(s-1)%4 >= swap_count); });
    
      
  }
  //return;
  if(prin) {
  std::cout<<"post ";
  for(auto x:tiling)
    std::cout<<+x<<" ";
  std::cout<<std::endl;
  }
}
/*
  for(uint8_t index=0;index<tiling.size();++index) {
   
    const uint8_t current_tile = tiling[index];
    if(current_tile && (swapped_already.count(std::make_pair(current_tile-(current_tile-1)%4,swap_count))+swapped_already.count(std::make_pair(swap_count,current_tile-(current_tile-1)%4)))==0) {

      swapped_already.insert(std::make_pair(current_tile-(current_tile-1)%4,swap_count));
      
      for(uint8_t cyclic=0;cyclic<4;++cyclic) { //
	uint8_t pre_swap=(current_tile-(current_tile-1)%4)+((current_tile-1)%4+cyclic)%4;
	if((pre_swap-1)/4==(swap_count-1)/4 && cyclic>1) {
	  swap_count+=2;
	  break;
	}
	/*
	if(std::find(tiling.begin(),tiling.end(),pre_swap)==tiling.end()) {
	  continue;
	  ++swap_count;
	}
	if(pre_swap>swap_count && (std::find(tiling.begin(),tiling.end(),pre_swap)<std::find(tiling.begin(),tiling.end(),swap_count))) {
	  std::cout<<"tripped"<<std::endl;
	  //std::cout<<(std::find(tiling.begin(),tiling.end(),pre_swap)<std::find(tiling.begin(),tiling.end(),swap_count))<<std::endl;
	  ++swap_count;
	  continue;
	}
	
	if(prin)
	  std::cout<<"swapping "<<+pre_swap<<" and "<<+swap_count<<std::endl;
	std::replace(tiling.begin(),tiling.end(),swap_count,uint8_t(255));
       	std::replace(tiling.begin(),tiling.end(),pre_swap,swap_count);
	std::replace(tiling.begin(),tiling.end(),uint8_t(255),pre_swap);
	++swap_count;
	//break;
      }
      
      //xfer[(t-(t-1)%4)]=xfer.size()+1;
    }
  }
  if(prin) {
    std::cout<<"post ";
for(auto x:tiling)
      std::cout<<+x<<" ";
    std::cout<<std::endl;
  }
}
*/

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

void InterfaceStrengths(std::vector<interface_type>& interfaces, std::vector<uint32_t>& strengths) {
  for(std::vector<interface_type>::const_iterator outer_face=interfaces.begin(); outer_face!=interfaces.end(); ++outer_face) {
    ++strengths[static_cast<uint8_t>(1.5*model_params::interface_size)+1-interface_model::SammingDistance(*outer_face,*outer_face)/2];
    for(std::vector<interface_type>::const_iterator inner_face=outer_face+1; inner_face!=interfaces.end(); ++inner_face) 
      ++strengths[model_params::interface_size-interface_model::SammingDistance(*outer_face,*inner_face)];
  }
}




