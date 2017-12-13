#include "interface_model.hpp"




namespace model_params
{
  //HARD CODED TO MATCH typedef//
  uint8_t interface_size=16;
  
  double temperature=1,mu_prob=0.2,unbound_factor=2,misbinding_rate=0,fitness_factor=2;
 
  std::binomial_distribution<uint8_t> b_dist(interface_size,mu_prob);
  std::uniform_real_distribution<double> real_dist(0, 1);

}

namespace interface_model
{
  xorshift RNG_Engine;  
  
  inline interface_type reverse_bits(interface_type v) {
    interface_type s = sizeof(v) * 8; // bit size; must be power of 2
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
  

  double SammingDistance(interface_type face1,interface_type face2) {
    //uint8_t x =model_params::interface_size-__builtin_popcount(face1 ^ reverse_bits(face2));
    //uint8_t y =__builtin_popcount(~face1 ^ reverse_bits(face2));
    //std::cout<<+x<<" vs "<<std::endl;
    return static_cast<double>(ArbitraryPopcount(face1 ^ reverse_bits(~face2)))/model_params::interface_size;
    
    //return 16-__builtin_popcount(face1 ^ reverse_bits(face2));
  }

  double SymmetryFactor(interface_type face1) {
    return static_cast<double>(ArbitraryPopcount((face1) ^ reverse_bits(face1)))/model_params::interface_size;
  }

  void MutateInterfaces(std::vector<interface_type>& binary_genome) {
    std::vector<uint8_t> interface_indices{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    for(interface_type& base : binary_genome) {
      std::shuffle(interface_indices.begin(), interface_indices.end(), RNG_Engine);
      for(uint8_t nth=0;nth<model_params::b_dist(RNG_Engine);++nth)
        base ^= (1U << interface_indices[nth]);
    }
  }


  double ProteinAssemblyOutcome(std::vector<interface_type> binary_genome,uint8_t N_repeats,PhenotypeTable* pt) {
    uint8_t dx,dy;
    std::vector<int8_t> assembly_information;//=AssembleProtein(binary_genome);//,assembly_information_prime;
    std::vector<uint8_t> spatial_information;//spatial_information_prime;
    std::vector<std::pair<uint8_t,uint32_t> > phenotype_IDs;phenotype_IDs.reserve(N_repeats);

    //spatial_information=SpatialGrid(assembly_information,dx,dy);
    //phenotype_IDs.emplace_back(std::accumulate(spatial_information.begin(),spatial_information.end(),0),pt->PhenotypeCheck(spatial_information,dx,dy));
    
    for(uint8_t nth=0;nth<N_repeats;++nth) {
      assembly_information=AssembleProtein(binary_genome);
      spatial_information=SpatialGrid(assembly_information,dx,dy);
      phenotype_IDs.emplace_back(std::accumulate(spatial_information.begin(),spatial_information.end(),0),pt->PhenotypeCheck(spatial_information,dx,dy));
    }
    return pt->GenotypeFitness(phenotype_IDs);
    
  }
  

  std::vector<int8_t> AssembleProtein(const std::vector<interface_type>& binary_genome) {
    std::vector<uint8_t> tile_types(binary_genome.size()/(4.));
    std::iota(tile_types.begin(), tile_types.end(), 0);
    std::vector<int8_t> placed_tiles{0,0}; //(x,y,tile type, orientation)
    std::vector<int8_t> growing_perimeter; //(x,y,direction,tile type, orientation)
    PerimeterGrowth(0,0,0,-1,0,growing_perimeter,placed_tiles);
    int8_t current_orientation, current_tile, current_direction, current_x, current_y;
    std::vector<uint8_t> faces{0,1,2,3};
    
    while(!growing_perimeter.empty()) {
      current_orientation=growing_perimeter.back();growing_perimeter.pop_back();
      current_tile=growing_perimeter.back();growing_perimeter.pop_back();
      current_direction=growing_perimeter.back();growing_perimeter.pop_back();
      current_y=growing_perimeter.back();growing_perimeter.pop_back();
      current_x=growing_perimeter.back();growing_perimeter.pop_back();
      std::shuffle(faces.begin(), faces.end(), RNG_Engine);
      for(uint8_t tile : tile_types) {
        for(uint8_t face : faces) {
          if(model_params::real_dist(RNG_Engine)<std::exp(-model_params::misbinding_rate-1*SammingDistance(binary_genome[current_tile*4+current_orientation],binary_genome[tile*4+face])/model_params::temperature)) {
            placed_tiles.insert(placed_tiles.end(),{current_x,current_y});
            PerimeterGrowth(current_x,current_y,(4+current_direction-face)%4,current_direction,tile,growing_perimeter,placed_tiles);
            goto endplacing;
          }
        }
      }
      endplacing:      
      if(placed_tiles.size()>(model_params::unbound_factor*binary_genome.size()*binary_genome.size()/2.)) {
        placed_tiles.clear();
        break;
      }
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
      else {
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
  
}//end interface_model namespace


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

std::vector<double> InterfaceStrengths(std::vector<interface_model::interface_type>& interfaces) {
  std::vector<double> strengths;
  strengths.reserve(interfaces.size()*(interfaces.size()+1)/2);
  for(std::vector<interface_model::interface_type>::const_iterator outer_face=interfaces.begin(); outer_face!=interfaces.end(); ++outer_face) {
    for(std::vector<interface_model::interface_type>::const_iterator inner_face=outer_face; inner_face!=interfaces.end(); ++inner_face) {
      strengths.emplace_back(1-interface_model::SammingDistance(*outer_face,*inner_face));
      //std::cout<<+*outer_face<<" "<<+*inner_face<<" = "<< 1-static_cast<double>(interface_model::SammingDistance(*outer_face,*inner_face))/model_params::interface_size<<std::endl;
    }
  }
  return strengths;
}

