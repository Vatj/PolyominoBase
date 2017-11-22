#include "interface_model.hpp"

namespace interface_model
{
  //std::random_device rd;
  xorshift RNG_Engine(124124124);//rd());

  
  const static double FAIL_THRESHOLD=0.75;
  
  template<typename T> inline T reverse_bits(T v) {
    T s = sizeof(v) * 8; // bit size; must be power of 2
    T mask = ~0;         
    while ((s >>= 1) > 0) 
      {
        mask ^= (mask << s);
        v = ((v >> s) & mask) | ((v << s) & ~mask);
      }
    return v;
  }

  uint8_t SammingDistance(uint16_t face1,uint16_t face2) {
    return 16-__builtin_popcount(face1 ^ reverse_bits<uint16_t>(face2));
  }

  void MutateInterfaces(std::vector<uint16_t>& binary_genome,double mu_prob) {
    std::vector<uint8_t> interface_indices{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    std::binomial_distribution<uint8_t> b_dist(interface_indices.size(),mu_prob);
    for(uint16_t& base : binary_genome) {
      std::shuffle(interface_indices.begin(), interface_indices.end(), RNG_Engine);
      uint8_t num_mutations= b_dist(RNG_Engine);
      for(uint8_t nth=0;nth<num_mutations;++nth)
        base ^= (1U << interface_indices[nth]);
    }
  }


  int ProteinAssemblyOutcome(std::vector<uint16_t> binary_genome,uint8_t N_repeats) {
    uint8_t dx,dy,dx_prime,dy_prime;
    double temperature=1;
    std::vector<int8_t> assembly_information=AssembleProtein(binary_genome,temperature);
    std::cout<<assembly_information.size()<<std::endl;
    for(int8_t a:assembly_information)
      std::cout<<+a<<" ";
    std::cout<<std::endl;
    
    
    std::vector<int8_t> assembly_information_prime;
    uint8_t UND_failures=assembly_information.size()==0 ? 1 : 0;
    std::vector<uint8_t> spatial_information;
    std::vector<uint8_t> spatial_information_prime;
    if(!assembly_information.empty()) {
    spatial_information=SpatialGrid(assembly_information,dx,dy);

    std::cout<<"size "<<  assembly_information.size()<<"dx "<<+dx<<", dy "<<+dy<<std::endl;
    for(uint8_t y=0;y<dy;++y) {
      for(uint8_t x=0;x<dx;++x) {
        std::cout<<+spatial_information[y*dx+x]<<" ";
      }
      std::cout<<std::endl;
    }
    std::cout<<std::endl;
    }
    

    for(uint8_t nth=0;nth<N_repeats;++nth) {
      assembly_information_prime=AssembleProtein(binary_genome,temperature);
      std::cout<<assembly_information_prime.size()<<std::endl;
      if(assembly_information_prime.size()==0) {
        ++UND_failures;
        std::cout<<"skipping"<<std::endl;
        continue;
      }
      else {
      spatial_information_prime=SpatialGrid(assembly_information_prime,dx_prime,dy_prime);
      std::cout<<"size "<<  assembly_information_prime.size()<<"dx "<<+dx_prime<<", dy "<<+dy_prime<<std::endl;
      for(uint8_t y=0;y<dy_prime;++y) {
        for(uint8_t x=0;x<dx_prime;++x) {
          std::cout<<+spatial_information_prime[y*dx_prime+x]<<" ";
        }
        std::cout<<std::endl;
      }
      std::cout<<std::endl;
    
      //if(!ComparePolyominoes(spatial_information,dx,dy,spatial_information_prime,dx_prime,dy_prime))
      //  ++UND_failures;      
      }
      if(UND_failures*1./N_repeats>FAIL_THRESHOLD)
        return -1;
      spatial_information=spatial_information_prime;
      dx=dx_prime;
      dy=dy_prime;
    }
    return UND_failures*1./N_repeats>FAIL_THRESHOLD ? -1 : 1;
    
  }
  

  std::vector<int8_t> AssembleProtein(const std::vector<uint16_t>& binary_genome,double temperature) {
    std::vector<uint8_t> faces{0,1,2,3};
    std::vector<uint8_t> tile_types(binary_genome.size()/(4.));
    std::iota(tile_types.begin(), tile_types.end(), 0);
    std::vector<int8_t> placed_tiles{0,0}; //(x,y,tile type, orientation)
    std::vector<int8_t> growing_perimeter; //(x,y,direction,tile type, orientation)
    PerimeterGrowth(0,0,0,-1,0,growing_perimeter,placed_tiles);   
    while(!growing_perimeter.empty()) {
      
      PlaceNextUnit(binary_genome,tile_types,faces,growing_perimeter,placed_tiles,temperature);  
      if(placed_tiles.size()>(binary_genome.size()*binary_genome.size()/4.)) {
        placed_tiles.clear();
        break;
      }
    }
    return placed_tiles;   
  }
    
  void PlaceNextUnit(const std::vector<uint16_t>& binary_genome,std::vector<uint8_t>& tile_types,std::vector<uint8_t>& faces,std::vector<int8_t>& growing_perimeter,std::vector<int8_t>& placed_tiles,double temperature) {
    uint8_t current_orientation=growing_perimeter.back();growing_perimeter.pop_back();
    uint8_t current_tile=growing_perimeter.back();growing_perimeter.pop_back();
    uint8_t current_direction=growing_perimeter.back();growing_perimeter.pop_back();
    int8_t current_y=growing_perimeter.back();growing_perimeter.pop_back();
    int8_t current_x=growing_perimeter.back();growing_perimeter.pop_back();
    std::shuffle(tile_types.begin(), tile_types.end(), RNG_Engine);
    std::shuffle(faces.begin(), faces.end(), RNG_Engine);
    std::uniform_real_distribution<double> uniform(0, 1);

    //std::cout<<"looking to place on "<<+current_x<<","<<+current_y<<std::endl;
    for(uint8_t tile : tile_types) {
      for(uint8_t face : faces) {
        uint8_t samming_energy=SammingDistance(binary_genome[current_tile*4+current_orientation],binary_genome[tile*4+face]);
        //std::cout<<+samming_energy<<", "<<+binary_genome[current_tile*4+current_orientation]<<" ^ "<<+binary_genome[tile*4+face]<<std::endl;
        //double d=uniform(RNG_Engine);
        if(uniform(RNG_Engine)<std::exp(-1*samming_energy/temperature)) {
          
          //if(d>=std::exp(-1*samming_energy/temperature))
          if(samming_energy!=0)
            std::cout<<std::exp(-1*samming_energy/temperature)<<std::endl;
          //std::cout<<"placed one "<<+current_x<<","<<+current_y<<","<<+tile<<std::endl;
          placed_tiles.insert(placed_tiles.end(),{current_x,current_y});//,(4+current_direction-face)%4});
          PerimeterGrowth(current_x,current_y,(4+current_direction-face)%4,current_direction,tile,growing_perimeter,placed_tiles);
          return;
        }

      }
    }
  }

  void PerimeterGrowth(int8_t x,int8_t y,uint8_t theta,uint8_t direction, uint8_t tile_type,std::vector<int8_t>& growing_perimeter,std::vector<int8_t>& placed_tiles) {
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
      if(!occupied_site) {
        //std::cout<<"adding tile for "<<+static_cast<int8_t>(x+dx)<<","<<+static_cast<int8_t>(y+dy)<<std::endl;
        growing_perimeter.insert(growing_perimeter.begin()+5*index_randomizer(RNG_Engine),{static_cast<int8_t>(x+dx),static_cast<int8_t>(y+dy),static_cast<int8_t>((f+2)%4),static_cast<int8_t>(tile_type),static_cast<int8_t>((4+f-theta)%4)});
      }
      
    }
  }
  
  std::vector<uint8_t> SpatialGrid(std::vector<int8_t>& placed_tiles, uint8_t& dx,uint8_t& dy) {
    std::vector<int8_t> x_locs, y_locs;
    x_locs.reserve(placed_tiles.size()/2);
    y_locs.reserve(placed_tiles.size()/2);
    for(std::vector<int8_t>::iterator check_iter = placed_tiles.begin();check_iter!=placed_tiles.end();check_iter+=2) {
      x_locs.emplace_back(*check_iter);
      y_locs.emplace_back(*(check_iter+1));
    }
    
    std::vector<int8_t>::iterator LEFT_X_Check,RIGHT_X_Check,TOP_Y_Check,BOTTOM_Y_Check;
    std::tie(LEFT_X_Check,RIGHT_X_Check)=std::minmax_element(x_locs.begin(),x_locs.end());
    std::tie(BOTTOM_Y_Check,TOP_Y_Check)=std::minmax_element(y_locs.begin(),y_locs.end());
    dx=*RIGHT_X_Check-*LEFT_X_Check+1;
    dy=*TOP_Y_Check-*BOTTOM_Y_Check+1;
    std::vector<uint8_t> spatial_grid(dx*dy); 
    for(uint16_t tileIndex=0;tileIndex<x_locs.size();++tileIndex)
      spatial_grid[(*TOP_Y_Check-y_locs[tileIndex])*dx + (x_locs[tileIndex]-*LEFT_X_Check)]=1;
    return spatial_grid;
  }

  bool ComparePolyominoes(std::vector<uint8_t>& Spatial_Occupation_Check,uint8_t& Delta_X_Check,uint8_t& Delta_Y_Check, std::vector<uint8_t>& Spatial_Occupation_Compare,uint8_t& Delta_X_Compare,uint8_t& Delta_Y_Compare) {
    if(std::accumulate(Spatial_Occupation_Check.begin(),Spatial_Occupation_Check.end(),0)!=std::accumulate(Spatial_Occupation_Compare.begin(),Spatial_Occupation_Compare.end(),0))
      return false;
    if(Delta_X_Check==Delta_X_Compare && Delta_Y_Check==Delta_Y_Compare && Delta_X_Check==Delta_Y_Check) { //bounding boxes match, symmetric
      if(Spatial_Occupation_Check==Spatial_Occupation_Compare) 
        return true;
      else {
        for(int rotation=0;rotation<3;++rotation) {
          ClockwiseRotation(Spatial_Occupation_Check,Delta_X_Check,Delta_Y_Check);
          if(Spatial_Occupation_Check==Spatial_Occupation_Compare)
            return true;
        }
        return false;
      }
    }
    if(Delta_X_Check==Delta_X_Compare && Delta_Y_Check==Delta_Y_Compare)  //bounding boxes match, asymmetric
      return Spatial_Occupation_Check==Spatial_Occupation_Compare ? true : Spatial_Occupation_Check==ClockwisePiRotation(Spatial_Occupation_Compare);
    if(Delta_X_Check==Delta_Y_Compare && Delta_Y_Check==Delta_X_Compare) { //bounding boxes pi/2 off, asymmetric
      ClockwiseRotation(Spatial_Occupation_Check,Delta_X_Check,Delta_Y_Check);
      return Spatial_Occupation_Check==Spatial_Occupation_Compare ? true : Spatial_Occupation_Check==ClockwisePiRotation(Spatial_Occupation_Compare);
    }
    return false;
  }
  void ClockwiseRotation(std::vector<uint8_t>& Spatial_Occupation,uint8_t& dx,uint8_t& dy) {
    std::vector<uint8_t> swapper;
    swapper.reserve(Spatial_Occupation.size());
    for(uint8_t column=0;column<dx;++column) 
      for(uint8_t row=dy;row!=0;--row) 
        swapper.emplace_back(Spatial_Occupation[(row-1)*dx+column]);
    std::swap(dx,dy);
    Spatial_Occupation=swapper;
 }
  std::vector<uint8_t> ClockwisePiRotation(std::vector<uint8_t>& spatial_occupation) {
    std::reverse(spatial_occupation.begin(),spatial_occupation.end());
    return spatial_occupation;
  }
  
}//end interface_model namespace

int main(int argc, char* argv[]) {
  std::cout<<"hi"<<std::endl;
  std::vector<uint16_t> g{0,0,0,0, 65535,31,3,31, 31,31,100,16383, 31,31,31,55807};
  // 0 binds with 65535
  // 3 binds with 16383
  // 100 binds with 55807
  // 31 is fairly neutral
  std::cout<<interface_model::ProteinAssemblyOutcome(g,10)<<std::endl;
}
