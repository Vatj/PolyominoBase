#include <stochastic_model.hpp>

namespace Stochastic
{
  std::mt19937 RNG_Generator(std::random_device{}());
  
  Phenotype_ID Analyse_Genotype_Outcome(std::vector<uint8_t> genome, uint8_t N_Repeated_Checks, StochasticPhenotypeTable* pt,uint8_t seed) {
    const uint8_t THRESHOLD_SIZE=(genome.size()*genome.size())/4;
    std::vector<int8_t> Placed_Tiles_Check=Stochastic_Polyomino_Builder(genome,THRESHOLD_SIZE,seed),Placed_Tiles_Compare;
    
    if(Placed_Tiles_Check.empty())
      return std::make_pair(0,0);
    if(Placed_Tiles_Check.size()/4 > THRESHOLD_SIZE)
      return std::make_pair(255,0);
    
    Phenotype phen1=Generate_Spatial_Occupancy(Placed_Tiles_Check,3),phen2;
    
    for(uint8_t nth_repeat=1;nth_repeat<N_Repeated_Checks;++nth_repeat) {
      Placed_Tiles_Compare=Stochastic_Polyomino_Builder(genome,THRESHOLD_SIZE,seed);
      if(Placed_Tiles_Compare.size()/4 > THRESHOLD_SIZE)
        return std::make_pair(255,0);
      if(Placed_Tiles_Compare.empty() || Placed_Tiles_Check.size()!=Placed_Tiles_Compare.size())
        return std::make_pair(0,0);
      phen2=Generate_Spatial_Occupancy(Placed_Tiles_Compare,3);
      if(!ComparePolyominoes(phen1,phen2))
        return std::make_pair(0,0);
      if(N_Repeated_Checks-nth_repeat>1) {
	phen1=phen2;
        Placed_Tiles_Check=Placed_Tiles_Compare;
      }
    }
    Phenotype_ID result;
    #pragma omp critical (table_lookup)
    {
      result=pt->GetPhenotypeID(phen2); //Placed_Tiles_Check.size()/4;
    }
    return result;
  }
  
  std::vector<int8_t> Stochastic_Polyomino_Builder(std::vector<uint8_t> genome, uint8_t THRESHOLD_SIZE, uint8_t initial_Tile) {
    std::vector<int8_t> Placed_Tiles{0,0,static_cast<int8_t>(initial_Tile),0},Interacting_Faces; //DEFINED AS (X,Y,Tile Type Number, Tile Rotation[in CW rotation])
    for(uint8_t face=0; face<4; ++face)
      if(genome[initial_Tile*4+face]!=0)
        Interacting_Adjacency(Interacting_Faces,genome[initial_Tile*4+face],face,0,0);
        
    while(!Interacting_Faces.empty()) {  
      for(int number_possibilities=(Interacting_Faces.size()/4)-1; number_possibilities >=0; --number_possibilities) {
        bool collision_Detected=false;
        for(std::vector<int8_t>::iterator occupied_iter = Placed_Tiles.begin(); occupied_iter!=Placed_Tiles.end(); occupied_iter+=4) {
          if(Interacting_Faces[number_possibilities*4]==*occupied_iter && Interacting_Faces[number_possibilities*4+1]==*(occupied_iter+1)) {
            if(Interacting_Faces[number_possibilities*4+2]!=genome[*(occupied_iter+2)*4+(Interacting_Faces[number_possibilities*4+3]-*(occupied_iter+3)+4)%4])
              return {}; //Steric mismatch, reject empty vector as sign
            Interacting_Faces.erase(Interacting_Faces.end()-4,Interacting_Faces.end()); 
            collision_Detected=true;
            break;
          }
        }
        if(collision_Detected)
          continue;
        else {
          uint8_t conjugate_count=std::count(genome.begin(),genome.end(),Interacting_Faces[number_possibilities*4+2]);
          std::uniform_int_distribution<uint8_t> Random_Count(0,conjugate_count-1);
          int nth_conjugate=Random_Count(RNG_Generator);
          std::vector<uint8_t>::iterator current_conjugate=std::find(genome.begin(),genome.end(),Interacting_Faces[number_possibilities*4+2]);
          for(int conj_cnt=1;conj_cnt<=nth_conjugate;++conj_cnt) {
            current_conjugate=std::find(current_conjugate+1,genome.end(),Interacting_Faces[number_possibilities*4+2]);
          }
          uint8_t new_Tile=(current_conjugate-genome.begin())/4,new_Face=(current_conjugate-genome.begin())%4,rotation= (Interacting_Faces[number_possibilities*4+3]-new_Face+4)%4;

          Placed_Tiles.insert(Placed_Tiles.end(),{Interacting_Faces[number_possibilities*4],Interacting_Faces[number_possibilities*4+1],static_cast<int8_t>(new_Tile),static_cast<int8_t>(rotation)});
          if(Placed_Tiles.size()/4 > THRESHOLD_SIZE) 
            return Placed_Tiles;
          uint8_t placed_X=Interacting_Faces[number_possibilities*4],placed_Y=Interacting_Faces[number_possibilities*4+1];
          Interacting_Faces.erase(Interacting_Faces.begin()+number_possibilities*4,Interacting_Faces.begin()+number_possibilities*4+4);
          for(uint8_t face=1;face<4;++face) {
            uint8_t temp_Face=genome[new_Tile*4+(new_Face+face)%4];
            if(temp_Face!=0)
              Interacting_Adjacency(Interacting_Faces,temp_Face,(new_Face+face+rotation)%4,placed_X,placed_Y);
          }
          break;
        }
      }
    }
    return Placed_Tiles;
  }

  void Interacting_Adjacency(std::vector<int8_t>& Interacting_Faces, uint8_t interacting_Face, uint8_t face_index, int8_t X, int8_t Y) {
    int8_t X_OFFSET=0,Y_OFFSET=0;
    switch(face_index) {
    case 0:
      X_OFFSET=0;Y_OFFSET=1;
      break;
    case 1:
      X_OFFSET=1;Y_OFFSET=0;
      break;
    case 2:
      X_OFFSET=0;Y_OFFSET=-1;
      break;
    case 3:
      X_OFFSET=-1;Y_OFFSET=0;
      break;
    }
    uint8_t conjugate_Face=Interaction_Matrix(interacting_Face);
    Interacting_Faces.insert(Interacting_Faces.begin()+std::uniform_int_distribution<uint8_t>{0,Interacting_Faces.size()/4}(RNG_Generator)*4,{static_cast<int8_t>(X+X_OFFSET),static_cast<int8_t>(Y+Y_OFFSET),static_cast<int8_t>(conjugate_Face),static_cast<int8_t>((face_index+2)%4)});
  }

  
  Phenotype Generate_Spatial_Occupancy(std::vector<int8_t>& Placed_Tiles_Check, uint8_t generate_mode) {
    std::vector<int8_t> X_Locs_Check, Y_Locs_Check;
    std::vector<uint8_t>  Tile_Type_Check,Tile_Orientation_Check;
    for(std::vector<int8_t>::iterator check_iter = Placed_Tiles_Check.begin();check_iter!=Placed_Tiles_Check.end();check_iter+=4) {
      X_Locs_Check.emplace_back(*check_iter);
      Y_Locs_Check.emplace_back(*(check_iter+1));
      Tile_Type_Check.emplace_back(*(check_iter+2));
      Tile_Orientation_Check.emplace_back(*(check_iter+3));
    }
    std::vector<int8_t>::iterator lx,rx,ty,by;
    std::tie(lx,rx)=std::minmax_element(X_Locs_Check.begin(),X_Locs_Check.end());
    std::tie(by,ty)=std::minmax_element(Y_Locs_Check.begin(),Y_Locs_Check.end());
    uint8_t dx=*rx-*lx+1;
    uint8_t dy=*ty-*by+1;
    std::vector<uint8_t> Spatial_Occupancy_Check(dx*dy);

    for(uint8_t tileIndex=0;tileIndex<X_Locs_Check.size();++tileIndex) {
      if(generate_mode==0) 
        Spatial_Occupancy_Check[(*ty-Y_Locs_Check[tileIndex])*dx + (X_Locs_Check[tileIndex]-*lx)]=1;
      if(generate_mode==3) 
        Spatial_Occupancy_Check[(*ty-Y_Locs_Check[tileIndex])*dx + (X_Locs_Check[tileIndex]-*lx)]=1+Tile_Type_Check[tileIndex]*4+Tile_Orientation_Check[tileIndex];
    }
    Phenotype phen{dx,dy,Spatial_Occupancy_Check};
    if(dy>dx)
      ClockwiseRotation(phen);
    MinimalTilingRepresentation(phen.tiling);
    return phen;
  }

}

